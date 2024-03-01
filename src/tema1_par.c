// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Strcture in order to keep the arguments of a thread
typedef struct {
    long id;
    ppm_image *image;
    ppm_image *scaled_image;
    unsigned char **grid;
    ppm_image **contour_map; 
    int number_of_threads;
    pthread_barrier_t *bar;
} thread_arguments;

// Function to determinates the minimum of two elements
int min(int a, int b) {
    if(a < b) {
        return a;
    } else {
        return b;
    }
}

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Initialise the grid for the image, using the coorditanes from rescaled one
unsigned char **init_grid(int p, int q) {
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    return grid;
}

// Initialise the memory for the rescaled image just if necessary
ppm_image *init_rescale(ppm_image *image) {

    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    return new_image;
}

// Rescale the image just in case the size will be decreased
void rescale_image(thread_arguments *arg) {

    if (arg->image->x == arg->scaled_image->x && arg->image->y == arg->scaled_image->y) {
        return;
    }

    uint8_t sample[3];

    int start = arg->id * (double)arg->scaled_image->x / arg->number_of_threads;
    int end = min((arg->id + 1) * (double)arg->scaled_image->x / arg->number_of_threads, arg->scaled_image->x);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < arg->scaled_image->y; j++) {
            float u = (float)i / (float)(arg->scaled_image->x - 1);
            float v = (float)j / (float)(arg->scaled_image->y - 1);
            sample_bicubic(arg->image, u, v, sample);

            arg->scaled_image->data[i * arg->scaled_image->y + j].red = sample[0];
            arg->scaled_image->data[i * arg->scaled_image->y + j].green = sample[1];
            arg->scaled_image->data[i * arg->scaled_image->y + j].blue = sample[2];
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample(thread_arguments *arg) {
    
    int step_x = STEP;
    int step_y = STEP;
    
    int p = arg->scaled_image->x / step_x;
    int q = arg->scaled_image->y / step_y;

    int start_p = arg->id * ((double)p / arg->number_of_threads);
    int end_p = min(p, (arg->id + 1) *(double)p / arg->number_of_threads);

    for (int i = start_p; i < end_p; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = arg->scaled_image->data[i * step_x * arg->scaled_image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                arg->grid[i][j] = 0;
            } else {
                arg->grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    // just the last thread has to do this operation
    if(arg->id == arg->number_of_threads - 1) {
        for (int i = 0; i < p; i++) {
            ppm_pixel curr_pixel = arg->scaled_image->data[i * step_x * arg->scaled_image->y + arg->scaled_image->x - 1];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                arg->grid[i][q] = 0;
            } else {
                arg->grid[i][q] = 1;
            }
        }

        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = arg->scaled_image->data[(arg->scaled_image->x - 1) * arg->scaled_image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                arg->grid[p][j] = 0;
            } else {
                arg->grid[p][j] = 1;
            }
        }
    }
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(thread_arguments *arg) {

    int step_x = STEP;
    int step_y = STEP;
    
    int p = arg->scaled_image->x / step_x;
    int q = arg->scaled_image->y / step_y;

    int start = arg->id * ((double)p / arg->number_of_threads);
    int end = min(p, (arg->id + 1) *(double)p / arg->number_of_threads);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * arg->grid[i][j] + 4 * arg->grid[i][j + 1] + 2 * arg->grid[i + 1][j + 1] + 1 * arg->grid[i + 1][j];
            
            int x = i * step_x;
            int y = j * step_y;
           
            for (int i = 0; i < arg->contour_map[k]->x; i++) {
                for (int j = 0; j < arg->contour_map[k]->y; j++) {
                    int contour_pixel_index = arg->contour_map[k]->x * i + j;
                    int image_pixel_index = (x + i) * arg->scaled_image->y + y + j;

                    arg->scaled_image->data[image_pixel_index].red = arg->contour_map[k]->data[contour_pixel_index].red;
                    arg->scaled_image->data[image_pixel_index].green = arg->contour_map[k]->data[contour_pixel_index].green;
                    arg->scaled_image->data[image_pixel_index].blue = arg->contour_map[k]->data[contour_pixel_index].blue;
                }
            }
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);

}

// Thread function
void *thread_function(void *arguments) {

    // Take the arguments
    thread_arguments* arg = (thread_arguments*) arguments;

    // 1. Rescale the image if necessary
    rescale_image(arg);

    // Wait all the threads to finish the rescale
    pthread_barrier_wait(arg->bar);

    // 2. Sample the grid
    sample(arg);

    // Wait all the threads to finish the grid
    pthread_barrier_wait(arg->bar);

    // 3. March the squares
    march(arg);

    // Exit the function
    pthread_exit(NULL);
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    // 0. Initialize and allocate the memory
    ppm_image *image = read_ppm(argv[1]);
    ppm_image **contour_map = init_contour_map();
    ppm_image *scaled_image = init_rescale(image);
    unsigned char **grid = init_grid(scaled_image->x / STEP, scaled_image->y / STEP);
    
    // The variables, the barrier and the array of structures for thread arguments
    int num_threads = atoi(argv[3]);
    pthread_t threads[num_threads];
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, num_threads);
    int r;
    long id;
    void *status;
    thread_arguments* args = malloc(num_threads * sizeof(thread_arguments));
    
    // Arguments of the thread
    for(id = 0; id < num_threads; id++) {
        args[id].id = (long) id;
        args[id].image = (ppm_image *) image;
        args[id].scaled_image = (ppm_image *) scaled_image;
        args[id].grid = (unsigned char **) grid;
        args[id].contour_map  = (ppm_image **) contour_map;
        args[id].number_of_threads = (int) num_threads;
        args[id].bar = (pthread_barrier_t *) &barrier;
    }

    // Create the threads
   for (id = 0; id < num_threads; id++) {       
        r = pthread_create(&threads[id], NULL, (void *) thread_function, &args[id]);
 
        if (r) {
            printf("Unable to create thread %ld\n", id);
            exit(-1);
        }
    }
 
    // Join the threads 
    for (id = 0; id < num_threads; id++) {
        r = pthread_join(threads[id], &status);
 
        if (r) {
            printf("Unable to wait the thread %ld\n", id);
            exit(-1);
        }
    }

    // 4. Write output
    write_ppm(scaled_image, argv[2]);

    // Free the resources
    free_resources(scaled_image, contour_map, grid, STEP);

    // Destroy the barrier
    pthread_barrier_destroy(&barrier);

    return 0;
}
