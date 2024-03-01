# Marching Squares, APD


### Structural Components

This code presents an implementation of the Marching Squares algorithm with parallelization. It involves the following structural components:

- Structure : \
``` thread_arguments ``` (This data structure is used to store thread arguments for the purpose of facilitating parsing.)

- Functions : \
``` int min(int a, int b)``` ( This function computes the minimum value of two input integers.) \
```ppm_image **init_contour_map()``` (Initializes memory for the contour map.) \
```unsigned char **init_grid(int p, int q)``` (Allocates memory for the grid.) \
``` ppm_image *init_rescale(ppm_image *image)``` (Allocates memory for the rescaled image.) \
``` void rescale_image(thread_arguments *arg)``` (Rescales an image using multiple threads.)\
```void sample(thread_arguments *arg)``` (Samples the grid using multiple threads.) \
```void march(thread_arguments *arg)``` (Applies the Marching Squares algorithm using multiple threads.) \
```void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x)``` (Releases the allocated memory resources.) \
```void *thread_function(void *arguments)``` ( This function invokes the necessary operations for each thread, using barriers to synchronize steps.) 


### Execution Workflow

The objective of this code is to parallelize an initial serial implementation of the Marching Squares algorithm to enhance its performance. The execution workflow involves the following steps: 

- Structural Organization: To enable parallelization, a structured approach was introduced, using a thread_arguments data structure to manage the thread arguments.

- Function Segmentation: The original functions were divided into two distinct parts: memory allocation and concurrent operations that can be executed by multiple threads.
Within each function, thread identifiers were employed to determine the portion of the image that each thread is responsible for.
- Sample Function Optimization: Special consideration was given to the sample function. The last column/row of the image is modified only at the very end of the operation, and as a result, these modifications are allocated to the last thread to ensure efficient resource utilization.

- Thread Synchronization: The thread_function plays an important role in arranging the execution of operations. It creates interdependencies between operations, being necessary the use of barriers to ensure that all threads complete a given step before progressing to the next one.

- Main Function Flow: In the main function, memory allocations are performed. After that, the creation and joining of threads optimize the overall code execution.


### Maslaev Andra, 333CD, AC-CTI