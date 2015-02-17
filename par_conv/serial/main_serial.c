/* 
 * File:   main_serial.c
 * Author: John
 *
 * Created on January 21, 2015, 11:22 AM
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>

#include <unistd.h>
#include <sys/times.h>

#include "../settings.h"
#include "../common/2d_malloc.h"
#include "../common/file_io.h"
#include "../common/filter.h"

/*
 * 
 */
int main_serial(int argc, char** argv)
{
    if (argc != 3)
    {
        printf("Usage: %s <iterations> <convergence>\n", argv[0]);
        return (EXIT_FAILURE);
    }

    int iterations = atoi(argv[1]);
    int convergence = atoi(argv[2]);

    printf("main_serial()\n");
    printf("Iterations: %d, Convergence: %d\n", iterations, convergence);
    printf("Threads: %d\n", 1);

    bool ok = true;

    /* Read input file into buffer. */

    unsigned char (**image_buffer)[CHANNELS];

    if (ok)
        ok = read_image((unsigned char ***) &image_buffer);

    /* Allocate memory for image data. */

    float (**prev_image)[CHANNELS];
    float (**curr_image)[CHANNELS];

    if (ok)
        ok = alloc_float_array((float ***) &prev_image, B + HEIGHT + B, B + WIDTH + B, CHANNELS);

    if (ok)
        ok = alloc_float_array((float ***) &curr_image, B + HEIGHT + B, B + WIDTH + B, CHANNELS);

    /* Convert input. */

    unsigned int i, j, c;

    if (ok)
    {
        for (i = 0; i < HEIGHT; i++)
            for (j = 0; j < WIDTH; j++)
                for (c = 0; c < CHANNELS; c++)
                    curr_image[i + B][j + B][c] = (float) image_buffer[i][j][c];
    }

    /* Start timing. */

    double t1, t2, real_time;
    struct tms tb1, tb2;
    double tickspersec = (double) sysconf(_SC_CLK_TCK);

    t1 = (double) times(&tb1);

    /* Apply filter. */

    unsigned int n;

    if (ok)
    {
        for (n = 0; (iterations == 0 || n < iterations) && (iterations != 0 || convergence != 0); n++)
        {
            /* Fill borders with outer image data. */

            // south
            for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
                for (j = B; j < B + WIDTH; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B + HEIGHT - 1][j][c];

            // north
            for (i = 0; i < B; i++) // north
                for (j = B; j < B + WIDTH; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B][j][c];

            // east
            for (i = B; i < B + HEIGHT; i++)
                for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[i][B + WIDTH - 1][c];

            // west
            for (i = B; i < B + HEIGHT; i++)
                for (j = 0; j < B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[i][B][c];

            // se
            for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
                for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B + HEIGHT - 1][B + WIDTH - 1][c];

            // nw
            for (i = 0; i < B; i++)
                for (j = 0; j < B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B][B][c];

            // sw
            for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
                for (j = 0; j < B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B + HEIGHT - 1][B][c];

            // ne
            for (i = 0; i < B; i++)
                for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
                    for (c = 0; c < CHANNELS; c++)
                        curr_image[i][j][c] = curr_image[B][B + WIDTH - 1][c];

            /* Apply filter. */

            apply_inner_filter(prev_image, curr_image, B + HEIGHT + B, B + WIDTH + B);

            apply_outer_filter(prev_image, curr_image, B + HEIGHT + B, B + WIDTH + B);

            /* Switch current / previous image buffers. */

            float (**temp)[CHANNELS];
            temp = prev_image;
            prev_image = curr_image;
            curr_image = temp;

            /* Check for convergence. */

            if (convergence > 0 && n % convergence == 0)
            {
                if (images_identical(curr_image, prev_image, B + HEIGHT + B, B + WIDTH + B))
                {
                    printf("Filter has converged after %d iterations.\n", n);
                    break;
                }
            }
        }
    }

    /* Stop time measurement, print time. */

    t2 = (double) times(&tb2);

    real_time = (double) (t2 - t1) / tickspersec;
    printf("Completed in %.3f sec\n", real_time);

    /* Convert output. */

    if (ok)
    {
        for (i = 0; i < HEIGHT; i++)
            for (j = 0; j < WIDTH; j++)
                for (c = 0; c < CHANNELS; c++)
                    image_buffer[i][j][c] = (unsigned char) curr_image[i + B][j + B][c];
    }

    /* Create output files, one for each channel. */

    write_channels(image_buffer, HEIGHT, WIDTH);

    /* Free allocated memory. */

    dealloc_uchar_array((unsigned char ***) &image_buffer);
    dealloc_float_array((float ***) &curr_image);
    dealloc_float_array((float ***) &prev_image);

    return (EXIT_SUCCESS);
}
