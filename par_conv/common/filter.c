/* 
 * File:   filter.c
 * Author: John
 *
 * Created on January 21, 2015, 10:53 AM
 */

#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include "filter.h"

void apply_inner_filter(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width)
{
    unsigned int i, j, c;
    int p, q;

    for (i = 2 * B; i < height - 2 * B; i++)
    {
        for (j = 2 * B; j < width - 2 * B; j++)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }
}

void apply_outer_filter(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width)
{
    unsigned int i, j, c;
    int p, q;

    /* Left border column. */

    for (i = B; i < height - B; i++)
    {
        for (j = B; j < 2 * B; j++)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }

    /* Right border column. */

    for (i = B; i < height - B; i++)
    {
        for (j = width - 1 - B; j > width - 1 - 2 * B; j--)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }

    /* Top border row, avoid recalculating overlap with columns. */

    for (j = 2 * B; j < width - 2 * B; j++)
    {
        for (i = B; i < 2 * B; i++)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }

    /* Bottom border row, avoid recalculating overlap with columns. */

    for (j = 2 * B; j < width - 2 * B; j++)
    {
        for (i = height - 1 - B; i > height - 1 - 2 * B; i--)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }
}

void apply_inner_filter_openmp(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width)
{
    unsigned int i, j, c;
    int p, q;

#pragma omp for
    for (i = 2 * B; i < height - 2 * B; i++)
    {
        for (j = 2 * B; j < width - 2 * B; j++)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                float value = 0.0f;

                for (p = -B; p <= B; p++)
                    for (q = -B; q <= B; q++)
                        value += input_image[i - p][j - q][c] * filter[p + B][q + B];

                output_image[i][j][c] = value;
            }
        }
    }
}

bool images_identical(float (**first_image)[CHANNELS], float (**second_image)[CHANNELS], int height, int width)
{
    unsigned int i, j, c;

    bool identical = true;

    for (i = 2 * B; i < height - 2 * B; i++)
    {
        for (j = 2 * B; j < width - 2 * B; j++)
        {
            for (c = 0; c < CHANNELS; c++)
            {
                if ((int) first_image[i][j][c] != (int) second_image[i][j][c])
                {
                    identical = false;
                    break;
                }
            }
        }
    }

    return identical;
}
