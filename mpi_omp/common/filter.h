/* 
 * File:   filter.h
 * Author: John
 *
 * Created on January 21, 2015, 10:54 AM
 */

#ifndef FILTER_H
#define	FILTER_H

#include <stdbool.h>

#include "../settings.h"

#define B 1

static const float filter[3][3] = {
    {0.0625f, 0.125f, 0.0625f},
    {0.1250f, 0.250f, 0.1250f},
    {0.0625f, 0.125f, 0.0625f},
};

void apply_inner_filter(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width);
void apply_outer_filter(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width);

void apply_inner_filter_openmp(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width);

bool images_identical(float (**output_image)[CHANNELS], float (**input_image)[CHANNELS], int height, int width);

#endif	/* FILTER_H */
