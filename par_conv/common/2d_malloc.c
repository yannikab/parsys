/* 
 * File:   2d_malloc.c
 * Author: John
 *
 * Created on 21 Ιανουάριος 2015, 10:21 πμ
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>

#include "2d_malloc.h"

bool alloc_uchar_array(unsigned char ***array, int rows, int columns, int channels)
{
	unsigned char *p;
	p = malloc(rows * columns * channels * sizeof (unsigned char));
	if (p == NULL)
	{
		perror("malloc");
		return false;
	}

	(*array) = malloc(rows * sizeof (unsigned char *));
	if ((*array) == NULL)
	{
		perror("malloc");
		free(p);
		return false;
	}

	int i;
	for (i = 0; i < rows; i++)
		(*array)[i] = &(p[i * columns * channels]);

	return true;
}

void dealloc_uchar_array(unsigned char ***array)
{
	free(&((*array)[0][0]));
	free(*array);
	*array = NULL;
}

bool alloc_float_array(float ***array, int rows, int columns, int channels)
{
	float *p;
	p = malloc(rows * columns * channels * sizeof (float));
	if (p == NULL)
	{
		perror("malloc");
		return false;
	}

	(*array) = malloc(rows * sizeof (float *));
	if ((*array) == NULL)
	{
		perror("malloc");
		free(p);
		return false;
	}

	int i;
	for (i = 0; i < rows; i++)
		(*array)[i] = &(p[i * columns * channels]);

	return true;
}

void dealloc_float_array(float ***array)
{
	free(&((*array)[0][0]));
	free(*array);
	*array = NULL;
}
