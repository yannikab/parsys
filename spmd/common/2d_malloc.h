/* 
 * File:   2d_malloc.h
 * Author: John
 *
 * Created on January 21, 2015, 10:48 AM
 */

#ifndef TWOD_MALLOC_H
#define	TWOD_MALLOC_H

bool alloc_uchar_array(unsigned char ***array, int rows, int columns, int channels);
void dealloc_uchar_array(unsigned char ***array);
bool alloc_float_array(float ***array, int rows, int columns, int channels);
void dealloc_float_array(float ***array);

#endif	/* TWOD_MALLOC_H */
