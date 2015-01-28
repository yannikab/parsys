/* 
 * File:   file_io.h
 * Author: John
 *
 * Created on January 21, 2015, 10:34 AM
 */

#ifndef FILE_IO_H
#define	FILE_IO_H

#include "../settings.h"

typedef enum {
    INPUT,
    OUTPUT,
} file_type;

char *create_file_name(file_type type, int channel);
bool read_image(unsigned char ***file_buffer);
bool write_channels(unsigned char (**file_buffer)[CHANNELS], int height, int width);

#endif	/* FILE_IO_H */
