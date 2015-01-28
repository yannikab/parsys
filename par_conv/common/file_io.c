/* 
 * File:   file_io.c
 * Author: John
 *
 * Created on January 21, 2015, 10:33 AM
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>

#include <string.h>

#include "2d_malloc.h"

#include "file_io.h"

#define STRSIZE 512

char *create_file_name(file_type type, int channel)
{
    char *file_name = NULL;

    bool ok = true;

    char chan_str[10];

    chan_str[0] = '\0';
    if (channel >= 0 && channel < CHANNELS)
        sprintf(chan_str, "%d_", channel);

    unsigned int length = 1;
    length += strlen(type == INPUT ? INDIR : OUTDIR);
    length += strlen(chan_str);
    length += strlen(type == INPUT ? INFILENAME : OUTFILENAME);

    if (ok)
    {
        ok = (file_name = malloc(length)) != NULL;

        if (!ok)
            perror("malloc");
    }

    if (ok)
    {
        file_name[0] = '\0';

        strcpy(file_name, type == INPUT ? INDIR : OUTDIR);
        strcat(file_name, chan_str);
        strcat(file_name, type == INPUT ? INFILENAME : OUTFILENAME);
    }

    return file_name;
}

bool read_image(unsigned char ***file_buffer)
{
    bool ok = true;

    /* Open input file. */

    char *in_file = NULL;

    if (ok)
        ok = (in_file = create_file_name(INPUT, -1)) != NULL;

    FILE *in_fp = NULL;

    if (ok)
    {
        ok = (in_fp = fopen(in_file, "rb")) != NULL;

        if (!ok)
            perror(in_file);
        //		else
        //			printf("Input file opened.\n");
    }

    free(in_file);

    /* Allocate memory for file buffer. */

    if (ok)
        ok = alloc_uchar_array(file_buffer, HEIGHT, WIDTH, CHANNELS);

    /* Read image data one line at a time. */

    unsigned int i;

    for (i = 0; ok && i < HEIGHT; i++)
    {
        ok = fread((*file_buffer)[i], 1, WIDTH * CHANNELS, in_fp) == WIDTH * CHANNELS;

        if (!ok)
            perror("read_image");
    }

    //	if (ok)
    //		printf("Image read.\n");

    /* If an error occurs, free allocated memory. */

    if (!ok)
        dealloc_uchar_array(file_buffer);

    return ok;
}

bool write_channels(unsigned char (**file_buffer)[CHANNELS], int height, int width)
{
    bool ok = true;

    /* Create one output buffer per channel. */

    unsigned char **channel_buffer[CHANNELS];
    unsigned int c;

    for (c = 0; ok && c < CHANNELS; c++)
        ok = alloc_uchar_array(&(channel_buffer[c]), height, width, 1);

    /* Copy each channel from image, with conversion to byte. */

    unsigned int i, j;

    if (ok)
        for (i = 0; i < height; i++)
            for (j = 0; j < width; j++)
                for (c = 0; c < CHANNELS; c++)
                    channel_buffer[c][i][j] = file_buffer[i][j][c];

    /* Create filename for each output channel. */

    char *out_filename[CHANNELS];

    for (c = 0; c < CHANNELS; c++)
        out_filename[c] = NULL;

    for (c = 0; ok && c < CHANNELS; c++)
    {
        out_filename[c] = create_file_name(OUTPUT, c);
        ok = out_filename[c] != NULL;
    }

    /* Write out each channel to a separate raw file. */

    FILE * out_fp[CHANNELS];
    for (c = 0; c < CHANNELS; c++)
        out_fp[c] = NULL;

    for (c = 0; ok && c < CHANNELS; c++)
    {
        out_fp[c] = fopen(out_filename[c], "wb");
        ok = out_fp[c] != NULL;

        if (!ok)
            perror(out_filename[c]);
    }

    for (c = 0; ok && c < CHANNELS; c++)
        for (i = 0; ok && i < height; i++)
            ok = fwrite(channel_buffer[c][i], 1, width, out_fp[c]) == width;

    for (c = 0; c < CHANNELS; c++)
    {
        ok = fclose(out_fp[c]) == 0;

        if (!ok)
            fprintf(stderr, "Could not close file: %s\n", out_filename[c]);
    }

    /* Free memory allocated for channel buffers. */

    for (c = 0; c < CHANNELS; c++)
        dealloc_uchar_array(&(channel_buffer[c]));

    /* Calculate md5sums. */

    // printf("\n");
    char command[STRSIZE];
    for (c = 0; ok && c < CHANNELS; c++)
    {
        // sprintf(command, "md5sum %s %s.tiff", out_filename[c], out_filename[c]);
        sprintf(command, "md5sum %s", out_filename[c]);
        // printf("%s\n", command);
        ok = system(command) == 0;
    }

    //	/* Convert output files to tiff format (ImageMagick). */
    //
    //	// printf("\n");
    //	for (c = 0; ok && c < CHANNELS; c++)
    //	{
    //		// sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", height, width, out_filename[c], out_filename[c]);
    //		sprintf(command, "convert -depth 8 -size %dx%d gray:%s -compress lzw %s.tiff", width, height, out_filename[c], out_filename[c]);
    //		// printf("%s\n", command);
    //		ok = system(command) == 0;
    //	}
    //
    //	/* Merge individual channel tiffs to a single tiff (ImageMagick). */
    //
    //	sprintf(command, "convert");
    //	if (ok)
    //	{
    //		for (c = 0; c < CHANNELS; c++)
    //		{
    //			strcat(command, " ");
    //			strcat(command, out_filename[c]);
    //			strcat(command, ".tiff");
    //		}
    //
    //		strcat(command, " -combine ");
    //		strcat(command, OUTDIR);
    //		strcat(command, OUTFILENAME);
    //		strcat(command, ".tiff");
    //
    //		// printf("%s\n", command);
    //
    //		ok = system(command) == 0;
    //	}
    //
    //	// printf("\n");

    /* Free memory allocated for filenames. */

    for (c = 0; c < CHANNELS; c++)
    {
        free(out_filename[c]);
        out_filename[c] = NULL;
    }

    return ok;
}
