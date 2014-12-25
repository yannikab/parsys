/* 
 * File:   main.c
 * Author: John
 *
 * Created on 16 Δεκέμβριος 2014, 6:35 μμ
 */

#include <stdio.h>
#include <stdlib.h>

#include <unistd.h>
#include <sys/times.h>

#include <string.h>
#include <assert.h>

#define WIDTH 1920
#define HEIGHT 2520
#define INFILENAME "waterfall_1920_2520.raw"
#define OUTFILENAME "waterfall_1920_2520_out.raw"
#define CHANNELS 3
#define ITERATIONS 5
#define STRSIZE 255

static FILE *in_fp = NULL;
static FILE *out_fp[CHANNELS];

static char *out_filename[CHANNELS];

static unsigned char (*input_buffer)[WIDTH][CHANNELS];
static unsigned char (*output_buffer)[WIDTH][CHANNELS];

static float (*input_image)[WIDTH][CHANNELS];
static float (*output_image)[WIDTH][CHANNELS];

static unsigned char (*output_channel[CHANNELS])[WIDTH];

//void apply_filter(char *image_data_out, const char *image_data_in, unsigned int height, unsigned int width);
void convert_input();
void convert_output();
void apply_filter();
//void split_channels();
//float h(int p, int q);
//float in(unsigned int i, unsigned int j, unsigned int b);
//void out(unsigned int i, unsigned int j, unsigned int b, float value);

/*
 * 
 */
int main(int argc, char** argv)
{
	unsigned int i, j, c;

	/* Open files. */

	in_fp = fopen(INFILENAME, "rb");
	if (in_fp == NULL)
	{
		perror(INFILENAME);
		return (EXIT_FAILURE);
	}

	for (c = 0; c < CHANNELS; c++)
	{
		out_filename[c] = NULL;

		out_filename[c] = malloc(STRSIZE * sizeof (char));

		if (out_filename[c] == NULL)
		{
			perror("malloc");
			for (i = 0; i < c; i++)
				free(out_filename[c]);
			return (EXIT_FAILURE);
		}

		out_filename[c][0] = '\0';

		char s[STRSIZE];
		s[0] = '\0';
		sprintf(s, "%d", c);
		strcat(out_filename[c], s);

		strcat(out_filename[c], "_");

		strcat(out_filename[c], OUTFILENAME);
	}

	for (c = 0; c < CHANNELS; c++)
	{
		out_fp[c] = NULL;
		out_fp[c] = fopen(out_filename[c], "wb");
		if (out_fp[c] == NULL)
		{
			perror(out_filename[c]);
			return (EXIT_FAILURE);
		}
	}

	/* Start timing. */

	double t1, t2, real_time;
	struct tms tb1, tb2;
	double tickspersec = (double) sysconf(_SC_CLK_TCK);

	t1 = (double) times(&tb1);

	/* Allocate memory for file buffers. */

	input_buffer = malloc(HEIGHT * sizeof (*input_buffer));
	if (input_buffer == NULL)
	{
		perror("malloc");
		return (EXIT_FAILURE);
	}

	output_buffer = malloc(HEIGHT * sizeof (*output_buffer));
	if (output_buffer == NULL)
	{
		perror("malloc");
		free(input_buffer);
		return (EXIT_FAILURE);
	}

	/* Allocate memory for image data. */

	input_image = malloc(HEIGHT * sizeof (*input_image));
	if (input_image == NULL)
	{
		perror("malloc");
		free(input_buffer);
		free(output_buffer);
		return (EXIT_FAILURE);
	}

	output_image = malloc(HEIGHT * sizeof (*output_image));
	if (output_image == NULL)
	{
		perror("malloc");
		free(input_buffer);
		free(output_buffer);
		free(input_image);
		return (EXIT_FAILURE);
	}

	for (c = 0; c < CHANNELS; c++)
	{
		output_channel[c] = NULL;

		output_channel[c] = malloc(HEIGHT * sizeof (*output_channel[c]));

		if (output_channel[c] == NULL)
		{
			perror("malloc");

			for (i = 0; i < c; i++)
				free(output_channel[i]);

			free(output_buffer);
			free(input_buffer);

			free(output_image);
			free(input_image);

			return (EXIT_FAILURE);
		}
	}

	/* Read image data one row at a time. */

	for (i = 0; i < HEIGHT; i++)
		fread(input_buffer[i], 1, WIDTH * CHANNELS, in_fp);

	/* Convert input. */

	for (i = 0; i < HEIGHT; i++)
		for (j = 0; j < WIDTH; j++)
			for (c = 0; c < CHANNELS; c++)
				output_image[i][j][c] = (float) input_buffer[i][j][c];
	// *(output_image + i) = (float) *(input_image_data + i);

	/* Apply filter. */

	for (i = 0; i < ITERATIONS; i++)
	{
		/* Use previous output as new input. */

		float (*tmp)[WIDTH][CHANNELS];
		tmp = input_image;
		input_image = output_image;
		output_image = tmp;

		apply_filter();
	}

	/* Convert output. */

	for (i = 0; i < HEIGHT; i++)
		for (j = 0; j < WIDTH; j++)
			for (c = 0; c < CHANNELS; c++)
				output_channel[c][i][j] = (char) output_image[i][j][c];
	// output_buffer[i][j][c] = (float) output_image[i][j][c];
	// *(output_image_data + i) = (unsigned char) *(output_image + i);

	/* Split output image data to separate channels. */

	// split_channels();

	/* Write out channels to separate files. */

	for (c = 0; c < CHANNELS; c++)
		for (i = 0; i < HEIGHT; i++)
			fwrite(output_channel[c][i], 1, WIDTH, out_fp[c]);

	/* Close files. */
	
	if (fclose(in_fp))
	{
		fprintf(stderr, "Could not close file: %s", INFILENAME);
		return (EXIT_FAILURE);
	}

	for (c = 0; c < CHANNELS; c++)
	{
		if (fclose(out_fp[c]))
		{
			fprintf(stderr, "Could not close file: %s", out_filename[c]);
			return (EXIT_FAILURE);
		}
	}

	/* Stop time measurement, print time. */

	t2 = (double) times(&tb2);

	real_time = (double) (t2 - t1) / tickspersec;
	printf("Completed in %.3f sec\n", real_time);

	/* Convert output to tiff format. */

	char command[STRSIZE];
	for (c = 0; c < CHANNELS; c++)
	{
		command[0] = '\0';
		sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", HEIGHT, WIDTH, out_filename[c], out_filename[c]);
		system(command);
	}

	return (EXIT_SUCCESS);
}


static const float filter[3][3] = {0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625};

void apply_filter()
{
	unsigned int i, j, c;
	int p, q;

	/* Null filter, memory copy. */
	//	memcpy(output_image_data, input_image_data, BYTES * HEIGHT * WIDTH * sizeof (char));
	//	return;

	//	memcpy(output_image, input_image, BYTES * HEIGHT * WIDTH * sizeof (float));
	//	return;

	for (i = 0; i < HEIGHT; i++)
	{
		for (j = 0; j < WIDTH; j++)
		{
			for (c = 0; c < CHANNELS; c++)
			{
				float value = 0.0;

				for (p = -1; p <= 1; p++)
				{
					if (i - p < 0 || i - p > HEIGHT - 1)
						continue;

					for (q = -1; q <= 1; q++)
					{
						if (j - q < 0 || j - q > WIDTH - 1)
							continue;

						value += input_image[i - p][j - q][c] * filter[p + 1][q + 1];
					}
				}

				output_image[i][j][c] = value;
			}
		}
	}
}

//void split_channels()
//{
//	unsigned int i, j, c;
//
//	for (i = 0; i < HEIGHT; i++)
//		for (j = 0; j < WIDTH; j++)
//			for (c = 0; c < CHANNELS; c++)
//				output_channel[c][i][j] = output_buffer[i][j][c];
//	//				*(output_channel[c] + i * WIDTH + j) = *(output_image_data + CHANNELS * (i * WIDTH + j) + c);
//}
