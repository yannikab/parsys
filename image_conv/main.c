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
#define BYTES 3
#define STRSIZE 255

static FILE *in_fp = NULL;
static FILE *out_fp[BYTES];

static char *out_filename[BYTES];

static unsigned char *input_image_data = NULL;
static unsigned char *output_image_data = NULL;

static float *input_image = NULL;
static float *output_image = NULL;

static unsigned char *output_channel[BYTES];

//void apply_filter(char *image_data_out, const char *image_data_in, unsigned int height, unsigned int width);
void convert_input();
void convert_output();
void apply_filter();
void split_channels();
double h(int p, int q);
void out(unsigned int i, unsigned int j, double value);
double in(unsigned int i, unsigned int j);

/*
 * 
 */
int main(int argc, char** argv)
{
	unsigned int b, c;

	/* Open files. */
	in_fp = fopen(INFILENAME, "rb");
	if (in_fp == NULL)
	{
		perror(INFILENAME);
		return (EXIT_FAILURE);
	}

	for (b = 0; b < BYTES; b++)
	{
		out_filename[b] = NULL;

		out_filename[b] = malloc(STRSIZE * sizeof (char));

		if (out_filename[b] == NULL)
		{
			perror("malloc");
			for (c = 0; c < b; c++)
				free(out_filename[b]);
			return (EXIT_FAILURE);
		}

		out_filename[b][0] = '\0';

		char s[STRSIZE];
		s[0] = '\0';
		sprintf(s, "%d", b);
		strcat(out_filename[b], s);

		strcat(out_filename[b], "_");

		strcat(out_filename[b], OUTFILENAME);
	}

	for (b = 0; b < BYTES; b++)
	{
		out_fp[b] = NULL;
		out_fp[b] = fopen(out_filename[b], "wb");
		if (out_fp[b] == NULL)
		{
			perror(out_filename[b]);
			return (EXIT_FAILURE);
		}
	}

	/* Start timing. */
	double t1, t2, real_time;
	struct tms tb1, tb2;
	double tickspersec = (double) sysconf(_SC_CLK_TCK);

	t1 = (double) times(&tb1);

	/* Allocate memory for image data. */
	input_image_data = malloc(BYTES * HEIGHT * WIDTH * sizeof (unsigned char));
	if (input_image_data == NULL)
	{
		perror("malloc");
		return (EXIT_FAILURE);
	}

	output_image_data = malloc(BYTES * HEIGHT * WIDTH * sizeof (unsigned char));
	if (output_image_data == NULL)
	{
		perror("malloc");
		free(input_image_data);
		return (EXIT_FAILURE);
	}

	/* Allocate memory for image data. */
	input_image = malloc(BYTES * HEIGHT * WIDTH * sizeof (float));
	if (input_image == NULL)
	{
		perror("malloc");
		free(input_image_data);
		free(output_image_data);
		return (EXIT_FAILURE);
	}

	output_image = malloc(BYTES * HEIGHT * WIDTH * sizeof (float));
	if (output_image == NULL)
	{
		perror("malloc");
		free(input_image_data);
		free(output_image_data);
		free(input_image);
		return (EXIT_FAILURE);
	}

	for (b = 0; b < BYTES; b++)
	{
		output_channel[b] = NULL;

		output_channel[b] = malloc(HEIGHT * WIDTH * sizeof (unsigned char));

		if (output_channel[b] == NULL)
		{
			perror("malloc");

			for (c = 0; c < b; c++)
				free(output_channel[b]);

			free(output_image_data);

			free(input_image_data);

			return (EXIT_FAILURE);
		}
	}

	//	unsigned char c;
	//	double d;
	//	for (c = 0; c < 255; c++)
	//	{
	//		d = (double) c;
	//		printf("%c %lf ", c, d);
	//	}

	unsigned int i, j;

	/* Read image data one row at a time. */
	for (i = 0; i < HEIGHT; i++)
	{
		//		for (j = 0; j < WIDTH; j++)
		//		{
		//			fread(image_data[i], 1, 1, in_fp);
		//			fwrite(image_data[i], 1, 1, out_fp);
		//		}

		//		fread(input_image_data[i], 1, WIDTH, in_fp);

		// read one row at a time
		fread(input_image_data + i * BYTES * WIDTH, 1, BYTES * WIDTH, in_fp);
	}

	convert_input();

	/* Apply filter. */
	//		apply_filter(output_image_data, input_image_data, HEIGHT, WIDTH);
	apply_filter();

	//		fwrite(input_image_data[i], 1, WIDTH, out_fp);
	//		fwrite(output_image_data + i * WIDTH, 1, WIDTH, out_fp);

	convert_output();

	// split output image data to separate channels
	split_channels();

	// write out channels to separate files
	for (b = 0; b < BYTES; b++)
		for (i = 0; i < HEIGHT; i++)
			fwrite(output_channel[b] + i * WIDTH, 1, WIDTH, out_fp[b]);

	/* Close files. */
	if (fclose(in_fp))
	{
		fprintf(stderr, "Could not close file: %s", INFILENAME);
		return (EXIT_FAILURE);
	}

	for (b = 0; b < BYTES; b++)
	{
		if (fclose(out_fp[b]))
		{
			fprintf(stderr, "Could not close file: %s", out_filename[b]);
			return (EXIT_FAILURE);
		}
	}

	/* Stop time measurement, print time. */

	t2 = (double) times(&tb2);

	real_time = (double) (t2 - t1) / tickspersec;
	printf("Completed in %.3f sec\n", real_time);

	/* Convert output to tiff format. */

	char command[STRSIZE];
	for (b = 0; b < BYTES; b++)
	{
		command[0] = '\0';
		sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", HEIGHT, WIDTH, out_filename[b], out_filename[b]);
		system(command);
	}

	//	command[0] = '\0';
	//	sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", HEIGHT, WIDTH, INFILENAME, INFILENAME);
	//	system(command);

	return (EXIT_SUCCESS);
}

void apply_filter()
{
	unsigned int i, j;
	int p, q;

	//	for (i = 0; i < height * width; i++)
	//	{
	//		*(image_data_out + i) = *(image_data_in + i);
	//	}

	/* Null filter, memory copy. */
	//	memcpy(output_image_data, input_image_data, BYTES * HEIGHT * WIDTH * sizeof (char));
	//	return;

	memcpy(output_image, input_image, BYTES * HEIGHT * WIDTH * sizeof (float));
	return;

	//	for (i = 0; i < HEIGHT; i++)
	//	{
	//		if (i == 0 || i == HEIGHT - 1)
	//			continue;
	//
	//		for (j = 0; j < WIDTH; j++)
	//		{
	//			if (j == 0 || i == HEIGHT - 1)
	//				continue;
	//
	//			double value = 0.0;
	//
	//			for (p = -1; p <= 1; p++)
	//			{
	//				for (q = -1; q <= 1; q++)
	//				{
	//					//					*(output_image_data + WIDTH * i + j) = *(input_image_data + WIDTH * i + j);
	//					double inval = in(i - p, j - q);
	//					inval *= h(p, q);
	//					value += inval;
	//				}
	//			}
	//
	//			out(i, j, value);
	//		}
	//	}
}

void split_channels()
{
	unsigned int i, j, b;

	for (i = 0; i < HEIGHT; i++)
	{
		for (j = 0; j < WIDTH; j++)
		{
			for (b = 0; b < BYTES; b++)
				*(output_channel[b] + i * WIDTH + j) = *(output_image_data + BYTES * (i * WIDTH + j) + b);
		}
	}
}

static const double filterArray[3][3] = {0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625};

double h(int p, int q)
{
	assert(p >= -1 && p <= 1);
	assert(q >= -1 && q <= 1);

	return filterArray[p + 1][q + 1];
}

double in(unsigned int i, unsigned int j)
{
	assert(i < HEIGHT);
	assert(j < WIDTH);

	return (double) *(input_image_data + WIDTH * i + j);
}

void out(unsigned int i, unsigned int j, double value)
{
	assert(i < HEIGHT);
	assert(j < WIDTH);
	assert(value >= 0.0 && value <= 255.0);

	//	*(output_image_data + WIDTH * i + j) = (unsigned char) value;
}

void convert_input()
{
	unsigned int i;
	for (i = 0; i < BYTES * HEIGHT * WIDTH; i++)
		*(input_image + i) = (float) *(input_image_data + i);
}

void convert_output()
{
	unsigned int i;
	for (i = 0; i < BYTES * HEIGHT * WIDTH; i++)
		*(output_image_data + i) = (unsigned char) *(output_image + i);
}
