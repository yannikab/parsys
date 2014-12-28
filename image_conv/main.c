/* 
 * File:   main.c
 * Author: John
 *
 * Created on 16 Δεκέμβριος 2014, 6:35 μμ
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <unistd.h>
#include <sys/times.h>

#include <string.h>
#include <assert.h>

#define WIDTH 1920
#define HEIGHT 2520
#define CHANNELS 3

#define B 1
static const float filter[3][3] = {0.0625, 0.125, 0.0625, 0.125, 0.25, 0.125, 0.0625, 0.125, 0.0625};

#define ITERATIONS 5

#define INFILENAME "waterfall_1920_2520.raw"
#define OUTFILENAME "waterfall_1920_2520_out.raw"
#define STRSIZE 255

//void apply_filter(char *image_data_out, const char *image_data_in, unsigned int height, unsigned int width);
//void split_channels();
//float h(int p, int q);
//float in(unsigned int i, unsigned int j, unsigned int b);
//void out(unsigned int i, unsigned int j, unsigned int b, float value);

bool alloc_uchar_array(unsigned char ***array, int rows, int columns, int channels)
{
	char *p;
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

void dealloc_uchar_array(unsigned char ***array)
{
	assert(array != NULL);

	free(&((*array)[0][0]));
	free(*array);
	*array = NULL;
}

void dealloc_float_array(unsigned char ***array)
{
	assert(array != NULL);

	free(&((*array)[0][0]));
	free(*array);
	*array = NULL;
}

bool read_image(unsigned char ***input_buffer)
{
	FILE *in_fp;
	int i;

	bool ok = true;

	/* Open input file. */

	in_fp = NULL;
	if (ok && (in_fp = fopen(INFILENAME, "rb")) == NULL)
	{
		perror(INFILENAME);
		ok = false;
	}

	//	if (ok)
	//		printf("Input file opened.\n");

	/* Allocate memory for file buffer. */

	if (ok && !alloc_uchar_array(input_buffer, HEIGHT, WIDTH, CHANNELS))
		ok = false;

	//	if (ok)
	//		printf("Buffer allocated.\n");

	/* Read image data one row at a time. */

	for (i = 0; ok && i < HEIGHT; i++)
		if (fread((*input_buffer)[i], 1, WIDTH * CHANNELS, in_fp) != WIDTH * CHANNELS)
			ok = false;

	//	if (ok)
	//		printf("Image read.\n");

	/* If an error occurs, free allocated memory. */

	if (!ok)
		dealloc_uchar_array(input_buffer);

	return ok;
}

bool write_channels(float (**input_image_data)[CHANNELS], int height, int width)
{
	int i, j, c;

	char *out_filename[CHANNELS];
	FILE * out_fp[CHANNELS];

	for (c = 0; c < CHANNELS; c++)
	{
		out_filename[c] = NULL;
		out_fp[c] = NULL;
	}

	bool ok = true;

	/* Create one output buffer per channel. */

	unsigned char **output_buffer[CHANNELS];

	for (c = 0; ok && c < CHANNELS; c++)
		ok = alloc_uchar_array(&(output_buffer[c]), height, width, 1);

	/* Copy each channel from image, with conversion to byte. */

	for (c = 0; ok && c < CHANNELS; c++)
		for (i = 0; i < height; i++)
			for (j = 0; j < width; j++)
				output_buffer[c][i][j] = (char) input_image_data[i][j][c];

	/* Create filename for each output channel. */

	for (c = 0; ok && c < CHANNELS; c++)
	{
		out_filename[c] = malloc(STRSIZE * sizeof (char));

		if (out_filename[c] == NULL)
		{
			perror("malloc");
			ok = false;
		}

		if (ok)
		{
			out_filename[c][0] = '\0';

			char s[STRSIZE];

			s[0] = '\0';
			sprintf(s, "../serial/%d", c);
			strcat(out_filename[c], s);
			strcat(out_filename[c], "_");

			strcat(out_filename[c], OUTFILENAME);
		}
	}

	for (c = 0; ok && c < CHANNELS; c++)
	{
		out_fp[c] = fopen(out_filename[c], "wb");
		if (out_fp[c] == NULL)
		{
			perror(out_filename[c]);
			ok = false;
		}
	}

	for (c = 0; ok && c < CHANNELS; c++)
		for (i = 0; ok && i < height; i++)
			ok = fwrite(output_buffer[c][i], 1, width, out_fp[c]) == width;

	for (c = 0; c < CHANNELS; c++)
	{
		if (fclose(out_fp[c]))
		{
			fprintf(stderr, "Could not close file: %s", out_filename[c]);
			ok = false;
		}
	}

	/* Convert output files to tiff format. */

	char command[STRSIZE];
	for (c = 0; ok && c < CHANNELS; c++)
	{
		command[0] = '\0';
		sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", height, width, out_filename[c], out_filename[c]);
		printf("%s\n", command);
		system(command);
	}

	/* Free memory for filenames. */
	for (c = 0; c < CHANNELS; c++)
	{
		free(out_filename[c]);
		out_filename[c] = NULL;
	}

	return ok;
}

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

	for (c = 0; c < CHANNELS; c++)
	{
		/* Left border column. */

		for (i = B; i < height - B; i++)
		{
			for (j = B; j < 2 * B; j++)
			{
				float value = 0.0f;

				for (p = -B; p <= B; p++)
					for (q = -B; q <= B; q++)
						value += input_image[i - p][j - q][c] * filter[p + B][q + B];

				output_image[i][j][c] = value;
			}
		}

		/* Right border column. */

		for (i = B; i < height - B; i++)
		{
			for (j = width - 1 - B; j > width - 1 - 2 * B; j--)
			{
				float value = 0.0f;

				for (p = -B; p <= B; p++)
					for (q = -B; q <= B; q++)
						value += input_image[i - p][j - q][c] * filter[p + B][q + B];

				output_image[i][j][c] = value;
			}
		}

		/* Top border row, avoid recalculating overlap with columns. */

		for (j = 2 * B; j < width - 2 * B; j++)
		{
			for (i = B; i < 2 * B; i++)
			{
				float value = 0.0f;

				for (p = -B; p <= B; p++)
					for (q = -B; q <= B; q++)
						value += input_image[i - p][j - q][c] * filter[p + B][q + B];

				output_image[i][j][c] = value;
			}
		}

		/* Bottom border row, avoid recalculating overlap with columns. */

		for (j = 2 * B; j < width - 2 * B; j++)
		{
			for (i = height - 1 - B; i > height - 1 - 2 * B; i--)
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

/*
 * 
 */
int main(int argc, char** argv)
{
	unsigned int i, j, c, n;

	bool ok = true;

	/* Read input file into buffer. */

	char (**input_buffer)[CHANNELS];

	if (ok)
		ok = read_image((unsigned char ***) &input_buffer);

	/* Allocate memory for image data. */

	float (**input_image_data)[CHANNELS];

	if (ok)
		ok = alloc_float_array((float ***) &input_image_data, B + HEIGHT + B, B + WIDTH + B, CHANNELS);

	float (**output_image_data)[CHANNELS];

	if (ok)
		ok = alloc_float_array((float ***) &output_image_data, B + HEIGHT + B, B + WIDTH + B, CHANNELS);

	/* Convert input. */

	for (c = 0; ok && c < CHANNELS; c++)
		for (i = 0; i < HEIGHT; i++)
			for (j = 0; j < WIDTH; j++)
				output_image_data[i + B][j + B][c] = (float) input_buffer[i][j][c];

	/* Start timing. */

	double t1, t2, real_time;
	struct tms tb1, tb2;
	double tickspersec = (double) sysconf(_SC_CLK_TCK);

	t1 = (double) times(&tb1);

	/* Apply filter. */

	for (n = 0; n < ITERATIONS; n++)
	{
		/* Fill borders with outer image data. */

		// south
		for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
			for (j = B; j < B + WIDTH; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B + HEIGHT - 1][j][c];

		// north
		for (i = 0; i < B; i++) // north
			for (j = B; j < B + WIDTH; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B][j][c];

		// east
		for (i = B; i < B + HEIGHT; i++)
			for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[i][B + WIDTH - 1][c];

		// west
		for (i = B; i < B + HEIGHT; i++)
			for (j = 0; j < B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[i][B][c];

		// se
		for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
			for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B + HEIGHT - 1][B + WIDTH - 1][c];

		// nw
		for (i = 0; i < B; i++)
			for (j = 0; j < B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B][B][c];

		// sw
		for (i = HEIGHT + B; i < HEIGHT + 2 * B; i++)
			for (j = 0; j < B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B + HEIGHT - 1][B][c];

		// ne
		for (i = 0; i < B; i++)
			for (j = WIDTH + B; j < WIDTH + 2 * B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = output_image_data[B][B + WIDTH - 1][c];

		/* Use previous output as new input. */

		float (**tmp)[CHANNELS];
		tmp = input_image_data;
		input_image_data = output_image_data;
		output_image_data = tmp;

		apply_inner_filter(output_image_data, input_image_data, B + HEIGHT + B, B + WIDTH + B);

		apply_outer_filter(output_image_data, input_image_data, B + HEIGHT + B, B + WIDTH + B);
	}

	/* Stop time measurement, print time. */

	t2 = (double) times(&tb2);

	real_time = (double) (t2 - t1) / tickspersec;
	printf("Completed in %.3f sec\n", real_time);

	write_channels(output_image_data, B + HEIGHT + B, B + WIDTH + B);

	return (EXIT_SUCCESS);
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
