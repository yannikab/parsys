/* 
 * File:   main.c
 * Author: John
 *
 * Created on 24 Δεκέμβριος 2014, 12:02 πμ
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#define WIDTH 1920
#define HEIGHT 2520
#define CHANNELS 3

#define B 1
static const float filter[3][3] = {0.0625f, 0.125f, 0.0625f, 0.125f, 0.25f, 0.125f, 0.0625f, 0.125f, 0.0625f};

#define ITERATIONS 20

static FILE *out_fp[CHANNELS];

static char *out_filename[CHANNELS];

#define INFILENAME "waterfall_1920_2520.raw"
#define OUTFILENAME "waterfall_1920_2520_out.raw"
#define STRSIZE 255

static unsigned char **output_buffer[CHANNELS]; // three output buffers of single bytes

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

	if (ok)
		printf("Input file opened.\n");

	/* Allocate memory for file buffer. */

	if (ok && !alloc_uchar_array(input_buffer, HEIGHT, WIDTH, CHANNELS))
		ok = false;

	if (ok)
		printf("Buffer allocated.\n");

	/* Read image data one row at a time. */

	for (i = 0; ok && i < HEIGHT; i++)
		if (fread((*input_buffer)[i], 1, WIDTH * CHANNELS, in_fp) != WIDTH * CHANNELS)
			ok = false;

	if (ok)
		printf("Image read.\n");

	/* If an error occurs, free allocated memory. */

	if (!ok)
		dealloc_uchar_array(input_buffer);

	return ok;
}

bool open_output_files(int rank)
{
	int c, i;

	for (c = 0; c < CHANNELS; c++)
	{
		out_filename[c] = NULL;

		out_filename[c] = malloc(STRSIZE * sizeof (char));

		if (out_filename[c] == NULL)
		{
			perror("malloc");
			for (i = 0; i < c; i++)
				free(out_filename[c]);
			return false;
		}

		out_filename[c][0] = '\0';

		char s[STRSIZE];

		s[0] = '\0';
		sprintf(s, "../%d", c);
		strcat(out_filename[c], s);
		strcat(out_filename[c], "_");

		s[0] = '\0';
		sprintf(s, "%d", rank - 1);
		strcat(out_filename[c], s);
		strcat(out_filename[c], "_");

		strcat(out_filename[c], OUTFILENAME);

		printf("%s\n", out_filename[c]);
	}

	for (c = 0; c < CHANNELS; c++)
	{
		out_fp[c] = NULL;
		out_fp[c] = fopen(out_filename[c], "wb");
		if (out_fp[c] == NULL)
		{
			perror(out_filename[c]);

			return false;
		}
	}
}

bool close_output_files()
{
	int c;

	/* Close files. */

	for (c = 0; c < CHANNELS; c++)
	{
		if (fclose(out_fp[c]))
		{
			fprintf(stderr, "Could not close file: %s", out_filename[c]);

			return false;
		}
	}

	return true;
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
	int size, rank, master;
	int i, j, c, r;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	master = size - 1;

	int total_area = WIDTH * HEIGHT;
	int local_area = total_area / (size - 1);

	int s = 1;
	while (s * s * local_area != total_area)
		s++;

	int local_width = WIDTH / s;
	int local_height = HEIGHT / s;

	//	printf("%d %d\n", local_width, local_height);

	/* Define s x s cartesian topology of slaves, excluding master node which is last. */

	//	MPI_Comm comm_slaves;
	//	int ndims = 2;
	//	int dims[] = {s, s};
	//	int periods[] = {0, 0};
	//	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, true, &comm_slaves);

	if (rank == master)
	{
		unsigned char (**input_buffer)[CHANNELS];

		/* Read input image. */

		read_image((unsigned char ***) &input_buffer);

		/* Convert byte data to float for applying arithmetic operations. */

		float (**input_image_data)[CHANNELS];

		alloc_float_array((float ***) &input_image_data, HEIGHT, WIDTH, CHANNELS);

		for (c = 0; c < CHANNELS; c++)
			for (i = 0; i < HEIGHT; i++)
				for (j = 0; j < WIDTH; j++)
					input_image_data[i][j][c] = (float) input_buffer[i][j][c];

		/* Create "subarray" datatype. */

		MPI_Datatype local_image_t;
		MPI_Type_vector(local_height, local_width * CHANNELS, WIDTH * CHANNELS, MPI_FLOAT, &local_image_t);
		MPI_Type_commit(&local_image_t);

		/* Send each subarray to corresponding process. */

		printf("Master: %d \n", rank);

		for (r = 0; r < size; r++)
		{
			if (r == master)
				continue;

			//			int d;
			//			MPI_Send(&d, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
			int x = r % s;
			int y = r / s;
			MPI_Send(&(input_image_data[y * local_height][x * local_width][0]), 1, local_image_t, r, 0, MPI_COMM_WORLD);

			printf("Send returned: %d \n", rank);
		}
	} else
	{
		float (**input_image_data)[CHANNELS];
		alloc_float_array((float ***) &input_image_data, B + local_height + B, B + local_width + B, CHANNELS);

		float (**output_image_data)[CHANNELS];
		alloc_float_array((float ***) &output_image_data, B + local_height + B, B + local_width + B, CHANNELS);

		/* Create inner image datatype. */

		MPI_Datatype inner_image_t;
		MPI_Type_vector(local_height, local_width * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &inner_image_t);
		MPI_Type_commit(&inner_image_t);

		/* Receive inner image. */

		printf("Slave: %d \n", rank);
		MPI_Status status;
		//		int d;
		//		MPI_Recv(&d, 1, MPI_INT, master, 0, MPI_COMM_WORLD, &status);
		//		MPI_Recv(&(output_image_data[1][1][0]), local_height * local_width * CHANNELS, MPI_FLOAT, master, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&(output_image_data[B][B][0]), 1, inner_image_t, master, 0, MPI_COMM_WORLD, &status);
		printf("Received at %d.\n", rank);

		int received;
		MPI_Get_count(&status, inner_image_t, &received);
		//		printf("Count: %d.\n", received);

		/* Apply filter. */

		for (i = 0; i < ITERATIONS; i++)
		{
			/* Use previous output as new input. */

			float (**tmp)[CHANNELS];
			tmp = input_image_data;
			input_image_data = output_image_data;
			output_image_data = tmp;

			// TODO: send and receive neighbouring node data

			apply_inner_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);

			apply_outer_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);
		}

		unsigned char **output_buffer[CHANNELS];

		for (c = 0; c < CHANNELS; c++)
			alloc_uchar_array(&(output_buffer[c]), local_height, local_width, 1);

		for (c = 0; c < CHANNELS; c++)
			for (i = 0; i < local_height; i++)
				for (j = 0; j < local_width; j++)
					output_buffer[c][i][j] = (char) output_image_data[i + B][j + B][c];

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
			sprintf(s, "../%d", c);
			strcat(out_filename[c], s);
			strcat(out_filename[c], "_");

			s[0] = '\0';
			sprintf(s, "%d", rank);
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

		for (c = 0; c < CHANNELS; c++)
			for (i = 0; i < local_height; i++)
				fwrite(output_buffer[c][i], 1, local_width, out_fp[c]);

		for (c = 0; c < CHANNELS; c++)
		{
			if (fclose(out_fp[c]))
			{
				fprintf(stderr, "Could not close file: %s", out_filename[c]);
				return (EXIT_FAILURE);
			}
		}

		/* Convert output files to tiff format. */

		char command[STRSIZE];
		for (c = 0; c < CHANNELS; c++)
		{
			command[0] = '\0';
			sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", local_height, local_width, out_filename[c], out_filename[c]);
			printf("%s\n", command);
			system(command);
		}

		for (c = 0; c < CHANNELS; c++)
		{
			free(out_filename[c]);
			out_filename[c] = NULL;
		}
	}

	MPI_Finalize();

	return (EXIT_SUCCESS);
}
