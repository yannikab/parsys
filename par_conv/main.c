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

/*
 * 
 */
int main(int argc, char** argv)
{
	int size, rank;
	int i, j, c;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int total_area = WIDTH * HEIGHT;
	int local_area = total_area / (size - 1);

	int s = 1;
	while (s * s * local_area != total_area)
		s++;

	int local_width = WIDTH / s;
	int local_height = HEIGHT / s;

	//	printf("%d %d\n", local_width, local_height);

	if (rank == 0)
	{
		unsigned char (**input_buffer)[CHANNELS];

		/* Read input image. */

		read_image((unsigned char ***) &input_buffer);

		/* Create "subarray" datatype. */

		MPI_Datatype local_image_t;
		MPI_Type_vector(local_height, local_width * CHANNELS, WIDTH * CHANNELS, MPI_UNSIGNED_CHAR, &local_image_t);
		MPI_Type_commit(&local_image_t);

		/* Send each subarray to corresponding process. */

		printf("Master: %d \n", rank);

		for (i = 1; i < size; i++)
		{
			int x = (i - 1) % s;
			int y = (i - 1) / s;
			MPI_Send(&(input_buffer[y * local_height][x * local_width][0]), 1, local_image_t, i, 0, MPI_COMM_WORLD);
			//			int d;
			//			MPI_Send(&d, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	} else
	{
		unsigned char (**local_buffer)[CHANNELS];
		alloc_uchar_array((unsigned char ***) &local_buffer, local_height, local_width, CHANNELS);

		printf("Slave: %d \n", rank);
		MPI_Status status;
		MPI_Recv(&(local_buffer[0][0][0]), local_height * local_width * CHANNELS, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		//		int d;
		//		MPI_Recv(&d, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		printf("Received at %d.\n", rank);
		int received;
		MPI_Get_count(&status, MPI_FLOAT, &received);
		printf("Count: %d.\n", received);

		unsigned char **output_buffer[CHANNELS];

		for (c = 0; c < CHANNELS; c++)
			alloc_uchar_array(&(output_buffer[c]), local_height, local_width, 1);

		for (c = 0; c < CHANNELS; c++)
			for (i = 0; i < local_height; i++)
				for (j = 0; j < local_width; j++)
					output_buffer[c][i][j] = local_buffer[i][j][c];

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
			sprintf(s, "%d", rank - 1);
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
