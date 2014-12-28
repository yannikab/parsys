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

#define ITERATIONS 5

#define INFILENAME "waterfall_1920_2520.raw"
#define OUTFILENAME "waterfall_1920_2520_out.raw"
#define STRSIZE 255

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
			sprintf(s, "../parallel/%d", c);
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

void get_neighbors(MPI_Comm comm, int *n_p, int *s_p, int *e_p, int *w_p, int *nw_p, int *se_p, int *ne_p, int *sw_p)
{
	int north[2], south[2], east[2], west[2];
	int nw[2], se[2], ne[2], sw[2];

	/* Get north/south/east/west ranks. */

	MPI_Cart_shift(comm, 0, 1, n_p, s_p);
	MPI_Cart_shift(comm, 1, -1, e_p, w_p);

	/* Get cartesian coordinates for n/s/e/w. */

	if (*n_p != MPI_PROC_NULL)
		MPI_Cart_coords(comm, *n_p, 2, north);
	if (*s_p != MPI_PROC_NULL)
		MPI_Cart_coords(comm, *s_p, 2, south);
	if (*e_p != MPI_PROC_NULL)
		MPI_Cart_coords(comm, *e_p, 2, east);
	if (*w_p != MPI_PROC_NULL)
		MPI_Cart_coords(comm, *w_p, 2, west);

	/* Get diagonal ranks. */

	if (*n_p != MPI_PROC_NULL && *w_p != MPI_PROC_NULL) // north west
	{
		nw[0] = north[0];
		nw[1] = west[1];
		MPI_Cart_rank(comm, nw, nw_p);
	} else
	{
		*nw_p = MPI_PROC_NULL;
	}

	if (*s_p != MPI_PROC_NULL && *e_p != MPI_PROC_NULL) // south east
	{
		se[0] = south[0];
		se[1] = east[1];
		MPI_Cart_rank(comm, se, se_p);
	} else
	{
		*se_p = MPI_PROC_NULL;
	}

	if (*n_p != MPI_PROC_NULL && *e_p != MPI_PROC_NULL) // north east
	{
		ne[0] = north[0];
		ne[1] = east[1];
		MPI_Cart_rank(comm, ne, ne_p);
	} else
	{
		*ne_p = MPI_PROC_NULL;
	}

	if (*s_p != MPI_PROC_NULL && *w_p != MPI_PROC_NULL) // south west
	{
		sw[0] = south[0];
		sw[1] = west[1];
		MPI_Cart_rank(comm, sw, sw_p);
	} else
	{
		*sw_p = MPI_PROC_NULL;
	}
}

bool in_even_row(int rank, MPI_Comm comm)
{
	int coords[2];

	MPI_Cart_coords(comm, rank, 2, coords);

	return coords[0] % 2 == 0;
}

bool in_even_column(int rank, MPI_Comm comm)
{
	int coords[2];

	MPI_Cart_coords(comm, rank, 2, coords);

	return coords[1] % 2 == 0;
}

/*
 * 
 */
int main(int argc, char** argv)
{
	int size, rank, master;
	int i, j, c, r, n;

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

	/* Define s x s cartesian topology of slaves, excluding master node which is last in rank. */

	MPI_Comm comm_slaves;
	int ndims = 2;
	int dims[] = {s, s};
	int periods[] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, true, &comm_slaves);

	MPI_Status status;

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

		//		printf("Master: %d \n", rank);

		for (r = 0; r < size; r++)
		{
			if (r == master)
				continue;

			//			int d;
			//			MPI_Send(&d, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
			int x = r % s;
			int y = r / s;
			MPI_Send(&(input_image_data[y * local_height][x * local_width][0]), 1, local_image_t, r, 0, MPI_COMM_WORLD);

			//			printf("Send returned: %d \n", rank);
		}

		/* Receive processed output from slave nodes. */

		for (r = 0; r < size; r++)
		{
			if (r == master)
				continue;

			int x = r % s;
			int y = r / s;
			MPI_Recv(&(input_image_data[y * local_height][x * local_width][0]), 1, local_image_t, r, 0, MPI_COMM_WORLD, &status);
		}

		write_channels(input_image_data, HEIGHT, WIDTH);

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

		/* Create border datatypes. */

		MPI_Datatype row_t, column_t, corner_t;
		MPI_Type_vector(B, local_width * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &row_t);
		MPI_Type_vector(local_height, B * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &column_t);
		MPI_Type_vector(B, B * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &corner_t);
		MPI_Type_commit(&row_t);
		MPI_Type_commit(&column_t);
		MPI_Type_commit(&corner_t);

		/* Receive inner image. */

		//		printf("Slave: %d \n", rank);
		//		int d;
		//		MPI_Recv(&d, 1, MPI_INT, master, 0, MPI_COMM_WORLD, &status);
		//		MPI_Recv(&(output_image_data[1][1][0]), local_height * local_width * CHANNELS, MPI_FLOAT, master, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&(output_image_data[B][B][0]), 1, inner_image_t, master, 0, MPI_COMM_WORLD, &status);
		//		printf("Received at %d.\n", rank);
		// int received;
		// MPI_Get_count(&status, MPI_FLOAT, &received);
		// printf("Rank: %d, North: %d.\n", rank, received);

		/* Get neighbors. */

		int r_n, r_s, r_e, r_w;
		int r_ne, r_nw, r_se, r_sw;

		get_neighbors(comm_slaves, &r_n, &r_s, &r_e, &r_w, &r_nw, &r_se, &r_ne, &r_sw);

		//		printf("rank: %d\n", rank);
		//		printf("north: %d, south: %d, east: %d, west: %d\n", north, south, east, west);
		//		printf("nw: %d, se: %d, ne: %d, sw: %d\n", nw, se, ne, sw);

		bool even_row = in_even_row(rank, comm_slaves);
		bool even_column = in_even_column(rank, comm_slaves);

		double start, finish, elapsed, max_elapsed;

		MPI_Barrier(comm_slaves);
		start = MPI_Wtime();

		/* Apply filter. */

		for (n = 0; n < ITERATIONS; n++)
		{
			if (even_row)
			{
				if (r_s != MPI_PROC_NULL) // sendrecv south
				{
					MPI_Send(&(output_image_data[local_height][B][0]), 1, row_t, r_s, 0, comm_slaves);
					MPI_Recv(&(output_image_data[local_height + B][B][0]), 1, row_t, r_s, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = B; j < B + local_width; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
				}

				if (r_n != MPI_PROC_NULL) // sendrecv north
				{
					MPI_Send(&(output_image_data[B][B][0]), 1, row_t, r_n, 0, comm_slaves);
					MPI_Recv(&(output_image_data[0][B][0]), 1, row_t, r_n, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = B; j < B + local_width; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][j][c];
				}

				if (r_se != MPI_PROC_NULL) // sendrecv southeast
				{
					MPI_Send(&(output_image_data[local_height][local_width][0]), 1, corner_t, r_se, 0, comm_slaves);
					MPI_Recv(&(output_image_data[local_height + B][local_width + B][0]), 1, corner_t, r_se, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][B + local_width - 1][c];
				}

				if (r_nw != MPI_PROC_NULL) // sendrecv northwest
				{
					MPI_Send(&(output_image_data[B][B][0]), 1, corner_t, r_nw, 0, comm_slaves);
					MPI_Recv(&(output_image_data[0][0][0]), 1, corner_t, r_nw, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][B][c];
				}

				if (r_sw != MPI_PROC_NULL) // sendrecv southwest
				{
					MPI_Send(&(output_image_data[local_height][B][0]), 1, corner_t, r_sw, 0, comm_slaves);
					MPI_Recv(&(output_image_data[local_height + B][0][0]), 1, corner_t, r_sw, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][B][c];
				}

				if (r_ne != MPI_PROC_NULL) // sendrecv northeast
				{
					MPI_Send(&(output_image_data[B][local_width][0]), 1, corner_t, r_ne, 0, comm_slaves);
					MPI_Recv(&(output_image_data[0][local_width + B][0]), 1, corner_t, r_ne, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][B + local_width - 1][c];
				}
			} else // odd row
			{
				if (r_n != MPI_PROC_NULL) // sendrecv north
				{
					MPI_Recv(&(output_image_data[0][B][0]), 1, row_t, r_n, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][B][0]), 1, row_t, r_n, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = B; j < B + local_width; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][j][c];
				}

				if (r_s != MPI_PROC_NULL) // sendrecv south
				{
					MPI_Recv(&(output_image_data[local_height + B][B][0]), 1, row_t, r_s, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[local_height][B][0]), 1, row_t, r_s, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = B; j < B + local_width; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
				}

				if (r_nw != MPI_PROC_NULL) // sendrecv northwest
				{
					MPI_Recv(&(output_image_data[0][0][0]), 1, corner_t, r_nw, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][B][0]), 1, corner_t, r_nw, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][B][c];
				}

				if (r_se != MPI_PROC_NULL) // sendrecv southeast
				{
					MPI_Recv(&(output_image_data[local_height + B][local_width + B][0]), 1, corner_t, r_se, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[local_height][local_width][0]), 1, corner_t, r_se, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][B + local_width - 1][c];
				}

				if (r_ne != MPI_PROC_NULL) // sendrecv northeast
				{
					MPI_Recv(&(output_image_data[0][local_width + B][0]), 1, corner_t, r_ne, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][local_width][0]), 1, corner_t, r_ne, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B][B + local_width - 1][c];
				}

				if (r_sw != MPI_PROC_NULL) // sendrecv southwest
				{
					MPI_Recv(&(output_image_data[local_height + B][0][0]), 1, corner_t, r_sw, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[local_height][B][0]), 1, corner_t, r_sw, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[B + local_height - 1][B][c];
				}
			}

			if (even_column)
			{
				if (r_e != MPI_PROC_NULL) // sendrecv east
				{
					MPI_Send(&(output_image_data[B][local_width][0]), 1, column_t, r_e, 0, comm_slaves);
					MPI_Recv(&(output_image_data[B][local_width + B][0]), 1, column_t, r_e, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = B; i < B + local_height; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
				}

				if (r_w != MPI_PROC_NULL) // sendrecv west
				{
					MPI_Send(&(output_image_data[B][B][0]), 1, column_t, r_w, 0, comm_slaves);
					MPI_Recv(&(output_image_data[B][0][0]), 1, column_t, r_w, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = B; i < B + local_height; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[i][B][c];
				}

			} else // odd column
			{
				if (r_w != MPI_PROC_NULL) // sendrecv west
				{
					MPI_Recv(&(output_image_data[B][0][0]), 1, column_t, r_w, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][B][0]), 1, column_t, r_w, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = B; i < B + local_height; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[i][B][c];
				}

				if (r_e != MPI_PROC_NULL) // sendrecv east
				{
					MPI_Recv(&(output_image_data[B][local_width + B][0]), 1, column_t, r_e, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][local_width][0]), 1, column_t, r_e, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = B; i < B + local_height; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
				}
			}

			/* Use previous output as new input. */

			float (**tmp)[CHANNELS];
			tmp = input_image_data;
			input_image_data = output_image_data;
			output_image_data = tmp;

			apply_inner_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);

			apply_outer_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);
		}

		finish = MPI_Wtime();

		elapsed = finish - start;

		MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm_slaves);

		if (rank == 0)
			printf("Time elapsed: %lf seconds\n", max_elapsed);

		MPI_Send(&(output_image_data[B][B][0]), 1, inner_image_t, master, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return (EXIT_SUCCESS);
}
