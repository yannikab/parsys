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

#include <unistd.h>

#include <mpi.h>

#define WIDTH 1920
#define HEIGHT 2520

#define CHANNELS 3

#define B 1

static const float filter[3][3] = {
	{0.0625f, 0.125f, 0.0625f},
	{0.1250f, 0.250f, 0.1250f},
	{0.0625f, 0.125f, 0.0625f},
};

#define ITERATIONS 15

#define INFILENAME "waterfall_1920_2520.raw"
#define OUTFILENAME "waterfall_1920_2520_out.raw"

#define INDIR "/tmp/stud0781/input/1.00x/"
#define OUTDIR "/tmp/stud0781/parallel/"

#define STRSIZE 512

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
	assert(array != NULL);

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
	assert(array != NULL);

	free(&((*array)[0][0]));
	free(*array);
	*array = NULL;
}

typedef enum {
	INPUT,
	OUTPUT,
}
file_type;

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

	printf("\n");
	char command[STRSIZE];
	for (c = 0; ok && c < CHANNELS; c++)
	{
		// sprintf(command, "md5sum %s %s.tiff", out_filename[c], out_filename[c]);
		sprintf(command, "md5sum %s", out_filename[c]);
		// printf("%s\n", command);
		ok = system(command) == 0;
	}

	/* Convert output files to tiff format (ImageMagick). */

	printf("\n");
	for (c = 0; ok && c < CHANNELS; c++)
	{
		// sprintf(command, "raw2tiff -l %d -w %d %s %s.tiff", height, width, out_filename[c], out_filename[c]);
		sprintf(command, "convert -depth 8 -size %dx%d gray:%s -compress lzw %s.tiff", width, height, out_filename[c], out_filename[c]);
		printf("%s\n", command);
		ok = system(command) == 0;
	}

	/* Merge individual channel tiffs to a single tiff (ImageMagick). */

	sprintf(command, "convert");
	if (ok)
	{
		for (c = 0; c < CHANNELS; c++)
		{
			strcat(command, " ");
			strcat(command, out_filename[c]);
			strcat(command, ".tiff");
		}

		strcat(command, " -combine ");
		strcat(command, OUTDIR);
		strcat(command, OUTFILENAME);
		strcat(command, ".tiff");

		printf("%s\n", command);

		ok = system(command) == 0;
	}

	printf("\n");

	/* Free memory allocated for filenames. */

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

void get_neighbours(MPI_Comm comm, int *n_p, int *s_p, int *e_p, int *w_p, int *nw_p, int *se_p, int *ne_p, int *sw_p)
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

	int size, rank;

	int i, j, c;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	printf("%d %d\n", size, rank);

	if (argc != 3)
		return (EXIT_FAILURE);

	int rows = atoi(argv[1]);
	int columns = atoi(argv[2]);

	if (rows * columns != size - 1)
		return (EXIT_FAILURE);

	//	printf("%d %d\n", rows, columns);

	int local_width = WIDTH / columns;
	int local_height = HEIGHT / rows;

	/* Create slaves communicator excluding rank 0 (master). */

	MPI_Group MPI_GROUP_WORLD, group_excl;
	MPI_Comm comm_excl;

	MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
	int ranks[] = {0};
	MPI_Group_excl(MPI_GROUP_WORLD, 1, ranks, &group_excl);
	MPI_Comm_create(MPI_COMM_WORLD, group_excl, &comm_excl);

	MPI_Status status;

	if (rank == 0) // master
	{
		//		char cwd[1024];
		//		if (getcwd(cwd, sizeof (cwd)) != NULL)
		//			fprintf(stdout, "Current working dir: %s\n", cwd);
		//		else
		//			perror("getcwd() error");

		printf("\nlocal_width: %d, local_height: %d\n", local_width, local_height);

		unsigned char (**image_buffer)[CHANNELS];

		/* Read input image. */

		read_image((unsigned char ***) &image_buffer);

		/* Allocate memory for coordinates map. */

		int (*coords)[2];

		coords = malloc(rows * columns * sizeof (*coords));

		/* Receive and store cartesian coordinates from each slave process. */

		unsigned int r;

		printf("\n");
		for (r = 0; r < rows * columns; r++)
		{
			// printf("size:%d, rank:%d\n", size, r);
			MPI_Recv(coords[r], 2, MPI_INT, r + 1, 0, MPI_COMM_WORLD, &status);
			printf("rank %d row:%d col:%d\n", r + 1, coords[r][0], coords[r][1]);
		}
		printf("\n");

		/* Create "subarray" datatype. */

		MPI_Datatype local_image_t;
		MPI_Type_vector(local_height, local_width * CHANNELS, WIDTH * CHANNELS, MPI_UNSIGNED_CHAR, &local_image_t);
		MPI_Type_commit(&local_image_t);

		/* Send each process its corresponding subarray. */

		for (r = 0; r < rows * columns; r++)
			MPI_Send(&(image_buffer[coords[r][0] * local_height][coords[r][1] * local_width][0]), 1, local_image_t, r + 1, 0, MPI_COMM_WORLD);

		/* Receive processed output from slave processes. */

		for (r = 0; r < rows * columns; r++)
			MPI_Recv(&(image_buffer[coords[r][0] * local_height][coords[r][1] * local_width][0]), 1, local_image_t, r + 1, 0, MPI_COMM_WORLD, &status);

		write_channels(image_buffer, HEIGHT, WIDTH);

	} else // slaves
	{
		/* Create cartesian coordinates communicator for slaves. */

		MPI_Comm comm_slaves;

		int ndims = 2;
		int dims[] = {rows, columns};
		int periods[] = {0, 0};

		// printf("rank:%d, comm_slaves:%d, comm_slaves_cart:%d\n", rank, comm_slaves, comm_slaves_cart);
		MPI_Cart_create(comm_excl, ndims, dims, periods, true, &comm_slaves);

		/* Determine rank in slaves communicator. */

		int slave_rank;
		MPI_Comm_rank(comm_slaves, &slave_rank);
		// printf("Rank: %d, Slave rank: %d\n", rank, slave_rank);

		/* Send cartesian coordinates to master. */

		int master = 0;

		int cart_coords[2];

		MPI_Cart_coords(comm_slaves, slave_rank, 2, cart_coords);
		MPI_Send(cart_coords, 2, MPI_INT, master, 0, MPI_COMM_WORLD);

		/* Allocate memory for local buffer. */

		unsigned char (**local_buffer)[CHANNELS];

		alloc_uchar_array((unsigned char ***) &local_buffer, B + local_height + B, B + local_width + B, CHANNELS);

		/* Create local buffer datatype for master receive/send. */

		MPI_Datatype local_buffer_t;
		MPI_Type_vector(local_height, local_width * CHANNELS, (B + local_width + B) * CHANNELS, MPI_UNSIGNED_CHAR, &local_buffer_t);
		MPI_Type_commit(&local_buffer_t);

		/* Receive local image data from master. */

		// int received;
		// printf("Slave: %d \n", rank);
		MPI_Recv(&(local_buffer[B][B][0]), 1, local_buffer_t, master, 0, MPI_COMM_WORLD, &status);
		// MPI_Get_count(&status, MPI_FLOAT, &received);
		// printf("Rank: %d, Received: %d.\n", rank, received);

		/* Allocate two float arrays for image processing (input, output). */

		float (**input_image_data)[CHANNELS];
		float (**output_image_data)[CHANNELS];

		alloc_float_array((float ***) &input_image_data, B + local_height + B, B + local_width + B, CHANNELS);

		alloc_float_array((float ***) &output_image_data, B + local_height + B, B + local_width + B, CHANNELS);

		/* Copy receive/send buffer data, converting to float for applying arithmetic operations. */

		for (i = 0; i < B + local_height + B; i++)
			for (j = 0; j < B + local_width + B; j++)
				for (c = 0; c < CHANNELS; c++)
					output_image_data[i][j][c] = (float) local_buffer[i][j][c];

		/* Get neighbouring process ranks. */

		int r_n, r_s, r_e, r_w;
		int r_ne, r_nw, r_se, r_sw;

		get_neighbours(comm_slaves, &r_n, &r_s, &r_e, &r_w, &r_nw, &r_se, &r_ne, &r_sw);

		//		if (rank == 3)
		//		{
		//			printf("rank: %d\n", rank);
		//			printf("north: %d, south: %d, east: %d, west: %d\n", r_n, r_s, r_e, r_w);
		//			printf("nw: %d, se: %d, ne: %d, sw: %d\n", r_nw, r_se, r_ne, r_sw);
		//		}

		/* Determine if process is in odd or even row/column. */
		
		bool even_row = in_even_row(slave_rank, comm_slaves);
		bool even_column = in_even_column(slave_rank, comm_slaves);

		/* Create border datatypes for communication between slaves. */

		MPI_Datatype row_t, column_t, corner_t;

		MPI_Type_vector(B, local_width * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &row_t);
		MPI_Type_commit(&row_t);

		MPI_Type_vector(local_height, B * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &column_t);
		MPI_Type_commit(&column_t);

		MPI_Type_vector(B, B * CHANNELS, (B + local_width + B) * CHANNELS, MPI_FLOAT, &corner_t);
		MPI_Type_commit(&corner_t);

		/* Set up timing. */

		double start, finish, elapsed, min_elapsed, max_elapsed, avg_elapsed;

		MPI_Barrier(comm_slaves);
		start = MPI_Wtime();

		/* Apply filter. */

		unsigned int n;

		for (n = 0; n < ITERATIONS; n++)
		{
			/* Send / receive vertical data. */

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
			}

			/* Send / receive horizontal data. */

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

			/* Send / receive diagonal data. */

			if (even_row)
			{
				if (r_se != MPI_PROC_NULL) // sendrecv southeast
				{
					MPI_Send(&(output_image_data[local_height][local_width][0]), 1, corner_t, r_se, 0, comm_slaves);
					MPI_Recv(&(output_image_data[local_height + B][local_width + B][0]), 1, corner_t, r_se, 0, comm_slaves, &status);
				} else // fill border buffer with edge image data
				{
					for (i = local_height + B; i < local_height + 2 * B; i++)
						for (j = local_width + B; j < local_width + 2 * B; j++)
							for (c = 0; c < CHANNELS; c++)
								if (r_s != MPI_PROC_NULL) // get data from south (received)
									output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
								else if (r_e != MPI_PROC_NULL) // get data from east (received)
									output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
								else // use corner data
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
								if (r_n != MPI_PROC_NULL) // get data from north (received)
									output_image_data[i][j][c] = output_image_data[i][B][c];
								else if (r_w != MPI_PROC_NULL) // get data from west (received)
									output_image_data[i][j][c] = output_image_data[B][j][c];
								else // use corner data
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
								if (r_s != MPI_PROC_NULL) // get data from south (received)
									output_image_data[i][j][c] = output_image_data[i][B][c];
								else if (r_w != MPI_PROC_NULL) // get data from west (received)
									output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
								else // use corner data
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
								if (r_n != MPI_PROC_NULL) // get data from north (received)
									output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
								else if (r_e != MPI_PROC_NULL) // get data from east (received)
									output_image_data[i][j][c] = output_image_data[B][j][c];
								else // use corner data
									output_image_data[i][j][c] = output_image_data[B][B + local_width - 1][c];
				}
			} else // odd row
			{
				if (r_nw != MPI_PROC_NULL) // sendrecv northwest
				{
					MPI_Recv(&(output_image_data[0][0][0]), 1, corner_t, r_nw, 0, comm_slaves, &status);
					MPI_Send(&(output_image_data[B][B][0]), 1, corner_t, r_nw, 0, comm_slaves);
				} else // fill border buffer with edge image data
				{
					for (i = 0; i < B; i++)
						for (j = 0; j < B; j++)
							for (c = 0; c < CHANNELS; c++)
								if (r_n != MPI_PROC_NULL) // get data from north (received)
									output_image_data[i][j][c] = output_image_data[i][B][c];
								else if (r_w != MPI_PROC_NULL) // get data from west (received)
									output_image_data[i][j][c] = output_image_data[B][j][c];
								else // use corner data
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
								if (r_s != MPI_PROC_NULL) // get data from south (received)
									output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
								else if (r_e != MPI_PROC_NULL) // get data from east (received)
									output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
								else // use corner data
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
								if (r_n != MPI_PROC_NULL) // get data from north (received)
									output_image_data[i][j][c] = output_image_data[i][B + local_width - 1][c];
								else if (r_e != MPI_PROC_NULL) // get data from east (received)
									output_image_data[i][j][c] = output_image_data[B][j][c];
								else // use corner data
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
								if (r_s != MPI_PROC_NULL) // get data from south (received)
									output_image_data[i][j][c] = output_image_data[i][B][c];
								else if (r_w != MPI_PROC_NULL) // get data from west (received)
									output_image_data[i][j][c] = output_image_data[B + local_height - 1][j][c];
								else // use corner data
									output_image_data[i][j][c] = output_image_data[B + local_height - 1][B][c];
				}
			}

			/* Use previous output as new input. */

			float (**tmp)[CHANNELS];
			tmp = input_image_data;
			input_image_data = output_image_data;
			output_image_data = tmp;

			/* Apply inner filter, does not require having border data available. */

			apply_inner_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);

			/* Apply outer filter, requires having border data available. */

			apply_outer_filter(output_image_data, input_image_data, B + local_height + B, B + local_width + B);
		}

		finish = MPI_Wtime();

		elapsed = finish - start;

		MPI_Reduce(&elapsed, &min_elapsed, 1, MPI_DOUBLE, MPI_MIN, 0, comm_slaves);
		MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm_slaves);
		MPI_Reduce(&elapsed, &avg_elapsed, 1, MPI_DOUBLE, MPI_SUM, 0, comm_slaves);
		avg_elapsed /= rows * columns;

		printf("Rank %d time elapsed: %lf seconds\n", rank, elapsed);

		MPI_Barrier(comm_slaves);

		if (slave_rank == 0)
			printf("Min: %lf, Max: %lf, Avg: %lf seconds\n", min_elapsed, max_elapsed, avg_elapsed);

		/* Convert float data back to byte for sending to master process. */

		for (i = 0; i < B + local_height + B; i++)
			for (j = 0; j < B + local_width + B; j++)
				for (c = 0; c < CHANNELS; c++)
					local_buffer[i][j][c] = (unsigned char) output_image_data[i][j][c];

		/* Send results back to master. */
		MPI_Send(&(local_buffer[B][B][0]), 1, local_buffer_t, master, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return (EXIT_SUCCESS);
}
