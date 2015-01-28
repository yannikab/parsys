/* 
 * File:   main_async.c
 * Author: John
 *
 * Created on January 21, 2015, 10:45 AM
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <assert.h>

#include <mpi.h>
#include <omp.h>

#include "settings.h"
#include "2d_malloc.h"
#include "file_io.h"
#include "topology.h"
#include "filter.h"

/*
 * 
 */
int main_async_omp(int argc, char** argv)
{
	int size, rank;

	int i, j, c;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//	printf("%d %d\n", size, rank);

	if (argc != 5)
	{
		if (rank == 0)
			printf("Usage: %s <iterations> <convergence> <rows> <columns>\n", argv[0]);
		return (EXIT_FAILURE);
	}

	int iterations = atoi(argv[1]);
	int convergence = atoi(argv[2]);
	int rows = atoi(argv[3]);
	int columns = atoi(argv[4]);

	if (rows * columns != size - 1)
		return (EXIT_FAILURE);

	//	printf("%d %d\n", rows, columns);

	int width = WIDTH / columns;
	int height = HEIGHT / rows;

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
		printf("main_async_omp()\n");
		printf("Iterations: %d, Convergence: %d\n", iterations, convergence);

		// printf("\nwidth: %d, height: %d\n", width, height);

		unsigned char (**image_buffer)[CHANNELS];

		/* Read input image. */

		read_image((unsigned char ***) &image_buffer);

		/* Allocate memory for coordinates map. */

		int (*coords)[2];

		coords = malloc(rows * columns * sizeof (*coords));

		/* Receive and store cartesian coordinates from each slave process. */

		unsigned int r;

		// printf("\n");
		for (r = 0; r < rows * columns; r++)
		{
			MPI_Recv(coords[r], 2, MPI_INT, r + 1, 0, MPI_COMM_WORLD, &status);
			// printf("rank %d row:%d col:%d\n", r + 1, coords[r][0], coords[r][1]);
		}
		// printf("\n");

		/* Create "subarray" datatype. */

		MPI_Datatype local_image_t;
		MPI_Type_vector(height, width * CHANNELS, WIDTH * CHANNELS, MPI_UNSIGNED_CHAR, &local_image_t);
		MPI_Type_commit(&local_image_t);

		/* Send each process its corresponding subarray. */

		for (r = 0; r < rows * columns; r++)
			MPI_Send(&(image_buffer[coords[r][0] * height][coords[r][1] * width][0]), 1, local_image_t, r + 1, 0, MPI_COMM_WORLD);

		/* Receive processed output from slave processes. */

		for (r = 0; r < rows * columns; r++)
			MPI_Recv(&(image_buffer[coords[r][0] * height][coords[r][1] * width][0]), 1, local_image_t, r + 1, 0, MPI_COMM_WORLD, &status);

		write_channels(image_buffer, HEIGHT, WIDTH);

	} else // slaves
	{
		/* Create cartesian coordinates communicator for slaves. */

		MPI_Comm comm_slaves;

		int ndims = 2;
		int dims[] = {rows, columns};
		int periods[] = {0, 0};

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

		alloc_uchar_array((unsigned char ***) &local_buffer, B + height + B, B + width + B, CHANNELS);

		/* Create local buffer datatype for master receive/send. */

		MPI_Datatype local_buffer_t;
		MPI_Type_vector(height, width * CHANNELS, (B + width + B) * CHANNELS, MPI_UNSIGNED_CHAR, &local_buffer_t);
		MPI_Type_commit(&local_buffer_t);

		/* Receive local image data from master. */

		// int received;
		MPI_Recv(&(local_buffer[B][B][0]), 1, local_buffer_t, master, 0, MPI_COMM_WORLD, &status);
		// MPI_Get_count(&status, MPI_FLOAT, &received);
		// printf("Rank: %d, Received: %d.\n", rank, received);

		/* Allocate two float arrays for image processing (input, output). */

		float (**prev_image)[CHANNELS];
		float (**curr_image)[CHANNELS];

		alloc_float_array((float ***) &prev_image, B + height + B, B + width + B, CHANNELS);

		alloc_float_array((float ***) &curr_image, B + height + B, B + width + B, CHANNELS);

		/* Copy receive/send buffer data, converting to float for applying arithmetic operations. */

		for (i = 0; i < B + height + B; i++)
			for (j = 0; j < B + width + B; j++)
				for (c = 0; c < CHANNELS; c++)
					curr_image[i][j][c] = (float) local_buffer[i][j][c];

		/* Get neighbouring process ranks. */

		int r_n, r_s, r_e, r_w;
		int r_ne, r_nw, r_se, r_sw;

		get_neighbours(comm_slaves, &r_n, &r_s, &r_e, &r_w, &r_nw, &r_se, &r_ne, &r_sw);

		/* Determine if process is in odd or even row/column. */

		bool even_row = in_even_row(slave_rank, comm_slaves);
		bool even_column = in_even_column(slave_rank, comm_slaves);

		/* Create border datatypes for communication between slaves. */

		MPI_Datatype row_t, column_t, corner_t;

		MPI_Type_vector(B, width * CHANNELS, (B + width + B) * CHANNELS, MPI_FLOAT, &row_t);
		MPI_Type_commit(&row_t);

		MPI_Type_vector(height, B * CHANNELS, (B + width + B) * CHANNELS, MPI_FLOAT, &column_t);
		MPI_Type_commit(&column_t);

		MPI_Type_vector(B, B * CHANNELS, (B + width + B) * CHANNELS, MPI_FLOAT, &corner_t);
		MPI_Type_commit(&corner_t);

		/* Set up timing. */

		double start, finish, elapsed, min_elapsed, max_elapsed, avg_elapsed;

		MPI_Barrier(comm_slaves);
		start = MPI_Wtime();

		/* Apply filter. */

		MPI_Request sends[8];
		MPI_Request recvs[8];
		MPI_Status send_status[8];
		MPI_Status recv_status[8];
		unsigned int p, q;

		bool converged = false;

#pragma omp parallel private (i,j,c)
		{
#pragma omp single
			if (slave_rank == 0)
				printf("Threads: %d\n", omp_get_num_threads());

			unsigned int n;

			for (n = 0; !converged && (iterations == 0 || n < iterations); n++)
			{
#pragma omp master // master thread handles MPI messaging
				{
					/* Reset send/recv indexes. */

					p = 0;
					q = 0;

					/* Send / receive vertical data. */

					if (r_s != MPI_PROC_NULL) // sendrecv south
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[height][B][0]), 1, row_t, r_s, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[height + B][B][0]), 1, row_t, r_s, 0, comm_slaves, &recvs[q++]);
						} else // odd row
						{
							MPI_Irecv(&(curr_image[height + B][B][0]), 1, row_t, r_s, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[height][B][0]), 1, row_t, r_s, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						for (i = height + B; i < height + 2 * B; i++)
							for (j = B; j < B + width; j++)
								for (c = 0; c < CHANNELS; c++)
									curr_image[i][j][c] = curr_image[B + height - 1][j][c];
					}

					if (r_n != MPI_PROC_NULL) // sendrecv north
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[B][B][0]), 1, row_t, r_n, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[0][B][0]), 1, row_t, r_n, 0, comm_slaves, &recvs[q++]);
						} else // odd row
						{
							MPI_Irecv(&(curr_image[0][B][0]), 1, row_t, r_n, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[B][B][0]), 1, row_t, r_n, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						for (i = 0; i < B; i++)
							for (j = B; j < B + width; j++)
								for (c = 0; c < CHANNELS; c++)
									curr_image[i][j][c] = curr_image[B][j][c];
					}

					/* Send / receive horizontal data. */

					if (r_e != MPI_PROC_NULL) // sendrecv east
					{
						if (even_column)
						{
							MPI_Isend(&(curr_image[B][width][0]), 1, column_t, r_e, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[B][width + B][0]), 1, column_t, r_e, 0, comm_slaves, &recvs[q++]);
						} else // odd column
						{
							MPI_Irecv(&(curr_image[B][width + B][0]), 1, column_t, r_e, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[B][width][0]), 1, column_t, r_e, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						for (i = B; i < B + height; i++)
							for (j = width + B; j < width + 2 * B; j++)
								for (c = 0; c < CHANNELS; c++)
									curr_image[i][j][c] = curr_image[i][B + width - 1][c];
					}

					if (r_w != MPI_PROC_NULL) // sendrecv west
					{
						if (even_column)
						{
							MPI_Isend(&(curr_image[B][B][0]), 1, column_t, r_w, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[B][0][0]), 1, column_t, r_w, 0, comm_slaves, &recvs[q++]);
						} else // odd column
						{
							MPI_Irecv(&(curr_image[B][0][0]), 1, column_t, r_w, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[B][B][0]), 1, column_t, r_w, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						for (i = B; i < B + height; i++)
							for (j = 0; j < B; j++)
								for (c = 0; c < CHANNELS; c++)
									curr_image[i][j][c] = curr_image[i][B][c];
					}

					/* Send / receive diagonal data. */

					if (r_se != MPI_PROC_NULL) // sendrecv southeast
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[height][width][0]), 1, corner_t, r_se, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[height + B][width + B][0]), 1, corner_t, r_se, 0, comm_slaves, &recvs[q++]);
						} else
						{ // odd row
							MPI_Irecv(&(curr_image[height + B][width + B][0]), 1, corner_t, r_se, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[height][width][0]), 1, corner_t, r_se, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						if (r_s == MPI_PROC_NULL && r_e == MPI_PROC_NULL) // use corner data
							for (i = height + B; i < height + 2 * B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B + height - 1][B + width - 1][c];
					}

					if (r_nw != MPI_PROC_NULL) // sendrecv northwest
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[B][B][0]), 1, corner_t, r_nw, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[0][0][0]), 1, corner_t, r_nw, 0, comm_slaves, &recvs[q++]);
						} else
						{ // odd row
							MPI_Irecv(&(curr_image[0][0][0]), 1, corner_t, r_nw, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[B][B][0]), 1, corner_t, r_nw, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						if (r_n == MPI_PROC_NULL && r_w == MPI_PROC_NULL) // use corner data
							for (i = 0; i < B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B][B][c];
					}

					if (r_sw != MPI_PROC_NULL) // sendrecv southwest
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[height][B][0]), 1, corner_t, r_sw, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[height + B][0][0]), 1, corner_t, r_sw, 0, comm_slaves, &recvs[q++]);
						} else
						{ // odd row
							MPI_Irecv(&(curr_image[height + B][0][0]), 1, corner_t, r_sw, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[height][B][0]), 1, corner_t, r_sw, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						if (r_s == MPI_PROC_NULL && r_w == MPI_PROC_NULL) // use corner data
							for (i = height + B; i < height + 2 * B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B + height - 1][B][c];
					}

					if (r_ne != MPI_PROC_NULL) // sendrecv northeast
					{
						if (even_row)
						{
							MPI_Isend(&(curr_image[B][width][0]), 1, corner_t, r_ne, 0, comm_slaves, &sends[p++]);
							MPI_Irecv(&(curr_image[0][width + B][0]), 1, corner_t, r_ne, 0, comm_slaves, &recvs[q++]);
						} else
						{ // odd row
							MPI_Irecv(&(curr_image[0][width + B][0]), 1, corner_t, r_ne, 0, comm_slaves, &recvs[q++]);
							MPI_Isend(&(curr_image[B][width][0]), 1, corner_t, r_ne, 0, comm_slaves, &sends[p++]);
						}
					} else // fill border buffer with edge image data
					{
						if (r_n == MPI_PROC_NULL && r_e == MPI_PROC_NULL) // use corner data
							for (i = 0; i < B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B][B + width - 1][c];
					}

					assert(p == q); // number of sends should match number of recvs
				}

				/* Apply inner filter using omp for, does not require having border data available. */

				apply_inner_filter_openmp(prev_image, curr_image, B + height + B, B + width + B);

				/* Wait for recvs, master thread handles MPI messaging. */
#pragma omp master
				MPI_Waitall(q, recvs, recv_status);

				/* Handle diagonal border data that require recvs to have been completed, master thread can safely do it. */
#pragma omp master
				{
					if (r_se == MPI_PROC_NULL) // southeast
					{
						if (r_s != MPI_PROC_NULL) // get data from south (received)
							for (i = height + B; i < height + 2 * B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[i][B + width - 1][c];
						else if (r_e != MPI_PROC_NULL) // get data from east (received)
							for (i = height + B; i < height + 2 * B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B + height - 1][j][c];
					}

					if (r_nw == MPI_PROC_NULL) // northwest
					{
						if (r_n != MPI_PROC_NULL) // get data from north (received)
							for (i = 0; i < B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[i][B][c];
						else if (r_w != MPI_PROC_NULL) // get data from west (received)
							for (i = 0; i < B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B][j][c];
					}

					if (r_sw == MPI_PROC_NULL) // southwest
					{
						if (r_s != MPI_PROC_NULL) // get data from south (received)
							for (i = height + B; i < height + 2 * B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[i][B][c];
						else if (r_w != MPI_PROC_NULL) // get data from west (received)
							for (i = height + B; i < height + 2 * B; i++)
								for (j = 0; j < B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B + height - 1][j][c];
					}

					if (r_ne == MPI_PROC_NULL) // northeast
					{
						if (r_n != MPI_PROC_NULL) // get data from north (received)
							for (i = 0; i < B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[i][B + width - 1][c];
						else if (r_e != MPI_PROC_NULL) // get data from east (received)
							for (i = 0; i < B; i++)
								for (j = width + B; j < width + 2 * B; j++)
									for (c = 0; c < CHANNELS; c++)
										curr_image[i][j][c] = curr_image[B][j][c];
					}
				}

				/* Apply outer filter, requires having all border data available. Master thread can safely do it. */
#pragma omp master
				apply_outer_filter(prev_image, curr_image, B + height + B, B + width + B);

				/* Wait for sends before we switch buffers. Master thread handles MPI messaging. */
#pragma omp master
				MPI_Waitall(p, sends, send_status);

				/* Switch current / previous image buffers. Master thread does it after sends have completed. */
#pragma omp master
				{
					float (**tmp)[CHANNELS];
					tmp = prev_image;
					prev_image = curr_image;
					curr_image = tmp;
				}

				/* Check for convergence. */
#pragma omp master
				if (convergence > 0 && n % convergence == 0)
				{
					int identical = images_identical(curr_image, prev_image, B + height + B, B + width + B) ? 1 : 0;
					int all_identical = 0;

					MPI_Allreduce(&identical, &all_identical, 1, MPI_INT, MPI_LAND, comm_slaves);

					if (all_identical)
					{
						if (slave_rank == 0)
							printf("Filter has converged after %d iterations.\n", n);

						converged = true; // break not allowed from OpenMP structured block
					}
				}

				/* Threads should not move to next iteration until master thread has finished switching buffers. */
#pragma omp barrier
			}
		}

		finish = MPI_Wtime();

		elapsed = finish - start;

		MPI_Reduce(&elapsed, &min_elapsed, 1, MPI_DOUBLE, MPI_MIN, 0, comm_slaves);
		MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm_slaves);
		MPI_Reduce(&elapsed, &avg_elapsed, 1, MPI_DOUBLE, MPI_SUM, 0, comm_slaves);
		avg_elapsed /= rows * columns;

		// printf("Rank %d time elapsed: %lf seconds\n", rank, elapsed);
		// MPI_Barrier(comm_slaves);

		if (slave_rank == 0)
			printf("Min: %lf, Max: %lf, Avg: %lf seconds\n", min_elapsed, max_elapsed, avg_elapsed);

		/* Convert float data back to byte for sending to master process. */

		for (i = 0; i < B + height + B; i++)
			for (j = 0; j < B + width + B; j++)
				for (c = 0; c < CHANNELS; c++)
					local_buffer[i][j][c] = (unsigned char) curr_image[i][j][c];

		/* Send results back to master. */

		MPI_Send(&(local_buffer[B][B][0]), 1, local_buffer_t, master, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return (EXIT_SUCCESS);
}