/* 
 * File:   topology.c
 * Author: John
 *
 * Created on January 21, 2015, 10:40 AM
 */

#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>

#include "topology.h"

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
