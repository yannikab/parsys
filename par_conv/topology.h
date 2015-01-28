/* 
 * File:   topology.h
 * Author: John
 *
 * Created on January 21, 2015, 10:40 AM
 */

#ifndef TOPOLOGY_H
#define	TOPOLOGY_H

#include <mpi.h>

void get_neighbors(MPI_Comm comm, int *n_p, int *s_p, int *e_p, int *w_p, int *nw_p, int *se_p, int *ne_p, int *sw_p);
bool in_even_row(int rank, MPI_Comm comm);
bool in_even_column(int rank, MPI_Comm comm);

#endif	/* TOPOLOGY_H */
