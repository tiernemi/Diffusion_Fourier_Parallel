/*
 * =====================================================================================
 *
 *       Filename:  mpi_utils.c
 *
 *    Description: Source file for mpi_utils 
 *
 *        Version:  1.0
 *        Created:  04/14/2016 09:43:31 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <mpi.h>
#include "mpi_utils.h"

void initMPI(int argc, char * argv[]) {
	MPI_Init(&argc,&argv) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
}

