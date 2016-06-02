#ifndef MPI_UTILS_H_8BRBW6XS
#define MPI_UTILS_H_8BRBW6XS

/*
 * =====================================================================================
 *
 *       Filename:  mpi_utils.h
 *
 *    Description:  Global MPI Utils  that contains useful information.
 *
 *        Version:  1.0
 *        Created:  04/14/2016 07:20:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

int rank ;
int nproc ;

void initMPI(int argc, char * argv[]) ;

#endif /* end of include guard: MPI_UTILS_H_8BRBW6XS */
