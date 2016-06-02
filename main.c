/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  Main function for simulating diffusion equation using FFTW. This
 *                  version of the program tracks the value of point A as a function 
 *                  of time. A is the bottom right corner of the grid such that
 *                  x = 3/4 gridDimX and y = gridDimY/4.
 *
 *        Version:  1.0
 *        Created:  04/14/2016 11:15:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include <stdio.h>
#include "getopt.h"
#include "stdlib.h"
#include "grid.h"
#include "mpi_utils.h"
#include "mpi.h"
#include "time_utils.h"

int main(int argc, char *argv[]) {
	
	initMPI(argc, argv) ;
	int gridDimX = 8 ; // Grid dimension X must be divisable by 4.
	int gridDimY = 8 ; // Grid dimension Y must be divisable by 4.
	int choice;
	while (1) {
		choice = getopt( argc, argv, "n:m:pi:t:b");
		if (choice == -1)
			break;
		switch(choice) {
			case 'n':
				gridDimY = atoi(optarg) ;
				break;
			case 'm':
				gridDimX = atoi(optarg) ;
				break;
			default:
				// Not sure how to get here... 
				printf("USAGE ./diff -n YDIM -m XDIM \n") ; 
				return EXIT_FAILURE;
		}
	}
	// Ensures grid dimensions are multiples of 4. //
	if (gridDimY%4 != 0 && gridDimY != 0) {
		fprintf(stderr, "Dimension Y isn't a multiple of 4\n");
		exit(-1) ;
	}
	if (gridDimX%4 != 0 && gridDimX != 0) {
		fprintf(stderr, "Dimension X isn't a multiple of 4\n");
		exit(-1) ;
	}
	
	Grid * grid = makeNewGrid(gridDimX,gridDimY) ;
	initialiseGrid(grid) ;
	simulateGrid(grid) ;
	freeGrid(grid) ;
	MPI_Finalize() ;

	return EXIT_SUCCESS ;
}



