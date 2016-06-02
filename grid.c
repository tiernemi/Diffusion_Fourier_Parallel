/*
 * =====================================================================================
 *
 *       Filename:  grid.c
 *
 *    Description:  Source file for grid object. 
 *
 *        Version:  1.0
 *        Created:  04/14/2016 11:24:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "grid.h"
#include "stdlib.h"
#include "stdio.h"
#include "fftw3-mpi.h"
#include "string.h"
#include "math.h"
#include "mpi_utils.h"

#define PI 3.14159265359

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  makeNewGrid
 *    Arguments:  int gridDim - Dimensions dimXdim of new grid.
 *      Returns:  Pointer to new grid object.
 *  Description:  Makes new grid object and sets all data points to 0.
 * =====================================================================================
 */
Grid * makeNewGrid(int gridDimX, int gridDimY) {
	fftw_mpi_init() ;
	Grid * newGrid = malloc(sizeof(Grid)) ;
	ptrdiff_t localDimY, globalOffset ;
	newGrid->allocScheme = fftw_mpi_local_size_2d(gridDimY, gridDimX, MPI_COMM_WORLD, &localDimY, &globalOffset);
	newGrid->localGridDimY = localDimY ;
	newGrid->localGridDimX = gridDimX ;
	newGrid->globalGridDimY = gridDimY ;
	newGrid->globalGridDimX = gridDimX ;
	newGrid->globalOffset = globalOffset ;
	newGrid->gridPoints = fftw_alloc_real(newGrid->allocScheme) ;
	return newGrid ;
}		/* -----  end of function makeNewGrid  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  initialiseGrid
 *    Arguments:  Grid * grid - Grid to initialise.
 *  Description:  Sets the top left quadrant of the grid to 1.
 * =====================================================================================
 */

void initialiseGrid(Grid * grid) {
	int localDimX = grid->localGridDimX ;
	int localDimY = grid->localGridDimY ;
	for (int x = 0 ; x < localDimX/4 ; ++x) {
		for (int y = 0 ; y < localDimY ; ++y) {
			if ((rank*localDimY+y)/(localDimX/4) == 3) {
				grid->gridPoints[y*localDimX+x] = 1 ;
			}
		}
	}
}		/* -----  end of function initialiseGrid  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printGrid
 *    Arguments:  Grid * grid - Grid to output to file.
 *                FILE * output - Output stream.
 *  Description:  Outputs grid to output.
 * =====================================================================================
 */

void printGrid(Grid * grid, FILE * output) {
	int localDimX = grid->localGridDimX ;
	int localDimY = grid->localGridDimY ;
	for (int y = localDimY-1 ; y >= 0 ; --y) {
		for (int x = 0 ; x < localDimX ; ++x) {
			fprintf(output, "%lf ", grid->gridPoints[y*localDimX+x]);
		}
		fprintf(output, "\n");
	}
}		/* -----  end of function printGrid  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printSplotData
 *    Arguments:  Grid * grid - Grid to output to file.
 *                FILE * output - Output stream.
 *  Description:  Outputs splot format of grid to output.
 * =====================================================================================
 */

void printSplotData(Grid * grid, FILE * output) {
	int localDimX = grid->localGridDimX ;
	int localDimY = grid->localGridDimY ;
	int globalOffset = grid->globalOffset ;
	for (int y = localDimY-1 ; y >= 0 ; --y) {
		for (int x = 0 ; x < localDimX ; ++x) {
			fprintf(output, "%d %d %lf\n", x, y+globalOffset, grid->gridPoints[y*localDimX+x]);
		}
	}
}		/* -----  end of function printSplotData  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  simulateGrid
 *    Arguments:  Grid * grid - Grid to simulate.
 *                bool mutliFlag - Turns on printing.
 *  Description:  Simulates grid described by diffusion equation using a parallel fft.
 *                Outputs the grid to a file in each timestep if printflag passed/
 * =====================================================================================
 */

void simulateGrid(Grid * grid) {
	int globalGridDimX = grid->globalGridDimX ;
	int globalGridDimY = grid->globalGridDimY ;
	int localDimX = grid->localGridDimX ;
	int localDimY = grid->localGridDimY ;
	int globalOffset = grid->globalOffset ;
	double * gridPointsTrTime0 = fftw_alloc_real(grid->allocScheme) ;
	double * gridPointsTrTimeT = fftw_alloc_real(grid->allocScheme) ;
	ptrdiff_t y,x,kx,ky ;

	fftw_plan planFor ;
	fftw_plan planBack ;

 	planFor =  fftw_mpi_plan_r2r_2d(globalGridDimY, globalGridDimX, grid->gridPoints, gridPointsTrTime0, MPI_COMM_WORLD, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE) ;
	planBack = fftw_mpi_plan_r2r_2d(globalGridDimY, globalGridDimX, gridPointsTrTimeT, grid->gridPoints, MPI_COMM_WORLD, FFTW_REDFT00 , FFTW_REDFT00, FFTW_ESTIMATE) ;

	// Transform phi_0
	fftw_execute(planFor) ;

	int numTimeSteps = 5E3 ;
	double deltaT = 0.002 ;
	double timeT = 0 ;
	int found = 0 ;
	FILE * out = fopen("valA.dat","w") ;
	for (int i = 0 ; i < numTimeSteps ; ++i) {
		timeT += deltaT ;
		// Advance transformed F(phi_0) to F(phi_t_i)
		for (kx = 0 ; kx < localDimX ; ++kx) {
			for (ky = 0 ; ky < localDimY ; ++ky) {
				ptrdiff_t globalky = ky+globalOffset ;
				gridPointsTrTimeT[ky*localDimX+kx] = gridPointsTrTime0[ky*localDimX+kx]*exp(-PI*PI*timeT*(kx*kx + globalky*globalky)) ;
			}
		}
		// Transform back from F(phi_t) to phi_t
		fftw_execute(planBack) ;
		// Multiply by normalisation constant 1/(2(N-1)) * (1/(2(N-1))). //
		for (x = 0 ; x < localDimX ; ++x) {
			for (y = 0 ; y < localDimY ; ++y) {
				grid->gridPoints[y*localDimX+x] /= (double)(4*(globalGridDimX-1)*(globalGridDimY-1)) ;
				if (globalOffset+y == globalGridDimY/4 && x == 3*(globalGridDimX/4)) {
					fprintf(out,  "%lf %lf\n", grid->gridPoints[y*localDimX+x], timeT) ;
					if (grid->gridPoints[y*localDimX+x] > 0.05 && found == 0) {
						found = 1 ;
						printf("%lf\n", timeT) ;
					}
				}
			}
		}
	}
	fclose(out) ;

	// Clean up resources. //
	fftw_destroy_plan(planFor);  
	fftw_destroy_plan(planBack);  
	fftw_free(gridPointsTrTime0) ;
	fftw_free(gridPointsTrTimeT) ;
}		/* -----  end of function simulateGrid  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  freeGrid
 *    Arguments:  Grid * grid
 *  Description:  Frees memory associated with grid.
 * =====================================================================================
 */

void freeGrid(Grid * grid) {
	fftw_free(grid->gridPoints) ;
	free(grid) ;
}		/* -----  end of function freeGrid  ----- */
