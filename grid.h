#ifndef GRID_H_6JLKZ8VN
#define GRID_H_6JLKZ8VN

/*
 * =====================================================================================
 *
 *       Filename:  grid.h
 *
 *    Description:  Header file for grid object. The grid simulates a diffusion PDE 
 *                  using VN boundary conditions.
 *
 *        Version:  1.0
 *        Created:  04/14/2016 11:20:25 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Michael Tierney (MT), tiernemi@tcd.ie
 *
 * =====================================================================================
 */

#include "stdio.h"
#include "fftw3-mpi.h"

/* 
 * ===  STRUCT  ========================================================================
 *         Name:  Grid
 *       Fields:  
 *  Description:  Grid object used to store values of Phi which obey a diffusion equation.
 *                The system is propogated in time using a FFT based solver.
 * =====================================================================================
 */

typedef struct Grid {
	double * gridPoints ;
	ptrdiff_t localGridDimY ;
	ptrdiff_t  localGridDimX ;
	ptrdiff_t globalGridDimY ;
	ptrdiff_t globalGridDimX ;
	ptrdiff_t allocScheme ;
	ptrdiff_t globalOffset ;
} Grid ;		/* -----  end of struct Grid  ----- */

Grid * makeNewGrid(int,int) ; 
void initialiseGrid(Grid *) ;
void printGrid(Grid *, FILE *) ;
void simulateGrid(Grid *,int, double, int) ;
void freeGrid(Grid *) ;

#endif /* end of include guard: GRID_H_6JLKZ8VN */
