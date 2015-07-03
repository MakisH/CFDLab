#ifndef _MPIHELPER_H_
#define _MPIHELPER_H_

#include "mpi.h"
#include "LBDefinitions.h"

// Initializes the MPI session
void initializeMPI( 
	int *rank,              /* Current process id */
	int *number_of_ranks,   /* Total number of processes */
	int argc,               /* Number of arguments */
	char ** argv            /* Array of arguments */
);

// Terminates the MPI session
void finalizeMPI();

// Transfer the overlapping cells
void swap(
	double * const * const sendBuffer,
	double * const * const readBuffer,
	const int * const sizeBuffer,
	const int direction,
	const int * const neighbor
);

typedef struct  {
	int x_start;
	int x_end;
	int y_start;
	int y_end;
	int z_start;
	int z_end;
}side ;

// Do the extraction
void extraction(
	double * const collideField,
	const int * const cpuDomain,
	double * const * const sendBuffer,
	const int direction,
	const side * const Bsides
);

// Do the injection
void injection(
	double * const collideField,
	const int * const cpuDomain,
	double * const * const readBuffer,
	const int direction,
	const side * const Bsides
);

#endif
