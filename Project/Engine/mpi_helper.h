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
	double * const sendBuffer,
	double * const readBuffer,
	const int rank
);

// Do the extraction
void extraction(
	double * const collideField,
	const int * const cpuDomain,
	double * const sendBuffer,
	const int direction,
	const int x_start,
	const int y_start,
	const int z_start,
	const int x_end,
	const int y_end,
	const int z_end
);

// Do the injection
void injection(
	double * const collideField,
	const int * const cpuDomain,
	double * const readBuffer,
	const int direction,
	const int x_start,
	const int y_start,
	const int z_start,
	const int x_end,
	const int y_end,
	const int z_end
);

#endif
