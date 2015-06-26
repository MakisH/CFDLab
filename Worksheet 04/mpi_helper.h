#ifndef _MPIHELPER_H_
#define _MPIHELPER_H_

#include "mpi.h"

// Initializes the MPI session
void initializeMPI( 
      int *rank,              /* Current process id */
      int *number_of_ranks,   /* Total number of processes */
      int argc,               /* Number of arguments */
      char *argv[]            /* Array of arguments */
);

// Terminates the MPI session
void finalizeMPI();

// Transfer the overlapping cells
void swap(
      double **sendBuffer, 
      double **readBuffer, 
      int *sizeBuffer, 
      int direction,
      int iProc, 
      int kProc, 
      int jProc, 
      int rank,
			int *neighbor
);

// Do the extraction
void extraction(
      double *collideField, 
      int *flagField, 
      int *xlength, 
      double **sendBuffer, 
      int boundary,
			int rank
);

// Do the injection
void injection(
      double *collideField, 
      int *flagField, 
      int *xlength, 
      double **readBuffer, 
      int boundary,
			int rank
);

#endif
