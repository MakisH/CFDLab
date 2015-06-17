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

// Do the extraction
void extraction(
      double *collideField, 
      int *flagField, 
      int *xlength, 
      double *sendBuffer[], 
      int boundary
);

// Do the injection
void injection(
      double *collideField, 
      int *flagField, 
      int *xlength, 
      double *readBuffer[], 
      int boundary
);

#endif