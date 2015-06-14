#include "stdlib.h"
#include "stdio.h"
#include "mpi_helper.h"

void initializeMPI( int *rank, int *number_of_ranks, int argc, char *argv[] ) {
  
  // Initialize n processes
  MPI_Init( &argc, &argv );
  
  // Ask for the size of the communicator (number_of_ranks)
  MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
  
  // Ask for the local process id (rank)
  MPI_Comm_rank( MPI_COMM_WORLD, rank );
  
  printf("Hello! I am rank %d of %d.\n", *rank, *number_of_ranks);
  
}

void finalizeMPI() {
  
  // Write all the output
  fflush(stdout);
  fflush(stderr);
  
  // Synchronize all processes
  MPI_Barrier( MPI_COMM_WORLD );
  
  // Terminate the MPI session
  MPI_Finalize();
}
