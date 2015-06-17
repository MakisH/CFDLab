#include "stdlib.h"
#include "stdio.h"
#include "mpi_helper.h"
#include "LBDefinitions.h"

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



void extraction(double *collideField, int *flagField, int *xlength, double *sendBuffer[], int boundary) {
  
  int each[5]; // trick to implement a "foreach" loop, for every i direction to be transfered
  int x_start, x_end, y_start, y_end, z_start, z_end;
  int currentCell;
  
  int SizeX = xlength[0] + 2; // Size of the extended domain in each direction
  int SizeY = xlength[1] + 2;
  int SizeZ = xlength[2] + 2;
  int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain
  
  switch (boundary) {

    // x- direction (left)
    case 0 :
      x_start = 0;          x_end = 0;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 1;
      each[1] = 5;
      each[2] = 8;
      each[3] = 11;
      each[4] = 15;
      break;

    // x+ direction (right)
    case 1 :
      x_start = SizeX - 1;  x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 3;
      each[1] = 7;
      each[2] = 10;
      each[3] = 13;
      each[4] = 17;
      break;

    // z+ direction (top)
    case 2 :
      x_start = 0;          x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = SizeZ - 1;  z_end = SizeZ - 1;
      
      each[0] = 14;
      each[1] = 15;
      each[2] = 16;
      each[3] = 17;
      each[4] = 18;
      break;      
      
    case 3:
      // z- direction (down)
      x_start = 0;          x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = 0;
      
      each[0] = 0;
      each[1] = 1;
      each[2] = 2;
      each[3] = 3;
      each[4] = 4;
      break;

    // y+ direction (front)
    case 4 :
      x_start = 0;          x_end = SizeX - 1;
      y_start = SizeY - 1;  y_end = SizeY - 1;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 4;
      each[1] = 11;
      each[2] = 12;
      each[3] = 13;
      each[4] = 18;
      break;      
      
    // y- direction (back)
    case 5 :
      x_start = 0;        x_end = SizeX - 1;
      y_start = 0;        y_end = 0;
      z_start = 0;        z_end = SizeZ - 1;
      
      each[0] = 0;
      each[1] = 5;
      each[2] = 6;
      each[3] = 7;
      each[4] = 14;
      break;
    
    default :
      x_start = 0;        x_end = -1;
      y_start = 0;        y_end = -1;
      z_start = 0;        z_end = -1;
      
  }
  
  int cell = 0;
  
  for (int x = x_start; x <= x_end; ++x) {
    for (int y = y_start; y <= y_end; ++y) {
      for (int z = z_start; z <= z_end; ++z) {
        
        cell++;
        currentCell = x + y*SizeX + z*SizeXY;
        
        for (int dir = 0; dir < 5; ++dir) {
          sendBuffer[boundary][5*cell + dir] = collideField[Q_NUMBER*currentCell + each[dir]];
        }
        
      }
    }
  }
  
}

void injection(double *collideField, int *flagField, int *xlength, double *readBuffer[], int boundary) {
  
  int each[5]; // trick to implement a "foreach" loop, for every i direction to be transfered
  int x_start, x_end, y_start, y_end, z_start, z_end;
  int currentCell;
  
  int SizeX = xlength[0] + 2; // Size of the extended domain in each direction
  int SizeY = xlength[1] + 2;
  int SizeZ = xlength[2] + 2;
  int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain
  
  switch (boundary) {

    // x- direction (left)
    case 0 :
      x_start = SizeX - 1;  x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 1;
      each[1] = 5;
      each[2] = 8;
      each[3] = 11;
      each[4] = 15;
      break;

    // x+ direction (right)
    case 1 :
      x_start = 0;          x_end = 0;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 3;
      each[1] = 7;
      each[2] = 10;
      each[3] = 13;
      each[4] = 17;
      break;

    // z+ direction (up)
    case 2 :
      x_start = 0;          x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = 0;          z_end = 0;
      
      each[0] = 14;
      each[1] = 15;
      each[2] = 16;
      each[3] = 17;
      each[4] = 18;
      break;      
      
    case 3:
      // z- direction (down) (in the case of injection the limits are opposite)
      x_start = 0;          x_end = SizeX - 1;
      y_start = 0;          y_end = SizeY - 1;
      z_start = SizeZ - 1;  z_end = SizeZ - 1;
      
      each[0] = 0;
      each[1] = 1;
      each[2] = 2;
      each[3] = 3;
      each[4] = 4;
      break;

    // y+ direction (front)
    case 4 :
      x_start = 0;          x_end = SizeX - 1;
      y_start = 0;          y_end = 0;
      z_start = 0;          z_end = SizeZ - 1;
      
      each[0] = 4;
      each[1] = 11;
      each[2] = 12;
      each[3] = 13;
      each[4] = 18;
      break;

    // y- direction (back)
    case 5 :
      x_start = 0;         x_end = SizeX - 1;
      y_start = SizeY - 1; y_end = SizeY - 1;
      z_start = 0;         z_end = SizeZ - 1;
      
      each[0] = 0;
      each[1] = 5;
      each[2] = 6;
      each[3] = 7;
      each[4] = 14;
      break;

    default :
      x_start = 0;        x_end = -1;
      y_start = 0;        y_end = -1;
      z_start = 0;        z_end = -1;
      
  }
  
  int cell = 0;
  
  for (int x = x_start; x <= x_end; ++x) {
    for (int y = y_start; y <= y_end; ++y) {
      for (int z = z_start; z <= z_end; ++z) {
        
        cell++;
        currentCell = x + y*SizeX + z*SizeXY;
        
        for (int dir = 0; dir < 5; ++dir) {
          collideField[Q_NUMBER*currentCell + each[dir]] = readBuffer[boundary][5*cell + dir];
        }
        
      }
    }
  }
  
}