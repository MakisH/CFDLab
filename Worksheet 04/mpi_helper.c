#include "stdlib.h"
#include "stdio.h"
#include "mpi_helper.h"
#include "LBDefinitions.h"

void initializeMPI( int *rank, int *number_of_ranks, int argc, char *argv[] ) {
    printf("hi from mpi\n!");

  // Initialize n processes
  MPI_Init( &argc, &argv );
  printf("test init\n!");
  // Ask for the size of the communicator (number_of_ranks)
  MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
  printf("mpi_comm_size\n!");
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

void swap(double **sendBuffer, double **readBuffer, int *sizeBuffer, int direction, int boundary, int iProc, int kProc, int jProc, int rank) {
  
  int neighbor_distance = 0;
  int neighborId_send = MPI_PROC_NULL;
  int neighborId_recv = MPI_PROC_NULL;
  MPI_Status status;
  
  switch (direction) {
    // x- direction (right-to-left)
    case DIRECTION_RL :
      neighbor_distance = -1;
      if ( rank % iProc == 0 ) { 
        // left boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;
      } else if ( rank % iProc == iProc - 1 ) { 
        // right boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    // x+ direction (left-to-right)
    case DIRECTION_LR :
      neighbor_distance = 1;
      if ( rank % iProc == 0 ) { 
        // left boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;
      } else if ( rank % iProc == iProc - 1 ) { 
        // right boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    // z+ direction (down-to-top)
    case DIRECTION_DT :
      neighbor_distance = iProc;
      if ( rank % (iProc*kProc) < iProc ) { 
        // bottom boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;
      } else if ( rank % (iProc*kProc) >= iProc*(kProc - 1) ) { 
        // top boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    // z- direction (top-to-down)
    case DIRECTION_TD :
      neighbor_distance = -iProc;
      if ( rank % (iProc*kProc) < iProc ) { 
        // bottom boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;
      } else if ( rank % (iProc*kProc) >= iProc*(kProc - 1) ) { 
        // top boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    // y+ direction (back-to-front)
    case DIRECTION_BF :
      neighbor_distance = -iProc * kProc;
      if ( rank >= iProc*(jProc-1)*kProc ) { 
        // back boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;
      } else if ( rank <= iProc*kProc - 1 ) { 
        // front boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    // y- direction (front-to-back)
    case DIRECTION_FB :
      neighbor_distance = iProc * kProc;
      if ( rank >= iProc*(jProc-1)*kProc ) { 
        // back boundary
        neighborId_send = MPI_PROC_NULL;
        neighborId_recv = rank - neighbor_distance;
      } else if ( rank <= iProc*kProc - 1 ) { 
        // front boundary
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = MPI_PROC_NULL;        
      } else { 
        // inner subdomain
        neighborId_send = rank + neighbor_distance;
        neighborId_recv = rank - neighbor_distance;        
      }
      break;
      
    default :
      neighborId_send = MPI_PROC_NULL;
      neighborId_recv = MPI_PROC_NULL;      
			printf("should never be here! Error in switch condition\n");
			return;
      break;
  }
  
  MPI_Send(sendBuffer[boundary], sizeBuffer[boundary], MPI_DOUBLE, neighborId_send, 1, MPI_COMM_WORLD);
  MPI_Recv(readBuffer[boundary], sizeBuffer[boundary], MPI_DOUBLE, neighborId_recv, 1, MPI_COMM_WORLD, &status);
  
}

void extraction(double *collideField, int *flagField, int *xlength, double **sendBuffer, int boundary) {
  
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
      x_start = 1;          x_end = 1;
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
      x_start = SizeX - 2;  x_end = SizeX - 2;
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
      z_start = SizeZ - 2;  z_end = SizeZ - 2;
      
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
      z_start = 1;          z_end = 1;
      
      each[0] = 0;
      each[1] = 1;
      each[2] = 2;
      each[3] = 3;
      each[4] = 4;
      break;

    // y+ direction (front)
    case 4 :
      x_start = 0;          x_end = SizeX - 1;
      y_start = SizeY - 2;  y_end = SizeY - 2;
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
      y_start = 1;        y_end = 1;
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
      break;
      
  }
  
  int cell = -1;
  
  for (int z = z_start; z <= z_end; ++z) {
    for (int y = y_start; y <= y_end; ++y) {
      for (int x = x_start; x <= x_end; ++x) {
        
        cell++;
        currentCell = x + y*SizeX + z*SizeXY;
        
        for (int dir = 0; dir < 5; ++dir) {
          sendBuffer[boundary][5*cell + dir] = collideField[Q_NUMBER*currentCell + each[dir]];
        }
        
      }
    }
  }
  
}

void injection(double *collideField, int *flagField, int *xlength, double **readBuffer, int boundary) {
  
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
      break;
      
  }
  
  int cell = -1;
  
  for (int z = z_start; z <= z_end; ++z) {
    for (int y = y_start; y <= y_end; ++y) {
      for (int x = x_start; x <= x_end; ++x) {
        
        cell++;
        currentCell = x + y*SizeX + z*SizeXY;
        
        for (int dir = 0; dir < 5; ++dir) {
          collideField[Q_NUMBER*currentCell + each[dir]] = readBuffer[boundary][5*cell + dir];
        }
        
      }
    }
  }
  
}
