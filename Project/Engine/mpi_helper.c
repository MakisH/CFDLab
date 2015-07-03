
//#include "stdlib.h" // not needed?

#include "mpi_helper.h"

// stdout
// stderr
#include <stdio.h>

// Q_NUMBER
#include "LBDefinitions.h"
#include "mpi.h"

void initializeMPI( int *rank, int *number_of_ranks, int argc, char ** argv ) {
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
	MPI_Comm_rank( MPI_COMM_WORLD, rank );
}

void finalizeMPI() {
	printf("bye mpi!\n");
	// Write all the output
	fflush(stdout);
	fflush(stderr);
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Finalize();
}

//void swap(double **sendBuffer, double **readBuffer, int *sizeBuffer, int direction, int boundary, int iProc, int kProc, int jProc, int rank) {
//	
//	int neighbor_distance = 0;
//	int neighborId_send = MPI_PROC_NULL;
//	int neighborId_recv = MPI_PROC_NULL;
//	MPI_Status status;
//	
//	switch (direction) {
//		// x- direction (right-to-left)
//		case DIRECTION_RL :
//			neighbor_distance = -1;
//			if ( rank % iProc == 0 ) { 
//				// left boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;
//			}
//			if ( rank % iProc == iProc - 1 ) { 
//				// right boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		// x+ direction (left-to-right)
//		case DIRECTION_LR :
//			neighbor_distance = 1;
//			if ( rank % iProc == 0 ) { 
//				// left boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;
//			} else if ( rank % iProc == iProc - 1 ) { 
//				// right boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		// z+ direction (down-to-top)
//		case DIRECTION_DT :
//			neighbor_distance = iProc;
//			if ( rank % (iProc*kProc) < iProc ) { 
//				// bottom boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;
//			} else if ( rank % (iProc*kProc) >= iProc*(kProc - 1) ) { 
//				// top boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		// z- direction (top-to-down)
//		case DIRECTION_TD :
//			neighbor_distance = -iProc;
//			if ( rank % (iProc*kProc) < iProc ) { 
//				// bottom boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;
//			} else if ( rank % (iProc*kProc) >= iProc*(kProc - 1) ) { 
//				// top boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		// y+ direction (back-to-front)
//		case DIRECTION_BF :
//			neighbor_distance = -iProc * kProc;
//			if ( rank >= iProc*(jProc-1)*kProc ) { 
//				// back boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;
//			} else if ( rank <= iProc*kProc - 1 ) { 
//				// front boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		// y- direction (front-to-back)
//		case DIRECTION_FB :
//			neighbor_distance = iProc * kProc;
//			if ( rank >= iProc*(jProc-1)*kProc ) { 
//				// back boundary
//				neighborId_send = MPI_PROC_NULL;
//				neighborId_recv = rank - neighbor_distance;
//			} else if ( rank <= iProc*kProc - 1 ) { 
//				// front boundary
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = MPI_PROC_NULL;        
//			} else { 
//				// inner subdomain
//				neighborId_send = rank + neighbor_distance;
//				neighborId_recv = rank - neighbor_distance;        
//			}
//			break;
//			
//		default :
//			neighborId_send = MPI_PROC_NULL;
//			neighborId_recv = MPI_PROC_NULL;      
//			printf("should never be here! Error in switch condition\n");
//			return;
//			break;
//	}
//		printf("MPI send & recv neighbors: %d %d\n",neighborId_send, neighborId_recv );
//	MPI_Send(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId_send, 1, MPI_COMM_WORLD);
//	MPI_Recv(readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId_recv, 1, MPI_COMM_WORLD, &status);
//	
//}
/// old swap
void swap(double * const * const sendBuffer, double * const * const readBuffer, const int * const sizeBuffer, const int direction, const int * const neighbor) {
	// MPI_Status status; // waiting for status makes it slower ?
	MPI_Request send_request;
	// // version 4
	//MPI_Sendrecv(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[inv_dir], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/// old
	//// version 3
	if(neighbor[direction] != MPI_PROC_NULL)	MPI_Isend(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, MPI_COMM_WORLD, &send_request);
		// opposite direction is : direction + 1  - direction % 2 * 2
	int inv_dir = direction + 1  - direction % 2 * 2;
	if(neighbor[inv_dir] != MPI_PROC_NULL){
		MPI_Recv(readBuffer[inv_dir], sizeBuffer[inv_dir], MPI_DOUBLE, neighbor[inv_dir], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if(neighbor[direction] != MPI_PROC_NULL)	MPI_Wait(&send_request,MPI_STATUS_IGNORE);
	// // version 2
		//MPI_Isend(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId, 0, MPI_COMM_WORLD, &send_request);
		//MPI_Recv(readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId, 0, MPI_COMM_WORLD, &status);
		//MPI_Wait(&send_request,MPI_STATUS_IGNORE);
	// // version 1
	//MPI_Sendrecv(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId, rank, readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborId, neighborId, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//for(int i = 0; i < sizeBuffer[direction]*5;++i){
	//	printf("%f ",readBuffer[direction][i]);
	//}
	//	printf("swap rank %d direction %d\n",rank, direction);
}

// before calling - check for NULL neighbor!
// DIR_R means that Right side of domain is streamed to the Right(3,7,10,13,17)
void extraction(double * const collideField, const int * const cpuDomain, double * const * const sendBuffer, const int direction, const side * const Bsides) {
	int currentCell;
	int cell = -1; // buffer cell
	for (int z = Bsides[direction].z_start; z <= Bsides[direction].z_end; ++z) {
		for (int y = Bsides[direction].y_start; y <= Bsides[direction].y_end; ++y) {
			for (int x = Bsides[direction].x_start; x <= Bsides[direction].x_end; ++x) {
				++cell; // buffer index - destination
				currentCell = x + y * (cpuDomain[0] + 2) + z * (cpuDomain[0] + 2) * (cpuDomain[1] + 2); // cpuDomain cell index - source
				for (int dir = 0; dir < 5; ++dir) { // loop over all velocities to be extracted
					sendBuffer[direction][5 * cell + dir] = collideField[Q_NUMBER * currentCell + each[direction][dir]];
				}
			}
		}
	}
}

// before calling - check for NULL neighbor!
// DIR_R means that we inject from Left to Right  on the Left side(x-)(3,7,10,13,17) 
void injection(double * const collideField, const int * const cpuDomain, double * const * const readBuffer, const int direction, const side * const Bsides) {
	int currentCell;
	int cell = -1;
	int inv_dir = direction + 1  - direction % 2 * 2;
	for (int z = Bsides[direction].z_start; z <= Bsides[direction].z_end; ++z) {
		for (int y = Bsides[direction].y_start; y <= Bsides[direction].y_end; ++y) {
			for (int x = Bsides[direction].x_start; x <= Bsides[direction].x_end; ++x) {
				++cell;
				currentCell = x + y * (cpuDomain[0] + 2) + z * (cpuDomain[0] + 2) * (cpuDomain[1] + 2);
				for (int vel_dir = 0; vel_dir < 5; ++vel_dir) {
					collideField[Q_NUMBER * currentCell + each[direction][vel_dir]] = readBuffer[inv_dir][5 * cell + vel_dir];
				}
			}
		}
	}
}
