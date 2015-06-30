#include "stdlib.h"
#include "stdio.h"
#include "mpi_helper.h"
#include "LBDefinitions.h"
#include "mpi.h"

void initializeMPI( int *rank, int *number_of_ranks, int argc, char *argv[] ) {
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
void swap(double **sendBuffer, double **readBuffer, int *sizeBuffer, int direction, int iProc, int kProc, int jProc, int rank, int *neighbor) {
	//printf("rank %d\n",rank);
	MPI_Status status;
	// // version 4
		MPI_Sendrecv(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction + 1  - direction % 2 * 2], MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		/// old
	//// version 3
		//MPI_Isend(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighborSendId, 0, MPI_COMM_WORLD, &send_request);
		//// opposite direction is : direction + 1  - direction % 2 * 2
		//MPI_Recv(readBuffer[direction + 1  - direction % 2 * 2], sizeBuffer[direction + 1  - direction % 2 * 2], MPI_DOUBLE, neighborRecvId, 0, MPI_COMM_WORLD, &status);
		//MPI_Wait(&send_request,MPI_STATUS_IGNORE);
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
void extraction(double *collideField, int *xlength/* == cpuDomain */, double **sendBuffer, int boundary, const side *Bsides) {
	int currentCell;
	int cell = -1; // buffer cell
	for (int z = Bsides[boundary].z_start; z <= Bsides[boundary].z_end; ++z) {
		for (int y = Bsides[boundary].y_start; y <= Bsides[boundary].y_end; ++y) {
			for (int x = Bsides[boundary].x_start; x <= Bsides[boundary].x_end; ++x) {
				++cell; // buffer index - destination
				currentCell = x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2); // cpuDomain cell index - source
				for (int dir = 0; dir < 5; ++dir) { // loop over all velocities to be extracted
					sendBuffer[boundary][5 * cell + dir] = collideField[Q_NUMBER * currentCell + each[boundary][dir]];
				}
			}
		}
	}
}

// before calling - check for NULL neighbor!
void injection(double *collideField, int *xlength, double **readBuffer, int boundary, const side *Bsides) {
	int currentCell;
	int cell = -1;
	for (int z = Bsides[boundary].z_start; z <= Bsides[boundary].z_end; ++z) {
		for (int y = Bsides[boundary].y_start; y <= Bsides[boundary].y_end; ++y) {
			for (int x = Bsides[boundary].x_start; x <= Bsides[boundary].x_end; ++x) {
				++cell;
				currentCell = x + y * (xlength[0] + 2) + z * (xlength[0] + 2) * (xlength[1] + 2);
				for (int dir = 0; dir < 5; ++dir) {
					collideField[Q_NUMBER * currentCell + each[boundary][dir]] = readBuffer[boundary][5 * cell + dir];
				}
			}
		}
	}
}
