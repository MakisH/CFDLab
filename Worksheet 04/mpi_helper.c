#include "stdlib.h"
#include "stdio.h"
#include "mpi_helper.h"
#include "LBDefinitions.h"
#include "mpi.h"

void initializeMPI( int *rank, int *number_of_ranks, int argc, char *argv[] ) {
	// Initialize n processes
	MPI_Init( &argc, &argv );
	// Ask for the size of the communicator (number_of_ranks)
	MPI_Comm_size( MPI_COMM_WORLD, number_of_ranks );
	// Ask for the local process id (rank)
	MPI_Comm_rank( MPI_COMM_WORLD, rank );
	//printf("hi mpi! rank %d of %d.\n", *rank, *number_of_ranks);
}

void finalizeMPI() {
	printf("bye mpi!\n");
	// Write all the output
	fflush(stdout);
	fflush(stderr);
	
	// Synchronize all processes
	MPI_Barrier( MPI_COMM_WORLD );
	// Terminate the MPI session
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
	//int neighborSendId = MPI_PROC_NULL;	// default action - do nothing if on the boundary
	//int neighborRecvId = MPI_PROC_NULL;	// keep it NULL, otherwise if-s underneath are wrong!
	MPI_Status status;//, st1;
	//MPI_Request send_request;
	
	//switch (direction) {
	//	case DIRECTION_LR :// x+ direction (left-to-right)
	//		if(rank % iProc != iProc - 1 ) // NOT right boundary => communicate with right
	//			neighborSendId = rank + 1; // send & receive are the same = rank + 1
	//		if(rank % iProc != 0)
	//			neighborRecvId = rank - 1;
	//		break;

	//	case DIRECTION_RL :// x- direction (right-to-left)
	//		if(rank % iProc != 0 ) // NOT left boundary => communicate with left
	//			neighborSendId = rank - 1;
	//		if(rank % iProc != iProc - 1 )
	//			neighborRecvId = rank + 1;
	//		break;

	//	case DIRECTION_DT :// z+ direction (down-to-top)
	//		if(rank / iProc / jProc != kProc - 1 ) // NOT top boundary => exch with top
	//			neighborSendId = rank + iProc * jProc;
	//		if(rank / iProc / jProc != 0 ) // NOT bottom boundary => exch with bottom
	//			neighborRecvId = rank - iProc * jProc;
	//		break;

	//	case DIRECTION_TD :// z- direction (top-to-down)
	//		if(rank / iProc / jProc != 0 ) // NOT bottom boundary => exch with bottom
	//			neighborSendId = rank - iProc * jProc;
	//		if(rank / iProc / jProc != kProc - 1 ) // NOT top boundary => exch with top
	//			neighborRecvId = rank + iProc * jProc;
	//		break;

	//	case DIRECTION_FB :// y+ direction (front-to-back)
	//		if(rank / iProc % jProc != jProc - 1) // NOT back boundary => exch with back
	//			neighborSendId = rank + iProc;
	//		if(rank / iProc % jProc != 0) // NOT front boundary => exch with front
	//			neighborRecvId = rank - iProc;
	//		break;

	//	case DIRECTION_BF :// y- direction (back-to-front)
	//		if(rank / iProc % jProc != 0) // NOT front boundary => exch with front
	//			neighborSendId = rank - iProc;
	//		if(rank / iProc % jProc != jProc - 1) // NOT back boundary => exch with back
	//			neighborRecvId = rank + iProc;
	//		break;

	//	default :
	//		neighborSendId = MPI_PROC_NULL;
	//		neighborRecvId = MPI_PROC_NULL;
	//		printf("should never be here! Error in switch condition\n");
	//		return;
	//		break;
	//}
	//printf("MPI send & recv neighbors in direction %d: %d\n", direction, neighborId);

	// // version 4
		MPI_Sendrecv(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, readBuffer[direction + 1  - direction % 2 * 2], sizeBuffer[direction + 1  - direction % 2 * 2], MPI_DOUBLE, neighbor[direction + 1  - direction % 2 * 2], MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	// // version 3

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
		printf("swap rank %d direction %d\n",rank, direction);
}

void extraction(double *collideField, int *flagField, int *xlength, double **sendBuffer, int boundary, int rank) {
	
	int each[5]; // trick to implement a "foreach" loop, for every i direction to be transfered
	int x_start, x_end, y_start, y_end, z_start, z_end;
	int currentCell;
	
	int SizeX = xlength[0] + 2; // Size of the extended domain in each direction
	int SizeY = xlength[1] + 2;
	int SizeZ = xlength[2] + 2;
	int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain
	
	switch (boundary) {

		// x- direction (left)
		case DIRECTION_RL :
			x_start = 1;          x_end = 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = SizeZ - 1;
			
			each[0] = 1;
			each[1] = 5;
			each[2] = 8;
			each[3] = 11;
			each[4] = 15;
			//printf("RL Left\n");
			break;

		// x+ direction (right)
		case DIRECTION_LR :
			x_start = SizeX - 2;  x_end = SizeX - 2;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = SizeZ - 1;
			
			each[0] = 3;
			each[1] = 7;
			each[2] = 10;
			each[3] = 13;
			each[4] = 17;
			//printf("LR Right\n");
			break;

		// z+ direction (top)
		case DIRECTION_DT :
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = SizeZ - 2;  z_end = SizeZ - 2;
			
			each[0] = 14;
			each[1] = 15;
			each[2] = 16;
			each[3] = 17;
			each[4] = 18;
			//printf("Top\n");
			break;      
			
		case DIRECTION_TD:
			// z- direction (down)
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 1;          z_end = 1;
			
			each[0] = 0;
			each[1] = 1;
			each[2] = 2;
			each[3] = 3;
			each[4] = 4;
			//printf("Down\n");
			break;

		// y+ direction (front)
		case DIRECTION_BF :
			x_start = 0;          x_end = SizeX - 1;
			y_start = SizeY - 2;  y_end = SizeY - 2;
			z_start = 0;          z_end = SizeZ - 1;
			
			each[0] = 0;
			each[1] = 5;
			each[2] = 6;
			each[3] = 7;
			each[4] = 14;
			//printf("Front\n");
			break;      
			
		// y- direction (back)
		case DIRECTION_FB :
			x_start = 0;        x_end = SizeX - 1;
			y_start = 1;        y_end = 1;
			z_start = 0;        z_end = SizeZ - 1;
			
			each[0] = 4;
			each[1] = 11;
			each[2] = 12;
			each[3] = 13;
			each[4] = 18;
			//printf("Back\n");
			break;
		
		default :
			printf("You should not be here! Are you in 4D ???\n");
			x_start = 0;        x_end = -1;
			y_start = 0;        y_end = -1;
			z_start = 0;        z_end = -1;
			break;
			
	}
		//printf("extr rank %d direction %d\n",rank,boundary);

	// buffer cell
	int cell = -1;
	// loop over the whole 3D wall in corresponding direction
	for (int z = z_start; z <= z_end; ++z) {
		for (int y = y_start; y <= y_end; ++y) {
			for (int x = x_start; x <= x_end; ++x) {
				// buffer index - destination
				cell++;
				// cpuDomain cell index - source
				currentCell = x + y * SizeX + z * SizeXY;
				// loop over all velocities to be extracted
				for (int dir = 0; dir < 5; ++dir) {
					sendBuffer[boundary][5 * cell + dir] = collideField[Q_NUMBER * currentCell + each[dir]];
					//if(rank == 4)
					//printf("direction %d - extract %f from coll %d into buff %d\n", boundary, collideField[Q_NUMBER * currentCell + each[dir]], Q_NUMBER * currentCell + each[dir], 5 * cell + dir );
				}
			}
		}
	}
	
}

void injection(double *collideField, int *flagField, int *xlength, double **readBuffer, int boundary, int rank) {
	
	int each[5]; // trick to implement a "foreach" loop, for every i direction to be transfered
	int x_start, x_end, y_start, y_end, z_start, z_end;
	int currentCell;
	
	int SizeX = xlength[0] + 2; // Size of the extended domain in each direction
	int SizeY = xlength[1] + 2;
	int SizeZ = xlength[2] + 2;
	int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain
	
	switch (boundary) {

		// x- direction (left)
		case DIRECTION_RL :
			x_start = 0;          x_end = 0;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = SizeZ - 1;
			
			each[0] = 3;
			each[1] = 7;
			each[2] = 10;
			each[3] = 13;
			each[4] = 17;
			break;

		// x+ direction (right)
		case DIRECTION_LR :
			x_start = SizeX - 1;  x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = SizeZ - 1;
			each[0] = 1;
			each[1] = 5;
			each[2] = 8;
			each[3] = 11;
			each[4] = 15;

			break;

		// z+ direction (up)
		case DIRECTION_DT :
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = SizeZ - 1;  z_end = SizeZ - 1;
			
			each[0] = 0;
			each[1] = 1;
			each[2] = 2;
			each[3] = 3;
			each[4] = 4;
			break;      
			
		case DIRECTION_TD:
			// z- direction (down) (in the case of injection the limits are opposite)
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = 0;

			each[0] = 14;
			each[1] = 15;
			each[2] = 16;
			each[3] = 17;
			each[4] = 18;
			break;

		// y+ direction (front)
		case DIRECTION_BF :
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
		case DIRECTION_FB :
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
			printf("Error, non-existiog direction!\n");
			x_start = 0;        x_end = -1;
			y_start = 0;        y_end = -1;
			z_start = 0;        z_end = -1;
			break;
	}
	//printf("inje rank %d direction %d\n",rank,boundary);
	int cell = -1;
	
	for (int z = z_start; z <= z_end; ++z) {
		for (int y = y_start; y <= y_end; ++y) {
			for (int x = x_start; x <= x_end; ++x) {
				
				cell++;
				currentCell = x + y * SizeX + z * SizeXY;
				
				for (int dir = 0; dir < 5; ++dir) {
					collideField[Q_NUMBER * currentCell + each[dir]] = readBuffer[boundary][5 * cell + dir];
					//if(rank == 4)
					//printf("rank %d, direction %d - inject %f from buff %d into coll %d\n", rank, boundary, readBuffer[boundary][5 * cell + dir], 5 * cell + dir, Q_NUMBER * currentCell + each[dir]);
				}
				
			}
		}
	}
	
}
