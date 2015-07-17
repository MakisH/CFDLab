
//#include "stdlib.h" // not needed?

#include "mpi_helper.h"
#include <unistd.h>

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

void swap(double * const * const sendBuffer, double * const *const readBuffer, const int rank) {
	// MPI_Status status; // waiting for status makes it slower ?
	MPI_Request * send_request = (MPI_Request *) malloc(neighbours_count[rank] * sizeof(MPI_Request));



	//version 5
	// send asynchronously to everyone
	for(int i = 0; i < neighbours_count[rank]; ++i){
		MPI_Isend(sendBuffer[i], neighbours_local_buffer_size[rank][i], MPI_DOUBLE, neighbours_procid[rank][i], neighbours_tag[rank][i], MPI_COMM_WORLD, &send_request[i]);

	}
	// receive from all neighbours asynchronously
	for(int i = 0; i < neighbours_count[rank]; ++i){
		//int inv_dir = neighbours_dir[rank][i]+ 1  - neighbours_dir[rank][i] % 2 * 2;
		MPI_Recv(readBuffer[i], neighbours_local_buffer_size[rank][i], MPI_DOUBLE, neighbours_procid[rank][i], neighbours_tag[rank][i], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("I am rank %d. I receive from cpu # %d, tag %d, value: %f\n", rank, neighbours_procid[rank][i], neighbours_tag[rank][i], *(readBuffer[i]+ 2 + neighbours_local_buffer_size[rank][i]/3));
		//printf("rank %d, I send: %f, to: %d -> tag: %d\n\n", rank, *(sendBuffer[i] + 2 + neighbours_local_buffer_size[rank][i]/3), neighbours_procid[rank][i], neighbours_tag[rank][i]);
		MPI_Wait(&send_request[i],MPI_STATUS_IGNORE);
	}
	free(send_request);
	// // version 4
	//MPI_Sendrecv(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, readBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[inv_dir], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/// old
	//// version 3

	//MPI_Request send_request ;
	//if(neighbor[direction] != MPI_PROC_NULL)	MPI_Isend(sendBuffer[direction], sizeBuffer[direction], MPI_DOUBLE, neighbor[direction], 0, MPI_COMM_WORLD, &send_request);
	//	// opposite direction is : direction + 1  - direction % 2 * 2
	//int inv_dir = direction + 1  - direction % 2 * 2;
	//if(neighbor[inv_dir] != MPI_PROC_NULL){
	//	MPI_Recv(readBuffer[inv_dir], sizeBuffer[inv_dir], MPI_DOUBLE, neighbor[inv_dir], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//}
	//if(neighbor[direction] != MPI_PROC_NULL)	MPI_Wait(&send_request,MPI_STATUS_IGNORE);


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
void extraction(double * const collideField, const int * const cpuDomain, double * const sendBuffer, 	const int direction,	const int x_start,	const int y_start,	const int z_start,	const int x_end, const int y_end,	const int z_end) {
	int currentCell;
	int cell = -1; // buffer cell
	for (int z = z_start; z <= z_end; ++z) {
		for (int y = y_start; y <= y_end; ++y) {
			for (int x = x_start; x <= x_end; ++x) {
				++cell; // buffer index - destination
				currentCell = x + y * (cpuDomain[0] + 2) + z * (cpuDomain[0] + 2) * (cpuDomain[1] + 2); // cpuDomain cell index - source
				for (int dir = 0; dir < 5; ++dir) { // loop over all velocities to be extracted
					sendBuffer[5 * cell + dir] = collideField[Q_NUMBER * currentCell + each[direction][dir]];

					//printf("%f ", sendBuffer[5*cell+dir]);
				}
			}
		}
	}
	//sleep(1000);
}

// before calling - check for NULL neighbor!
// DIR_R means that we inject from Left to Right  on the Left side(x-)(3,7,10,13,17)
void injection(double * const collideField, const int * const cpuDomain, double * const readBuffer, const int direction,	const int x_start,	const int y_start,	const int z_start,	const int x_end, const int y_end,	const int z_end ) {
	int currentCell;
	int cell = -1;
	int inv_dir = direction + 1  - direction % 2 * 2;
	for (int z = z_start; z <= z_end; ++z) {
		for (int y = y_start; y <= y_end; ++y) {
			for (int x = x_start; x <= x_end; ++x) {
				++cell;
				currentCell = x + y * (cpuDomain[0] + 2) + z * (cpuDomain[0] + 2) * (cpuDomain[1] + 2);
			//	printf ("I write to cell # %d\n", currentCell);
				for (int vel_dir = 0; vel_dir < 5; ++vel_dir) {
					collideField[Q_NUMBER * currentCell + each[inv_dir][vel_dir]] = readBuffer[5 * cell + vel_dir];
				}
			}
		}
	}
//	sleep(1000);
}
