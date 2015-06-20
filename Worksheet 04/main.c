#ifndef _MAIN_C_
#define _MAIN_C_

#include "mpi_helper.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){
	int rank = 5;
	int number_of_ranks;

		// Start MPI
	initializeMPI( &rank, &number_of_ranks, argc, argv );

	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;

	int xlength[3];
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	int iProc, jProc, kProc;

	// send and read buffers for all possible directions :
	// [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
	double *sendBuffer[6];
	double *readBuffer[6];
	int sizeBuffer[6];
	int cpuDomain[3];
	int cpuDomain_size;

	if(0 == rank){
		// Read the config file using only one thread
		if(readParameters( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv ) == 1) return 1; // reading parameters from the file.

		// check if iProc and xlength are legitimate
		if(0 == iProc || 0 == jProc || 0 == kProc){
			printf("Invalid number of processors in some dimension(0)!\n");
			return 1;
		}
		else if(xlength[0] % iProc || xlength[1] % jProc || xlength[2] % kProc){
			printf("Non-integer ratio xlength/ijkProc - Invalid !\n");
			return 1;
		}
		cpuDomain[0] = xlength[0] / iProc;
		cpuDomain[1] = xlength[1] / jProc;
		cpuDomain[2] = xlength[2] / kProc;

		// Allocating the main three arrays.
		cpuDomain_size = (cpuDomain[0] + 2) * (cpuDomain[1] + 2) * (cpuDomain[2] + 2);
	}
	printf("before bcas\n");
	MPI_Bcast( xlength, 3, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &iProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &jProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &kProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( cpuDomain, 3, MPI_INT, 0, MPI_COMM_WORLD );
	printf("before bcas\n");
	MPI_Bcast( &cpuDomain_size, 1, MPI_INT, 0, MPI_COMM_WORLD );

		collideField = (double *) malloc(Q_NUMBER * cpuDomain_size * sizeof(double));
		streamField = (double *) malloc(Q_NUMBER * cpuDomain_size * sizeof(double));
		flagField = (int *) malloc(cpuDomain_size * sizeof(int));


		printf("before Init fields\n");
		printf("values:\n xlength: %d %d %d\n Proc ijk %d %d %d\n cpuDomain: %d %d %d\n ",xlength[0],xlength[1],xlength[2],iProc,jProc,kProc,cpuDomain[0],cpuDomain[1],cpuDomain[2]);
	// Init the main three arrays.
	initialiseFields( collideField, streamField, flagField, cpuDomain, iProc, jProc, kProc, rank);

	//// check for insufficient number of processors
	//if(number_of_ranks < iProc * jProc * kProc) {
	//	printf("There are not enough processors for this simulation, at least %d needed!", iProc * jProc * kProc);
	//	return 1;
	//}

	// Each processor responsible for its own domain.

	// allocate the buffers
	initialiseBuffers(sendBuffer, readBuffer, cpuDomain, sizeBuffer);
	printf("before t loop\n");
	for(int t = 0; t <= timesteps; t++){
		double *tmp = NULL;

		// TODO: maybe move all these to a separate function?
		// Do extraction, swap, injection for x+ (left to right)
		printf("before extr\n");
		if ( rank % iProc != iProc - 1 ) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_LR );
		printf("before swap\n");
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_LR, iProc, kProc, jProc, rank);
		printf("before injection\n");
		if ( rank % iProc != iProc - 1 ) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_LR );

		printf("before extr2\n");
		// Do extraction, swap, injection for x- (right to left)
		if ( rank % iProc != 0 ) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_RL );
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_RL, iProc, kProc, jProc, rank);
		if ( rank % iProc != 0 ) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_RL );
		printf("before extr3\n");

 //   // Do extraction, swap, injection for y+ (back to forth)
		if ( rank /(iProc * jProc ) != 0) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_BF);
		printf("before sw3\n");
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_BF, iProc, kProc, jProc, rank);
		printf("before inj3\n");
		if ( rank /(iProc * jProc ) != 0) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_BF );
		printf("middle streaming transport\n");
	 // Do extraction, swap, injection for y- (forth to back)
		if ( rank /(iProc * jProc ) != kProc - 1) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_FB );
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_FB, iProc, kProc, jProc, rank);
		if ( rank /(iProc * jProc ) != kProc - 1) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_FB );

	 // Do extraction, swap, injection for z+ (down to up)
		if ( (rank % kProc) / iProc != jProc - 1 ) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_DT );
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_DT, iProc, kProc, jProc, rank);
		if ( (rank % kProc) / iProc != jProc - 1 ) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_DT );
		
		// Do extraction, swap, injection for z- (up to down)
		if ( (rank % kProc) / iProc != 0 ) extraction( collideField, flagField, cpuDomain, sendBuffer, DIRECTION_TD );
		swap( sendBuffer, readBuffer, sizeBuffer, DIRECTION_TD, iProc, kProc, jProc, rank);
		if ( (rank % kProc) / iProc != 0 ) injection( collideField, flagField, cpuDomain, readBuffer, DIRECTION_TD );
		printf("before streaming\n");
		doStreaming( collideField, streamField, flagField, cpuDomain );

		tmp = collideField;
		collideField = streamField;
		streamField = tmp;
		printf("before collision\n");
		doCollision( collideField, flagField, &tau, cpuDomain );
		printf("before boundary\n");
		treatBoundary( collideField, flagField, velocityWall, cpuDomain );
		printf("%d time\n",t);
		if ( t % timestepsPerPlotting == 0 ) {
				printf("%d cpu x\n", cpuDomain[0]);
				printf("%d cpu y\n", cpuDomain[1]);
				printf("%d cpu z\n", cpuDomain[2]);
				printf("%d rank\n", rank);
				printf("%d xlength\n", xlength[0]);
				printf("%d ylength\n", xlength[1]);
				printf("%d zlength\n", xlength[2]);
				printf("%d iproc\n", iProc);
				printf("%d jproc\n", jProc);
				printf("%d kproc\n\n", kProc);
				printf("Writing the vtk file for timestep # %d \n", t);
			writeVtkOutput( collideField, flagField, "pics/simLB", t, cpuDomain, rank, xlength, iProc, jProc, kProc );
		}
	}
	free(collideField);
	free(streamField);
	free(flagField);

	free(readBuffer[0]);
	free(readBuffer[1]);
	free(readBuffer[2]);
	free(readBuffer[3]);
	free(readBuffer[4]);
	free(readBuffer[5]);
	
	free(sendBuffer[0]);
	free(sendBuffer[1]);
	free(sendBuffer[2]);
	free(sendBuffer[3]);
	free(sendBuffer[4]);
	free(sendBuffer[5]);
	// Terminate MPI
	finalizeMPI();
	
	return 0;
}
#endif

