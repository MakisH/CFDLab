#ifndef _MAIN_C_
#define _MAIN_C_

#include "mpi_helper.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
//sleep()
#include <unistd.h>

// free complains without this library ???
#include "stdlib.h"
int main(int argc, char *argv[]){
	printf("hi!\n");
	int rank;
	int np;

	
	// Start MPI
	// slow for more processors ???
	initializeMPI( &rank, &np, argc, argv);
	int xlength[3];
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	int iProc, jProc, kProc;
	int error_code;

	int cpuDomain[3];
	// Read the config file using only one thread
	if(0 == rank){
		error_code = readParameters( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv);
// Error checking
		if(error_code) return error_code;
		else if(0 == iProc || 0 == jProc || 0 == kProc){
			printf("Invalid number of processors in some dimension(cannot be 0)!\n");
			return 3;
		}
		else if(np < iProc * jProc * kProc) {
			printf("not enough processors - %d needed, but only %d are present! (adjust -np or ijkProc parameters)", iProc * jProc * kProc, np);
			return 4;
		}
		else if(xlength[0] % iProc || xlength[1] % jProc || xlength[2] % kProc){
			printf("Non-integer ratio xlength/ijkProc - Invalid !\n");
			return 5;
		}
		cpuDomain[0] = xlength[0] / iProc;
		cpuDomain[1] = xlength[1] / jProc;
		cpuDomain[2] = xlength[2] / kProc;
	}
	MPI_Bcast( xlength, 3, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &iProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &jProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &kProc, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( cpuDomain, 3, MPI_INT, 0, MPI_COMM_WORLD );
	/*
//// test for correct data reading
//	printf("\
//xlength x y z		= %d %d %d\n \
//tau			= %f\n \
//timesteps		= %d\n \
//timestepsPerPlotting	= %d\n \
//velocityWall x y z	= %f %f %f\n \
//iProc			= %d\n \
//jProc			= %d\n \
//kProc			= %d\n \
//cpuDomain x y z	= %d %d %d\n", xlength[0], xlength[1], xlength[2], tau, timesteps, timestepsPerPlotting, velocityWall[0], velocityWall[1], velocityWall[2], iProc, jProc, kProc, cpuDomain[0], cpuDomain[1], cpuDomain[2]);
//	sleep(1000);
*/

	//printf("values:\n xlength: %d %d %d\n Proc ijk %d %d %d\n cpuDomain: %d %d %d\n ",xlength[0],xlength[1],xlength[2],iProc,jProc,kProc,cpuDomain[0],cpuDomain[1],cpuDomain[2]);
	// Init the main three arrays.

	double *collideField = (double *) malloc(Q_NUMBER * (cpuDomain[0] + 2) * (cpuDomain[1] + 2) * (cpuDomain[2] + 2) * sizeof(double));
	double *streamField = (double *) malloc(Q_NUMBER * (cpuDomain[0] + 2) * (cpuDomain[1] + 2) * (cpuDomain[2] + 2) * sizeof(double));
	int *flagField = (int *) malloc((cpuDomain[0] + 2) * (cpuDomain[1] + 2) * (cpuDomain[2] + 2) * sizeof(int));

	int neighbor[6];	// default action - do nothing if on the boundary
	initialiseFields( collideField, streamField, flagField, cpuDomain, iProc, jProc, kProc, rank, neighbor);
	//printf("neighbors:\n%d %d %d %d %d %d\n",neighbor[0], neighbor[1], neighbor[2], neighbor[3], neighbor[4], neighbor[5]);
	// precomputed values for extraction and injection bounds depending on the direction
	const side Bsides_ext[6] = {{cpuDomain[0], cpuDomain[0], 0, cpuDomain[1] + 1, 0, cpuDomain[2] + 1},
															{1, 1,											 0, cpuDomain[1] + 1, 0, cpuDomain[2] + 1},
															{0, cpuDomain[0] + 1, 0, cpuDomain[1] + 1, cpuDomain[2], cpuDomain[2]},
															{0, cpuDomain[0] + 1, 0, cpuDomain[1] + 1, 1, 1},
															{0, cpuDomain[0] + 1, cpuDomain[1], cpuDomain[1], 0, cpuDomain[2] + 1},
															{0, cpuDomain[0] + 1, 1, 1,												0, cpuDomain[2] + 1}};
	const side Bsides_inj[6] = {{0, 0, 0, cpuDomain[1] + 1, 0, cpuDomain[2] + 1},
															{cpuDomain[0] + 1, cpuDomain[0] + 1, 0, cpuDomain[1] + 1, 0, cpuDomain[2] + 1},
															{0, cpuDomain[0] + 1, 0, cpuDomain[1] + 1, 0, 0},
															{0, cpuDomain[0] + 1, 0, cpuDomain[1] + 1, cpuDomain[2] + 1, cpuDomain[2] + 1},
															{0, cpuDomain[0] + 1, 0, 0, 0, cpuDomain[2] + 1},
															{0, cpuDomain[0] + 1, cpuDomain[1] + 1, cpuDomain[1] + 1, 0, cpuDomain[2] + 1}};

// send and read buffers for all possible directions :
	// [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
	double *sendBuffer[6] = {NULL, NULL, NULL, NULL, NULL, NULL}; // needed for free, otherwise gives a stupid error.
	double *readBuffer[6] = {NULL, NULL, NULL, NULL, NULL, NULL};;
	int sizeBuffer[6];
	initialiseBuffers(sendBuffer, readBuffer, cpuDomain, sizeBuffer, neighbor);
	double *tmp = NULL;
	for(int t = 0; t <= timesteps; t++){
		// TODO: maybe move all these to a separate function?
		//if(!rank)
		//	printf("t = %d, rank = %d\n",t,rank);
		// Do extraction, swap, injection for

		// swap needs two checks so they're in the function

				// x+ (left to right)
		if (neighbor[DIR_R] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_R, Bsides_ext);
		//printf("after extr\n");
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_R, neighbor);
		//printf("after swap\n");
		if (neighbor[DIR_L] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer, DIR_R, Bsides_inj);
		//printf("after inj\n");
		// x- (right to left)
		if (neighbor[DIR_L] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_L, Bsides_ext);
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_L, neighbor);
		if (neighbor[DIR_R] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer, DIR_L, Bsides_inj);
		//printf("after 2\n");
		// z+ (down to up)
		if (neighbor[DIR_T] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_T, Bsides_ext);
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_T, neighbor);
		if (neighbor[DIR_D] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer, DIR_T, Bsides_inj);
		//printf("after 3\n");
		// z- (up to down)
		if (neighbor[DIR_D] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_D, Bsides_ext);
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_D, neighbor);
		if (neighbor[DIR_T] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer, DIR_D, Bsides_inj);
		//printf("after 4\n");
		// y- (forth to back)
		if (neighbor[DIR_B] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_B, Bsides_ext);
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_B, neighbor);
		if (neighbor[DIR_F] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer, DIR_B, Bsides_inj);
		//printf("after 5\n");
		// y+ (back to forth)
		if (neighbor[DIR_F] != MPI_PROC_NULL)	extraction( collideField, cpuDomain, sendBuffer, DIR_F, Bsides_ext);
		swap( sendBuffer, readBuffer, sizeBuffer, DIR_F, neighbor);
		if (neighbor[DIR_B] != MPI_PROC_NULL)	injection( collideField, cpuDomain, readBuffer,DIR_F, Bsides_inj);

		//printf("before boundary\n");
		treatBoundary( collideField, flagField, velocityWall, cpuDomain);
		//printf("before streaming\n");
		doStreaming( collideField, streamField, flagField, cpuDomain);

		tmp = collideField;
		collideField = streamField;
		streamField = tmp;
		//printf("before collision\n");
		doCollision( collideField, flagField, &tau, cpuDomain );
		//printf("after collision\n");
		
		//printf("%d time\n",t);
		if ( t % timestepsPerPlotting == 0 ) {
				//printf("%d cpu x\n", cpuDomain[0]);
				//printf("%d cpu y\n", cpuDomain[1]);
				//printf("%d cpu z\n", cpuDomain[2]);
				//printf("%d rank\n", rank);
				//printf("%d xlength\n", xlength[0]);
				//printf("%d ylength\n", xlength[1]);
				//printf("%d zlength\n", xlength[2]);
				//printf("%d iproc\n", iProc);
				//printf("%d jproc\n", jProc);
				//printf("%d kproc\n\n", kProc);
			if(!rank)
				printf("Write vtk for time # %d \n", t);
			writeVtkOutput( collideField, flagField, "pics/simLB", t, cpuDomain, rank, xlength, iProc, jProc, kProc );
		}
	}
	free((void *)collideField);
	free((void *)streamField);
	free((void *)flagField);

	free((void *)readBuffer[0]);
	free((void *)readBuffer[1]);
	free((void *)readBuffer[2]);
	free((void *)readBuffer[3]);
	free((void *)readBuffer[4]);
	free((void *)readBuffer[5]);
	
	free((void *)sendBuffer[0]);
	free((void *)sendBuffer[1]);
	free((void *)sendBuffer[2]);
	free((void *)sendBuffer[3]);
	free((void *)sendBuffer[4]);
	free((void *)sendBuffer[5]);
	finalizeMPI();
	
	return 0;
}
#endif