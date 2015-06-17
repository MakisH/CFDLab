#ifndef _MAIN_C_
#define _MAIN_C_

#include "mpi_helper.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"

int main(int argc, char *argv[]){

	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength[3]; // TO-DO: expand xlength in xlength[3] and correct readParameters !!!!!!!!
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
  int iProc, jProc, kProc;
  int rank = 5;
  int number_of_ranks;

  // send and read buffers for all possible directions :
  // [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
	double *sendBuffer[6];
	double *readBuffer[6];
	int sizeBuffer[6];

  // Read the config file
	readParameters( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv ); // reading parameters from the file.

  // Start MPI
  initializeMPI( &rank, &number_of_ranks, argc, argv );

  // Read the config file
	readParameters( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv ); // reading parameters from the file.

	// Each CPU is going to work in its own subdomain.
	int cpuDomain[3];
	int cpuDomain_size;
	cpuDomain[0] = xlength[0]/iProc;
	cpuDomain[1] = xlength[1]/jProc;
	cpuDomain[2] = xlength[2]/kProc;
	cpuDomain_size = (cpuDomain[0] + 2) * (cpuDomain[1] + 2) * (cpuDomain[2] + 2);

	// Allocating the main three arrays.
	collideField = (double *) malloc(Q_NUMBER * cpuDomain_size * sizeof(double));
	streamField = (double *) malloc(Q_NUMBER * cpuDomain_size * sizeof(double));
	flagField = (int *) malloc(cpuDomain_size * sizeof(int));

	// Init the main three arrays.
	initialiseFields( collideField, streamField, flagField, cpuDomain, iProc, jProc, kProc, rank);

	printf("%d\n", iProc);

	// allocate the buffers
	initialiseBuffers(sendBuffer, readBuffer, cpuDomain, sizeBuffer);
	printf("haha\n");
	int t = 0;
	//for(int t = 0; t <= timesteps; t++){
	//	double *swap = NULL;
 //   
 //   // TODO: maybe move all these to a separate function?
 //   // Do extraction, swap, injection for x+ (left to right)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 1 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 1 );
 //   
 //   // Do extraction, swap, injection for x- (right to left)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 0 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 0 );
 //   
 //   // Do extraction, swap, injection for y+ (back to forth)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 4 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 4 );
 //   
 //   // Do extraction, swap, injection for y- (forth to back)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 5 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 5 );
 //   
 //   // Do extraction, swap, injection for z+ (down to up)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 2 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 2 );
 //   
 //   // Do extraction, swap, injection for z- (up to down)
 //   extraction( collideField, flagField, cpuDomain, sendBuffer, 3 );
 //   // TODO: swap
 //   injection( collideField, flagField, cpuDomain, readBuffer, 3 );
 //   
	//	doStreaming( collideField, streamField, flagField, cpuDomain );

	//	swap = collideField;
	//	collideField = streamField;
	//	streamField = swap;

	//	doCollision( collideField, flagField, &tau, cpuDomain );

	//	treatBoundary( collideField, flagField, velocityWall, cpuDomain );
	//	printf("haha\n");

	//	if ( t % timestepsPerPlotting == 0 ) {
      printf("Writing the vtk file for timestep # %d \n", t);
      writeVtkOutput( collideField, flagField, "pics/simLB", t, cpuDomain, rank, xlength, iProc, jProc, kProc );
 //   }
 //   
	//}

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

