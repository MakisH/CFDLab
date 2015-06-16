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
  int rank;
  int number_of_ranks;

  // send and read buffers for all possible directions :
  // [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
   double *sendBuffer[6];
   double *readBuffer[6];

  // Start MPI
  initializeMPI( &rank, &number_of_ranks, argc, argv );

  // Read the config file
	readParameters( &xlength[0], &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv ); // reading parameters from the file.

	// Each CPU is going to work in its own subdomain.
	int innerDomain[3];
	innerDomain[0] = xlength[0]/iProc + 2;
	innerDomain[1] = xlength[1]/jProc + 2;
	innerDomain[2] = xlength[2]/kProc + 2;
	innerDomain_size = innerDomain[0] * innerDomain[1] * innerDomain[2];

	// Allocating the main three arrays.
	collideField = (double *) malloc(Q_NUMBER * innerDomain_size * sizeof(double));
	streamField = (double *) malloc(Q_NUMBER * innerDomain_size * sizeof(double));
	flagField = (int *) malloc(innerDomain_size * sizeof(int));

	// Init the main three arrays.
	initialiseFields( collideField, streamField, flagField, innerDomain );

	// allocate the buffers
	initialiseBuffers(sendBuffer, readBuffer, innerDomain);

	for(int t = 0; t <= timesteps; t++){
		double *swap = NULL;
		doStreaming( collideField, streamField, flagField, innerDomain );

		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision( collideField, flagField, &tau, innerDomain );

		treatBoundary( collideField, flagField, velocityWall, innerDomain );

		if ( t % timestepsPerPlotting == 0 ) {
      printf("Writing the vtk file for timestep # %d \n", t);
      writeVtkOutput( collideField, flagField, "pics/simLB", t, innerDomain );
    }
    
	}

	free(collideField);
	free(streamField);
	free(flagField);
  
  // Terminate MPI
  finalizeMPI();
  
	return 0;
}

#endif

