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
	int xlength;
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
  int iProc, jProc, kProc;
  int rank;
  int number_of_ranks;
  // send and read buffers for all possible directions :
  // [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
//   double *sendBuffer[6];
//   double *readBuffer[6];
  
  // Start MPI
  initializeMPI( &rank, &number_of_ranks, argc, argv );

  // Read the config file
	readParameters( &xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, &iProc, &jProc, &kProc, argc, argv ); // reading parameters from the file.

	// Allocating the main three arrays.
	int domain = (xlength + 2) * (xlength + 2) * (xlength + 2);
	collideField = (double *) malloc(Q_NUMBER * domain * sizeof(double));
	streamField = (double *) malloc(Q_NUMBER * domain * sizeof(double));
	flagField = (int *) malloc(domain * sizeof(int));

	// Init the main three arrays.
	initialiseFields( collideField, streamField, flagField, xlength );

	for(int t = 0; t <= timesteps; t++){
		double *swap = NULL;
		doStreaming( collideField, streamField, flagField, xlength );

		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision( collideField, flagField, &tau, xlength );

		treatBoundary( collideField, flagField, velocityWall, xlength );

		if ( t % timestepsPerPlotting == 0 ) {
      printf("Writing the vtk file for timestep # %d \n", t);
      writeVtkOutput( collideField, flagField, "pics/simLB", t, xlength );
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

