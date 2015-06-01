#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include <string.h>

int main(int argc, char *argv[]){
	// declare all input variables
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	char problem[32]; // could be smaller ?
	int xlength[3];
	double tau;
	int timesteps;
	int timestepsPerPlotting;
	double velocityIn;
	double densityIn;
	double densityRef;
	int		initxyzXYZ[6];
	double velocityWall[3];
	// those parameters are becoming too many ... but do we care and can we do something about it ?
	// reading parameters from the file.
	if( readParameters( problem, xlength, &tau, &timesteps, &timestepsPerPlotting, &velocityIn, &densityIn, &densityRef, velocityWall, initxyzXYZ, argc, argv )) return 1; 
	char problem_path[80];
	sprintf(problem_path,"results/%s/%s",problem,problem);
	// Allocating the 3 main arrays.
	int domain = (xlength[0] + 2) * (xlength[1] + 2) * (xlength[2] + 2);
	collideField = (double *) malloc(Q_NUMBER * domain * sizeof(double));
	streamField = (double *) malloc(Q_NUMBER * domain * sizeof(double));
	flagField = (int *) malloc(domain * sizeof(int));
	printf("before init \n");
	// Init the 3 main arrays.
	if ( initialiseFields( collideField, streamField, flagField, xlength, problem, initxyzXYZ)) return 1;
	for(int t = 0; t < timesteps; t++){
		double *swap = NULL;
		printf("martosss is done\n");
		doStreaming( collideField, streamField, flagField, xlength );

		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision( collideField, flagField, &tau, xlength );

		treatBoundary( collideField, flagField, velocityWall, xlength, &densityRef, &velocityIn, &densityIn );

		if ( t % timestepsPerPlotting == 0 ) {
			printf("Writing the vtk file for timestep # %d \n", t);
			writeVtkOutput( collideField, flagField, problem_path, t, xlength ); // needs fixing
		}
	
	}

	free(collideField);
	free(streamField);
	free(flagField);
	return 0;
}

#endif

