#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  
		READ_INT( szFileName, *xlength );
		READ_DOUBLE( szFileName, *tau );
		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );
		read_double( szFileName, "velocityWall1", &velocityWall[0] );
		read_double( szFileName, "velocityWall2", &velocityWall[1] );
		read_double( szFileName, "velocityWall3", &velocityWall[2] );
	}
	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	// might be faster with if-else, but it is insignificant compared to wtkoutput & co
	int x, y, z, i;
	int xlen2 = xlength + 2;
	int xlen2sq = xlen2 * xlen2;

	/* stream & collide Fields initialization. */
	for (z = 0; z < xlen2; ++z){
		for (y = 0; y < xlen2; ++y){
			for (x = 0; x < xlen2; ++x){
				for (i = 0; i < Q_NUMBER; ++i){
					streamField[Q_NUMBER * (x + y * xlen2 + z * xlen2sq) + i] = LATTICEWEIGHTS[i];
					collideField[Q_NUMBER * (x + y * xlen2 + z * xlen2sq) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}
	/* flagField init: boundary vs. fluid. */
	/* Boundary initialization: (5 walls: no-slip; 1 wall: moving wall). Fluid: inner part. */

	// Fluid init (inner part of flagField).
	for (z = 1; z <= xlength; ++z) {
		for (y = 1; y <= xlength; ++y) {
			for (x= 1; x <= xlength; ++x) {
				flagField[x + y * xlen2 + z * xlen2sq] = FLUID;
			}
		}
	}

	// Boundary init.
	z = xlength + 1;
	for (x = 0; x < xlen2; ++x){ // treat x, y and z as pure iterators, the dimension depends on where we have 0 or xlength+1 and iteration happens.
		for (y = 0; y < xlen2; ++y){
			// add all other walls in the same loop, simply switch the indices
			// - x and y are iterators, use x with xlensq for memory locality
			// z is "0" or "xlength + 1"
			// the only question is whether this is OK with the cache memory
			// in fact we can move these iterations above at line 19(before the third loop from 0 to xlength + 1)
			flagField[	  y * xlen2 + x * xlen2sq	] = NO_SLIP;		// x- dimension	
			flagField[y				+ x * xlen2sq	] = NO_SLIP;		// y- dimension
			flagField[y + x * xlen2					] = NO_SLIP;		// z- dimension

			flagField[z + y * xlen2 + x * xlen2sq	] = NO_SLIP;		// x+ dimension
			flagField[y + z * xlen2 + x * xlen2sq	] = NO_SLIP;		// y+ dimension
			flagField[y + x * xlen2 + z * xlen2sq	] = MOVING_WALL;	// z+ dimension
		}
	}
}
