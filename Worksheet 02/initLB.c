#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  
		READ_INT( szFileName, xlength[0] );
		READ_INT( szFileName, xlength[1] );
		READ_INT( szFileName, xlength[2] );
		READ_DOUBLE( szFileName, *tau );

		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );
		read_double( szFileName, "velocityWall1", &velocityWall[0] );
		read_double( szFileName, "velocityWall2", &velocityWall[1] );
		read_double( szFileName, "velocityWall3", &velocityWall[2] );
	}
	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength){
	// might be faster with if-else, but it is insignificant compared to wtkoutput & co
	int x, y, z, i;
	int xlen2 = xlength[0] + 2;
	int ylen2 = xlength[1] + 2;
	int zlen2 = xlength[2] + 2;
	//int xlen2sq = xlen2 * xlen2;
	//int ylen2sq = ylen2 * ylen2;
	int xylen2 = xlen2 * ylen2;

	/* stream & collide Fields initialization. */
	for (z = 0; z < zlen2; ++z){
		for (y = 0; y < ylen2; ++y){
			for (x = 0; x < xlen2; ++x){
				for (i = 0; i < Q_NUMBER; ++i){
					streamField[Q_NUMBER * (x + y * ylen2 + z * xylen2) + i] = LATTICEWEIGHTS[i];
					collideField[Q_NUMBER * (x + y * ylen2 + z * xylen2) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}
	/* flagField init: boundary vs. fluid. */
	/* Boundary initialization: (5 walls: no-slip; 1 wall: moving wall). Fluid: inner part. */

	// Fluid init (inner part of flagField).
	for (z = 1; z <= xlength[2]; ++z) {
		for (y = 1; y <= xlength[1]; ++y) {
			for (x = 1; x <= xlength[0]; ++x) {
				flagField[x + y * ylen2 + z * xylen2] = FLUID;
			}
		}
	}

	// Boundary init.
	for (z = 0; z < zlen2; ++z){ // treat x, y and z as pure iterators, the dimension depends on where we have 0 or xlength+1 and iteration happens.
		for (y = 0; y < ylen2; ++y){
			for (x = 0; x < xlen2; ++x){
				// add all other walls in the same loop, simply switch the indices
				// - x and y are iterators, use x with xlensq for memory locality
				// z is "0" or "xlength + 1"
				// the only question is whether this is OK with the cache memory
				// in fact we can move these iterations above at line 19(before the third loop from 0 to xlength + 1)
				flagField[	  y * ylen2 + z * xylen2	] = NO_SLIP;		// x- dimension	
				flagField[x				+ z * xylen2	] = NO_SLIP;		// y- dimension
				flagField[x + y * ylen2					] = NO_SLIP;		// z- dimension

				flagField[xlen2 - 1 +	y * ylen2 + z * xylen2	] = NO_SLIP;		// x+ dimension
				flagField[x + (ylen2 - 1) * ylen2 + z * xylen2	] = NO_SLIP;		// y+ dimension
				flagField[x + y * ylen2 + (zlen2 - 1) * xylen2	] = MOVING_WALL;	// z+ dimension
			}
		}
	}
}
