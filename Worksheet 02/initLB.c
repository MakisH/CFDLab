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
	for (x = 0; x < xlen2; ++x){
		for (y = 0; y < xlen2; ++y){
			for (z = 0; z < xlen2; ++z){
				for (i = 0; i < Q_NUMBER; ++i){
					streamField[Q_NUMBER * (z + y * xlen2 + x * xlen2sq) + i] = LATTICEWEIGHTS[i];
					collideField[Q_NUMBER * (z + y * xlen2 + x * xlen2sq) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}
	/* flagField init: boundary vs. fluid. */
	/* Boundary initialization: (5 walls: no-slip; 1 wall: moving wall). Fluid: inner part. */
	/* Why that many for statements? We used the loop unrolling approach to get rid of the if-statements, which would be present in the third for loop. */

	// Fluid init (inner part of flagField).
	for (x = 1; x <= xlength; ++x) {
		for (y = 1; y <= xlength; ++y) {
			for (z = 1; z <= xlength; ++z) {
				flagField[z + y * xlen2 + x * xlen2sq] = FLUID;
			}
		}
	}

	//unsigned int iter1,iter2, fixed_iter_xy; <- Why do we need that?
	z = xlength + 1;
	for (x = 0; x < xlen2; ++x){
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

/*
	// No slip. For every y, z we set x = 0.
	x = 0;
	for (y=0; y<xlength+2; y++){
		for (z=0; z<xlength+2; z++){
			flagField[x + y*xlength + z*xlength*xlength] = NO_SLIP;
		}
	}

	// No slip. For every x, z we set y = 0.
	y = 0;
	for (x = 0; x < xlength+2; x++){
		for (z = 0; z < xlength+2; z++){
			flagField[x + y*xlength + z*xlength*xlength] = NO_SLIP;
		}
	}

	// No slip. For every x, y we set z = 0.
	z = 0;
		for (x = 0; x < xlength+2; x++){
				for (y = 0; y < xlength+2; y++){
						flagField[x + y*xlength + z*xlength*xlength] = NO_SLIP;
				}
		}

	// No slip. For every y, z we set x = xlength+1.
	x=xlength+1;
		for (y = 0; y < xlength+2; y++){
				for (z = 0; z < xlength+2; z++){
						flagField[x + y*xlength + z*xlength*xlength] = NO_SLIP;
				}
		}

	// No slip. For every x, z we set y = xlength+1;
		y = xlength+1;
		for (x = 0; x < xlength+2; x++){
				for (z = 0; z < xlength+2; z++){
						flagField[x + y*xlength + z*xlength*xlength] = NO_SLIP;
				}
		}

*/

}
