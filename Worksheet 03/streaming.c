#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int *xlength){

	unsigned int xylen2 = (xlength[0] + 2) * (xlength[1] + 2);

	for(unsigned int z = 1; z <= xlength[2]; ++z){					// iterate through all FLUID cells
		unsigned int izz = xylen2 * z;											// improve code readability and optimize computations

		for(unsigned int y = 1; y <= xlength[1]; ++y){
			unsigned int jy = (xlength[1] + 2) * y;							// improve code readability and optimize computations

			for(unsigned int x = 1; x <= xlength[0]; ++x){
				unsigned int current_idx = (x + jy + izz) * Q_NUMBER;

				for(unsigned int Q_iter = 0; Q_iter < Q_NUMBER; ++Q_iter){	// iterate through all directions
					// this is very cache ineficient, because ~9 different cache lines are accessed for each cell update(3 per z-plane for each row)

					// - find neighboring cells and take their respective directions ... tricky...
					// - we need the opposite cell to the Q - direction(hence the "-" sign in front of the offset with Lattice velocities)
					// with the velocity in the same direction(hence the 2nd Q_iter without minus)
					streamField[current_idx + Q_iter] = 
						collideField[current_idx	- (LATTICEVELOCITIES[Q_iter][0]
																			+  LATTICEVELOCITIES[Q_iter][1] * (xlength[1] + 2)
																			+  LATTICEVELOCITIES[Q_iter][2] * xylen2) * Q_NUMBER + Q_iter] ;
				}
			}
		}
	}
}
