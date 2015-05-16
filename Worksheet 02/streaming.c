#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

	unsigned int xlength2 = (xlength + 2) * (xlength + 2);
	for(unsigned int i = 1; i <= xlength; ++i){								// iterate through all FLUID cells
		unsigned int ixx = xlength2 * i;									// improve code readability and optimize computations
		for(unsigned int j = 1; j <= xlength; ++j){
			unsigned int jx = (xlength + 2) * j;							// improve code readability and optimize computations
			for(unsigned int k = 1; k <= xlength; ++k){
				unsigned int current_idx = (k + jx + ixx) * Q_NUMBER;
				for(unsigned int Q_iter = 0; Q_iter < Q_NUMBER; ++Q_iter){	// iterate through all directions
					// - find neighboring cells and take their respective directions ... tricky...
					// - we need the opposite cell to the Q - direction(hence the "-" sign in front of the offset with Lattice velocities)
					// with the velocity in the same direction(hence the Q_iter without minus)
					streamField[current_idx + Q_iter] = 
						collideField[current_idx - (LATTICEVELOCITIES[Q_iter][0]
												  + LATTICEVELOCITIES[Q_iter][1] * (xlength + 2)
												  + LATTICEVELOCITIES[Q_iter][2] * xlength2) * Q_NUMBER + Q_iter] ;
				}
			}
		}
	}
}
