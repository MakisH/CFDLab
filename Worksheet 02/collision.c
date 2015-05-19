#include "collision.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
  for (int i = 0; i < Q_NUMBER; ++i) {
		//currentCell[i] -= ( currentCell[i] - feq[i]) / *tau;
		currentCell[i] = ( currentCell[i] * (*tau - 1) + feq[i]) / *tau; // 10-15 sec faster ...
  }
}

/* This is a dirty version of doCollision. However, it is implemented like that because it is faster! */
void doCollision(double *collideField, int *flagField, const double * const tau, int xlength){
	unsigned int xlength2 = (xlength + 2) * (xlength + 2) * Q_NUMBER;
	double density;
	double velocity[3];
	double feq[Q_NUMBER];

	// move the pointer colllideField instead of calculating its offset every time => 10% faster for 100 grid (15-20sec)
	double *stop1;
	double *stop2;
	double *stop3;
	unsigned int stop2_offset = Q_NUMBER * (xlength + 2) * xlength;
	unsigned int stop3_offset = Q_NUMBER * xlength;
	collideField += xlength2 + (xlength + 3) * Q_NUMBER;
	for (stop1 			= collideField + xlength2 * xlength;collideField < stop1; collideField += 2 * Q_NUMBER * (xlength + 2)){	// skip y-boundary rows at plane end
		for (stop2 		= collideField + stop2_offset;		collideField < stop2; collideField += 2 * Q_NUMBER){					// skip x-boundary cells at row end
			for (stop3 	= collideField + stop3_offset; 		collideField < stop3; collideField += Q_NUMBER){
				// locally update cells, can be easily parallelized
				computeDensity (collideField, &density);
				computeVelocity (collideField, &density, velocity);
				computeFeq (&density, velocity, feq);
				computePostCollisionDistributions (collideField, tau, feq);
			}
		}
	}
}

