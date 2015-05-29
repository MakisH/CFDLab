#include "collision.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	for (int i = 0; i < Q_NUMBER; ++i) {
		//currentCell[i] -= ( currentCell[i] - feq[i]) / *tau;
		currentCell[i] = ( currentCell[i] * (*tau - 1) + feq[i]) / *tau; // 10-15 sec faster ...
	}
}

// TO DO - better simply go through all cells and check whether FLUID or not - the only solution for an arbitrary geometry?
void doCollision(double *collideField, int *flagField, const double * const tau, int *xlength){
	unsigned int xylen2 = (xlength[0] + 2) * (xlength[1] + 2);
	double density;
	double velocity[3];
	double feq[Q_NUMBER];

	// move the pointer colllideField instead of calculating its offset every time => 10% faster for 100 grid (15-20sec)
	double *stopx;
	double *stopy;
	double *stopz;
	unsigned int stopx_offset = xlength[0] * Q_NUMBER; // process xlength[0] cells in x-dimension
	unsigned int stopy_offset = (xlength[0] + 2) * xlength[1] * Q_NUMBER; // process xlength[1] rows of cells in y-dimension
	collideField += (xylen2 + xlength[0] + 3) * Q_NUMBER; // xylen2 for z-ofset, xlength[0] + 2 for y-offset and 1 for x-offset
	for (stopz			= collideField + xylen2 * xlength[2] * Q_NUMBER;collideField < stopz; collideField += 2 * Q_NUMBER * (xlength[0] + 2)){ // skip y-boundary rows at plane end
		for (stopy		= collideField + stopy_offset;			collideField < stopy; collideField += 2 * Q_NUMBER){ // skip x-boundary cells at row end
			for (stopx	= collideField + stopx_offset;			collideField < stopx; collideField += Q_NUMBER){
				// locally update cells, can be easily parallelized
				computeDensity (collideField, &density);
				computeVelocity (collideField, &density, velocity);
				computeFeq (&density, velocity, feq);
				computePostCollisionDistributions (collideField, tau, feq);
			}
		}
	}
}

