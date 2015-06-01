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

	int x, y, z;

	double density;
	double velocity[3];
	double feq[Q_NUMBER];
	double *currentCell;

	for (z = 1; z <= xlength[2]; ++z){
		for (y = 1; y <= xlength[1]; ++y){
			for (x = 1; x <= xlength[0]; ++x){

				if (flagField[x + (xlength[1] + 2) * y + (xlength[2] + 2) * (xlength[2] + 2) * z] == FLUID)
				{
					// address to the -> first <- distribution function within the respective cell.
					currentCell = collideField + Q_NUMBER * (x + (xlength[1] + 2) * y + (xlength[2] + 2) * (xlength[2] + 2) * z);
					computeDensity (currentCell, &density);
					computeVelocity (currentCell, &density, velocity);
					computeFeq (&density, velocity, feq);
					computePostCollisionDistributions (currentCell, tau, feq);
				}
			}
		}
	}


// OPTIMIZED VERSION
/*	unsigned int xylen2 = (xlength[0] + 2) * (xlength[1] + 2);
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
				if (*(flagField + collideField) == FLUID)
				{
					// locally update cells, can be easily parallelized
					computeDensity (collideField, &density);
					computeVelocity (collideField, &density, velocity);
					computeFeq (&density, velocity, feq);
					computePostCollisionDistributions (collideField, tau, feq);
				}
			}
		}
	}
*/
}
