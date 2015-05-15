#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
/** computes the density from the pdf-s stored at currentCell.
 *  currentCell thus denotes the address of the first pdf of the
 *  respective cell. The result is stored in density.
 */
	*density = *currentCell;
	for(unsigned int i = 1; i < Q_NUMBER; ++i)
		*density += currentCell[i];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
  /** computes the velocity within currentCell and stores the result in velocity */
	
// if this is true, it can be heavily optimized!
// Table1(page3) - loop i => +/-/0 depend on currentCell ( use enums somehow ?)
	velocity[0]	= *currentCell * LATTICEVELOCITIES[0][0];
	velocity[1] = *currentCell * LATTICEVELOCITIES[0][1];
	velocity[2] = *currentCell * LATTICEVELOCITIES[0][2];
	for(unsigned int i = 1; i < Q_NUMBER; ++i){
		velocity[0] += currentCell[i] * LATTICEVELOCITIES[i][0];
		velocity[1] += currentCell[i] * LATTICEVELOCITIES[i][1];
		velocity[2] += currentCell[i] * LATTICEVELOCITIES[i][2];
	}
	velocity[0] /=*density;
	velocity[1] /=*density;
	velocity[2] /=*density;
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
 /* TODO */
}

