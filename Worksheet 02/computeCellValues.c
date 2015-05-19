#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
	*density = *currentCell;
	for(unsigned int i = 1; i < Q_NUMBER; ++i)
		*density += currentCell[i];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
// We can unroll loop => many multiplications saved, since Lattice velocities are 0 or +/-1
	velocity[0] = *currentCell * LATTICEVELOCITIES[0][0];
	velocity[1] = *currentCell * LATTICEVELOCITIES[0][1];
	velocity[2] = *currentCell * LATTICEVELOCITIES[0][2];
	for(unsigned int i = 1; i < Q_NUMBER; ++i){
		velocity[0] += currentCell[i] * LATTICEVELOCITIES[i][0];
		velocity[1] += currentCell[i] * LATTICEVELOCITIES[i][1];
		velocity[2] += currentCell[i] * LATTICEVELOCITIES[i][2];
	}
	velocity[0] /= *density;
	velocity[1] /= *density;
	velocity[2] /= *density;
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
	// precalculate inner product 1 - u*u / (2*c_s^2)
	double OneMinusu_u_2c2 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]) * 0.5 / C_S_sq;
	for (int i = 0; i < Q_NUMBER; ++i) {
		// precalculate inner product c*u/c_s^2
		double c_u_c2 = (LATTICEVELOCITIES[i][0] * velocity[0] + LATTICEVELOCITIES[i][1] * velocity[1] + LATTICEVELOCITIES[i][2] * velocity[2]) / C_S_sq;
		// calculate the equilibrium distribution element
		feq[i] = LATTICEWEIGHTS[i] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );
	}
}

