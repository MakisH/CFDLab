#include "collision.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	for (int i = 0; i < Q_NUMBER; ++i) {
		//currentCell[i] -= ( currentCell[i] - feq[i]) / *tau;
		currentCell[i] = ( currentCell[i] * (*tau - 1) + feq[i]) / *tau; // 10-15 sec faster ...
	}
}

void doCollision(double *collideField, int *flagField, const double * const tau, int *xlength){

	int x, y, z;

	double density;
	double velocity[3];
	double feq[Q_NUMBER];
	double *currentCell;

	for (z = 1; z <= xlength[2]; ++z){
		for (y = 1; y <= xlength[1]; ++y){
			for (x = 1; x <= xlength[0]; ++x){

				if (flagField[x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z] == FLUID)
				{
					// address to the -> first <- distribution function within the respective cell.
					currentCell = collideField + Q_NUMBER * (x + (xlength[0] + 2) * y + (xlength[0] + 2) * (xlength[1] + 2) * z);
					computeDensity (currentCell, &density);
					computeVelocity (currentCell, &density, velocity);
					computeFeq (&density, velocity, feq);
					computePostCollisionDistributions (currentCell, tau, feq);
				}
			}
		}
	}
}
