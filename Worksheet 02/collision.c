#include "collision.h"
#include "LBDefinitions.h"


void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
  /* TODO */ // Be careful! A C guy should check if what I do really makes sense!
  for (int i=0; i<Q_NUMBER; i++) {
        *(currentCell + i) -= ( 1 / tau ) * ( *(currentCell + i) - feq[i]; 
  }

}

void doCollision(double *collideField, int *flagField, const double * const tau, int xlength){

	int x, y, z;

	double density;
	double velocity[3];
	double feq[Q_NUMBER];
	double *currentCell;

	for (x = 1; x < xlength + 1; x++){
		for (y = 1; y < xlength + 1; y++){
			for (z = 1; z < xlength + 1; z++){

				// address to the -> first <- distribution function within the respective cell.
				currentCell = collideField + Q_NUMBER * (x + xlength*y + xlength*xlength*z);

				computeDensity (currentCell, &density);
				computeVelocity (currentCell, &density, velocity);
				computeFeq (&density, velocity, feq);
				computePostCollisionDistributions (currentCell, tau, feq);

			}
		}
	}


}

