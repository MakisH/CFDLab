#include "boundary.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  /* TODO */

  int inv_i;
  int neighborCell;
  double f_inv_i;
  double density;
  double c_uwall;

  // z = 0 plane 
  // i affected:             1,  2,  3,  4,  5
  // corresponding inv(i):  19, 18, 17, 16, 15
  // general rule: inv(i) = Q_NUMBER + 1 - i
  for (int x = 0; x <= xlength + 1; x++) {
    for (int y = 0; y <= ylength + 1; y++) {
      int currentCell = x + y*(xlength+2);
      for (int i = 1; i<=5; i++) {
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i);

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          density =  computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }
    }
  }


}

