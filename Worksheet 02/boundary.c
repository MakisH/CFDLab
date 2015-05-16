#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  /* TODO */

  int inv_i;
  int currentCell;
  int neighborCell;
  double f_inv_i;
  double density;
  double c_uwall;
  
  double * each;
  int each_size;

  // z-boundaries
  for (int x = 0; x <= xlength+1; x++) {
    for (int y = 0; y <= xlength+1; y++) {
      // z = -1 plane 
      // i affected:             1,  2,  3,  4,  5
      // corresponding inv(i):  19, 18, 17, 16, 15
      // general rule: inv(i) = Q_NUMBER + 1 - i
      currentCell = x + y*(xlength+2) -1*(xlength+2)*(xlength+2);
      for (int i = 1; i<=5; i++) {
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }

      // z = +1 plane
      // i affected:            15, 16, 17, 18, 19
      // corresponding inv(i):   5,  4,  3,  2,  1
      // general rule: inv(i) = Q_NUMBER + 1 - i
      currentCell = x + y*(xlength+2) +1*(xlength+2)*(xlength+2);
      for (int i = 15; i<=19; i++) {
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }
    }
  }

  // y-boundaries
  for (int x=0; x<=xlength+1; x++) {
    for (int z=0; z<=xlength+1; z++) {

      // y = -1 plane 
      // i affected:             1*,  6,  7,  8, 15*
      // corresponding inv(i):  19, 14, 13,  2,  5
      // general rule: inv(i) = Q_NUMBER + 1 - i
      // *: already touched
      currentCell = x -1*(xlength+2) + z*(xlength+2)*(xlength+2);
      each = {6, 7, 8}; each_size = 3;
      for (int e=0; e<each_size; e++) {
        i = each[e];
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }

      // y = +1 plane 
      // i affected:             5*,  12, 13, 14, 19*
      // corresponding inv(i):  15,  8,  7,  6,  1
      // general rule: inv(i) = Q_NUMBER + 1 - i
      // *: already touched
      currentCell = x +1*(xlength+2) + z*(xlength+2)*(xlength+2);
      each = {12, 13, 14}; each_size = 3;
      for (int e=0; e<each_size; e++) {
        i = each[e];
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }
    }
  }

  // x-boundaries
  for (int y=0; y<=xlength+1; y++) {
    for (int z=0; z<=xlength+1; z++) {
      
      // x = -1 plane 
      // i affected:             2*,  6*,  9, 12*, 16*
      // corresponding inv(i):  18, 14, 11, 18,  4
      // general rule: inv(i) = Q_NUMBER + 1 - i
      // *: already touched
      currentCell = -1 + y*(xlength+2) + z*(xlength+2)*(xlength+2);
      each = {9}; each_size = 1;
      for (int e=0; e<each_size; e++) {
        i = each[e];
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }

      // x = +1 plane 
      // i affected:             4*,  8*, 11, 14*, 18*
      // corresponding inv(i):  16, 12,  9,  6,  2
      // general rule: inv(i) = Q_NUMBER + 1 - i
      // *: already touched
      currentCell = +1 + y*(xlength+2) + z*(xlength+2)*(xlength+2);
      each = {11}; each_size = 1;
      for (int e=0; e<each_size; e++) {
        i = each[e];
        inv_i = Q_NUMBER + 1 - i;
        neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1]*xlength + LATTICEVELOCITIES[inv_i][2]*xlength*xlength + LATTICEVELOCITIES[inv_i][3]*xlength*xlength*xlength;
        f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

        if ( flagField[currentCell] == NO_SLIP ) {

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

        } else if (flagField[currentCell] == MOVING_WALL) {

          computeDensity(collideField+neighborCell, &density);

          c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

          collideField[ Q_NUMBER*currentCell + i] = f_inv_i
            + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

        }
      }
    }
  }



}

