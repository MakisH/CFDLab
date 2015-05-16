#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

  int i, inv_i, currentCell, neighborCell;
  int x_start, x_end, y_start, y_end, z_start, z_end;
  double f_inv_i, density, c_uwall;
  
  double each[5]; // trick to implement a "foreach" loop, for every i to be touched
  int each_size;

  double SizeX = (xlength+2);
  double SizeY = (xlength+2);
  double SizeZ = (xlength+2);
  double SizeXY = SizeX * SizeY;

  for (int boundary = 1; boundary <= 6; boundary++) {

    // which boundary am I processing now?
    switch (boundary) {
      // z = 0 (no-slip by default)
      case 1 :
        x_start = 0;      x_end = SizeX;
        y_start = 0;      y_end = SizeY;
        z_start = 0;      z_end = 0;
        each = {1, 2, 3, 4, 5}; eachSize = 5;
        break;

      // z = SizeZ (moving wall by default)
      case 2 :
        x_start = 0;      x_end = SizeX;
        y_start = 0;      y_end = SizeY;
        z_start = SizeZ;  z_end = SizeZ;
        each = {15, 16, 17, 18, 19}; eachSize = 5;
        break;

      // y = 0 (no-slip by default)
      case 3 :
        x_start = 0;      x_end = SizeX;
        y_start = 0;      y_end = 0;
        z_start = 0;      z_end = SizeZ;
        each = {6, 7, 8}; eachSize = 3;
        break;

      // y = SizeY (no-slip by default)
      case 4 :
        x_start = 0;      x_end = SizeX;
        y_start = SizeY;  y_end = SizeY;
        z_start = 0;      z_end = SizeZ;
        each = {12, 13, 14}; eachSize = 3;
        break;

      // x = 0 (no-slip by default)
      case 5 :
        x_start = 0;      x_end = 0;
        y_start = 0;      y_end = SizeY;
        z_start = 0;      z_end = SizeZ;
        each = {9}; eachSize = 1;
        break;

      // x = SizeX (no-slip by default)
      case 6 :
        x_start = SizeX;  x_end = SizeX;
        y_start = 0;      y_end = SizeY;
        z_start = 0;      z_end = SizeZ;
        each = {11}; eachSize = 1;
        break;
    }

    // TODO: check correctness of index limits
    for (int x = x_start; x <= x_end; x++) {
      for (int y = y_start; y <= y_end; y++) {
        for (int z = z_start; z <= z_end; z++) {
       
          // Index of the current cell on the 3D grid (e.g. of flagField). Q not counted.
          currentCell = x + y*SizeX + z*SizeXY;

          for (int e=0; e<each_size; e++) { 
            // "Foreach" i on the boundary. We may kill prefetching but there is no foreach in C...
            i = each[e];
            // inv(i) - inverse direction of i
            inv_i = Q_NUMBER + 1 - i;
            // Neighbor cell of current cell in inv(i) direction
            neighborCell = currentCell + LATTICEVELOCITIES[inv_i][1] + LATTICEVELOCITIES[inv_i][2]*SizeX + LATTICEVELOCITIES[inv_i][3]*SizeXY;
            // We use f*_inv(i) in both cases (no-slip and moving wall)
            f_inv_i = collideField[ Q_NUMBER*neighborCell + inv_i];

            // What type of boundary condition do we have?
            if ( flagField[currentCell] == NO_SLIP ) {

              // update the boundary
              collideField[ Q_NUMBER*currentCell + i] = f_inv_i;

            } else if (flagField[currentCell] == MOVING_WALL) {

              // density in the neighbor cell
              computeDensity(collideField+neighborCell, &density);

              // vector product c_i * u_wall
              c_uwall = LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2] + LATTICEVELOCITIES[i][3] * wallVelocity[3];

              // update the boundary
              collideField[ Q_NUMBER*currentCell + i] = f_inv_i + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;

            } // if flag
          } // for each    

        } // for z
      } // for y
    } // for x

  } // for boundary

} // function

