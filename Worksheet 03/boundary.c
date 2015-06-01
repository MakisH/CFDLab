#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int *xlength, const double * const ref_density, const double * const ref_velocity,
		   const double * const density_in){

  int i, inv_i, mirror_i, currentCell, neighborCell;
  int neighborX, neighborY, neighborZ;
  double f_inv_i, density, c_uwall, velocity;

  double feq[Q_NUMBER];

  int SizeX = (xlength + 2); // Size of the extended domain in each direction
  int SizeY = (xlength + 2);
  int SizeZ = (xlength + 2);
  int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain

  int affected[5], mirror[5]; // Arrays to implement a "foreach" structure for the free-slip condition
  
  // Traverse all the (extended) domain (we need this as we don't have information about the obstacles)
  // Note: we could improve performance by providing a "hint" about in which regions boundaries exist.
  //       Then, we could process the outer boundaries like in Worksheet2 and process just the obstacle regions separately. 
  for (int x = 0; x <= SizeX-1; ++x) {
    for (int y = 0; y <= SizeY-1; ++y) {
      for (int z = 0; z <= SizeZ-1; ++z) {

        // Index of the current cell on the 3D grid (e.g. of flagField). Q not counted.
        currentCell = x + y*SizeX + z*SizeXY; // current boundary cell

        // What kind of (boundary) cell do we process now?
        switch (flagField[currentCell]) {
          
          //----- FLUID -------------------------------------------------------------------------//
          case FLUID :
            break;
          
            
          //----- NO_SLIP -----------------------------------------------------------------------//
          case NO_SLIP :
            // For each direction in the current cell
            for (int i = 0; i < Q_NUMBER; ++i) {

              // Neighbor cell of current cell in i-direction
              neighborX = x + LATTICEVELOCITIES[i][0];
              neighborY = y + LATTICEVELOCITIES[i][1];
              neighborZ = z + LATTICEVELOCITIES[i][2];

              // Check if the neighbor cell coordinates are valid (and not on or outside the limits, that is, outter boundaries)
              if ( neighborX > 0 && neighborX < SizeX-1 && neighborY > 0 && neighborY < SizeY-1 && neighborZ > 0 && neighborZ < SizeZ-1 ) {

                // Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
                neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;

                // inv(i) - inverse direction of i
                inv_i = Q_NUMBER - i - 1;
                
                // Index of the inverse direction of the neighbor cell.
                f_inv_i = collideField[ Q_NUMBER * neighborCell + inv_i ];

                // update the boundary
                collideField[ Q_NUMBER*currentCell + i ] = f_inv_i;
                
              } // if neighbor
            } // for each direction
            break;
            
            
          //----- MOVING_WALL -------------------------------------------------------------------//
          case MOVING_WALL :
            // For each direction in the current cell
            for (int i = 0; i < Q_NUMBER; ++i) {

              // Neighbor cell of current cell in i-direction
              neighborX = x + LATTICEVELOCITIES[i][0];
              neighborY = y + LATTICEVELOCITIES[i][1];
              neighborZ = z + LATTICEVELOCITIES[i][2];

              // Check if the neighbor cell coordinates are valid (and not on or outside the limits, that is, outter boundaries)
              if ( neighborX > 0 && neighborX < SizeX-1 && neighborY > 0 && neighborY < SizeY-1 && neighborZ > 0 && neighborZ < SizeZ-1 ) {

                // Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
                neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;

                // inv(i) - inverse direction of i
                inv_i = Q_NUMBER - i - 1;
                
                // Index of the inverse direction of the neighbor cell.
                f_inv_i = collideField[ Q_NUMBER * neighborCell + inv_i ];

                // density in the neighbor cell
                computeDensity(collideField+Q_NUMBER*neighborCell, &density);

                // vector product c_i * u_wall
                c_uwall = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];
                
                // update the boundary
                collideField[ Q_NUMBER * currentCell + i] = f_inv_i + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;
                
              } // if neighbor
            } // for each direction 
            break;
            
            
          //----- FREE_SLIP ---------------------------------------------------------------------//           
          case FREE_SLIP :
            // For each direction in the current cell
            for (int i = 0; i < Q_NUMBER; ++i) {

              // We consider as neighbors only the cells that share a face.
              // So, we need only the directions perpendicular to the faces.
              if ( i==2 || i==6 || i==8 || i==10 || i==12 || i==16 ) {
                
                // Neighbor cell of current cell in i-direction
                neighborX = x + LATTICEVELOCITIES[i][0];
                neighborY = y + LATTICEVELOCITIES[i][1];
                neighborZ = z + LATTICEVELOCITIES[i][2];

                // Check if the neighbor cell coordinates are valid (and not on or outside the limits, that is, outter boundaries)
                if ( neighborX > 0 && neighborX < SizeX-1 && neighborY > 0 && neighborY < SizeY-1 && neighborZ > 0 && neighborZ < SizeZ-1 ) {

                  // Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
                  neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;

                  // affected: array of affected directions in the current cell
                  // mirror: array of mirrored directions of affected directions
                  // The mirroring depends on the mirroring plane, that is perpendicular to i
                  switch (i) {
                    case 2: // down face [ 0  1  2  3  4 ] --> [ 14 15 16 17 18 ]
                      affected[0] = 0;  mirror[0] = 14;
                      affected[1] = 1;  mirror[1] = 15;
                      affected[2] = 2;  mirror[2] = 16;
                      affected[3] = 3;  mirror[3] = 17;
                      affected[4] = 4;  mirror[4] = 18;
                      break;
                      
                    case 6: // foreground face [ 0  5  6  7 14 ] --> [ 4 11 12 13 18 ]
                      affected[0] = 0;  mirror[0] = 4;
                      affected[1] = 5;  mirror[1] = 11;
                      affected[2] = 6;  mirror[2] = 12;
                      affected[3] = 7;  mirror[3] = 13;
                      affected[4] = 14; mirror[4] = 18;
                      break;
                      
                    case 8: // left face [ 1 11  8  5 15 ] --> [ 3 13 10  7 17 ]
                      affected[0] = 1;  mirror[0] = 3;
                      affected[1] = 11; mirror[1] = 13;
                      affected[2] = 8;  mirror[2] = 10;
                      affected[3] = 5;  mirror[3] = 7;
                      affected[4] = 15; mirror[4] = 17;
                      break;
                      
                    case 10: // right face [ 3 13 10  7 17 ] --> [ 1 11  8  5 14 ]
                      affected[0] = 3;  mirror[0] = 1;
                      affected[1] = 13; mirror[1] = 11;
                      affected[2] = 10; mirror[2] = 8;
                      affected[3] = 7;  mirror[3] = 5;
                      affected[4] = 17; mirror[4] = 14;
                      break;
                      
                    case 12: // background face [ 4 11 12 13 18 ] --> [ 0  5  6  7 14 ]
                      affected[0] = 4;  mirror[0] = 0;
                      affected[1] = 11; mirror[1] = 5;
                      affected[2] = 12; mirror[2] = 6;
                      affected[3] = 13; mirror[3] = 7;
                      affected[4] = 18; mirror[4] = 14;
                      break;
                      
                    case 16: // up face [ 14 15 16 17 18 ] --> [ 0  1  2  3  4 ]
                      affected[0] = 14; mirror[0] = 0;
                      affected[1] = 15; mirror[1] = 1;
                      affected[2] = 16; mirror[2] = 2;
                      affected[3] = 17; mirror[3] = 3;
                      affected[4] = 18; mirror[4] = 4;
                      break;
                  } // switch i                  
                  
                  // foreach affected direction
                  for (e=0; e<=4; ++e) {
                    // Index of the mirrored direction of the neighbor cell.
                    f_mirror = collideField[ Q_NUMBER * neighborCell + mirror[e] ];

                    // update the boundary
                    collideField[ Q_NUMBER*currentCell + affected[e] ] = f_mirror;
                  } // foreach

                } // if neighbor
              } // if i perpendicular
            } // for each direction 
            break;           
            
            
          //----- INFLOW ------------------------------------------------------------------------//  
          case INFLOW :
	    computeFeq(*ref_density, *ref_velocity, feq);
	    collideField[Q_NUMBER * currentCell + i] = feq[i];
            break;
            
            
          //----- OUTFLOW -----------------------------------------------------------------------//           
          case OUTFLOW :
	  computeDensity(currentCell, &density);
	  computeVelocity(currentCell, density, &velocity);

	  computeFeq(*ref_density, velocity, feq);
	  collideField[Q_NUMBER * currentCell + i] = feq[i] + feq[inv_i] - collideField[Q_NUMBER * currentCell + inv_i];
            break;
            
            
          //----- PRESSURE_IN -------------------------------------------------------------------//
          case PRESSURE_IN :
	  computeDensity(currentCell, &density);
	  computeVelocity(currentCell, density, &velocity);

	  computeFeq(*density_in, velocity, feq);
	  collideField[Q_NUMBER * currentCell + i] = feq[i] + feq[inv_i] - collideField[Q_NUMBER * currentCell + inv_i];
            break;
            
            
        } // switch flagField
        
      } // for z
    } // for y
  } // for x

} // function




/*
SPECIAL BOUNDARY CASES TO BE ADDED TO THE MAIN CODE.

// Variables to add to the begining of this file.
double feq[Q_NUMBER]; 
read the ref_velocity;
read ref_density;
read density_in;
double density;
double velocity;


if ( flagField[currentCell] == INFLOW ){
	computeFeq(ref_density, ref_velocity, feq);
	collideField[Q_NUMBER * currentCell + i] = feq[i];
}

if ( flagField[currentCell] == OUTFLOW ){
	computeDensity(currentCell, &density);
	computeVelocity(currentCell, density, &velocity);

	computeFeq(ref_density, velocity, feq);
	collideField[Q_NUMBER * currentCell + i] = feq[i] + feq[inv_i] - collideField[Q_NUMBER * currentCell + inv_i];
}

if ( flagField[currentCell] == INPRESSURE ){
	computeDensity(currentCell, &density);
	computeVelocity(currentCell, density, &velocity);

	computeFeq(density_in, velocity, feq);
	collideField[Q_NUMBER * currentCell + i] = feq[i] + feq[inv_i] - collideField[Q_NUMBER * currentCell + inv_i];
}

*/
