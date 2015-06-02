#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField,
									 int* flagField,
									 const double * const wallVelocity,
									 int *xlength,
									 const double * const ref_density,
									 const double * const velocityIn,
									 const double * const density_in){

	int i, inv_i, currentCell, neighborCell;
	int neighborX, neighborY, neighborZ;
	double f_inv_i, density, c_uwall, velocity;
	double feq[Q_NUMBER];

	int SizeX = xlength[0] + 2; // Size of the extended domain in each direction
	int SizeY = xlength[1] + 2;
	int SizeZ = xlength[2] + 2;
	int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain

	// Arrays to implement a "foreach" structure for the free-slip condition
	int perpendicular[6] = {2, 6, 8, 10, 12, 16};
	int affected[5], mirror[5]; 
	
	// Traverse all the (extended) domain (we need this as we don't have information about the obstacles)
	// Note: we could improve performance by providing a "hint" about in which regions boundaries exist.
	//       Then, we could process the outer boundaries like in Worksheet2 and process just the obstacle regions separately. 
	//printf("Sizes X Y Z XY %d %d %d %d\n", SizeX, SizeY, SizeZ, SizeXY);
	
	for (int z = 0; z < SizeZ; ++z) {
		for (int y = 0; y < SizeY; ++y) {
			for (int x = 0; x < SizeX; ++x) {
				// Index of the current cell on the 3D grid (e.g. of flagField). Q not counted.
				currentCell = x + y*SizeX + z*SizeXY; // current boundary cell
			//	if(flagField[currentCell] != FLUID) printf("%d %d start \n",currentCell,flagField[currentCell]);
			//	if(flagField[currentCell] == 1) printf("%d is NO SLIP\n",currentCell);
			//	if( 636 <  currentCell) printf("%d\n", flagField[currentCell]);
				// What kind of (boundary) cell do we process now?
				switch (flagField[currentCell]) {
					
					//----- FLUID -------------------------------------------------------------------------//
					case FLUID :
					//	printf("%d \n",currentCell);
						break;
					//----- NO_SLIP -----------------------------------------------------------------------//
					case NO_SLIP :
					//	printf("Enter No-Hell %d is NO SLIP\n",currentCell);
						// For each direction in the current cell
						for (int i = 0; i < Q_NUMBER; ++i) {

							// Neighbor cell of current cell in i-direction
							neighborX = x + LATTICEVELOCITIES[i][0];
							neighborY = y + LATTICEVELOCITIES[i][1];
							neighborZ = z + LATTICEVELOCITIES[i][2];
								
							// Check if the the coordinates of the neighbor cell are valid
							if ( neighborX >= 0 && neighborX <= SizeX-1 && neighborY >= 0 && neighborY <= SizeY-1 && neighborZ >= 0 && neighborZ <= SizeZ-1 ) {
							
								// Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
								neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;
								
								// Check if the neighbor cell is fluid
								if ( flagField[neighborCell] == FLUID ) {
									
									// inv(i) - inverse direction of i
									inv_i = Q_NUMBER - i - 1;
									
									// Index of the inverse direction of the neighbor cell.
									f_inv_i = collideField[ Q_NUMBER * neighborCell + inv_i ];

									// update the boundary
									collideField[ Q_NUMBER * currentCell + i ] = f_inv_i;
								//  printf("boundary temp %d\n", inv_i);

								} // if neighbor is fluid
							} // if neighbor coordinates
						} // for each direction
				//		printf("%d \n",currentCell);
						break;
						
					//----- MOVING_WALL -------------------------------------------------------------------//
					case MOVING_WALL :
						// For each direction in the current cell
						for (int i = 0; i < Q_NUMBER; ++i) {

							// Neighbor cell of current cell in i-direction
							neighborX = x + LATTICEVELOCITIES[i][0];
							neighborY = y + LATTICEVELOCITIES[i][1];
							neighborZ = z + LATTICEVELOCITIES[i][2];
								
							// Check if the the coordinates of the neighbor cell are valid
							if ( neighborX >= 0 && neighborX <= SizeX-1 && neighborY >= 0 && neighborY <= SizeY-1 && neighborZ >= 0 && neighborZ <= SizeZ-1 ) {

								// Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
								neighborCell = neighborX + neighborY * SizeX + neighborZ * SizeXY;

								// Check if the neighbor cell is fluid
								if ( flagField[neighborCell] == FLUID ) {

									// inv(i) - inverse direction of i
									inv_i = Q_NUMBER - i - 1;
									
									// Index of the inverse direction of the neighbor cell.
									f_inv_i = collideField[ Q_NUMBER * neighborCell + inv_i ];

									// density in the neighbor cell
									computeDensity(collideField + Q_NUMBER * neighborCell, &density);

									// vector product c_i * u_wall
									c_uwall = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];
									
									// update the boundary
									collideField[ Q_NUMBER * currentCell + i] = f_inv_i + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;
									
								} // if neighbor is fluid
							} // if neighbor coordinates
						} // for each direction 
			//			printf("%d \n",currentCell);
						break;
						
						
					//----- FREE_SLIP ---------------------------------------------------------------------//           
					case FREE_SLIP :
						// For each direction in the current cell that is perpendicular to a face
						for (int p = 0; p < 6; ++p) {

							// We want i to take only the values-directions that are perpendicular to a face,
							// as we assume that neighbors for the free-slip are only the cells that share an interface.
							i = perpendicular[p];

							// Neighbor cell of current cell in i-direction
							neighborX = x + LATTICEVELOCITIES[i][0];
							neighborY = y + LATTICEVELOCITIES[i][1];
							neighborZ = z + LATTICEVELOCITIES[i][2];
							
							// Check if the the coordinates of the neighbor cell are valid
							if ( neighborX >= 0 && neighborX <= SizeX-1 && neighborY >= 0 && neighborY <= SizeY-1 && neighborZ >= 0 && neighborZ <= SizeZ-1 ) {

								// Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
								neighborCell = neighborX + neighborY * SizeX + neighborZ * SizeXY;
							
								// Check if the neighbor cell is fluid
								if ( flagField[neighborCell] == FLUID ) {
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
											
										case 10: // right face [ 3 13 10  7 17 ] --> [ 1 11  8  5 15 ]
											affected[0] = 3;  mirror[0] = 1;
											affected[1] = 13; mirror[1] = 11;
											affected[2] = 10; mirror[2] = 8;
											affected[3] = 7;  mirror[3] = 5;
											affected[4] = 17; mirror[4] = 15;
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

										default : printf("Error we should never be here!\n"); break;
									} // switch i                  
									
									// foreach affected direction
									for (int e = 0; e < 5; ++e) {
										
										// update the boundary														// Index of the mirrored direction of the neighbor cell.
										collideField[ Q_NUMBER * currentCell + affected[e] ] = collideField[ Q_NUMBER * neighborCell + mirror[e] ];
									} // foreach affected
									
								} // if neighbor is fluid
							} // if neighbor coordinates
						} // for each direction 
			//			printf("%d \n",currentCell);
						break;           
						
					//----- INFLOW ------------------------------------------------------------------------//  
					case INFLOW :
						// Compute the equilibrium distribution for the reference density and velocity
						computeFeq(ref_density, velocityIn, &collideField[Q_NUMBER * currentCell]);
		//				printf("%d \n",currentCell);
						break;
					//----- OUTFLOW -----------------------------------------------------------------------//           
					case OUTFLOW :

						// For each direction in the current cell
						for (int i = 0; i < Q_NUMBER; ++i) { 
							// Neighbor cell of current cell in i-direction
							neighborX = x + LATTICEVELOCITIES[i][0];
							neighborY = y + LATTICEVELOCITIES[i][1];
							neighborZ = z + LATTICEVELOCITIES[i][2];
						
							// Check if the the coordinates of the neighbor cell are valid
							if ( neighborX >= 0 && neighborX <= SizeX-1 && neighborY >= 0 && neighborY <= SizeY-1 && neighborZ >= 0 && neighborZ <= SizeZ-1 ) {

                // Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
                neighborCell = neighborX + neighborY * SizeX + neighborZ * SizeXY;
                
                // Check if the neighbor cell is fluid                
								if ( flagField[neighborCell] == FLUID ) {
                  computeDensity(collideField + neighborCell * Q_NUMBER, &density);
                  computeVelocity(collideField + neighborCell * Q_NUMBER, &density, &velocity);
                  computeFeq(ref_density, &velocity, feq);

									inv_i = Q_NUMBER - i - 1;
									collideField[Q_NUMBER * currentCell + i] = feq[inv_i] + feq[i] - collideField[Q_NUMBER * neighborCell + inv_i];
								} // if neighbor is fluid
							} // if neighbor coordinates
						} // for each direction
						break;
					//----- PRESSURE_IN -------------------------------------------------------------------//
					case PRESSURE_IN :

						// For each direction in the current cell
						for (int i = 0; i < Q_NUMBER; ++i) { 
							// Neighbor cell of current cell in i-direction
							neighborX = x + LATTICEVELOCITIES[i][0];
							neighborY = y + LATTICEVELOCITIES[i][1];
							neighborZ = z + LATTICEVELOCITIES[i][2];
							
							// Check if the the coordinates of the neighbor cell are valid
							if ( neighborX >= 0 && neighborX <= SizeX-1 && neighborY >= 0 && neighborY <= SizeY-1 && neighborZ >= 0 && neighborZ <= SizeZ-1 ) {

                // Index of the neighbor cell on the 3D grid (e.g. of flagField). Q not counted.
                neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;
                
                // Check if the neighbor cell is fluid
								if ( flagField[neighborCell] == FLUID ) {
									computeDensity(collideField + neighborCell * Q_NUMBER, &density);
									computeVelocity(collideField + neighborCell * Q_NUMBER, &density, &velocity);
									computeFeq(density_in, &velocity, feq);
									inv_i = Q_NUMBER - i - 1;
									collideField[Q_NUMBER * currentCell + i] = feq[inv_i] + feq[i] - collideField[Q_NUMBER * neighborCell + inv_i];
								} // if neighbor is fluid
								
							} // if neighbor coordinates
						} // for each direction
						break;
						
				} // switch flagField
		 //  printf("%d done \n",currentCell);
			} // for x
		} // for y
	} // for z

} // function
