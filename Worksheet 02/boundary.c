#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

	int i, inv_i, currentCell, neighborCell;
	int neighborX, neighborY, neighborZ;
	int x_start, x_end, y_start, y_end, z_start, z_end;
	double f_inv_i, density, c_uwall;

	// each is an array that stores each of the i-directions of the boundary walls that are going
	// to be updated each time.
	int each[5]; // trick to implement a "foreach" loop, for every i direction to be touched

	int SizeX = (xlength + 2); // Size of the extended domain in each direction
	int SizeY = (xlength + 2);
	int SizeZ = (xlength + 2);
	int SizeXY = SizeX * SizeY; // Size of the XY plane of the extended domain

	for (int boundary = 1; boundary <= 6; boundary++) {

		// which boundary am I processing now?
		// Start by traversing each whole 2D boundary plane.
		// Select the indexes each time in a way that no edges/corners are touched two/three times.
		switch (boundary) {
		// z = 0 (no-slip by default)
		case 1 :
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = 0;          z_end = 0;
			each[0] = 19; //1
			each[1] = 18; //2
			each[2] = 17; //3
			each[3] = 16; //4
			each[4] = 15; //5
			break;

		// z = SizeZ (moving wall by default)
		case 2 :
			x_start = 0;          x_end = SizeX - 1;
			y_start = 0;          y_end = SizeY - 1;
			z_start = SizeZ - 1;  z_end = SizeZ - 1;
			each[0] = 5; //15
			each[1] = 4; //16
			each[2] = 3; //17
			each[3] = 2; //18
			each[4] = 1; //19
			break;

		// y = 0 (no-slip by default)
		case 3 :
			x_start = 0;        x_end = SizeX - 1;
			y_start = 0;        y_end = 0;
			z_start = 1;        z_end = SizeZ - 2;
			each[0] = 1; //19
			each[1] = 6; //14
			each[2] = 7; //13
			each[3] = 8; //12
			each[4] = 5; //15
			break;

		// y = SizeY (no-slip by default)
		case 4 :
			x_start = 0;          x_end = SizeX - 1;
			y_start = SizeY - 1;  y_end = SizeY - 1;
			z_start = 1;          z_end = SizeZ - 2;
			each[0] = 15; //5
			each[1] = 8; //12
			each[2] = 7; //13
			each[3] = 6; //14
			each[4] = 1; //19
			break;

		// x = 0 (no-slip by default)
		case 5 :
			x_start = 0;          x_end = 0;
			y_start = 1;          y_end = SizeY - 2;
			z_start = 1;          z_end = SizeZ - 2;
			each[0] = 18; //2
			each[1] = 14; //6
			each[2] = 11; //9
			each[3] = 8; //12
			each[4] = 4; //16
			break;

		// x = SizeX (no-slip by default)
		case 6 :
			x_start = SizeX - 1;  x_end = SizeX - 1;
			y_start = 1;          y_end = SizeY - 2;
			z_start = 1;          z_end = SizeZ - 2;
			each[0] = 16; //4
			each[1] = 12; //8
			each[2] = 9; //11
			each[3] = 6; //14
			each[4] = 2; //18
			break;
		}

		for (int x = x_start; x <= x_end; ++x) {
			for (int y = y_start; y <= y_end; ++y) {
				for (int z = z_start; z <= z_end; ++z) {

					// Index of the current cell on the 3D grid (e.g. of flagField). Q not counted.
					currentCell = x + y*SizeX + z*SizeXY; // current boundary cell

					for (int e = 0; e < 5; ++e) {
						// "Foreach" i direction. We may kill prefetching but there is no foreach in C...
						// The index of the array is in [0,18]
						i = each[e] - 1;

						// inv(i) - inverse direction of i
						inv_i = Q_NUMBER - i - 1;

						// Neighbor cell of current cell in i-direction
						neighborX = x + LATTICEVELOCITIES[i][0];
						neighborY = y + LATTICEVELOCITIES[i][1];
						neighborZ = z + LATTICEVELOCITIES[i][2];

						// Check if the neighbor cell coordinates are valid (and not on or outside the limits, that is, boundaries)
						if ( neighborX > 0 && neighborX < SizeX-1 && neighborY > 0 && neighborY < SizeY-1 && neighborZ > 0 && neighborZ < SizeZ-1 ) {

							neighborCell = neighborX + neighborY*SizeX + neighborZ*SizeXY;

							// We use f*_inv(i) in both cases (no-slip and moving wall)
							f_inv_i = collideField[ Q_NUMBER * neighborCell + inv_i];

							// What type of boundary condition do we have? We could avoid this if by hard-coding different loops for different
							// kinds of boundary conditions, by we would decrease generality.
							if ( flagField[currentCell] == NO_SLIP ) {

								// update the boundary
								collideField[ Q_NUMBER*currentCell + i] = f_inv_i;
							} else if (flagField[currentCell] == MOVING_WALL) {

								// density in the neighbor cell
								computeDensity(collideField+Q_NUMBER*neighborCell, &density);

								// vector product c_i * u_wall
								c_uwall = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];

								// update the boundary
								collideField[ Q_NUMBER * currentCell + i] = f_inv_i + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;
							} // if flag
						} // if neighbor
					} // for each

				} // for z
			} // for y
		} // for x

	} // for boundary

} // function
