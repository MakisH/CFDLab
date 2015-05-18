#include "LBDefinitions.h"
#include "boundary.h"
#include "computeCellValues.h"
#include <stdio.h>

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
    int x,nx,y,ny,z,nz,i,step=xlength+2;
    double density,dotProd;
    /* NOTE:you can improve the performance of this function by looping over only the boundaries,
     * it will save a number memory access calls, comparisons and control jumps.
     * 
     * However, for our case the performance gain was so insignificant that we decided to
     * stick to this implementation which is more clear and has all the formulas in one place
     * which reduces the probability of error. 
     * 
     * If you are interested in a more efficient implementation, please check 
     * https://github.com/POWER-Morzh/CFDLab02/blob/master/boundary.c
     * where we implemented this function in two ways (second one in comments) and
     * we will gladly substitute current implementation with the mentioned one.
     *  */
    for(x=0;x<step;x++){
        for(y=0;y<step;y++){
            for(z=0;z<step;z++){
                if(flagField[x+y*step+z*step*step]!=FLUID){
                    for(i=0;i<Q_NUMBER;i++){
                        nx=x+LATTICEVELOCITIES[i][0];
                        ny=y+LATTICEVELOCITIES[i][1];
                        nz=z+LATTICEVELOCITIES[i][2];

                        /* We don't need the values outside of our extended domain */
                        if(0<nx && nx<step-1 && 0<ny && ny<step-1 && 0<nz && nz<step-1){
                            if (flagField[x+y*step+z*step*step]==MOVING_WALL){
                                /* Compute density in the neighbour cell */
                                computeDensity(&collideField[Q_NUMBER*(nx+ny*step+nz*step*step)],&density);
                                /* Compute dot product */
                                dotProd=LATTICEVELOCITIES[i][0]*wallVelocity[0]+
                                        LATTICEVELOCITIES[i][1]*wallVelocity[1]+
                                        LATTICEVELOCITIES[i][2]*wallVelocity[2];
                                /* Assign the boudary cell value */
                                collideField[Q_NUMBER*(x+y*step+z*step*step)+i]=	
                                        collideField[Q_NUMBER*(nx+ny*step+nz*step*step)+Q_NUMBER-1-i]+
                                        2*LATTICEWEIGHTS[i]*density*3.0*dotProd;
                            }else if(flagField[x+y*step+z*step*step]==NO_SLIP){
                                collideField[Q_NUMBER*(x+y*step+z*step*step)+i]=
                                        collideField[Q_NUMBER*(nx+ny*step+nz*step*step)+Q_NUMBER-1-i];
                            }
                        }
                    }
                }
            }
        }
    }
}


//void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
//	// idea - on each baoundary we're going to have 5 boundary conditions ... loop over them with 5 double loops and do what you have to do ... 6*x*x < x*x*x ?
//	// what are we doing in the corner between moving wall and no_slip ???
//	int i, inv_i, currentCell, neighborCell;
//	int x_start, x_end, y_start, y_end, z_start, z_end;
//	double f_inv_i, density, c_uwall;
//
//	int each[5]; // trick to implement a "foreach" loop, for every i to be touched
//	int eachSize;
//
//	double SizeX = (xlength + 2);
//	double SizeY = (xlength + 2);
//	double SizeZ = (xlength + 2);
//	double SizeXY = SizeX * SizeY;
//
//	for (int boundary = 1; boundary <= 6; boundary++) {
//		//printf("boundary = %d \n", boundary);
//		// which boundary am I processing now?
//		switch (boundary) {
//		// z = 0 (no-slip by default)
//		case 1 :
//			x_start = 0;        x_end = SizeX - 1;
//			y_start = 0;        y_end = SizeY - 1;
//			z_start = 0;        z_end = 0;
//			each[0] = 1;
//			each[1] = 2;
//			each[2] = 3;
//			each[3] = 4;
//			each[4] = 5;
//			eachSize = 5;
//			break;
//
//		// z = SizeZ (moving wall by default)
//		case 2 :
//			x_start = 0;        x_end = SizeX - 1;
//			y_start = 0;        y_end = SizeY - 1;
//			z_start = SizeZ - 1;  z_end = SizeZ - 1;
//			each[0] = 15;
//			each[1] = 16;
//			each[2] = 17;
//			each[3] = 18;
//			each[4] = 19;
//			eachSize = 5;
//			break;
//
//		// y = 0 (no-slip by default)
//		case 3 :
//			x_start = 0;        x_end = SizeX - 1;
//			y_start = 0;        y_end = 0;
//			z_start = 1;        z_end = SizeZ - 2;
//			each[0] = 6;
//			each[1] = 7;
//			each[2] = 8;
//			each[3] = 0;
//			each[4] = 0;
//			eachSize = 3;
//			break;
//
//		// y = SizeY (no-slip by default)
//		case 4 :
//			x_start = 0;        x_end = SizeX - 1;
//			y_start = SizeY - 1;  y_end = SizeY - 1;
//			z_start = 1;        z_end = SizeZ - 2;
//			each[0] = 12;
//			each[1] = 13;
//			each[2] = 14;
//			each[3] = 0;
//			each[4] = 0;
//			eachSize = 3;
//			break;
//
//		// x = 0 (no-slip by default)
//		case 5 :
//			x_start = 0;        x_end = 0;
//			y_start = 1;        y_end = SizeY - 2;
//			z_start = 1;        z_end = SizeZ - 2;
//			each[0] = 9;
//			each[1] = 0;
//			each[2] = 0;
//			each[3] = 0;
//			each[4] = 0;
//			eachSize = 1;
//			break;
//
//		// x = SizeX (no-slip by default)
//		case 6 :
//			x_start = SizeX - 1;  x_end = SizeX - 1;
//			y_start = 1;        y_end = SizeY - 2;
//			z_start = 1;        z_end = SizeZ - 2;
//			each[0] = 11;
//			each[1] = 0;
//			each[2] = 0;
//			each[3] = 0;
//			each[4] = 0;
//			eachSize = 1;
//			break;
//		}
//
//		// TODO: check correctness of index limits
//		for (int x = x_start; x <= x_end; ++x) {
//			for (int y = y_start; y <= y_end; ++y) {
//				for (int z = z_start; z <= z_end; ++z) {
//					// printf(" iteration %d %d %d Started boundary %d!\n",x,y,z,boundary);
//					// Index of the current cell on the 3D grid (e.g. of flagField). Q not counted.
//					currentCell = x + y*SizeX + z*SizeXY; // current cell from the boundary
//
//					for (int e = 0; e < eachSize; ++e) { 
//						// "Foreach" i on the boundary. We may kill prefetching but there is no foreach in C...
//						i = each[e];
//						// inv(i) - inverse direction of i
//						inv_i = Q_NUMBER - i;
//						// Neighbor cell of current cell in inv(i) direction
//						neighborCell = currentCell + LATTICEVELOCITIES[inv_i][0] + LATTICEVELOCITIES[inv_i][1] * SizeX + LATTICEVELOCITIES[inv_i][2] * SizeXY;
//						// We use f*_inv(i) in both cases (no-slip and moving wall)
//						f_inv_i = collideField[ Q_NUMBER * neighborCell + i];
//
//						// What type of boundary condition do we have?
//						if ( flagField[currentCell] == NO_SLIP ) {
//
//							// update the boundary
//							collideField[ Q_NUMBER*currentCell + i] = f_inv_i;
//
//						} else if (flagField[currentCell] == MOVING_WALL) {
//
//							// density in the neighbor cell
//							computeDensity(collideField+neighborCell, &density);
//
//							// vector product c_i * u_wall
//							c_uwall = LATTICEVELOCITIES[i][0] * wallVelocity[0] + LATTICEVELOCITIES[i][1] * wallVelocity[1] + LATTICEVELOCITIES[i][2] * wallVelocity[2];
//
//							// update the boundary
//							collideField[ Q_NUMBER * currentCell + i] = f_inv_i + 2 * LATTICEWEIGHTS[i] * density * c_uwall / C_S_sq;
//
//						} // if flag
//					} // for each
//
//				} // for z
//			} // for y
//		} // for x
//
//	} // for boundary
//
//} // function

