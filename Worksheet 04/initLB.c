#include "initLB.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int *iProc, int *jProc, int *kProc, int argc, char *argv[]){
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  
		read_int( szFileName, "xlength", xlength );
		read_int( szFileName, "ylength", xlength+1 );
		read_int( szFileName, "zlength", xlength+2 );
		READ_DOUBLE( szFileName, *tau );
		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );
		read_double( szFileName, "velocityWall1", &velocityWall[0] );
		read_double( szFileName, "velocityWall2", &velocityWall[1] );
		read_double( szFileName, "velocityWall3", &velocityWall[2] );
		READ_INT( szFileName, *iProc );
		READ_INT( szFileName, *jProc );
		READ_INT( szFileName, *kProc );
	}
	return 0;
}

void initialiseBuffers(double **sendBuffer, double **readBuffer,  int *xlength, int *sizeBuffer){

	int x = xlength[0]+2;
	int y = xlength[1]+2;
	int z = xlength[2]+2;

	int domain = 5; // because we have 5 possible directions to be extracted to buffer

	// We should substitute the sizes in malloc, but we don't have time now.
	sizeBuffer[0] = y * z * domain;
	sizeBuffer[1] = y * z * domain;
	sizeBuffer[2] = x * y * domain;
	sizeBuffer[3] = x * y * domain;
	sizeBuffer[4] = x * z * domain;
	sizeBuffer[5] = x * z * domain;

	// We initilise 6 different buffers.
	// sendBuffer planes[0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]

	sendBuffer[0] = (double *) malloc(y * z * domain * sizeof(double)); // left plane
	sendBuffer[1] = (double *) malloc(y * z * domain * sizeof(double)); // right plane
	sendBuffer[2] = (double *) malloc(x * y * domain * sizeof(double)); // top plane
	sendBuffer[3] = (double *) malloc(x * y * domain * sizeof(double)); // bottom plane
	sendBuffer[4] = (double *) malloc(x * z * domain * sizeof(double)); // front plane
	sendBuffer[5] = (double *) malloc(x * z * domain * sizeof(double)); // back plane

	// readBuffer planes[0:right sendBuffer, 1:left sendBuffer, 2:bottom sendBuffer, 3:top sendBuffer, 4:back sendBuffer, 5:front sendBuffer]
	readBuffer[0] = (double *) malloc(y * z * domain * sizeof(double)); // left plane
	readBuffer[1] = (double *) malloc(y * z * domain * sizeof(double)); // right plane
	readBuffer[2] = (double *) malloc(x * y * domain * sizeof(double)); // top plane
	readBuffer[3] = (double *) malloc(x * y * domain * sizeof(double)); // bottom plane
	readBuffer[4] = (double *) malloc(x * z * domain * sizeof(double)); // front plane
	readBuffer[5] = (double *) malloc(x * z * domain * sizeof(double)); // back plane

}

void initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength, int iProc, int jProc, int kProc, int rank){


	int x, y, z, i;
	int xlen2 = xlength[0]+2;
	int ylen2 = xlength[1]+2;
	int zlen2 = xlength[2]+2;

	// Global domain, CPU order: (iProc = x_axis, jProc = y_axis, kProc = z_axis)

	// Example (CUBE iProc*jProc*kProc = 4 * 3 * 2: -> for columns, three rows, two slices
	// first x,z plane
	// 8 9 10 11
	// 4 5 6 7
	// 1 2 3 4

	// second x,z plane
	// 20 21 22 23
	// 16 17 18 19
	// 12 13 14 15


	// Boundary init.
	// Left boundary. If true, then we pick the left plane A=A(x=0, y, z) of this process and define it as no-slip.
	if (rank % iProc == 0){
		for (z = 0; z < zlen2; z++) {
			for (y = 0; y < ylen2; y++) {
				flagField[y * xlen2 + z * xlen2*ylen2] = NO_SLIP;
			}
		}
	} else {

		for (z = 0; z < zlen2; z++) {
			for (y = 0; y < ylen2; y++) {
				flagField[y * xlen2 + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
			}
		}
	}


	// Right boundary. If true, then we pick the right plane A=A(x=xlen, y, z) of this process and define it as no-slip.
	if (rank % iProc == iProc - 1){
		for (z = 0; z < zlen2; z++) {
			for (y = 0; y < ylen2; y++) {
				flagField[(xlen2 - 1) + y * xlen2 + z * xlen2*ylen2] = NO_SLIP;
			}
		}
	} else {
			for (z = 0; z < zlen2; z++) {
				for (y = 0; y < ylen2; y++) {
				flagField[(xlen2 - 1) + y * xlen2 + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
			}
		}
}


	// Front boundary. If true, then we pick the front plane A=A(x,y=0,z) of this process and define it as no-slip.
	if (rank <= iProc*kProc - 1){
		for (z = 0; z < zlen2; z++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + z * xlen2*ylen2] = NO_SLIP;
			}
		}
	} else {
			for (z = 0; z < zlen2; z++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
		}

	// Back boundary. If true, then we pick the back plane A=A(x,y=ylen,z) of this process and define it as no-slip.
	if (rank >= iProc*(jProc-1)*kProc){
		for (z = 0; z < zlen2; z++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + (ylen2 - 1)*xlen2 + z * xlen2*ylen2] = NO_SLIP;
			}
		}
	} else {

			for (z = 0; z < zlen2; z++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + (ylen2 - 1)*xlen2 + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
	}

	// Top boundary. If true, then we pick the top plane A=A(x,y,z=zlen2) of this process and define it as moving-boundary(!).
	if (rank % iProc*jProc <= iProc*kProc - 1 && rank % iProc*jProc >= iProc*(kProc - 1)) {
		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2 + (zlen2 - 1) * xlen2*ylen2] = MOVING_WALL;
			}
		}
	} else {

			for (y = 0; y < ylen2; y++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + y*xlen2 + (zlen2 - 1) * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
		}

	// Bottom boundary. If true, then we pick the bottom plane A=A(x,y,z=0) of this process and define it as no-slip.
	if (rank % iProc*jProc < iProc - 1) {
		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2] = NO_SLIP;
			}
		}
	} else {

		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2] = PARALLEL_BOUNDARY;
			}
		}
	}


	/* stream & collide Fields initialization. */
	for (z = 0; z < zlen2; ++z){
		for (y = 0; y < ylen2; ++y){
			for (x = 0; x < xlen2; ++x){
				for (i = 0; i < Q_NUMBER; ++i){
					streamField[Q_NUMBER * (x + y * xlen2 + z * xlen2*ylen2) + i] = LATTICEWEIGHTS[i];
					collideField[Q_NUMBER * (x + y * xlen2 + z * xlen2*ylen2) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}

	// Fluid init (inner part of flagField).
	for (z = 1; z <= zlen2-2; ++z) {
		for (y = 1; y <= ylen2-2; ++y) {
			for (x= 1; x <= xlen2-2; ++x) {
				flagField[x + y * xlen2 + z * xlen2*ylen2] = FLUID;
			}
		}
	}

}
