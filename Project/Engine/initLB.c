#include "initLB.h"
#include "mpi.h"
#include "helper.h"
#include "LBDefinitions.h"

int readParameters(int * const xlength, double * const tau, double * const velocityWall, int * const timesteps, int * const timestepsPerPlotting, double_3d * const inflow, double * const pressure_in, int * const iProc, int * const jProc, int * const kProc,  int argc,  char *  *  argv){
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  

		read_int( szFileName, "xlength", xlength );
		read_int( szFileName, "ylength", xlength + 1 );
		read_int( szFileName, "zlength", xlength + 2 );
		if(xlength[0] < 2 || xlength[1] < 2 || xlength[2] < 2){
			printf("Dimensions xyzlength must be > 1, please fix the geometry!\n");
			return 2;
		}

		READ_DOUBLE( szFileName, *tau );
		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );

		// who needs loops? :)
		read_double( szFileName, "inflow_0_v_x", &inflow[0].x);
		read_double( szFileName, "inflow_0_v_y", &inflow[0].y);
		read_double( szFileName, "inflow_0_v_z", &inflow[0].z);

		read_double( szFileName, "inflow_1_v_x", &inflow[1].x);
		read_double( szFileName, "inflow_1_v_y", &inflow[1].y);
		read_double( szFileName, "inflow_1_v_z", &inflow[1].z);

		read_double( szFileName, "inflow_2_v_x", &inflow[2].x);
		read_double( szFileName, "inflow_2_v_y", &inflow[2].y);
		read_double( szFileName, "inflow_2_v_z", &inflow[2].z);

		read_double( szFileName, "inflow_3_v_x", &inflow[3].x);
		read_double( szFileName, "inflow_3_v_y", &inflow[3].y);
		read_double( szFileName, "inflow_3_v_z", &inflow[3].z);

		read_double( szFileName, "inflow_4_v_x", &inflow[4].x);
		read_double( szFileName, "inflow_4_v_y", &inflow[4].y);
		read_double( szFileName, "inflow_4_v_z", &inflow[4].z);

		read_double( szFileName, "inflow_5_v_x", &inflow[5].x);
		read_double( szFileName, "inflow_5_v_y", &inflow[5].y);
		read_double( szFileName, "inflow_5_v_z", &inflow[5].z);

		read_double( szFileName, "pressure_in_0_d", &pressure_in[0]);
		read_double( szFileName, "pressure_in_1_d", &pressure_in[1]);
		read_double( szFileName, "pressure_in_2_d", &pressure_in[2]);
		read_double( szFileName, "pressure_in_3_d", &pressure_in[3]);
		read_double( szFileName, "pressure_in_4_d", &pressure_in[4]);
		read_double( szFileName, "pressure_in_5_d", &pressure_in[5]);

		read_double( szFileName, "velocityWall1", velocityWall );
		read_double( szFileName, "velocityWall2", velocityWall + 1 );
		read_double( szFileName, "velocityWall3", velocityWall + 2 );
		READ_INT( szFileName, *iProc );
		READ_INT( szFileName, *jProc );
		READ_INT( szFileName, *kProc );
		//for( int i = 0; i < INFLOW_COUNT; ++i){
		//	printf("inflow[%d].x = %f\n", i, inflow[i].x);
		//}
	}
	return 0;
}

void initialiseBuffers(double **sendBuffer, double **readBuffer, const int * const cpuDomain, int * sizeBuffer, const int * const neighbor){

	int xlen2 = cpuDomain[0] + 2;
	int ylen2 = cpuDomain[1] + 2;
	int zlen2 = cpuDomain[2] + 2;

	int domain = 5; // because we have 5 possible directions to be extracted to buffer

	// We should substitute the sizes in malloc, but we don't have time now.
	sizeBuffer[0] = ylen2 * zlen2 * domain;
	sizeBuffer[1] = ylen2 * zlen2 * domain; // =sizeBuffer[0]
	sizeBuffer[2] = xlen2 * ylen2 * domain;
	sizeBuffer[3] = xlen2 * ylen2 * domain; // =sizeBuffer[2]
	sizeBuffer[4] = xlen2 * zlen2 * domain;
	sizeBuffer[5] = xlen2 * zlen2 * domain; // =sizeBuffer[4]

	// We initilise 6 different buffers.
	// sendBuffer planes[0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]

	if(neighbor[0] != MPI_PROC_NULL) sendBuffer[0] = (double *) malloc(sizeBuffer[0] * sizeof(double)); // left plane
	if(neighbor[1] != MPI_PROC_NULL) sendBuffer[1] = (double *) malloc(sizeBuffer[1] * sizeof(double)); // right plane
	if(neighbor[2] != MPI_PROC_NULL) sendBuffer[2] = (double *) malloc(sizeBuffer[2] * sizeof(double)); // top plane
	if(neighbor[3] != MPI_PROC_NULL) sendBuffer[3] = (double *) malloc(sizeBuffer[3] * sizeof(double)); // bottom plane
	if(neighbor[4] != MPI_PROC_NULL) sendBuffer[4] = (double *) malloc(sizeBuffer[4] * sizeof(double)); // front plane
	if(neighbor[5] != MPI_PROC_NULL) sendBuffer[5] = (double *) malloc(sizeBuffer[5] * sizeof(double)); // back plane

	// readBuffer planes[0:right sendBuffer, 1:left sendBuffer, 2:bottom sendBuffer, 3:top sendBuffer, 4:back sendBuffer, 5:front sendBuffer]
	if(neighbor[0] != MPI_PROC_NULL) readBuffer[0] = (double *) malloc(sizeBuffer[0] * sizeof(double)); // left plane
	if(neighbor[1] != MPI_PROC_NULL) readBuffer[1] = (double *) malloc(sizeBuffer[1] * sizeof(double)); // right plane
	if(neighbor[2] != MPI_PROC_NULL) readBuffer[2] = (double *) malloc(sizeBuffer[2] * sizeof(double)); // top plane
	if(neighbor[3] != MPI_PROC_NULL) readBuffer[3] = (double *) malloc(sizeBuffer[3] * sizeof(double)); // bottom plane
	if(neighbor[4] != MPI_PROC_NULL) readBuffer[4] = (double *) malloc(sizeBuffer[4] * sizeof(double)); // front plane
	if(neighbor[5] != MPI_PROC_NULL) readBuffer[5] = (double *) malloc(sizeBuffer[5] * sizeof(double)); // back plane

}

void initialiseFields(double *collideField, double *streamField, int *flagField, const int * const cpuDomain, const int iProc, const int jProc, const int kProc, const int rank, int * const neighbor){
	// local domain is altogether Dlength + 2, where the first and last cells are either buffer(parallel boundary) or global domain(no slip)

	int x, y, z, i;
	int xlen2 = cpuDomain[0] + 2;
	int ylen2 = cpuDomain[1] + 2;
	int zlen2 = cpuDomain[2] + 2;
	//int xyzlen2 = xlen2 * ylen2 * zlen2;

	

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
				flagField[y * xlen2 + z * xlen2 * ylen2] = NO_SLIP;
			}
		}
		neighbor[DIR_L] = MPI_PROC_NULL;
	} else {

		for (z = 0; z < zlen2; z++) {
			for (y = 0; y < ylen2; y++) {
				flagField[y * xlen2 + z * xlen2 * ylen2] = PARALLEL_BOUNDARY;
			}
		}
		neighbor[DIR_L] = rank - 1;
	}

	// Right boundary. If true, then we pick the right plane A=A(x=xlen, y, z) of this process and define it as no-slip.
	if (rank % iProc == iProc - 1){
		for (z = 0; z < zlen2; z++) {
			for (y = 0; y < ylen2; y++) {
				flagField[(xlen2 - 1) + y * xlen2 + z * xlen2*ylen2] = NO_SLIP;
			}
		}
		neighbor[DIR_R] = MPI_PROC_NULL;
	} else {
			for (z = 0; z < zlen2; z++) {
				for (y = 0; y < ylen2; y++) {
				flagField[(xlen2 - 1) + y * xlen2 + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
			}
		}
		neighbor[DIR_R] = rank + 1;
}

	// Front boundary. If true, then we pick the front plane A=A(x,y=0,z) of this process and define it as no-slip.
	if (rank / iProc %jProc == 0){
		for (z = 0; z < zlen2; z++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + z * xlen2*ylen2] = NO_SLIP;
			}
		}
		neighbor[DIR_F] = MPI_PROC_NULL;
	} else {
			for (z = 0; z < zlen2; z++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
			neighbor[DIR_F] = rank - iProc;
		}

	// Back boundary. If true, then we pick the back plane A=A(x,y=ylen,z) of this process and define it as no-slip.
	if (rank / iProc % jProc == jProc -1){
		for (z = 0; z < zlen2; z++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + (ylen2 - 1)*xlen2 + z * xlen2*ylen2] = NO_SLIP;
			}
		}
		neighbor[DIR_B] = MPI_PROC_NULL;
	} else {
			for (z = 0; z < zlen2; z++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + (ylen2 - 1)*xlen2 + z * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
			neighbor[DIR_B] = rank + iProc;
	}

	// Top boundary. If true, then we pick the top plane A=A(x,y,z=zlen2) of this process and define it as moving-boundary(!).
	if (rank / iProc / jProc == kProc - 1) {
		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2 + (zlen2 - 1) * xlen2*ylen2] = MOVING_WALL;
			}
		}
		neighbor[DIR_T] = MPI_PROC_NULL;
	} else {

			for (y = 0; y < ylen2; y++) {
				for (x = 0; x < xlen2; x++) {
					flagField[x + y*xlen2 + (zlen2 - 1) * xlen2*ylen2] = PARALLEL_BOUNDARY;
				}
			}
			neighbor[DIR_T] = rank + iProc * jProc;
		}

	// Bottom boundary. If true, then we pick the bottom plane A=A(x,y,z=0) of this process and define it as no-slip.
	if (rank / iProc / jProc == 0) {
		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2] = NO_SLIP;
			}
		}
		neighbor[DIR_D] = MPI_PROC_NULL;
	} else {

		for (y = 0; y < ylen2; y++) {
			for (x = 0; x < xlen2; x++) {
				flagField[x + y*xlen2] = PARALLEL_BOUNDARY;
			}
		}
		neighbor[DIR_D] = rank - iProc * jProc;
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
				flagField[x + y * xlen2 + z * xlen2 * ylen2] = FLUID;
			}
		}
	}

	// print flagfield initialization for debug
	printf("domain %d %d %d\n",xlen2,ylen2,zlen2);
	for(z = zlen2 - 1;z >= 0; --z){
		for(y = ylen2 - 1;y >= 0; --y){
			for(x = 0;x < xlen2; ++x){
				printf("%d ",flagField[x + y * xlen2 + z * xlen2 * ylen2]);
			}
		printf("plane %d rank _%d\n",z, rank);
		}
		printf("\n");
	}
	//printf("exit initLB \n");

}
