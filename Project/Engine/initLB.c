#include "initLB.h"
#include "mpi.h"
#include "helper.h"
#include "LBDefinitions.h"

int read_assign_PGM (int *flagField, char *fileName, int *cpuDomain)
{
	/*
		One CPU, one flagField, streamField, collideField. Martin gives pointers to proper sections of this arrays
		We initilize flagField with values of PGM file.
		We (just in case, if we'll need that somewhere) also copy the dimensions of PGM from PGM file.
	*/

	/* READ THE PGM FILE, STORE IN FLAGFIELD */

	FILE *file = NULL;
	char line[1024];
	int xsize, ysize;

	int xlen2 = cpuDomain[0];
	int ylen2 = cpuDomain[1];
	int xylen2 = xlen2 * ylen2;

	if ((file=fopen(fileName,"rb"))==0){
		char szBuff[80];
				sprintf( szBuff, "Can not read file %s !!!", fileName );
				//ERROR( szBuff );
				return 1;
	}

	/* check for the right "magic number" */
	if ( fread(line,1,3,file)!=3 )
	{
					fclose(file);
					//ERROR("Error Wrong Magic field!");
					return 1;
	}

	do
	if(fgets(line,sizeof line,file));
	while(*line=='#');

	/* read the width and height */
	sscanf(line,"%d %d\n",&xsize,&ysize);

	/* Rows are x dimension, columns are y dimension */
	int z = 0;
	for (int y = 0; y < ysize; y++){
		for (int x = 0; x < xsize; x++){
			int byte;
			if(fscanf(file, "%d", &byte));
			printf("%d ", byte);

			flagField[x + y * xlen2 + z * xylen2] = byte;

			if (byte==EOF){
				fclose(file);
				//ERROR("read failed");
				return 1;
			}
		}
	printf("\n");
	}
	fclose(file);

	return 1;
}


int readParameters(int * const xlength, double * const tau, double * const velocityWall, int * const timesteps, int * const timestepsPerPlotting, double_3d * const inflow, double * const pressure_in, double * const ref_density, int * const iProc, int * const jProc, int * const kProc,  int argc,  char *  *  argv){
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

		read_double( szFileName, "densityRef", ref_density);

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


	// ALERT!! -> How is cpuDomain defined?? cpuDomain+2 might be too much!
	int x, y, z, i;
	int xlen2 = cpuDomain[0] + 2;
	int ylen2 = cpuDomain[1] + 2;
	int zlen2 = cpuDomain[2] + 2;
	//int xyzlen2 = xlen2 * ylen2 * zlen2;

	int xylen2 = xlen2 * ylen2;
	// Now apply free slip to every x,y slice on z=1 to the end.
	for (int z = 1; z < cpuDomain[2]; z++){
		for (int y = 0; y < cpuDomain[1]; y++){
			for (int x = 0; x < cpuDomain[0]; x++){
				flagField[x + y * xlen2 + z * xylen2] = FREE_SLIP;
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
