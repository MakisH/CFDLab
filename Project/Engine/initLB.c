#include "initLB.h"
#include "mpi.h"
#include "helper.h"
#include "LBDefinitions.h"

int read_assign_PGM (int *flagField, char *fileName, const int * const cpuDomain)
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
	int z = 1; // The middle layer is PGM, the 0-th and 2-nd layer are free slip!
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


int readParameters(double * const tau, double * const velocityWall, int * const timesteps, int * const timestepsPerPlotting, double_3d * const inflow, double * const pressure_in, double * const ref_density, int argc,  char *  *  argv){
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];

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

		//for( int i = 0; i < INFLOW_COUNT; ++i){
		//	printf("inflow[%d].x = %f\n", i, inflow[i].x);
		//}
	}
	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, const int * const cpuDomain){
	// local domain is altogether Dlength + 2, where the first and last cells are either buffer(parallel boundary) or global domain(no slip)

	// ALERT!! -> How is cpuDomain defined?? cpuDomain+2 might be too much!
	int x, y, z, i;
	int xlen2 = cpuDomain[0];
	int ylen2 = cpuDomain[1];
	int zlen2 = cpuDomain[2];
	//int xyzlen2 = xlen2 * ylen2 * zlen2;

	int xylen2 = xlen2 * ylen2;

	/* stream & collide Fields initialization. */
	for (z = 0; z < zlen2; ++z){
		for (y = 0; y < ylen2; ++y){
			for (x = 0; x < xlen2; ++x){
				for (i = 0; i < Q_NUMBER; ++i){
					streamField[Q_NUMBER * (x + y * xlen2 + z * xylen2) + i] = LATTICEWEIGHTS[i];
					collideField[Q_NUMBER * (x + y * xlen2 + z * xylen2) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}

	// init free slip.
	z = 0; // The 0-th plane.
	for (y = 0; y<ylen2; y++){
		for (x = 0; x<xlen2, x++){
			flagField[x + y*xlen2 + z*xylen2] = FREE_SLIP;
		}
	}

	// The 1-st plane has the PGM geometry inside.

	z = 2; // The 2-nd plane.
	for (y = 0; y<ylen2; y++){
		for (x = 0; x<xlen2, x++){
			flagField[x + y*xlen2 + z*xylen2] = FREE_SLIP;
		}
	}


	// print flagfield initialization for debug
	printf("domain %d %d %d\n",xlen2,ylen2,zlen2);
	//for(z = zlen2 - 1; z >= 0; --z){
	//	for(y = ylen2 - 1; y >= 0; --y){
	//		for(x = 0; x < xlen2; ++x){
	//			printf("%d ",flagField[x + y * xlen2 + z * xlen2 * ylen2]);
	//		}
	//	printf("plane %d rank _%d\n",z, rank);
	//	}
	//	printf("\n");
	//}
	//printf("exit initLB \n");

}
