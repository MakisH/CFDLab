#include "initLB.h"

///// we might need some of these for the initialiseFields() function - for the pgm part
#include <ctype.h>
#include <errno.h>
#include <stdio.h>

int readParameters(char		*problem,
									 int		*xlength,
									 double *tau,
									 int		*timesteps,
									 int		*timestepsPerPlotting,
									 double *velocityIn,
									 double *densityIn,
									 double *densityRef,
									 double *velocityWall,
									 int		*initxyzXYZ,
									 int		argc,
									 char		*argv[]){
	// we need to read the problem char variable and act accordingly to read the other variables(switch....) 
	// ... or read all and think later :)
	if ( argc != 2 ) {
		printf("Usage: ./sim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  

		READ_STRING(szFileName, problem);

		read_int( szFileName, "xlength", &xlength[0] );
		read_int( szFileName, "ylength", &xlength[1] );
		read_int( szFileName, "zlength", &xlength[2] );

		READ_DOUBLE( szFileName, *tau );

		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );

		read_double( szFileName, "velocityInX", &velocityIn[0] );
		read_double( szFileName, "velocityInY", &velocityIn[1] );
		read_double( szFileName, "velocityInZ", &velocityIn[2] );

		READ_DOUBLE(szFileName, *densityIn);
		READ_DOUBLE(szFileName, *densityRef);

		read_double( szFileName, "velocityWall1", &velocityWall[0] );
		read_double( szFileName, "velocityWall2", &velocityWall[1] );
		read_double( szFileName, "velocityWall3", &velocityWall[2] );
		
		// careful here! // domain starts from bottom left
		read_int( szFileName, "wallLeft",				&initxyzXYZ[2] ); // z-
		read_int( szFileName, "wallRight",			&initxyzXYZ[5] ); // z+
		read_int( szFileName, "wallTop",				&initxyzXYZ[3] ); // x-
		read_int( szFileName, "wallBottom",			&initxyzXYZ[0] ); // x+
		read_int( szFileName, "wallBackground",	&initxyzXYZ[4] ); // y+ 
		read_int( szFileName, "wallForeground", &initxyzXYZ[1] ); // y-

	}
	return 0;
}

int initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength, char *problem, int *initxyzXYZ){
	// TO DO - should be fixed - implement exit if scenario is wrong
	// TO DO - should be fixed - think about the direction of initiliasiation - all scenarios should start numbering from top left corner
	int x, y, z, i;
	int xlen2 = xlength[0] + 2;
	int ylen2 = xlength[1] + 2;
	int zlen2 = xlength[2] + 2;
	int xyz_field_domain = xlen2 * ylen2 * zlen2 * Q_NUMBER;
	int xylen2 = xlen2 * ylen2;
//	unsigned int X_min, Y_min, Z_min, X_max, Y_max, Z_max;

	if (!strcmp(problem,"lbm_tilted_plate"))
	{
		printf("sc1 %d\n", strcmp(problem,"lbm_tilted_plate"));
	// need to scale input file according to dimensions, if they're different we throw an error
	// this is a repetitions of the pgm_read() function,
	//		but without using the tricky array functions and scaling to support different x-,y-,z- dimensions.

		FILE *input = NULL;
		char line[1024];
		int levels;
		int xsize, ysize;
		int i1, j1;
//		int **pic = NULL;
		const char *filename = "./configs/lbm_tilted_plate.pgm";//lbm_tilted_plate.vtk";
		int scale_x;
		int scale_y;
		if ((input=fopen(filename,"rb"))==0)
		{
		   char szBuff[80];
		       sprintf( szBuff, "Can not read file %s !!!", filename );
		       ERROR( szBuff );
					 return 1;
		}

		/* check for the right "magic number" */
		if ( fread(line,1,3,input)!=3 )
		{
		        fclose(input);
		        ERROR("Error Wrong Magic field!");
						return 1;
		}

		/* skip the comments */
		do
		if(fgets(line,sizeof line,input));
		while(*line=='#');

		/* read the width and height */
		sscanf(line,"%d %d\n",&xsize,&ysize);

		//printf("Image size: %d x %d\n", xsize,ysize);

		/* read # of gray levels */
		if(fgets(line,sizeof line,input));
		sscanf(line,"%d\n",&levels);

		/* allocate memory for image */ // we don't need this, we have our own space !
		//pic = imatrix(0,xsize+2,0,ysize+2);
		//printf("Image initialised...\n");

		/* read pixel row by row */
		if( xsize == 0 ||  ysize == 0)
		{
			printf("Dimensions are invalid(0)!\n");
			return 1;
		}
		else if(xlength[2] % xsize != 0 || xlength[0] % ysize != 0)
		{
			printf("xlen2 %d\nxsize %d\nxlen0 %d\nysize %d \n",xlength[2],xsize,xlength[0],ysize);
			printf("Dimensions and xlength mismatch!\n");
			return 1;
		}
		// need independent dim-s, so xlength has to be multiple of corresponding(!) xsize and ysize
		scale_x = xlength[2] / xsize; // our implementation uses different naming
		scale_y = xlength[0] / ysize; // - xsize is xlength[2] and ysize is xlength[0]  
		for(j1 = 1; j1 <= ysize; ++j1){
			for (i1 = 1; i1 <= xsize; ++i1){
				int byte;
				if(fscanf(input, "%d", &byte)); // can someone explain the if-statement ?

				if (byte==EOF)
				{
					fclose(input);
					ERROR("read failed");
					return 1;
				}
				else
				{
			//		if(i1 < 150 && i1 > 100) printf("%d ",byte);
// every read initialises a "cube" from the 3D flagField, scaled according to x-,y-,z- length
					// could be done in different ways, e.g. save the input values and then traverse the flagField and initialize it sequentially
					for (z = (i1 - 1) * scale_x + 1; z <= i1 * scale_x + scale_x; ++z) {
						for (y = 1; y <= xlength[1]; ++y) { // height is always traversed fully
							for (x = (j1 - 1) * scale_y + 1; x <= j1 * scale_y + scale_y; ++x) {
								flagField[x + y * ylen2 + z * xylen2] = byte;
								//printf("%d ",x + y * ylen2 + z * xylen2);
							}
						}
					}
				}
			}
		}
		/* close file */
		fclose(input);

	}
	else if (!strcmp(problem,"lbm_plane_shear_flow"))
	{
		printf("sc2 %d\n", strcmp(problem,"lbm_plane_shear_flow"));
		// Fluid init (inner part of flagField).
		for (z = 1; z <= xlength[2]; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]; ++x) {
					flagField[x + y * ylen2 + z * xylen2] = FLUID; // we could use a single pointer instead of calculating the whole offset multiple times ... requires effort vs small performance since only 1 call
					printf("%d ",x + y * ylen2 + z * xylen2);
				}
			}
		}

	}
	else if (!strcmp(problem,"lbm_flow_over_step"))
	{
		printf("sc3 %d\n", strcmp(problem,"lbm_flow_over_step"));
		// Fluid init (inner part of flagField).
		for (z = 1; z <= xlength[2]; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]; ++x) {
					flagField[x + y * ylen2 + z * xylen2] = FLUID; // we could use a single pointer instead of calculating the whole offset multiple times ... requires effort vs small performance since only 1 call
				//	printf("%d ",x + y * ylen2 + z * xylen2);
				}
			}
		}

		// add "step" initialization - "cube" has z = x = xlen/2, y = ylen ... x-> dimension starts from top left!
		for (z = 1; z <= xlength[2]/2; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]/2; ++x) { //cube is in the 2nd half of the x dimension
					flagField[x + y * ylen2 + z * xylen2] = NO_SLIP;
					//printf("%d ",x + y * ylen2 + z * xylen2);
				}
			}
		}
		printf("init done!\n");
	}
	else
	{
				printf("sc4 %d\n", strcmp(problem,"lbm_flow_over_step"));

		printf("\n\nUnknown / Misspelled scenario name!\n\n");
		printf("Known scenarios are:\n\n");
		printf("A tilted plane\n");
		printf("Plane shear flow\n");
		printf("Flow over a step\n");
		return 1;
		// we could implement a loop to input the correct name, or 0 to exit

	}

	// identical for all scenarios => called once outside of if-statement.

	// check for forbidden boundary cells
	int counter;
	for (z = 1; z <= xlength[2]; ++z) {
    for (y = 1; y <= xlength[1]; ++y) {
      for (x = 1; x <= xlength[0]; ++x) {
        counter = 0;
        if (flagField[x + y * ylen2 + z * xylen2] != FLUID) {
          if (flagField[(x-1) + y * ylen2 + z * xylen2] == FLUID) counter++;
          if (flagField[(x+1) + y * ylen2 + z * xylen2] == FLUID) counter++;
          if (flagField[x + (y-1) * ylen2 + z * xylen2] == FLUID) counter++;
          if (flagField[x + (y+1) * ylen2 + z * xylen2] == FLUID) counter++;            
          if (flagField[x + y * ylen2 + (z-1) * xylen2] == FLUID) counter++;          
          if (flagField[x + y * ylen2 + (z+1) * xylen2] == FLUID) counter++;          
        }
        if (counter > 2) {
          printf("The domain contains forbidden boundary cells, that is, boundary cells with more than 2 fluid neighbors.");
          return 1;
        }
      }
    }
	}
	
	// stream & collide Fields initialization.
	for(x = 0; x < xyz_field_domain; x += Q_NUMBER){
		for (i = 0; i < Q_NUMBER; ++i){ // we need this array to loop over Latticeweights ... if they're hardcoded, we can skip it(less memory use) should not matter that much, since it's done only once
			streamField[i+x] = LATTICEWEIGHTS[i];
			collideField[i+x] = LATTICEWEIGHTS[i];
		}
	}
		/** flagField init: 
		* 0 - FLUID
		* 1 - NO_SLIP
		* 2 - MOVING_WALL
		* 3 - FREE-SLIP
		* 4 - INFLOW
		* 5 - OUTFLOW
		* 6 - PRESSURE_IN
	**/
	/* Boundary initialization: using input parameters */

		for (z = 0; z < zlen2; ++z){ // treat x, y and z as pure iterators, the dimension depends on where we have 0 or xlength+1 and iteration happens.
			for (y = 0; y < ylen2; ++y){
				for (x = 0; x < xlen2; ++x){
					// add all other walls in the same loop, simply switch the indices
					// - x and y are iterators, use x with xlensq for memory locality
					// z is "0" or "xlength + 1"
					// the only question is whether this is OK with the cache memory
					// in fact we can move these iterations above at line 19(before the third loop from 0 to xlength + 1)
					flagField[	  y * ylen2 + z * xylen2	] = initxyzXYZ[0];						// x- dimension	
										// printf("%d ",y * ylen2 + z * xylen2);
					flagField[x							+ z * xylen2	] = initxyzXYZ[1];						// y- dimension
					flagField[x + y * ylen2								] = initxyzXYZ[2];						// z- dimension

					flagField[xlen2 - 1 +	y * ylen2 + z * xylen2		] = initxyzXYZ[3];	// x+ dimension
					flagField[x + (ylen2 - 1) * ylen2 + z * xylen2	] = initxyzXYZ[4];	// y+ dimension
					flagField[x + y * ylen2 + (zlen2 - 1) * xylen2	] = initxyzXYZ[5];	// z+ dimension
				}
			}
		}
		printf("init done done\n");
		return 0;
}