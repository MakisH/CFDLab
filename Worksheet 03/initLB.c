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
			int xsize, ysize;

	if (!strcmp(problem,"lbm_tilted_plate"))
	{
		printf("sc1 %d\n", strcmp(problem,"lbm_tilted_plate"));
	// need to scale input file according to dimensions, if they're different we throw an error
	// this is a repetitions of the pgm_read() function,
	//		but without using the tricky array functions and scaling to support different x-,y-,z- dimensions.

		FILE *input = NULL;
		char line[1024];
		int levels;
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
		else if(xlength[0] % xsize != 0 || xlength[2] % ysize != 0)
		{
			printf("xlen2 %d\nxsize %d\nxlen0 %d\nysize %d \n",xlength[2],xsize,xlength[0],ysize);
			printf("Dimensions and xlength mismatch!\n");
			return 1;
		}
		// need independent dim-s, so xlength has to be multiple of corresponding(!) xsize and ysize
		scale_x = xlength[0] / xsize; // our implementation uses different naming
		scale_y = xlength[2] / ysize; // - xsize is xlength[2] and ysize is xlength[0]  
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
					//printf("%d",byte);
				//	if(byte == 1)printf("dim x y z %d %d %d .... x y %d %d\n",xlength[0],xlength[1],xlength[2],i1,j1);
// every read initialises a "cube" from the 3D flagField, scaled according to x-,y-,z- length
					// could be done in different ways, e.g. save the input values and then traverse the flagField and initialize it sequentially

					for (z = (j1 - 1) * scale_y + 1; z <= j1 * scale_y + scale_y; ++z) {
						for (y = 1; y <= xlength[1]; ++y) { // height is always traversed fully
							for (x = (i1 - 1) * scale_x + 1; x <= i1 * scale_x + scale_x; ++x) {
								flagField[x + y * xlen2 + z * xylen2] = byte;
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
					//printf("%d ",x + y * ylen2 + z * xylen2);
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

		//for (z = 0; z < zlen2; ++z){ // treat x, y and z as pure iterators, the dimension depends on where we have 0 or xlength+1 and iteration happens.
		//	for (y = 0; y < ylen2; ++y){
		//		for (x = 0; x < xlen2; ++x){
		//			// add all other walls in the same loop, simply switch the indices
		//			// - x and y are iterators, use x with xlensq for memory locality
		//			// z is "0" or "xlength + 1"
		//			// the only question is whether this is OK with the cache memory
		//			// in fact we can move these iterations above at line 19(before the third loop from 0 to xlength + 1)
		//			flagField[	  y * xlen2 + z * xylen2	] = initxyzXYZ[0];						// x- dimension	
		//								
		//			flagField[x							+ z * xylen2	] = initxyzXYZ[1];						// y- dimension
		//			flagField[x + y * xlen2								] = initxyzXYZ[2];						// z- dimension
		//			 printf("%d ",x + y * xlen2	);
		//			flagField[xlen2 - 1 +	y * xlen2 + z * xylen2		] = initxyzXYZ[3];	// x+ dimension
		//			flagField[x + (ylen2 - 1) * xlen2 + z * xylen2	] = initxyzXYZ[4];	// y+ dimension
		//			flagField[x + y * xlen2 + (zlen2 - 1) * xylen2	] = initxyzXYZ[5];	// z+ dimension
		//		}
		//	}
		//}
	for(int y = 0; y < ylen2; ++y){
		for(x = 0; x < xlen2; ++x){
					flagField[x + y * xlen2													] = initxyzXYZ[2];	// z- dimension
					flagField[x + y * xlen2 + (zlen2 - 1) * xylen2	] = initxyzXYZ[5];	// z+ dimension
		}
	}
	for(int z = 0; z < zlen2; ++z){
		for(x = 0; x < xlen2; ++x){
					flagField[x												+ z * xylen2	] = initxyzXYZ[1];	// y- dimension
					flagField[x + (ylen2 - 1) * xlen2 + z * xylen2	] = initxyzXYZ[4];	// y+ dimension
		}
	}
	for(int z = 0; z < zlen2; ++z){
		for(y = 0; y < ylen2; ++y){
					flagField[						y * xlen2 + z * xylen2		] = initxyzXYZ[0];	// x- dimension	
					flagField[xlen2 - 1 +	y * xlen2 + z * xylen2		] = initxyzXYZ[3];	// x+ dimension
		}
	}


	// Check for forbidden boundary cells (no_slip or free_slip)
	// fluidNeighbors: counter of neighbor cells that are fluid
	// Although we have 3D problems, we assume that boundaries have no fluid neighbors
	// in the y-direction (foreground-background).
	// This check has a big optimization potential but it is executed only once at the start.
	int fluidNeighbors;
	for (z = 1; z <= xlength[2]; ++z) {
    for (y = 1; y <= xlength[1]; ++y) {
      for (x = 1; x <= xlength[0]; ++x) {
        
        // How many fluid neighbors does the current cell have?
        fluidNeighbors = 0;
        if (flagField[x + y * ylen2 + z * xylen2] != FLUID) {
          if (flagField[(x-1) + y * ylen2 + z * xylen2] == FLUID) ++fluidNeighbors;
          if (flagField[(x+1) + y * ylen2 + z * xylen2] == FLUID) ++fluidNeighbors;
          if (flagField[x + y * ylen2 + (z-1) * xylen2] == FLUID) ++fluidNeighbors;
          if (flagField[x + y * ylen2 + (z+1) * xylen2] == FLUID) ++fluidNeighbors;
        }

        // Depending on the number of neighbors, is it an inner cell, an edge, a valid corner or something forbidden?
        if (fluidNeighbors > 2) {
          printf("The domain contains at least one boundary cell with more than 2 fluid neighbors.\n");
          printf("Coordinates: x = %d, y = %d, z = %d \n", x, y, z);
          return 1;
        } else if (fluidNeighbors == 0) continue; // it is an inner obstacle cell
        else if (fluidNeighbors == 1) continue; // it is an edge
        else if (fluidNeighbors == 2) {
            // Check if the two fluid neighbors share a corner (check all the possibilities, remember we don't have general 3D cases)
            if (flagField[(x-1) + y * ylen2 + z * xylen2] == FLUID && flagField[x + y * ylen2 + (z+1) * xylen2] == FLUID) continue;
            if (flagField[(x+1) + y * ylen2 + z * xylen2] == FLUID && flagField[x + y * ylen2 + (z+1) * xylen2] == FLUID) continue;
            if (flagField[(x+1) + y * ylen2 + z * xylen2] == FLUID && flagField[x + y * ylen2 + (z-1) * xylen2] == FLUID) continue;
            if (flagField[(x-1) + y * ylen2 + z * xylen2] == FLUID && flagField[x + y * ylen2 + (z-1) * xylen2] == FLUID) continue;
            // The execution didn't move to the next iteration, so this is not a valid corner!
            printf("The domain contains at least one boundary cell with maximum 2 fluid neighbors that is forbidden.\n");
            printf("Coordinates: x = %d, y = %d, z = %d \n", x, y, z);
            return 1;
        }
      }
    }
	}
	//printf("init done done\n");
	//printf("domain size is %d %d %d with xsize and ysize %d %d\n",xlength[0], xlength[1], xlength[2], xsize, ysize);
	//for(int i = xlength[0] + 1; i >= 0 ; --i){
	//	for(int j = 0; j < xlength[2] + 2; ++j){
	//		printf("%d",flagField[i + xlen2 + xylen2 * j]);
	//	}
	//	printf("\n");
	//}
		return 0;
}