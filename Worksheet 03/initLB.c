#include "initLB.h"

///// we might need some of these for the initialiseFields() function - for the pgm part
//#include <ctype.h>
//#include <errno.h>
//#include <stdio.h>

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
	// we need to read the problem char variable and act accordingly to read the other variables(switch....)
	if ( argc != 2 ) {
		printf("Usage: ./lbsim input_file");
		return 1;
	}
	else {
		const char *szFileName = NULL;
		szFileName = argv[1];  
		READ_INT( szFileName, xlength[0] );
		READ_INT( szFileName, xlength[1] );
		READ_INT( szFileName, xlength[2] );
		READ_DOUBLE( szFileName, *tau );

		READ_INT( szFileName, *timesteps );
		READ_INT( szFileName, *timestepsPerPlotting );
		read_double( szFileName, "velocityWall1", &velocityWall[0] );
		read_double( szFileName, "velocityWall2", &velocityWall[1] );
		read_double( szFileName, "velocityWall3", &velocityWall[2] );
	}
	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength, const char *problem){
	// TO DO - implement exit if scenario is wrong
	// TO DO - think about the direction of initiliasiation - a) is initialized with x starting from upper left corner, whereas c is initialised with x starting from bottomleft corner
	int x, y, z, i;
	int xlen2 = xlength[0] + 2;
	int ylen2 = xlength[1] + 2;
	int zlen2 = xlength[2] + 2;
	int xyz_field_domain = xlen2 * ylen2 * zlen2 * Q_NUMBER;
	int xylen2 = xlen2 * ylen2;
	unsigned int X_min, Y_min, Z_min, X_max, Y_max, Z_max;

	/** flagField init: 
		* 0 - FLUID
		* 1 - NO_SLIP
		* 2 - MOVING_WALL
		* 3 - FREE-SLIP
		* 4 - INFLOW
		* 5 - OUTFLOW
		* 6 - PRESSURE_IN
	**/
	/* Boundary initialization: According to scenario */

	if (problem == "karman_vortex_street")
	{

		X_min = NO_SLIP;		// x- dimension	
		Y_min = FREE_SLIP;	// y- dimension
		Z_min = PRESSURE_IN;// z- dimension

		X_max = NO_SLIP;		// x+ dimension
		Y_max = FREE_SLIP;	// y+ dimension
		Z_max = OUTFLOW;		// z+ dimension

// hard part
// need to scale input file according to dimensions, if they're different we throw an error
// this is a repetitions of the pgm_read() function, but without using the tricky array functions and scaling to support different x-,y-,z- dimensions.

		FILE *input = NULL;
		char line[1024];
		int levels;
		int xsize, ysize;
		int i1, j1;
		int **pic = NULL;
		const char *filename = "./configs/lbm_tilted_plate.vtk";
		int scale_x;
		int scale_y;
		if ((input=fopen(filename,"rb"))==0)
		{
		   char szBuff[80];
		       sprintf( szBuff, "Can not read file %s !!!", filename );
		       ERROR( szBuff );
		}

		/* check for the right "magic number" */
		if ( fread(line,1,3,input)!=3 )
		{
		        fclose(input);
		        ERROR("Error Wrong Magic field!");
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
			return;
		}
		else if(xlength[2] % xsize != 0 || xlength[0] % ysize != 0)
		{
			printf("Dimensions and xlength mismatch!\n");
			return;
		}
		// need independent dim-s, so xlength has to be multiple of corresponding(!) xsize and ysize
		scale_x = xlength[2] / xsize; // our implementation uses different naming
		scale_y = xlength[0] / ysize; // - xsize is xlength[2] and ysize is xlength[0]  
		for(j1 = 1; j1 < ysize + 1; ++j1){
			for (i1 = 1; i1 < xsize+1; ++i1){
				int byte;
				if(fscanf(input, "%d", &byte));

				if (byte==EOF)
				{
					fclose(input);
					ERROR("read failed");
				}
				else
				{
// every read initialises a "cube" from the 3D flagField, scaled according to x-,y-,z- length
					for (z = xsize * scale_x; z < xsize * scale_x + scale_x; ++z) {
						for (y = 1; y <= xlength[1]; ++y) { // height is always traversed fully
							for (x = ysize * scale_y; x < ysize * scale_y + scale_y; ++x) {
								flagField[x + y * ylen2 + z * xylen2] = byte;
							}
						}
					}
				}
			}
		}
		/* close file */
		fclose(input);

	}
	else if (problem == "plane_shear_flow")
	{

		X_min = FREE_SLIP;	// x- dimension	
		Y_min = NO_SLIP;		// y- dimension
		Z_min = PRESSURE_IN;// z- dimension

		X_max = FREE_SLIP;	// x+ dimension
		Y_max = NO_SLIP;		// y+ dimension
		Z_max = OUTFLOW;		// z+ dimension

		// Fluid init (inner part of flagField).
		for (z = 1; z <= xlength[2]; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]; ++x) {
					flagField[x + y * ylen2 + z * xylen2] = FLUID; // we could use a single pointer instead of calculating the whole offset multiple times ... requires effort vs small performance since only 1 call
				}
			}
		}

	}
	else if (problem == "flow_over_step")
	{

		X_min = NO_SLIP;	// x- dimension	
		Y_min = NO_SLIP;	// y- dimension
		Z_min = INFLOW;		// z- dimension

		X_max = NO_SLIP;	// x+ dimension
		Y_max = NO_SLIP;	// y+ dimension
		Z_max = OUTFLOW;	// z+ dimension

		// Fluid init (inner part of flagField).
		for (z = 1; z <= xlength[2]; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]; ++x) {
					flagField[x + y * ylen2 + z * xylen2] = FLUID; // we could use a single pointer instead of calculating the whole offset multiple times ... requires effort vs small performance since only 1 call
				}
			}
		}

		// add "step" initialization - "cube" has z = x = xlen/2, y = ylen
		for (z = 1; z <= xlength[0]/2; ++z) {
			for (y = 1; y <= xlength[1]; ++y) {
				for (x = 1; x <= xlength[0]/2; ++x) {
					flagField[x + y * ylen2 + z * xylen2] = NO_SLIP;
				}
			}
		}

	}
	else
	{

		printf("\n\nUnknown / Misspelled scenario name!\n\n");
		printf("Known scenarios are:\n\n");
		printf("A tilted plane\n");
		printf("Plane shear flow\n");
		printf("Flow over a step\n");
		return;
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

		// Boundary initialization.
		for (z = 0; z < zlen2; ++z){ // treat x, y and z as pure iterators, the dimension depends on where we have 0 or xlength+1 and iteration happens.
			for (y = 0; y < ylen2; ++y){
				for (x = 0; x < xlen2; ++x){
					// add all other walls in the same loop, simply switch the indices
					// - x and y are iterators, use x with xlensq for memory locality
					// z is "0" or "xlength + 1"
					// the only question is whether this is OK with the cache memory
					// in fact we can move these iterations above at line 19(before the third loop from 0 to xlength + 1)
					flagField[	  y * ylen2 + z * xylen2	] = X_min;						// x- dimension	
					flagField[x							+ z * xylen2	] = Y_min;						// y- dimension
					flagField[x + y * ylen2								] = Z_min;						// z- dimension

					flagField[xlen2 - 1 +	y * ylen2 + z * xylen2		] = X_max;	// x+ dimension
					flagField[x + (ylen2 - 1) * ylen2 + z * xylen2	] = Y_max;	// y+ dimension
					flagField[x + y * ylen2 + (zlen2 - 1) * xylen2	] = Z_max;	// z+ dimension
				}
			}
		}

}