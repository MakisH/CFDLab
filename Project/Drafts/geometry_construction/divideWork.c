#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void readPGM_cpu(int *flagField, int *cpuDomain, int *flagField_cpu, int *xlength, int *scale){

	// INPUT: flagField ... The Big flagField
	// OUTPUT: flagField_cpu ... flagField belonging to current cpu.
	// cpuDomain .. the domain size of each CPU, without boundaries
	// flagField_cpu ... The Small flagField of each CPU, size with boundaries included, so cpuDomain+2!
	// rank ... rank of CPU
	// xlength ... dimensions of flagField (like mentioned, the flagField doesn't contain any boundaries
	// scale ... vector [a, b, c]. rank should be put to a component that we take an interval from.

	// Factors for looping through 1D array
	int xlen2_cpu = cpuDomain[0];
	int xylen2_cpu = xlen2_cpu  * (cpuDomain[1]);

	// Again factors.
	int xlen2 = xlength[0];
	int xylen2 = (xlength[0]) * (xlength[1]);

	int x, y, z; // iterators for flagField
	int x_cpu = 0, y_cpu = 0, z_cpu = 0; // iterators for flagField_cpu. We start from 1 because flagField_cpu contains boundary and flagField doesn't. We copy from whole flagField to inner part of flagField_cpu.

	// We do a mapping from global coordinate system (flagField) to local coordinate system (flagField_cpu)
	for (z = scale[2]*cpuDomain[2]; z < (scale[2] + 1) * (cpuDomain[2]); z++){
		for (y = scale[1]*cpuDomain[1]; y < (scale[1] + 1) * (cpuDomain[1]); y++){
			for (x = scale[0]*cpuDomain[0]; x < (scale[0] + 1) * (cpuDomain[0]); x++){
				flagField_cpu[x_cpu + xlen2_cpu * y_cpu + xylen2_cpu * z_cpu] = flagField[x + xlen2 * y + xylen2 * z];
				x_cpu++;
			}
			x_cpu = 0;
			y_cpu++;
		}
		z_cpu++;
		y_cpu = 0;
	}
}
/* USELESS FOR NOW.
void cpuBoundary_connect(int *flagField_cpu, int cpu_connection, int *cpuDomain, int spaceFinger, int *xlength, char *buildingSide){
	// Input:
	// cpuDomain .. the domain size of each CPU, without boundaries
	// flagField_cpu ... The Small flagField of each CPU, size with boundaries included, so cpuDomain+2!

	// cpu_connection ... a flag that shall be put on the connection with finger. It represents the number of CPU (rank) to which we connect within this finger.
	// spaceFinger ... how many units(!) until the start of the finger.
	// buildingSide ... do you want to connect to finger on Lower or Upper side?

	int x, y, z;
	int interval[2];
	interval[0] = spaceFinger * xlength[0]/16; // 16 because we have 1:16 ratio between the size of the whole building and a finger.
	interval[1] = interval[0] + xlength[0]/16; // 1*cpu_domain[0] because finger has a ratio-style width of 1.

	int xlen2 = cpuDomain[0]+2;
	int xylen2 = (cpuDomain[0]+2) * (cpuDomain[1]+2);

	if (strcmp(buildingSide, "lower") == 0) y = cpuDomain[1]+1; // this is the last index of y, not cpuDomain[1]+2!!
	else if(strcmp(buildingSide, "upper") == 0) y = 0;

	for (z = 0; z < cpuDomain[2]+2; z++){
		for (x = interval[0]; x < interval[1]; x++){
			flagField_cpu[x + xlen2 * y + xylen2 * z] = cpu_connection;
		}
	}
}

void cpuBoundary(int *flagField_cpu, int typeofBoundary, char *plane, int *cpuDomain){

	// The function sets a given typeofBoundary for a give plane of the cube.
	// planes: xy, xz, yz, xy_max, xz_max, yz_max

	int x, y, z;
	int xlen2 = cpuDomain[0]+2;
	int xylen2 = xlen2 * (cpuDomain[1]+2);

	if(strcmp(plane, "xy") == 0) {
		z = 0;
		for (y = 0; y < cpuDomain[1]+2; y++){
			for (x = 0; x < cpuDomain[0]+2; x++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	} else if (strcmp(plane, "xy_max") == 0) {
		z = cpuDomain[2]+2;
		for (y = 0; y < cpuDomain[1]+2; y++){
			for (x = 0; x < cpuDomain[0]+2; x++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	} else if (strcmp(plane, "xz") == 0) {
		y = 0;
		for (z = 0; z < cpuDomain[2]+2; z++){
			for (x = 0; x < cpuDomain[0]+2; x++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	} else if (strcmp(plane, "xz_max") == 0) {
		y = cpuDomain[1] + 2;
		for (z = 0; z < cpuDomain[2]+2; z++){
			for (x = 0; x < cpuDomain[0]+2; x++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	} else if (strcmp(plane, "yz") == 0) {
		x = 0;
		for (z = 0; z < cpuDomain[2]+2; z++){
			for (y = 0; y < cpuDomain[1]+2; y++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	} else if (strcmp(plane, "yz_max") == 0) {
		x = cpuDomain[0] + 2;
		for (z = 0; z < cpuDomain[2]+2; z++){
			for (y = 0; y < cpuDomain[1]+2; y++){
				flagField_cpu[x + y*xlen2 + z*xylen2] = typeofBoundary;
			}
		}
	}
}
*/
