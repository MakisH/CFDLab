#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>

#define Q_NUMBER 19
#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2
// #define FREE_SLIP 3
// #define INFLOW 4
// #define OUTFLOW 5
// #define PRESSURE_IN 6
#define PARALLEL_BOUNDARY 7

	static const int LATTICEVELOCITIES[Q_NUMBER][3] = {{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}, {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
	static const double LATTICEWEIGHTS[Q_NUMBER] = {1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36};

	static const double C_S = 0.577350269189626;

	// Precalculations for performance
	static const double C_S_sq = 0.577350269189626 * 0.577350269189626;

#endif
