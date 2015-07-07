#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>

#define Q_NUMBER 19

#define FLUID 0
#define NO_SLIP 1
#define FREE_SLIP 2
#define OUTFLOW 3

#define PARALLEL_BOUNDARY 7

#define MOVING_WALL 10

#define INFLOW 30
#define INFLOW_1 31
#define INFLOW_2 32
#define INFLOW_3 33
#define INFLOW_4 34
#define INFLOW_5 35

#define PRESSURE_IN 50
#define PRESSURE_IN_1 51
#define PRESSURE_IN_2 52
#define PRESSURE_IN_3 53
#define PRESSURE_IN_4 54
#define PRESSURE_IN_5 55

#define DIR_R 0
#define DIR_L 1
#define DIR_T 2
#define DIR_D 3
#define DIR_B 4
#define DIR_F 5

// opposite direction : direction + 1  - direction % 2 * 2
static const int LATTICEVELOCITIES[Q_NUMBER][3] = {{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {0, -1, 0}, {1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}, {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
static const double LATTICEWEIGHTS[Q_NUMBER] = {1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36};

static const double C_S = 0.577350269189626;

// Precalculations for performance
static const double C_S_sq = 0.577350269189626 * 0.577350269189626;

const static int each[6][5] = {	{3, 7, 10,13,17},	// x+ R
																{1, 5, 8, 11,15},	// x- L
																{14,15,16,17,18},	// z+ T
																{0, 1, 2, 3, 4},	// z- D
																{4, 11,12,13,18},	// y+ B
																{0, 5, 6, 7, 14}};// y- F

#endif
