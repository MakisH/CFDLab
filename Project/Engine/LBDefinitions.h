#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>
#include "stdlib.h"

#define Q_NUMBER 19

#define FLUID 0
#define NO_SLIP 1
#define FREE_SLIP 2
#define OUTFLOW 3

#define PARALLEL_BOUNDARY 7

#define MOVING_WALL 10

#define INFLOW_COUNT 6
#define INFLOW 30
#define INFLOW_1 31
#define INFLOW_2 32
#define INFLOW_3 33
#define INFLOW_4 34
#define INFLOW_5 35

#define PRESSURE_IN_COUNT 6
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

static const int * const neighbours[8] =  {(int[]){1,4,6},(int[]){0,2,4,6},(int[]){1,3,4,7},(int[]){2,5,7},(int[]){0,1,2},(int[]){3},(int[]){0,1},(int[]){2,3}}; 
// those are the neighbouring processors for each processor, e.g. proc.0(first element of array) has neighbours processors 1,4 and 6
// int[] is needed to make a pointer to the inside element array
static const int neighbours_count[8] = {3,5,4,4,3,2,3,2};
static const int * const neighbours_dir[8] = {(int[]){DIR_R,DIR_F,DIR_B},
																							(int[]){DIR_L,DIR_R,DIR_F,DIR_B,DIR_B},
																							(int[]){DIR_L,DIR_R,DIR_F,DIR_B},
																							(int[]){DIR_L,DIR_F,DIR_F,DIR_B},
																							(int[]){DIR_B,DIR_B,DIR_B},
																							(int[]){DIR_B,DIR_B},
																							(int[]){DIR_F,DIR_F,DIR_F},
																							(int[]){DIR_F,DIR_F}}; // in one direction there may be several transactions(as in cpu 5 and 7)

static const int * const neighbours_procid[8] = {(int[]){1,4,6},
																								 (int[]){0,2,4,6,6},
																								 (int[]){1,3,4,7},
																								 (int[]){2,5,5,7},
																								 (int[]){0,1,2},
																								 (int[]){3,3},
																								 (int[]){0,1,1},
																								 (int[]){2,3}};
// this makes the proccess coupling in the communication step - labeling is done from 1 to 13, starting from process 0 and going in "send" direction up to process 3, labeling only new connections
static const int * const neighbours_tag[8] =		{(int[]){1,4,9},
																								 (int[]){1,2,5,10,11},
																								 (int[]){2,3,6,12},
																								 (int[]){3,7,8,13},
																								 (int[]){4,5,6},
																								 (int[]){7,8},
																								 (int[]){9,10,11},
																								 (int[]){12,13}};

// this is where the sevens are ... meaning it works for injection(because we inject in 7)
static const int * const neighbours_local_start_ext_x[8] = {(int[]){64,48,16},
																														(int[]){0,65,33,1,49},
																														(int[]){0,65,17,33},
																														(int[]){0,1,49,17},
																														(int[]){0,0,0},
																														(int[]){0,0},
																														(int[]){0,0,0},
																														(int[]){0,0}};

static const int * const neighbours_local_start_ext_y[8] = {(int[]){0,0,31},
																														(int[]){0,0,0,31,31},
																														(int[]){0,0,0,31},
																														(int[]){0,0,0,31},
																														(int[]){64,64,64},
																														(int[]){64,64},
																														(int[]){0,0,0},
																														(int[]){0,0}};

static const int * const neighbours_local_start_ext_z[8] = {(int[]){0,0,0},
																														(int[]){0,0,0,0,0},
																														(int[]){0,0,0,0},
																														(int[]){0,0,0,0},
																														(int[]){0,0,0},
																														(int[]){0,0},
																														(int[]){0,0,0},
																														(int[]){0,0}};

static const int * const neighbours_local_end_ext_x[8] = {(int[]){64,63,31},
																													(int[]){0,65,48,16,64},
																													(int[]){0,65,32,48},
																													(int[]){0,16,64,32},
																													(int[]){15,15,15},
																													(int[]){15,15},
																													(int[]){15,15,15},
																													(int[]){15,15}};

static const int * const neighbours_local_end_ext_y[8] = {(int[]){31,0,31},
																													(int[]){31,31,0,31,31},
																													(int[]){31,31,0,31},
																													(int[]){31,0,0,31},
																													(int[]){64,64,64},
																													(int[]){64,64},
																													(int[]){0,0,0},
																													(int[]){0,0}};

static const int * const neighbours_local_end_ext_z[8] = {(int[]){2,2,2},
																													(int[]){2,2,2,2,2},
																													(int[]){2,2,2,2},
																													(int[]){2,2,2,2},
																													(int[]){2,2,2},
																													(int[]){2,2},
																													(int[]){2,2,2},
																													(int[]){2,2}};



static const int * const neighbours_local_start_inj_x[8] = {(int[]){63,48,16},
																														(int[]){1,64,33,1,49},
																														(int[]){1,64,17,33},
																														(int[]){1,1,49,17},
																														(int[]){0,0,0},
																														(int[]){0,0},
																														(int[]){0,0,0},
																														(int[]){0,0}};

static const int * const neighbours_local_start_inj_y[8] = {(int[]){0,1,30},
																														(int[]){0,0,1,30,30},
																														(int[]){0,0,1,30},
																														(int[]){0,1,1,30},
																														(int[]){63,63,63},
																														(int[]){63,63},
																														(int[]){1,1,1},
																														(int[]){1,1}};

static const int * const neighbours_local_start_inj_z[8] = {(int[]){0,0,0},
																														(int[]){0,0,0,0,0},
																														(int[]){0,0,0,0},
																														(int[]){0,0,0,0},
																														(int[]){0,0,0},
																														(int[]){0,0},
																														(int[]){0,0,0},
																														(int[]){0,0}};

static const int * const neighbours_local_end_inj_x[8] = {(int[]){63,63,31},
																													(int[]){1,64,48,16,64},
																													(int[]){1,64,32,48},
																													(int[]){1,16,64,32},
																													(int[]){15,15,15},
																													(int[]){15,15},
																													(int[]){15,15,15},
																													(int[]){15,15}};

static const int * const neighbours_local_end_inj_y[8] = {(int[]){31,1,30},
																													(int[]){31,31,1,30,30},
																													(int[]){31,31,1,30},
																													(int[]){31,1,1,30},
																													(int[]){63,63,63},
																													(int[]){63,63},
																													(int[]){1,1,1},
																													(int[]){1,1}};

static const int * const neighbours_local_end_inj_z[8] = {(int[]){2,2,2},
																													(int[]){2,2,2,2,2},
																													(int[]){2,2,2,2},
																													(int[]){2,2,2,2},
																													(int[]){2,2,2},
																													(int[]){2,2},
																													(int[]){2,2,2},
																													(int[]){2,2}};

static const int * const neighbours_local_buffer_size[8] = {(int[]){480,240,240},
																														(int[]){480,480,240,240,240},
																														(int[]){480,480,240,240},
																														(int[]){480,240,240,240},
																														(int[]){240,240,240},
																														(int[]){240,240},
																														(int[]){240,240,240},
																														(int[]){240,240}};

static const int chunk_count[8] = {1,1,1,1,3,2,3,2};

static const int * const chunk_id[8] = {(int[]){0},
																				(int[]){1},
																				(int[]){2},
																				(int[]){3},
																				(int[]){4,5,6},
																				(int[]){7,8},
																				(int[]){9,10,11},
																				(int[]){12,13}};

static const int * const chunk_begin_offset[8] = {(int[]){0},
																									(int[]){0},
																									(int[]){0},
																									(int[]){0},
																									(int[]){0,3120,6240},
																									(int[]){0,3120},
																									(int[]){0,3120,6240},
																									(int[]){0,3120}};

static const int * const chunk_size[8] = {(int[]){6240},
																					(int[]){6240},
																					(int[]){6240},
																					(int[]){6240},
																					(int[]){3120,3120,3120},
																					(int[]){3120,3120},
																					(int[]){3120,3120,3120},
																					(int[]){3120,3120}}; // actually this can be computed from cpuDomain values
typedef int int3d[3];
static const int3d * const cpuDomain[8] = {(const int3d[1]){{63,30,1}},
																					 (const int3d[1]){{63,30,1}},
																					 (const int3d[1]){{63,30,1}},
																					 (const int3d[1]){{63,30,1}},
																					 (const int3d[3]){{63,14,1},{63,14,1},{63,14,1}},
																					 (const int3d[2]){{63,14,1},{63,14,1}},
																					 (const int3d[3]){{63,14,1},{63,14,1},{63,14,1}},
																					 (const int3d[2]){{63,14,1},{63,14,1}}};

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
