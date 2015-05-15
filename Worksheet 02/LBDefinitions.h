#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>

  static const int LATTICEVELOCITIES[19][3] = {{0, -1, -1}, {-1, 0, -1}, {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {-1, -1, 0}, {-1, 0, 0}, {0, 0, 0}, {1, 0, 0}, {-1, 1, 0}, {0, 1, 0}, {1, 1, 0}, {0, -1, 1}, {-1, 0, 1}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
  static const double LATTICEWEIGHTS[19] = {1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 12.0/36, 2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36};

  static const double C_S = 0.577350269189626;

  static const unsigned int Q_NUMBER = 19;
  static const unsigned int FLUID = 0;
  static const unsigned int NO_SLIP = 1;
  static const unsigned int MOVING_WALL = 2;

  // Precalculations for performance
  static const double C_S_sq = 0.577350269189626 * 0.577350269189626;

#endif

