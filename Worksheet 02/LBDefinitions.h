#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#include <math.h>

  static const int LATTICEVELOCITIES[19][3];
  static const double LATTICEWEIGHTS[19];
  static const double C_S = 1.0 / sqrt(3.0);

  static const unsigned int Q_NUMBER = 19;
  static const unsigned int FLUID = 0;
  static const unsigned int NO_SLIP = 1;
  static const unsigned int MOVING_WALL = 2;

#endif

