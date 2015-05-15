#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q_NUMBER 19

  static const int LATTICEVELOCITIES[Q_NUMBER][3];
  static const double LATTICEWEIGHTS[Q_NUMBER];
  static const double C_S = 1.0/sqrt(3.0);

  static const unsigned int FLUID = 0;
  static const unsigned int NO_SLIP = 1;
  static const unsigned int MOVING_WALL = 2;

  // Precalculations for performance
  static const double C_S_sq = C_S * C_S;

#endif

