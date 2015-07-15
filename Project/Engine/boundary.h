#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include "initLB.h"

/** handles the boundaries in our simulation setup */
void treatBoundary(double *collideField,
									 int* flagField,
									 const double * const wallVelocity,
									 const double * const ref_density,
                   const int * const cpuDomain,
                   double_3d *velocityIn,
                   double *density_in);

#endif
