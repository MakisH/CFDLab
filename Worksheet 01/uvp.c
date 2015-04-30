#include "uvp.h"
#include "helper.h"
#include <math.h>

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){

	// Let's find the max abs values of U and V.
	double maxU=0, maxV=0; // we can afford setting the init values to 0, since we deal with abs values.

	for (int i=0; i<imax; ++i){
		for (int j=0; j<jmax; ++j){
			maxU=fmax(maxU, abs(U[i][j]));
			maxV=fmax(maxV, abs(V[i][j]));
		}
	}

	// Finding the min of the prescribed three values.
	double min = fmin(dx/maxU, dy/maxV);
	min = fmin(min, 1/2 * Re * pow(1/(dx*dy) + 1/(dy*dy), -1));

	// Finally, the calc of dt.
	*dt = tau * min;
	// BE CAREFUL! IF TAU<0, THEN dt FROM THE FILE (VIA read_parameters(..)) SHOULD BE USED. THIS HAS TO BE IMPLEMENTED IN MAIN FILE.
}

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
){

}

