#include "uvp.h"
#include "helper.h"
#include <math.h>

void calculate_dt (double Re,
	      double tau,
	      double *dt,
	      double dx,
	      double dy, int imax, int jmax, double **U, double **V)
{

   /* if tau negative we don't do anything in this function, therefore the dt from file is used */
  if (tau >= 0){

  // Let's find the max abs values of U and V.
  double maxU = 0, maxV = 0;	// we can afford setting the init values to 0, since we deal with abs values.

  for (int i = 1; i <= imax; ++i)
    {
      for (int j = 1; j <= jmax; ++j)
	{
	  maxU = fmax (maxU, fabs (U[i][j]));
	  maxV = fmax (maxV, fabs (V[i][j]));
	}
    }

  // Finding the min of the prescribed three values.
  double min = fmin (dx / maxU, dy / maxV);
  min = fmin (min, 1.0 / 2 * Re * pow (1 / (dx * dy) + 1 / (dy * dy), -1));
  //printf("%f   %f   %f   %f \n", dx, maxU, dy, maxV);
  //printf("%f", 1.0 / 2 * Re * pow (1 / (dx * dy) + 1 / (dy * dy), -1));
  //printf("%f   %f \n",tau, min);
  // Finally, the calc of dt.
  *dt = tau * min;
 }
}

void calculate_fg (double Re,
	      double GX,
	      double GY,
	      double alpha,
	      double dt,
	      double dx,
	      double dy,
	      int imax,
	      int jmax, double **U, double **V, double **F, double **G)
{

  double u, u1, u2, u3, u4, u5, v, v1, v2, v3, v4, v6; // temp variables for easier construction of derivatives
  double u2x, uv_y, uxx, uyy, v2y, uv_x, vxx, vyy; // derivatives


  for (int i = 1; i <= imax - 1; ++i)
    {
      for (int j = 1; j <= jmax; ++j)
	{
	  u = U[i][j];
	  u1 = U[i + 1][j];
	  u2 = U[i][j + 1];
	  u3 = U[i - 1][j];
	  u4 = U[i][j - 1];
	  u5 = U[i - 1][j + 1];

	  v = V[i][j];
	  v1 = V[i + 1][j];
	  v2 = V[i][j + 1];
	  v3 = V[i - 1][j];
	  v4 = V[i][j - 1];
	  v6 = V[i + 1][j - 1];

	  u2x =
	    1 / dx * (pow (1 / 2 * (u + u1), 2) - pow (1 / 2 * (u3 + u), 2)) +
	    alpha / (2 * dx) * (abs (u + u1) * (u - u1) -
				abs (u3 + u) * (u3 - u));
	  uv_y =
	    1 / (2 * dy) * ((v + v1) * (u + u2) - (v4 + v6) * (u4 + u)) +
	    alpha / (2 * dy) * (abs (v + v1) * (u - u2) -
				abs (v4 + v6) * (u4 - u));

	  uxx = 1 / (dx * dx) * (u1 - 2 * u + u3);
	  uyy = 1 / (dy * dy) * (u2 - 2 * u + u4);

	  F[i][j] = u + dt * (1 / Re * (uxx + uyy) - u2x - uv_y + GX); // calculation of F
	}
    }



  for (int i = 1; i <= imax; ++i)
    {
      for (int j = 1; j <= jmax - 1; ++j)
	{
	  u = U[i][j];
	  u1 = U[i + 1][j];
	  u2 = U[i][j + 1];
	  u3 = U[i - 1][j];
	  u4 = U[i][j - 1];
	  u5 = U[i - 1][j + 1];

	  v = V[i][j];
	  v1 = V[i + 1][j];
	  v2 = V[i][j + 1];
	  v3 = V[i - 1][j];
	  v4 = V[i][j - 1];
	  v6 = V[i + 1][j - 1];

	  v2y =
	    1 / dy * (pow (1 / 2 * (v + v2), 2) - pow (1 / 2 * (v4 + v), 2)) +
	    alpha / (2 * dy) * (abs (v + v2) * (v - v2) -
				abs (v4 + v) * (v4 - v));
	  uv_x =
	    1 / (2 * dx) * ((u + u2) * (v + v1) - (u3 + u5) * (v3 + v)) +
	    alpha / (2 * dx) * (abs (u + u2) * (v - v1) -
				abs (u3 + u5) * (v3 - v));

	  vxx = 1 / (dx * dx) * (v1 - 2 * v + v3);
	  vyy = 1 / (dy * dy) * (v2 - 2 * v + v4);

	  G[i][j] = v + dt * (1 / Re * (vxx + vyy) - v2y - uv_x + GY); // Calculation of G.
	}
    }


  /* Boundary conditions */

  for (int i = 1; i <= imax; ++i)
    {
      G[i][0] = V[i][0];
      G[i][jmax] = V[i][jmax];
    }

  for (int j = 1; j <= jmax; ++j)
    {
      F[0][j] = U[0][j];
      F[imax][j] = U[imax][j];
    }

}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
){
	double fracdxdt = 1 / (dx * dt);
	double fracdydt = 1 /  (dy * dt);
	double *RSi;
	double *Fi;
	double *Fiminus1;
	double *Gi;
    // calculate rs for the fluid cells
    for(unsigned int i = 1; i <= imax + 1; ++i){
    	// could pull R[i], G[i] and F[i] calculation here
 	RSi = RS[i];
 	Fi = F[i];
 	Fiminus1 = F[i - 1];
 	Gi = G[i];
        for(unsigned int j = 1; j <= jmax + 1; ++j){
		// calculate only interior cells, but use the boundary
		RSi[j] = (Fi[j] - Fiminus1[j]) * fracdxdt + (Gi[j] - Gi[j - 1]) * fracdydt;
        }
    }
}

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{

  // Pre-calculate the -dt/dx and -dt/dy factors
  double _dt_dx = - dt / dx;
  double _dt_dy = - dt / dy; 

  // Update the horizontal velocities (explicit Euler method)
  for (int i = 1; i < imax; ++i)
  {
    for (int j = 1; j < jmax + 1; ++j)
    {
      U[i][j] = F[i][j] + _dt_dx * ( P[i+1][j] - P[i][j] );
    }
  }

  // Update the vertival velocities (explicit Euler method)
  for (int i = 1; i < imax + 1; ++i)
  {
    for (int j = 1; j < jmax; ++j)
    {
      V[i][j] = G[i][j] + _dt_dy * ( P[i][j+1] - P[i][j] );
    }
  }

}
