#include "uvp.h"
#include "helper.h"
#include <math.h>

void calculate_dt (double Re,
	      double tau,
	      double *dt,
	      double dx,
	      double dy, int imax, int jmax, double **U, double **V)
{


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
  min = fmin ( min, 0.5 * Re * 1.0 / (1 / (dx * dx) + 1 / (dy * dy)) );

  // Finally, the calc of dt.
  *dt = tau * min;

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
	    1 / dx * ( ((u + u1)*(u + u1)/4) - ((u3 + u)*(u3 + u)/4) ) +
	    alpha / (4 * dx) * (fabs (u + u1) * (u - u1) -
				fabs (u3 + u) * (u3 - u));
	  uv_y =
	    1 / (4 * dy) * ((v + v1) * (u + u2) - (v4 + v6) * (u4 + u)) +
	    alpha / (4 * dy) * (fabs (v + v1) * (u - u2) -
				fabs (v4 + v6) * (u4 - u));

	  uxx = 1 / (dx * dx) * (u1 - 2 * u + u3); 
	  uyy = 1 / (dy * dy) * (u2 - 2 * u + u4); 

	  F[i][j] = u + dt * (1.0 / Re * (uxx + uyy) - u2x - uv_y + GX); // calculation of F
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
	    1 / dy * ( (0.5 * (v + v2))*(0.5 * (v + v2)) - (0.5 * (v4 + v))*(0.5 * (v4 + v)) ) +
	    alpha / (4 * dy) * (fabs (v + v2) * (v - v2) -
				fabs (v4 + v) * (v4 - v));
	  uv_x =
	    1 / (4 * dx) * ((u + u2) * (v + v1) - (u3 + u5) * (v3 + v)) +
	    alpha / (4 * dx) * (fabs (u + u2) * (v - v1) -
				fabs (u3 + u5) * (v3 - v));

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
    double dx_dt = 1 / (dx * dt);
    double dy_dt = 1 /  (dy * dt);
	
    for(unsigned int i = 1; i <= imax; ++i)
    {
      for(unsigned int j = 1; j <= jmax; ++j)
      {
	      RS[i][j] = (F[i][j] - F[i-1][j]) * dx_dt + (G[i][j] - G[i][j-1]) * dy_dt;
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

  for (int i = 1; i < imax; ++i)
  {
    for (int j = 1; j <= jmax; ++j)
    {
      // Update the horizontal velocities (explicit Euler method)
      U[i][j] = F[i][j] + _dt_dx * ( P[i+1][j] - P[i][j] );
      // Update the vertival velocities (explicit Euler method)
      V[i][j] = G[i][j] + _dt_dy * ( P[i][j+1] - P[i][j] );
    }
  }

}
