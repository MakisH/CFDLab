#include "helper.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){

  // Usage info
  if ( argn > 3 ) {
    printf("Usage: sim [input file] [problem name] \n");
    return 1;
  }

  // Input file
  char *filename = NULL;
  if ( argn == 2 ) {
    filename = args[1];
  } else {
    filename = "cavity100.dat";
  }

  // Name of the problem (to be used in the output filename)
  char *problem = NULL;
  if ( argn == 3) {
    problem = args[2];
  } else {
    problem = "cavity100";
  }

  // Parameters declaration
  // Geometry data
  double * xlength = malloc(sizeof(double));
  double * ylength = malloc(sizeof(double));
  int    * imax = malloc(sizeof(int)); 
  int    * jmax = malloc(sizeof(int));
  double * dx = malloc(sizeof(double));
  double * dy = malloc(sizeof(double));

  // Time-stepping data
  double * t_end = malloc(sizeof(double));
  double * dt = malloc(sizeof(double));
  double * tau = malloc(sizeof(double));
  double * dt_value = malloc(sizeof(double));
  
  // Pressure iteration data
  int    * itermax = malloc(sizeof(int));
  double * eps = malloc(sizeof(double));
  double * omg = malloc(sizeof(double));
  double * alpha = malloc(sizeof(double));

  // Problem-dependent quantities
  double * Re = malloc(sizeof(double));
  double * GX = malloc(sizeof(double)); 
  double * GY = malloc(sizeof(double));
  double * UI = malloc(sizeof(double));
  double * VI = malloc(sizeof(double));
  double * PI = malloc(sizeof(double));

  // Read the input file
  read_parameters(filename, Re, UI, VI, PI, GX, GY, t_end, 
                  xlength, ylength, dt, dx, dy, imax, jmax, 
                  alpha, omg, tau, itermax, eps, dt_value);

  // Setup arrays
  double **U  = (double**)matrix(0, *imax+1, 0, *jmax+1);
  double **V  = (double**)matrix(0, *imax+1, 0, *jmax+1);
  double **P  = (double**)matrix(0, *imax+1, 0, *jmax+1);
  double **RS = (double**)matrix(0, *imax+1, 0, *jmax+1);
  double **F  = (double**)matrix(0, *imax+1, 0, *jmax+1);
  double **G  = (double**)matrix(0, *imax+1, 0, *jmax+1);
  
  // Some help variables
  double  t = 0.0;
  int     n = 0;
  int     it = 0;
  double * res = malloc(sizeof(double));
  *res = 100 * *eps; // just larger than eps
  
  // Assign initial values
  init_uvp( *UI, *VI, *PI, *imax, *jmax, U, V, P );
  
  // Time loop
  while ( t <= *t_end )
  {
    // Calculate dt if read tau is not negative
    if (*tau > 0) calculate_dt(*Re, *tau, dt, *dx, *dy, *imax, *jmax, U, V);

    // Set the boundary values
    boundaryvalues(*imax, *jmax, U, V);
    
    // Calculate F and G terms of the pressure Poisson equation
    calculate_fg(*Re, *GX, *GY, *alpha, *dt, *dx, *dy, *imax, *jmax, U, V, F, G);

    // Calculate the right-hand side of the pressure Poisson equation
    calculate_rs(*dt, *dx, *dy, *imax, *jmax, F, G, RS);

    // SOR loop
    it = 0; // reset the number of sor iterations
    *res = 100 * *eps; // just larger than eps
    while ( it < *itermax && *res > *eps ) 
    {
      sor(*omg, *dx, *dy, *imax, *jmax, P, RS, res); // one SOR iteration
      it++;
    }

    // Update the velocities
    calculate_uv(*dt, *dx, *dy, *imax, *jmax, U, V, F, G, P);

    // Output of u, v, p for visualization
    if (n%100==0)
    write_vtkFile(problem, n, *xlength, *ylength, *imax, *jmax, *dx, *dy, U, V, P);

    // Update loop state
    t = t + *dt;
    n = n + 1;

  }

  //printf("result = %f \n",U[*imax/2][7 * *jmax/8]);
  
  // Free arrays
  free_matrix( U, 0, *imax+1, 0, *jmax+1);
  free_matrix( V, 0, *imax+1, 0, *jmax+1);
  free_matrix( P, 0, *imax+1, 0, *jmax+1);
  free_matrix( RS, 0, *imax+1, 0, *jmax+1);
  free_matrix( F, 0, *imax+1, 0, *jmax+1);
  free_matrix( G, 0, *imax+1, 0, *jmax+1);

  // Free other pointers
  free(xlength);
  free(ylength);
  free(imax);
  free(jmax);
  free(dx);
  free(dy);
  free(t_end);
  free(dt);
  free(tau);
  free(dt_value);
  free(itermax);
  free(eps);
  free(omg);
  free(alpha);
  free(res);

  return 0;

}
