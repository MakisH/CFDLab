void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
)
{


  for (int j=1; j<=jmax; ++j)
    {
      // Horizontal velocities on vertical boundaries
      U[0][j] = 0;
      U[imax][j] = 0;
      // Vertical boundaries on vertical boundaries
      V[0][j]       = - V[1][j];
      V[imax+1][j]  = - V[imax][j];
    }

  for (int i=1; i<=imax; ++i)
    {
      // Vertical velocities on horizontal boundaries
      V[i][0]       = 0;
      V[i][jmax]    = 0;

      // Horizontal velocities on horizontal boundaries
      U[i][0]       = - U[i][1];
      U[i][jmax+1]  = 2.0 - U[i][jmax];
    }

}
