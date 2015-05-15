#include "computeCellValues.h"

void computeDensity(const double * const currentCell, double *density){
  /* TODO */
}



void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
  /* TODO */
}



void computeFeq(const double * const density, const double * const velocity, double *feq){
  
  // precalculate inner product u*u
  double u_u = velocity[1] * velocity[1] + velocity[2] * velocity[2] + velocity[3] * velocity[3];

  for (int i=0; i<Q_NUMBER; i++) {
    // precalculate inner product c*u
    double c_u = LATTICEVELOCITIES[i][1] * velocity[1] +   LATTICEVELOCITIES[i][2] * velocity[2] +   LATTICEVELOCITIES[i][3] * velocity[3]
    // calculate the equilibrium distribution element
    feq[i] = LATTICEWEIGHTS[i] * density * (  1 
                                              + c_u / C_S_sq
                                              + c_u * c_u / ( 2 * C_S_sq * C_S_sq )
                                              - u_u / ( 2 * C_S_sq )   )           
  }
  

}

