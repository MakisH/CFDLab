#ifndef _VISUALLB_H_
#define _VISUALLB_H_
#include <stdio.h>


/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(const double * const collideField, const int * const flagField, char *filename, unsigned int t, const int * const xlength, const int part_id);

void write_vtkHeader( FILE *fp, const int * const xlength);

void write_vtkPointCoordinates( FILE *fp, int part_id);

#endif
