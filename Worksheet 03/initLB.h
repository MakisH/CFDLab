#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"
#include "LBDefinitions.h"

/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(char		*problem,							// scenario name - "flow_over_step", "karman_vortex_street" or "plane_shear_flow"
									 int		*xlength,							// reads domain size. Parameter name: "xlength" 
									 double *tau,									// relaxation parameter tau. Parameter name: "tau" 
									 int		*timesteps,						// number of timesteps. Parameter name: "timesteps" 
									 int		*timestepsPerPlotting,// timesteps between subsequent VTK plots. Parameter name: "vtkoutput" 
									 double *velocityIn,					// velocity for sc 1 and 3
									 double *densityIn,						// densityin for sc2
									 double *densityRef,					// densityref for sc2
									 double *velocityWall,				// velocity of the lid. Parameter name: "characteristicvelocity" 
									 int		*initxyzXYZ,					// boundary flags array - 1-6 for no_slip, moving_wall, free_slip, inflow, outflow, pressure_in
									 int		argc,									// number of arguments. Should equal 2 (program + name of config file 
									 char		*argv[]){							// argv[1] shall contain the path to the config file 

/* initialises the particle distribution functions and the flagfield */
int initialiseFields(double *collideField, double *streamField, int *flagField, int *xlength, char *problem, int *initxyzXYZ);

#endif

