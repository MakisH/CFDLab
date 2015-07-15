#ifndef _INITLB_H_
#define _INITLB_H_

typedef struct {
	double x;
	double y;
	double z;
} double_3d;

/* Reads PGM file and assigns it to flagField. */
int read_assign_PGM (int *flagField, char *fileName, const int * const cpuDomain);


/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
		double * const tau,                        /* relaxation parameter tau. Parameter name: "tau" */
		double * const velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
		int * const timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
		int * const timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
		double_3d * const inflow,											/* array of structures (AoS) holding inflow values for x, y and z directions */
		double * const pressure_in,								/* AoS holding pressure_in values for x, y and z directions */
		double * const ref_density,
		 int argc,                            /* number of arguments. Should equal 2 (program + name of config file */
		 char *  *  argv            /* argv[1] shall contain the path to the config file */
);

/* initialises the particle distribution functions and the flagfield */
void initialiseFields(
	double *collideField,
	double *streamField,
	const int * const cpuDomain
);

#endif
