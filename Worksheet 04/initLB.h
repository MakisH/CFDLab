#ifndef _INITLB_H_
#define _INITLB_H_

/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
		int * const xlength,                       /* reads domain size. Parameter name: "xlength" */
		double * const tau,                        /* relaxation parameter tau. Parameter name: "tau" */
		double * const velocityWall,               /* velocity of the lid. Parameter name: "characteristicvelocity" */
		int * const timesteps,                     /* number of timesteps. Parameter name: "timesteps" */
		int * const timestepsPerPlotting,          /* timesteps between subsequent VTK plots. Parameter name: "vtkoutput" */
		int * const iProc,                         /* Number of processes-subdomains per x-direction */
		int * const jProc,                         /* Number of processes-subdomains per y-direction */
		int * const kProc,                         /* Number of processes-subdomains per z-direction */
		 int argc,                            /* number of arguments. Should equal 2 (program + name of config file */
		 char *  *  argv            /* argv[1] shall contain the path to the config file */
);

/* initialises the particle distribution functions and the flagfield */
void initialiseFields(
	double *collideField,
	double *streamField,
	int *flagField,
	const int * const cpuDomain,
	const int iProc,
	const int jProc,
	const int kProc,
	const int rank,
	int * const neighbor
);

void initialiseBuffers(
	double **sendBuffer,
	double **readBuffer,
	const int * const cpuDomain,
	int *sizeBuffer,
	const int * const neighbor
);
#endif

