#ifndef _MAIN_C_
#define _MAIN_C_

#include "mpi_helper.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "random_wings.h"
#include <time.h>
//#include <unistd.h>

int main(int argc, char *argv[]){
	printf("hi!\n");
	int rank;
	int np;

	// Start MPI
	initializeMPI( &rank, &np, argc, argv);
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	double_3d * inflow = (double_3d *) malloc(INFLOW_COUNT * sizeof(double_3d));
	double * pressure_in = (double *) malloc(PRESSURE_IN_COUNT * sizeof(double));
	int error_code;
	double ref_density;

	// Read the config file using only one thread
	if(0 == rank){
		error_code = readParameters(&tau, velocityWall, &timesteps, &timestepsPerPlotting, inflow, pressure_in, &ref_density, argc, argv);
		// Error checking
		if(error_code) return error_code;
	}
	MPI_Bcast( &tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( &timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	for(int i = 0; i < INFLOW_COUNT; ++i){
		MPI_Bcast( &inflow[i].x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( &inflow[i].y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
		MPI_Bcast( &inflow[i].z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	}

	MPI_Bcast( pressure_in, PRESSURE_IN_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &ref_density, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	int domain_size = 0;
	for(int i = 0; i < chunk_count[rank];++i){
		domain_size += chunk_whole_flag_size[rank][i];
	}
	domain_size *= 3; // because we're in 3D
	double *collideField = (double *) malloc(Q_NUMBER * domain_size * sizeof(double));
	double *streamField = (double *) malloc(Q_NUMBER * domain_size * sizeof(double));
	int *flagField = (int *) malloc(domain_size * sizeof(int));
	char pgm_read_file[1024];

	for(int i = 0; i < chunk_count[rank]; ++i){
		initialiseFields( collideField + Q_NUMBER * chunk_begin_offset[rank][i], streamField + Q_NUMBER * chunk_begin_offset[rank][i], flagField + chunk_begin_offset[rank][i], cpuDomain[rank][i]); // collide and stream

		sprintf( pgm_read_file, "./pgm/cpu_%d.pgm",chunk_id[rank][i]);
		printf("%s rank = %d \n\n\n\n",pgm_read_file, rank);
		read_assign_PGM(flagField + chunk_begin_offset[rank][i],pgm_read_file,cpuDomain[rank][i]);

	/// Output FlagField for CPU0
		//if(rank == 5){
		//	int x, y, z;
		//	int xlen2 = cpuDomain[rank][i][0]+2;
		//	int ylen2 = cpuDomain[rank][i][1]+2;
		//	int zlen2 = cpuDomain[rank][i][2]+2;
		//			printf("domain %d %d %d\n",xlen2,ylen2,zlen2);
		//		for(z = zlen2 - 1; z >= 0; --z){
		//			for(y = ylen2 - 1; y >= 0; --y){
		//				for(x = 0; x < xlen2; ++x){
		//					printf("%d ",flagField[x + y * xlen2 + z * xlen2 * ylen2]);
		//				}
		//			printf("plane %d rank _%d\n",z, rank);
		//			}
		//			printf("\n");
		//		}
		//		printf("exit initLB \n");
		//}
		//		sleep(1000);
	}

	// send and read buffers for all possible directions :
	// [0:left, 1:right, 2:top, 3:bottom, 4:front, 5:back]
	double * *sendBuffer = (double * *) malloc(neighbours_count[rank] * sizeof(double *));
	double * *readBuffer = (double * *) malloc(neighbours_count[rank] * sizeof(double *));
	for(int i = 0; i < neighbours_count[rank]; ++i){
		sendBuffer[i] = (double *) malloc(neighbours_local_buffer_size[rank][i] * sizeof(double));
		readBuffer[i] = (double *) malloc(neighbours_local_buffer_size[rank][i] * sizeof(double));
	}

	#if defined(INCLUDE_RANDOM_DOORS)
	// INIT OF RANDOM DOOR CLOSENING / OPENING
	srand(time(NULL)*rank); // Initialise the random seed for each cpu
	int offset = 20; // For how long do we leave the doors closed?
	//int frequency_factor = 10; // bigger the f, more frequent we open/close the door

	int *random_timestep = NULL;
	int *wall_trigger = NULL;


	if (rank == 4 || rank == 6){
		random_timestep = malloc(3 * sizeof(int));
		wall_trigger = calloc(3, sizeof(int));

		random_timestep[0] = abs(rand()) % (offset);
		random_timestep[1] = abs(rand()) % (offset);
		random_timestep[2] = abs(rand()) % (offset);

		wall_trigger[0] = 1;
		wall_trigger[1] = 1;
		wall_trigger[2] = 1;

	} else if (rank == 5 || rank == 7) {
		wall_trigger = calloc(2, sizeof(int));
		random_timestep = malloc(2 * sizeof(int));

		random_timestep[0] = abs(rand()) % (offset);
		random_timestep[1] = abs(rand()) % (offset);

		wall_trigger[0] = 1;
		wall_trigger[1] = 1;

	}
	#endif

	double *tmp = NULL;
	for(int t = 0; t <= timesteps; t++){

		for(int i = 0; i < neighbours_count[rank]; ++i){
			extraction( collideField + Q_NUMBER * chunk_begin_offset[rank][neighbours_chunk_id[rank][i]],
									cpuDomain[rank][neighbours_chunk_id[rank][i]],
									sendBuffer[i],
									neighbours_dir[rank][i],
									neighbours_local_start_ext_x[rank][i],
									neighbours_local_start_ext_y[rank][i],
									neighbours_local_start_ext_z[rank][i],
									neighbours_local_end_ext_x[rank][i],
									neighbours_local_end_ext_y[rank][i],
									neighbours_local_end_ext_z[rank][i]);
		}
			swap( sendBuffer, readBuffer, rank);
		for(int i = 0; i < neighbours_count[rank]; ++i){
			injection( collideField + Q_NUMBER * chunk_begin_offset[rank][neighbours_chunk_id[rank][i]],
									cpuDomain[rank][neighbours_chunk_id[rank][i]],
									readBuffer[i],
									neighbours_dir[rank][i],
									neighbours_local_start_inj_x[rank][i],
									neighbours_local_start_inj_y[rank][i],
									neighbours_local_start_inj_z[rank][i],
									neighbours_local_end_inj_x[rank][i],
									neighbours_local_end_inj_y[rank][i],
									neighbours_local_end_inj_z[rank][i]);
		}
		for(int i = 0; i < chunk_count[rank]; ++i){
			treatBoundary( collideField + Q_NUMBER * chunk_begin_offset[rank][i], flagField + chunk_begin_offset[rank][i], velocityWall, &ref_density, cpuDomain[rank][i], inflow, pressure_in);
			doStreaming( collideField + Q_NUMBER * chunk_begin_offset[rank][i], streamField + Q_NUMBER * chunk_begin_offset[rank][i], flagField + chunk_begin_offset[rank][i], cpuDomain[rank][i]);
		}
		tmp = collideField;
		collideField = streamField;
		streamField = tmp;
		for(int i = 0; i < chunk_count[rank]; ++i){
			doCollision( collideField + Q_NUMBER * chunk_begin_offset[rank][i], flagField + chunk_begin_offset[rank][i], &tau, cpuDomain[rank][i] );
		}
		if ( t % timestepsPerPlotting == 0 ) {
			if(!rank){
				printf("Write vtk for time # %d \n", t);
			}
			for(int i = 0; i < chunk_count[rank]; ++i){

				writeVtkOutput( collideField + Q_NUMBER * chunk_begin_offset[rank][i], flagField + chunk_begin_offset[rank][i], "pics/simLB", t, cpuDomain[rank][i], chunk_id[rank][i] );
			}
		}

		#if defined(INCLUDE_RANDOM_DOORS)
		// Open and close fingers randomly.
		randomFingerOpenClose (t, rank, random_timestep, wall_trigger, flagField, offset);
		#endif
	}

	free((void *)collideField);
	free((void *)streamField);
	free((void *)flagField);

	free(inflow);
	free(pressure_in);

	for(int i = 0; i < neighbours_count[rank]; ++i){
		free(sendBuffer[i]);
		free(readBuffer[i]);
	}
	free(sendBuffer);
	free(readBuffer);

	finalizeMPI();

	return 0;
}
#endif
