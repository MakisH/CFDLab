#include "helper.h"
#include "mpi_helper.h"
#include "LBDefinitions.h"
#include "time.h"

void randomFingerOpenClose(int t, int rank, int *random_timestep, int *wall_trigger,
                           int *flagField, int offset);
