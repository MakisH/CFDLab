#include "random_wings.h"


void randomFingerOpenClose(int t, int rank, int *random_timestep, int *wall_trigger,
                           int *flagField, int offset)
{

  srand(time(NULL)*rank*100);

  // CLOSING / OPENING WALL
  // Only do that for rank \iselement [4, 7]
  if (rank >= 4 && rank <= 7){
      for (int finger = 0; finger < neighbours_count[rank]; finger++){
        printf ("RANK: %d, WALL AT TIMESTEP: %d \n", rank, random_timestep[finger]);
          if (t == random_timestep[finger] && wall_trigger[finger] == 1){

            int *address = flagField + chunk_begin_offset[rank][finger];
            int xlen2 = cpuDomain[rank][finger][0] + 2;
            int ylen2 = cpuDomain[rank][finger][1] + 2;
            int xylen2 = xlen2 * ylen2;

            int z = 1; // We change flagField on z=1 layer.
            int y;

            // Should I put wall on lower finger or the upper one?
            if (rank == 4 || rank == 5) y = ylen2 - 15;
            else y = 14;

              for (int x = 0; x < xlen2; x++){
                address[x + xlen2 * y + xylen2 * z] = NO_SLIP;
              }

            printf("I HAVE SET UP THE WALL! -> at timestep = %d\n", t);
            wall_trigger[finger] = 0;
        } else if (t == random_timestep[finger] + offset){

          int *address = flagField + chunk_begin_offset[rank][finger];
          int xlen2 = cpuDomain[rank][finger][0] + 2;
          int ylen2 = cpuDomain[rank][finger][1] + 2;
          int xylen2 = xlen2 * ylen2;

          int z = 1; // We change flagField on z=1 layer.
          int y;

            // Should I put wall on lower finger or the upper one?
            if (rank == 4 || rank == 5) y = ylen2 - 15;
            else y = 14;

            for (int x = 1; x < xlen2 - 1; x++){
              address[x + xlen2 * y + xylen2 * z] = FLUID;
            }
          wall_trigger[finger] = 1;
          random_timestep[finger] = abs((rand())) % (offset) + t; // recalculate the time
          printf("NEW RANDOM TIMESTEP IS %d\n", random_timestep[finger]);
      }
    }
  }
}
