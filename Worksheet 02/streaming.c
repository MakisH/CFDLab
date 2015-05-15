#include "streaming.h"
#include "LBDefinitions.h"

#include <stdio.h>
#include <string.h>
// used for memcpy ?

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){ 
  for(int z = 1; z <= xlength; z += 2 * xlength + 6){ // no FLUID cells on the boundary! => interval 1-xlength
	  unsigned int col_temp = xlength*xlength*Q_NUMBER; 
	  for(int y = 1; y <= xlength; y += 2){
		  // row-wise copy 
		memcpy(streamField + Q_NUMBER * (y * xlength + z * col_temp,collideField) + Q_NUMBER * (y * xlength + z * col_temp), 
		 // row-wise copy 
		sizeof(double) * xlength * Q_NUMBER*) ;
	  }
  }
}

