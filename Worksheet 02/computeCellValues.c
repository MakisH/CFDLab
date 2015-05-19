#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
	// 10% speedup(+15sec)
		*density =	currentCell[0] + currentCell[1] + currentCell[2] + currentCell[3] + currentCell[4] + currentCell[5] + currentCell[6] + currentCell[7] + 
					currentCell[8] + currentCell[9] + currentCell[10] + currentCell[11] + currentCell[12] + currentCell[13] + currentCell[14] + currentCell[15] + 
					currentCell[16] + currentCell[17] + currentCell[18];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){

	velocity[0] = (-currentCell[1]+currentCell[3]-currentCell[5]+currentCell[7]-currentCell[8]+currentCell[10]-currentCell[11]+currentCell[13]-currentCell[15]+currentCell[17])/ *density;
	velocity[1] = (-currentCell[0]+currentCell[4]-currentCell[5]-currentCell[6]-currentCell[7]+currentCell[11]+currentCell[12]+currentCell[13]-currentCell[14]+currentCell[18])/ *density;
	velocity[2] = (-currentCell[0]-currentCell[1]-currentCell[2]-currentCell[3]-currentCell[4]+currentCell[14]+currentCell[15]+currentCell[16]+currentCell[17]+currentCell[18])/ *density;

}

void computeFeq(const double * const density, const double * const velocity, double *feq){
	// precalculate inner product 1 - u*u / (2*c_s^2)
	double OneMinusu_u_2c2 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]) * 0.5 / C_S_sq;

		// precalculate inner product c*u/c_s^2
		double c_u_c2;
		c_u_c2 = ( -velocity[1] -velocity[2]) / C_S_sq;
		feq[0] = LATTICEWEIGHTS[0] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = ( -velocity[0] -velocity[2]) / C_S_sq;
		feq[1] = LATTICEWEIGHTS[1] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = (				-velocity[2]) / C_S_sq;
		feq[2] = LATTICEWEIGHTS[2] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = (  velocity[0] -velocity[2]) / C_S_sq;
		feq[3] = LATTICEWEIGHTS[3] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = (  velocity[1] -velocity[2]) / C_S_sq;
		feq[4] = LATTICEWEIGHTS[4] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );

		c_u_c2 = ( -velocity[0] -velocity[1]) / C_S_sq;
		feq[5] = LATTICEWEIGHTS[5] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = (				-velocity[1]) / C_S_sq;
		feq[6] = LATTICEWEIGHTS[6] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		c_u_c2 = (  velocity[0] -velocity[1]) / C_S_sq;
		feq[7] = LATTICEWEIGHTS[7] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );

		c_u_c2 = ( -velocity[0]				) / C_S_sq;
		feq[8] = LATTICEWEIGHTS[8] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 )	 );
		feq[9] = LATTICEWEIGHTS[9] * *density * (  OneMinusu_u_2c2									 );
		c_u_c2 = (  velocity[0]				) / C_S_sq;
		feq[10] = LATTICEWEIGHTS[10] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );
		c_u_c2 = ( -velocity[0] +velocity[1]) / C_S_sq;
		feq[11] = LATTICEWEIGHTS[11] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );
		c_u_c2 = (				 velocity[1]) / C_S_sq;
		feq[12] = LATTICEWEIGHTS[12] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );
		c_u_c2 = (  velocity[0] +velocity[1]) / C_S_sq;
		feq[13] = LATTICEWEIGHTS[13] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );

		c_u_c2 = ( -velocity[1] +velocity[2]) / C_S_sq;
		feq[14] = LATTICEWEIGHTS[14] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );

		c_u_c2 = ( -velocity[0] +velocity[2]) / C_S_sq;
		feq[15] = LATTICEWEIGHTS[15] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );

		c_u_c2 = (				 velocity[2]) / C_S_sq;
		feq[16] = LATTICEWEIGHTS[16] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );

		c_u_c2 = (  velocity[0] +velocity[2]) / C_S_sq;
		feq[17] = LATTICEWEIGHTS[17] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );

		c_u_c2 = (  velocity[1] +velocity[2]) / C_S_sq;
		feq[18] = LATTICEWEIGHTS[18] * *density * (  OneMinusu_u_2c2 + c_u_c2 * (1 + c_u_c2 * 0.5 ) );
}
