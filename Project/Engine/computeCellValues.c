#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
		*density = *currentCell;
		for(unsigned int i = 1; i < Q_NUMBER; ++i)
			*density += currentCell[i];
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
		// Update velocity in x,y,z direction as inner product density*

	velocity[0] = ( -currentCell[1] +currentCell[3] -currentCell[5] +currentCell[7] -currentCell[8]
					+currentCell[10]-currentCell[11]+currentCell[13]-currentCell[15]+currentCell[17])/ *density;
	velocity[1] = ( -currentCell[0] +currentCell[4] -currentCell[5] -currentCell[6] -currentCell[7]
					+currentCell[11]+currentCell[12]+currentCell[13]-currentCell[14]+currentCell[18])/ *density;
	velocity[2] = ( -currentCell[0] -currentCell[1] -currentCell[2] -currentCell[3] -currentCell[4]
					+currentCell[14]+currentCell[15]+currentCell[16]+currentCell[17]+currentCell[18])/ *density;

	//	// original solution ... it is 20% slower with -O3 in comparison to above
	//	// We can unroll loop => many multiplications saved, since Lattice velocities are 0 or +/-1 => less memory access and less FLOPS
	//	velocity[0] = LATTICEVELOCITIES[0][0] * currentCell[0];
	//	velocity[1] = LATTICEVELOCITIES[0][1] * currentCell[0];
	//	velocity[2] = LATTICEVELOCITIES[0][2] * currentCell[0];
	//for(unsigned int i = 1; i < Q_NUMBER; ++i){
	//	velocity[0] += LATTICEVELOCITIES[i][0] * currentCell[i];
	//	velocity[1] += LATTICEVELOCITIES[i][1] * currentCell[i];
	//	velocity[2] += LATTICEVELOCITIES[i][2] * currentCell[i];
	//}
	//	velocity[0] /= *density;
	//	velocity[1] /= *density;
	//	velocity[2] /= *density;
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
	// precalculate inner product 1 - u*u / (2*c_s^2)
		// precalculate inner product c*u/c_s^2
		double vel0 = velocity[0] / C_S_sq;
		double vel1 = velocity[1] / C_S_sq;
		double vel2 = velocity[2] / C_S_sq;
		double OneMinusu_u_2c2 = 1 - (velocity[0] * vel0 + velocity[1] * vel1 + velocity[2] * vel2) * 0.5;

		double local_density = *density / 36;
		double c_u_c2;
		c_u_c2 = ( -vel1 -vel2);
		feq[0] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = ( -vel0 -vel2);
		feq[1] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (		 -vel2);
		feq[2] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel0 -vel2);
		feq[3] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel1 -vel2);
		feq[4] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);

		c_u_c2 = ( -vel0 -vel1);
		feq[5] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (				-vel1)  ;
		feq[6] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel0 -vel1);
		feq[7] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = ( -vel0	);
		feq[8] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		feq[9] = 12 * local_density * (  OneMinusu_u_2c2							 );
		c_u_c2 = (  vel0	);
		feq[10] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = ( -vel0 +vel1);
		feq[11] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (		  vel1);
		feq[12] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel0 +vel1);
		feq[13] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);

		c_u_c2 = ( -vel1 +vel2);
		feq[14] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = ( -vel0 +vel2);
		feq[15] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (		  vel2);
		feq[16] = 2 * local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel0 +vel2);
		feq[17] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);
		c_u_c2 = (  vel1 +vel2);
		feq[18] = local_density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2*c_u_c2 * 0.5);

	//// original unrolled solution  ... 20% slower with -O3 in comparison to above
	//// precalculate inner product 1 - u*u / (2*c_s^2)
	//double OneMinusu_u_2c2 = 1 - (velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]) * 0.5 / C_S_sq;
	//	double vel0 = velocity[0] / C_S_sq;
	//	double vel1 = velocity[1] / C_S_sq;
	//	double vel2 = velocity[2] / C_S_sq;
	//for (int i = 0; i < Q_NUMBER; ++i) {
	//	// precalculate inner product c*u/c_s^2
	//		double c_u_c2;

	//	c_u_c2 = LATTICEVELOCITIES[i][0] * vel0 + LATTICEVELOCITIES[i][1] * vel1 + LATTICEVELOCITIES[i][2] * vel2;
	//	// calculate the equilibrium distribution element
	//	feq[i] = LATTICEWEIGHTS[i] * *density * (  OneMinusu_u_2c2 + c_u_c2 + c_u_c2 * c_u_c2 * 0.5 );
	//}
}
