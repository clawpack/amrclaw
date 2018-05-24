#ifndef __FUSED_RIEMANN_LIMITER_HEADERS_H__
#define __FUSED_RIEMANN_LIMITER_HEADERS_H__

#define EPSILON (real)0.00000001f	// flying constant, however this is our epsilon
#include "real.H"
#include "common.H"
#include <cmath>
#include <cassert>

using std::sqrt;

// The Waves and Wave Speeds lie in the shared memory.
// The concern at this stage is not coalescing and alignment (not as it would be in global)
// but bank conflicts, different schemes are likely to yield different performances
// Note however that the Riemann solver depends on the function below as they call them 
// with their proper arguments to access the correct location for the speeds and waves.
// The accessor method take all dimension information and use all but one. This is done
// to make code changes easy when a different accessing scheme is to be used.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////      Waves    //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline __device__ int general_indexing(int a, int sizeA, int b, int sizeB, int c, int sizeC, int d, int sizeD)
{
	return (a*sizeB*sizeC*sizeD + b*sizeC*sizeD + c*sizeD + d);
}
inline __device__ int getIndex_SharedWave(int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth)
{
	//return ( general_indexing(row, blockHeight, col, blockWidth, wave, numWaves, state, numStates) ); // 2.55, 22r 23r
	//return ( general_indexing(row, blockHeight, col, blockWidth, state, numStates, wave, numWaves) ); // 2.55, 22r 23r
	//return ( general_indexing(row, blockHeight, state, numStates, col, blockWidth, wave, numWaves) ); // 2.67
	//return ( general_indexing(row, blockHeight, state, numStates, wave, numWaves, col, blockWidth) ); // 2.88
	//return ( general_indexing(row, blockHeight, wave, numWaves, col, blockWidth, state, numStates) ); // 2.67
	return ( general_indexing(row, blockHeight, wave, numWaves, state, numStates, col, blockWidth) ); // 2.88, 24r 26r <------- this vs
	//return ( (((row*numWaves) + wave)*numStates + state)*blockWidth + col);
	//return ( row*blockWidth*numWaves*numStates + col*numWaves*numStates + wave*numStates + state);	//					<------- this
	//
	//return ( general_indexing(col, blockWidth, row, blockHeight, wave, numWaves, state, numStates) );
	//return ( general_indexing(col, blockWidth, row, blockHeight, state, numStates, wave, numWaves) );
	//return ( general_indexing(col, blockWidth, wave, numWaves, row, blockHeight, state, numStates) );
	//return ( general_indexing(col, blockWidth, wave, numWaves, state, numStates, row, blockHeight) );
	//return ( general_indexing(col, blockWidth, state, numStates, wave, numWaves, row, blockHeight) );
	//return ( general_indexing(col, blockWidth, state, numStates, row, blockHeight, wave, numWaves) );
	//
	//return ( general_indexing(wave, numWaves, row, blockHeight, col, blockWidth, state, numStates) );
	//return ( general_indexing(wave, numWaves, row, blockHeight, state, numStates, col, blockWidth) );
	//return ( general_indexing(wave, numWaves, col, blockWidth, row, blockHeight, state, numStates) );
	//return ( general_indexing(wave, numWaves, col, blockWidth, state, numStates, row, blockHeight) );
	//return ( general_indexing(wave, numWaves, state, numStates, col, blockWidth, row, blockHeight) );
	//return ( general_indexing(wave, numWaves, state, numStates, row, blockHeight, col, blockWidth) );

	//return ( general_indexing(state, numStates, col, blockWidth, row, blockHeight, wave, numWaves) );	// generally bad
	//return ( general_indexing(state, numStates, col, blockWidth, wave, numWaves, row, blockHeight) );
	//return ( general_indexing(state, numStates, row, blockHeight, col, blockWidth, wave, numWaves) );
	//return ( general_indexing(state, numStates, row, blockHeight, wave, numWaves, col, blockWidth) );
	//return ( general_indexing(state, numStates, wave, numWaves, row, blockHeight, col, blockWidth) );
	//return ( general_indexing(state, numStates, wave, numWaves, col, blockWidth, row, blockHeight) );
}
inline __device__ real &getSharedWave(real* sharedData, int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth)
{
	return sharedData[getIndex_SharedWave(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)];
}
inline __device__ void setSharedWave(real* sharedData, int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth, real newVal)
{
	sharedData[getIndex_SharedWave(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)] = newVal;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////  Wave Speeds  //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline __device__ int getIndex_SharedSpeed(int row, int col, int wave, int numWaves, int blockHeight, int blockWidth)
{
        // TODO: this arranging is weird?
	return (row*numWaves*blockWidth + wave*blockWidth + col);		// seems like the best
	//return ( row*blockWidth*numWaves + col*numWaves + wave);
	//return (col*numWaves*blockHeight + wave*blockHeight + row);
	//return (col*blockHeight*numWaves + row*numWaves + wave);
	//return (wave*blockHeight*blockWidth + row*blockWidth + col);
	//return (wave*blockWidth*blockHeight + col*blockHeight + row);
}
inline __device__ real &getSharedSpeed(real* sharedData, int row, int col, int wave, int numWaves, int blockHeight, int blockWidth)
{
	return sharedData[getIndex_SharedSpeed(row, col, wave, numWaves, blockHeight, blockWidth)];
}
inline __device__ void setSharedSpeed(real* sharedData, int row, int col, int wave, int numWaves, int blockHeight, int blockWidth, real newVal)
{
	sharedData[getIndex_SharedSpeed(row, col, wave, numWaves, blockHeight, blockWidth)] = newVal;
}


// As a general rule a Riemann solver must have the following arguments:
// Input:
// - 2 cells, one left one right (equivalently up and down)
// - 2 sets of coefficients one for the left cell the other for the right
// - The position (row, column) of the interface
// - The number of states
// - The number of waves
// Output:
// - a location for storing the set of waves, base address 
// - a location for storing the set of wave speeds
//
// The input will come from the registers (read from global into registers)
// and output will be to shared memory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////   Riemann Solvers   ////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct acoustics_homo_2d_horizontal
{
        // u_left and u_right are empty and not used here
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
                assert(NUMCOEFF == 0); // this is the case for homogeneous acoustics
		real c = sqrt(BULK/RHO);
		real z = c*RHO;
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c);

		real alpha1 = ( q_left[0] - q_right[0] + z*(q_right[1] - q_left[1])) / (2*z);
		real alpha2 = ( q_right[0] - q_left[0] + z*(q_right[1] - q_left[1])) / (2*z);
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
	}
};
struct acoustics_horizontal
{
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
		real rho_l = u_left[0];
		real bulk_l = u_left[1];

		real rho_r = u_right[0];
		real bulk_r = u_right[1];

		real c_l = sqrt(bulk_l/rho_l);	// sound speed
		real z_l = c_l*rho_l;			// impedance

		real c_r = sqrt(bulk_r/rho_r);
		real z_r = c_r*rho_r;
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c_l);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c_l);

		real alpha1 = ( q_left[0] - q_right[0] + z_r*(q_right[1] - q_left[1])) / (z_l+z_r);	// ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
		real alpha2 = ( q_right[0] - q_left[0] + z_l*(q_right[1] - q_left[1])) / (z_l+z_r);	// (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z_l);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z_r);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
	}
};
struct acoustics_homo_2d_vertical
{
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
                assert(NUMCOEFF == 0); // this is the case for homogeneous acoustics
		real c = sqrt(BULK/RHO);
		real z = c*RHO;

		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c);
		
		real alpha1 = ( q_left[0] - q_right[0] + z*(q_right[2] - q_left[2])) / (2*z);	// ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
		real alpha2 = ( q_right[0] - q_left[0] + z*(q_right[2] - q_left[2])) / (2*z);	// (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
	}
};
struct acoustics_vertical
{
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
		real rho_l = u_left[0];
		real bulk_l = u_left[1];

		real rho_r = u_right[0];
		real bulk_r = u_right[1];

		real c_l = sqrt(bulk_l/rho_l);	// sound speed
		real z_l = c_l*rho_l;			// impedance

		real c_r = sqrt(bulk_r/rho_r);
		real z_r = c_r*rho_r;
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c_l);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c_l);
		
		real alpha1 = ( q_left[0] - q_right[0] + z_r*(q_right[2] - q_left[2])) / (z_l+z_r);	// ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
		real alpha2 = ( q_right[0] - q_left[0] + z_l*(q_right[2] - q_left[2])) / (z_l+z_r);	// (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z_l);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z_r);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
	}
};
struct shallow_water_horizontal
{
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
		// try rearanging, or using ul and ur, vl and vr variables instead of dividing every time!
		real g = 9.8f;

		real h_left   = q_left[0];
		real hu_left  = q_left[1];
		real hv_left  = q_left[2];
		
		real h_right  = q_right[0];
		real hu_right = q_right[1];
		real hv_right = q_right[2];

		real h_bar = 0.5f*(h_left + h_right);
		real sqrt_h_left  = sqrt(h_left);
		real sqrt_h_right = sqrt(h_right);

		//real sum_sqrt_hleft_hright = sqrt_h_left + sqrt_h_right;

		real u_hat = ((hu_left/sqrt_h_left)+(hu_right/sqrt_h_right))/(sqrt_h_left + sqrt_h_right);
		real v_hat = ((hv_left/sqrt_h_left)+(hv_right/sqrt_h_right))/(sqrt_h_left + sqrt_h_right);

		//real a_left  = h_left*u_hat - hu_left;
		//real a_right = hu_right - h_right*u_hat;
		//real v_hat   = (a_left*(hv_left/h_left) + a_right*(hv_right/h_right)) / (a_left+a_right);

		real c_hat = sqrt(g*h_bar);
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, u_hat - c_hat);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x, u_hat);
		setSharedSpeed(waveSpeeds, row, col, 2, numWaves, blockDim.y, blockDim.x, u_hat + c_hat);
		
		real alpha1 = 0.5f*((u_hat + c_hat)*(h_right - h_left) - (hu_right - hu_left))/c_hat;
		real alpha2 = (hv_right-hv_left)-v_hat*(h_right - h_left);
		real alpha3 = 0.5f*((c_hat - u_hat)*(h_right - h_left) + (hu_right - hu_left))/c_hat;

		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha1*(u_hat - c_hat));
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha1*v_hat);
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
		
		setSharedWave(wave, row, col, 2, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha3);
		setSharedWave(wave, row, col, 2, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha3*(u_hat + c_hat));
		setSharedWave(wave, row, col, 2, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha3*v_hat);
	}
};
struct shallow_water_vertical
{
	__device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,			// input
								int row, int col, int numStates, int numWaves,
								real* wave, real* waveSpeeds)										// output
	{
		// try rearanging, or using ul and ur, vl and vr variables instead of dividing every time!
		real g = 9.8f;

		real h_left   = q_left[0];
		real hu_left  = q_left[1];
		real hv_left  = q_left[2];
		
		real h_right  = q_right[0];
		real hu_right = q_right[1];
		real hv_right = q_right[2];

		real h_bar = 0.5f*(h_left + h_right);
		real sqrt_h_left  = sqrt(h_left);
		real sqrt_h_right = sqrt(h_right);

		//real sum_sqrt_hleft_hright = sqrt_h_left + sqrt_h_right;

		//real u_hat = ((hu_left/sqrt_h_left)+(hu_right/sqrt_h_right))/(sqrt_h_left + sqrt_h_right);
		//real v_hat = ((hv_left/sqrt_h_left)+(hv_right/sqrt_h_right))/(sqrt_h_left + sqrt_h_right);
		//real a_left  = h_left*u_hat - hu_left;
		//real a_right = hu_right - h_right*u_hat;
		//real v_hat   = (a_left*(hv_left/h_left) + a_right*(hv_right/h_right)) / (a_left+a_right);

		real u_hat = (hu_left/sqrt_h_left + hu_right/sqrt_h_right)/(sqrt_h_left + sqrt_h_right);
		real v_hat = (hv_left/sqrt_h_left + hv_right/sqrt_h_right)/(sqrt_h_left + sqrt_h_right);


		real c_hat = sqrt(g*h_bar);

		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, v_hat - c_hat);
		setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x, v_hat);
		setSharedSpeed(waveSpeeds, row, col, 2, numWaves, blockDim.y, blockDim.x, v_hat + c_hat);

		real alpha1 = 0.5f*((v_hat + c_hat)*(h_right - h_left) - (hv_right - hv_left))/c_hat;
		real alpha2 = -(hu_right-hu_left) + u_hat*(h_right - h_left);
		real alpha3 = 0.5f*((c_hat - v_hat)*(h_right - h_left) + (hv_right - hv_left))/c_hat;

		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
		setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha1*u_hat);
		setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha1*(v_hat - c_hat));
		
		setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, -alpha2);
		setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, 0.0f);
		
		setSharedWave(wave, row, col, 2, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha3);
		setSharedWave(wave, row, col, 2, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha3*u_hat);
		setSharedWave(wave, row, col, 2, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha3*(v_hat + c_hat));
	}
};
struct Burger_horizontal
{
	inline __device__ void operator() (	real* q_left, real* q_right, real* u_left, real* u_right,			// input
										int row, int col, int numStates, int numWaves,
										real* wave, real* waveSpeeds)										// output
	{
		real theta = u_left[0];
		real a = 0.5f*cos(theta);
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, q_right[0] - q_left[0]);
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, a*(q_left[0] + q_right[0]));
	}
};
struct Burger_vertical
{
	inline __device__ void operator() (	real* q_left, real* q_right, real* u_left, real* u_right,			// input
										int row, int col, int numStates, int numWaves,
										real* wave, real* waveSpeeds)										// output
	{
		real theta = u_left[0];
		real a = 0.5f*sin(theta);
		
		setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, a*(q_left[0] + q_right[0]));
		
		setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, q_right[0] - q_left[0]);
	}
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////   Limiters   //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__device__ real limiting_shared_h_l (Limiter phi, real* waves, int row, int col, int waveNum)
{
	real main_wave = getSharedWave(waves, row, col,   waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real aux_wave  = getSharedWave(waves, row, col+1, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real main_wave_norm_square = main_wave*main_wave;
	real aux_wave_dot_main_wave = aux_wave*main_wave;
	#pragma unroll
	for (int i = 1; i < numStates; i++)
	{
		main_wave = getSharedWave(waves, row, col,   waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);
		aux_wave  = getSharedWave(waves, row, col+1, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);

		main_wave_norm_square += main_wave*main_wave;
		aux_wave_dot_main_wave += aux_wave*main_wave;
	}

	if (main_wave_norm_square < EPSILON)
		return (real)1.0f;
	return phi(aux_wave_dot_main_wave/main_wave_norm_square);
}
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__device__ real limiting_shared_h_r (Limiter phi, real* waves, int row, int col, int waveNum)
{
	real main_wave = getSharedWave(waves, row, col,   waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real aux_wave  = getSharedWave(waves, row, col-1, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real main_wave_norm_square = main_wave*main_wave;
	real aux_wave_dot_main_wave = aux_wave*main_wave;
	#pragma unroll
	for (int i = 1; i < numStates; i++)
	{
		main_wave = getSharedWave(waves, row, col,   waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);
		aux_wave  = getSharedWave(waves, row, col-1, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);

		main_wave_norm_square += main_wave*main_wave;
		aux_wave_dot_main_wave += aux_wave*main_wave;
	}

	if (main_wave_norm_square < EPSILON)
		return (real)1.0f;
	return phi(aux_wave_dot_main_wave/main_wave_norm_square);
}
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__device__ real limiting_shared_v_u (Limiter phi, real* waves, int row, int col, int waveNum)
{
	real main_wave = getSharedWave(waves, row,   col, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real aux_wave  = getSharedWave(waves, row-1, col, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real main_wave_norm_square = main_wave*main_wave;
	real aux_wave_dot_main_wave = aux_wave*main_wave;
	#pragma unroll
	for (int i = 1; i < numStates; i++)
	{
		main_wave = getSharedWave(waves, row,   col, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);
		aux_wave  = getSharedWave(waves, row-1, col, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);

		main_wave_norm_square += main_wave*main_wave;
		aux_wave_dot_main_wave += aux_wave*main_wave;
	}

	if (main_wave_norm_square < EPSILON)
		return (real)1.0f;
	return phi(aux_wave_dot_main_wave/main_wave_norm_square);
}
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__device__ real limiting_shared_v_d (Limiter phi, real* waves, int row, int col, int waveNum)
{
	real main_wave = getSharedWave(waves, row,   col, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real aux_wave  = getSharedWave(waves, row+1, col, waveNum, 0, numWaves, numStates, blockSizeY, blockSizeX);
	real main_wave_norm_square = main_wave*main_wave;
	real aux_wave_dot_main_wave = aux_wave*main_wave;
	#pragma unroll
	for (int i = 1; i < numStates; i++)
	{
		main_wave = getSharedWave(waves, row,   col, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);
		aux_wave  = getSharedWave(waves, row+1, col, waveNum, i, numWaves, numStates, blockSizeY, blockSizeX);

		main_wave_norm_square += main_wave*main_wave;
		aux_wave_dot_main_wave += aux_wave*main_wave;
	}

	if (main_wave_norm_square < EPSILON)
		return (real)1.0f;
	return phi(aux_wave_dot_main_wave/main_wave_norm_square);
}
////////////////////////////////////////////////////////////////////

struct limiter_none
{
	__device__ real operator() (real theta)
	{
		return (real)0.0f;
	}
};
struct limiter_LaxWendroff
{
	__device__ real operator() (real theta)
	{
		return (real)1.0f;
	}
};
struct limiter_MC
{
	__device__ real operator() (real theta)
	{
		real minimum = fmin((real)2.0f, (real)2.0f*theta);
		return fmax((real)0.0f, fmin(((real)1.0f+theta)/(real)2.0f, minimum));
	}
};
struct limiter_superbee
{
	__device__ real operator() (real theta)
	{
		real maximum = fmax((real)0.0f, fmin((real)1.0f,(real)2.0f*theta));
		return fmax(maximum, fmin((real)2.0f,theta));
	}
};
struct limiter_VanLeer
{
	__device__ real operator() (real theta)
	{
		real absTheta = fabs(theta);
		return (theta + absTheta) / ((real)1.0f + absTheta);
	}
};

#endif
