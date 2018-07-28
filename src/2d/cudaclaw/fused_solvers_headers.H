#ifndef __FUSED_RIEMANN_LIMITER_HEADERS_H__
#define __FUSED_RIEMANN_LIMITER_HEADERS_H__

#include "real.H"
#include "common.H"
#include "staticArray.H"
#include <cmath>
#include <cassert>

#define epsilon 0.00000001

using std::sqrt;

// constexpr real epsilon = 0.00000001;

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
inline __device__ 
int general_indexing(int a, int sizeA, int b, int sizeB, int c, int sizeC, int d, int sizeD)
{
    return (a*sizeB*sizeC*sizeD + b*sizeC*sizeD + c*sizeD + d);
}
inline __device__ 
int getIndex_SharedWave(int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth)
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

inline __device__ 
real &getSharedWave(real* sharedData, const int row, const int col, const int wave, const int state, const int numWaves, const int numStates, const int blockHeight, const int blockWidth)
{
    return sharedData[getIndex_SharedWave(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)];
}
inline __device__ 
void setSharedWave(real* sharedData, int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth, real newVal)
{
    sharedData[getIndex_SharedWave(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)] = newVal;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////  Wave Speeds  //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline __device__ 
int getIndex_SharedSpeed(int row, int col, int wave, int numWaves, int blockHeight, int blockWidth)
{
    // return (row*numWaves*blockWidth + wave*blockWidth + col);		// seems like the best
    //
    //return ( row*blockWidth*numWaves + col*numWaves + wave);
    //return (col*numWaves*blockHeight + wave*blockHeight + row);
    //return (col*blockHeight*numWaves + row*numWaves + wave);
    return (wave*blockHeight*blockWidth + row*blockWidth + col);
    //return (wave*blockWidth*blockHeight + col*blockHeight + row);
}
inline __device__ 
real &getSharedSpeed(real* sharedData, const int row, const int col, const int wave, const int numWaves, const int blockHeight, const int blockWidth)
{
    return sharedData[getIndex_SharedSpeed(row, col, wave, numWaves, blockHeight, blockWidth)];
}
inline __device__ 
void setSharedSpeed(real* sharedData, int row, int col, int wave, int numWaves, int blockHeight, int blockWidth, real newVal)
{
    sharedData[getIndex_SharedSpeed(row, col, wave, numWaves, blockHeight, blockWidth)] = newVal;
}


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

    if (main_wave_norm_square < epsilon)
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

    if (main_wave_norm_square < epsilon)
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

    if (main_wave_norm_square < epsilon)
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

    if (main_wave_norm_square < epsilon)
        return (real)1.0f;
    return phi(aux_wave_dot_main_wave/main_wave_norm_square);
}
////////////////////////////////////////////////////////////////////

struct limiter_none
{
    __device__ real operator() (real theta)
    {
        return claw_zero;
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
        return fmax(claw_zero, fmin(((real)1.0f+theta)/(real)2.0f, minimum));
    }
};
struct limiter_superbee
{
    __device__ real operator() (real theta)
    {
        real maximum = fmax(claw_zero, fmin((real)1.0f,(real)2.0f*theta));
        return fmax(maximum, fmin((real)2.0f,theta));
    }
};
struct limiter_VanLeer
{
    __forceinline__ __device__ 
    real operator() (real theta)
    {
        real absTheta = fabs(theta);
        return (theta + absTheta) / ((real)1.0f + absTheta);
    }
};


#endif
