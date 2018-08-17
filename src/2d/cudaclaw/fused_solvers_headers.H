#ifndef __FUSED_SOLVERS_HEADERS_H__
#define __FUSED_SOLVERS_HEADERS_H__

#include "real.H"
#include "common.H"
#include "arrayIndex.H"
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
real &getSharedWave(real* sharedData, const int row, const int col, const int wave, const int state, const int numWaves, const int numStates, const int blockHeight, const int blockWidth)
{
    return sharedData[GET_INDEX_SHARED_WAVE_0INDEX(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)];
}
inline __device__ 
void setSharedWave(real* sharedData, int row, int col, int wave, int state, int numWaves, int numStates, int blockHeight, int blockWidth, real newVal)
{
    sharedData[GET_INDEX_SHARED_WAVE_0INDEX(row, col, wave, state, numWaves, numStates, blockHeight, blockWidth)] = newVal;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////  Wave Speeds  //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline __device__ 
real &getSharedSpeed(real* sharedData, const int row, const int col, const int wave, const int numWaves, const int blockHeight, const int blockWidth)
{
    return sharedData[GET_INDEX_SHARED_SPEED_0INDEX(row, col, wave, numWaves, blockHeight, blockWidth)];
}
inline __device__ 
void setSharedSpeed(real* sharedData, int row, int col, int wave, int numWaves, int blockHeight, int blockWidth, real newVal)
{
    sharedData[GET_INDEX_SHARED_SPEED_0INDEX(row, col, wave, numWaves, blockHeight, blockWidth)] = newVal;
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
    __device__ 
    real operator() (real theta)
    {
        real absTheta = fabs(theta);
        return (theta + absTheta) / ((real)1.0f + absTheta);
    }
};


#endif
