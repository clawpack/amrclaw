#ifndef __BURGER_RIEMANN_SOLVER_H__
#define __BURGER_RIEMANN_SOLVER_H__

#include "real.H"
#include "common.H"
#include "staticArray.H"
#include <cmath>
#include <cassert>

#ifdef CUDA
// As a general rule a Riemann solver must have the following arguments:
// Input:
// - 2 cells, one left one right (equivalently up and down)
// - 2 sets of coefficients one for the left cell the other for the right, if needed
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
struct Burger_horizontal
{
    inline __device__ 
    void operator() (	real* q_left, real* q_right, real* u_left, real* u_right,			// input
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
    inline __device__ 
    void operator() (	real* q_left, real* q_right, real* u_left, real* u_right,			// input
            int row, int col, int numStates, int numWaves,
            real* wave, real* waveSpeeds)										// output
    {
        real theta = u_left[0];
        real a = 0.5f*sin(theta);

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, a*(q_left[0] + q_right[0]));

        setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, q_right[0] - q_left[0]);
    }
};

#endif
#endif
