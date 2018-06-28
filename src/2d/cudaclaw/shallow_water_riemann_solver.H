#ifndef __SHALLOW_WATER_RIEMANN_SOLVER_H__
#define __SHALLOW_WATER_RIEMANN_SOLVER_H__

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
struct shallow_water_horizontal
{
    __device__ void 
    operator() (real* q_left, real* q_right, real* u_left, real* u_right,	
            int row, int col, int numStates, int numWaves,
            real* wave, real* waveSpeeds)	
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

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha2);

        setSharedWave(wave, row, col, 2, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha3);
        setSharedWave(wave, row, col, 2, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha3*(u_hat + c_hat));
        setSharedWave(wave, row, col, 2, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha3*v_hat);
    }
};

struct shallow_water_vertical
{
    __device__ void 
    operator() (real* q_left, real* q_right, real* u_left, real* u_right,
            int row, int col, int numStates, int numWaves,
            real* wave, real* waveSpeeds)	
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

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, -alpha2);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);

        setSharedWave(wave, row, col, 2, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha3);
        setSharedWave(wave, row, col, 2, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha3*u_hat);
        setSharedWave(wave, row, col, 2, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha3*(v_hat + c_hat));
    }
};

#endif
#endif
