#ifndef __FUSED_SOLVERS_H__
#define __FUSED_SOLVERS_H__


#include "common.H"
#include "real.H"
#include "params.H"
#include "fused_solvers_headers.H"
#include "../clawpack_CUDA_helper.H"
#include "reduce_Max.H"
#include "clawpack_device.H"
#include "sharedMemory.H"
#include "staticArray.H"
#include "acoustics_riemann_solver.H"

#include <cmath>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <math_functions.h>
#include <math_constants.h>
#include <cstdio>

#if __CUDA_ARCH__ >= 600
    #define KERNEL_MIN_BLOCKS 2
#else
    #define KERNEL_MIN_BLOCKS 1
#endif

#define DEBUG

__device__ static double atomicMax(volatile double* address, double val)
{
    unsigned long long* address_as_ull = (unsigned long long*) address;
    unsigned long long old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_ull, assumed,
            __double_as_longlong(::fmax(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}


template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY>
inline __device__ void first_order_update(volatile real* waves, volatile real* waveSpeeds, real* amdq, real* apdq, const int row, const int col)
{
#pragma unroll
    for (int w = 0; w < numWaves; w++)
    {
        real waveSpeed =  getSharedSpeed(waveSpeeds, row, col, w, numWaves, blockSizeY, blockSizeX);
        if (waveSpeed < claw_zero)
        {
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                amdq[k] +=	waveSpeed * wave_state;
            }
        }
        else
        {
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                apdq[k] +=	waveSpeed * wave_state;
            }
        }
    }
}



/*
 * This regular wave approach does not handle capacity
 * In case where capacity variable exists, use fwave approach below.
 */
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__forceinline__ __device__ 
void first_second_order_update_horizontal(pdeParam &param, volatile real* waves, volatile real* waveSpeeds, real* amdq, real* apdq, const int row, const int col, const Limiter phi)
{
    real F_tilde[numStates];
    for (int k = 0; k < numStates; k++){
        F_tilde[k] = claw_zero;
    }
#pragma unroll
    for (int w = 0; w < numWaves; w++)
    {
        real waveSpeed =  getSharedSpeed(waveSpeeds, row, col, w, numWaves, blockSizeY, blockSizeX);

        if (waveSpeed < claw_zero)
        {
            real limiting_factor = limiting_shared_h_l<numStates, numWaves, blockSizeX, blockSizeY>(phi, waves, row, col, w);

#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                amdq[k] += waveSpeed * wave_state;
                F_tilde[k] += -0.5*waveSpeed*(1+(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
            }
        }
        else
        {
            real limiting_factor = limiting_shared_h_r<numStates, numWaves, blockSizeX, blockSizeY>(phi, waves, row, col, w);

#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                apdq[k] += waveSpeed * wave_state;
                F_tilde[k] += 0.5*waveSpeed*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
            }
        }
    }


#pragma unroll
    for (int k = 0; k < numStates; k++) {
        amdq[k] += F_tilde[k];
        apdq[k] -= F_tilde[k];
    }
}

/*
 * This regular wave approach does not handle capacity
 * In case where capacity variable exists, use fwave approach below.
 */
template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
__forceinline__ __device__ 
void first_second_order_update_vertical(pdeParam &param, volatile real* waves, volatile real* waveSpeeds, real* amdq, real* apdq, const int row, const int col, const Limiter phi)
{
    real F_tilde[numStates];
    for (int k = 0; k < numStates; k++){
        F_tilde[k] = claw_zero;
    }
#pragma unroll
    for (int w = 0; w < numWaves; w++)
    {
        real waveSpeed =  getSharedSpeed(waveSpeeds, row, col, w, numWaves, blockSizeY, blockSizeX);

        if (waveSpeed < claw_zero)
        {
            real limiting_factor = limiting_shared_v_d<numStates, numWaves, blockSizeX, blockSizeY>(phi, waves, row, col, w);
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                F_tilde[k] += -0.5*waveSpeed*(1+(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
                amdq[k] += waveSpeed * wave_state;
            }
        }
        else
        {
            real limiting_factor = limiting_shared_v_u<numStates, numWaves, blockSizeX, blockSizeY>(phi, waves, row, col, w);
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(waves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                F_tilde[k] += 0.5*waveSpeed*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
                apdq[k] += waveSpeed * wave_state;
            }
        }
    }

#pragma unroll
    for (int k = 0; k < numStates; k++) {
        amdq[k] += F_tilde[k];
        apdq[k] -= F_tilde[k];
    }
}

template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, const bool hasCapa,
    class Limiter>
__forceinline__ __device__ 
void first_second_order_update_fwave_horizontal(pdeParam &param, volatile real* fwaves, const volatile real* waveSpeeds, real* amdq, real* apdq, const int row, const int col, const Limiter phi,
        const real capa = claw_nan)
{
    real F_tilde[numStates];
#pragma unroll
    for (int k = 0; k < numStates; k++){
        F_tilde[k] = claw_zero;
    }
#pragma unroll
    for (int w = 0; w < numWaves; w++)
    {
        real waveSpeed =  getSharedSpeed(waveSpeeds, row, col, w, numWaves, blockSizeY, blockSizeX);

        if (waveSpeed < claw_zero)
        {
            real limiting_factor = limiting_shared_h_l<numStates, numWaves, blockSizeX, blockSizeY>(phi, fwaves, row, col, w);

            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                if (hasCapa)
                    F_tilde[k] += 
                        -0.5*(1+(param.dt/param.dx/capa)*waveSpeed)
                        *limiting_factor*wave_state;
                else
                    F_tilde[k] += 
                        -0.5*(1+(param.dt/param.dx)*waveSpeed)
                        *limiting_factor*wave_state;
                amdq[k] += wave_state;
            }
        }
        else if (waveSpeed > claw_zero)
        {
            real limiting_factor = limiting_shared_h_r<numStates, numWaves, blockSizeX, blockSizeY>(phi, fwaves, row, col, w);

            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                if (hasCapa)
                    F_tilde[k] += 
                        0.5*(1-(param.dt/param.dx/capa)*waveSpeed)*limiting_factor*wave_state;
                else
                    F_tilde[k] += 0.5*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
                apdq[k] += wave_state;
            }
        }
        else 
        {
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                amdq[k] += 0.5*wave_state;
                apdq[k] += 0.5*wave_state;
            }
        }
    }


#pragma unroll
    for (int k = 0; k < numStates; k++) {
        amdq[k] += F_tilde[k];
        apdq[k] -= F_tilde[k];
    }
}

template<const unsigned int numStates, const unsigned int numWaves, 
    const unsigned int blockSizeX, const unsigned int blockSizeY, const bool hasCapa, 
    class Limiter>
__forceinline__ __device__ 
void first_second_order_update_fwave_vertical(pdeParam &param, volatile real* fwaves, const volatile real* waveSpeeds, real* amdq, real* apdq, 
        const int row, const int col, const Limiter phi, const real capa = claw_nan)
{
    real F_tilde[numStates];
#pragma unroll
    for (int k = 0; k < numStates; k++){
        F_tilde[k] = claw_zero;
    }
#pragma unroll
    for (int w = 0; w < numWaves; w++)
    {
        real waveSpeed =  getSharedSpeed(waveSpeeds, row, col, w, numWaves, blockSizeY, blockSizeX);

        if (waveSpeed < claw_zero)
        {
            real limiting_factor = limiting_shared_v_d<numStates, numWaves, blockSizeX, blockSizeY>(phi, fwaves, row, col, w);
            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                if (hasCapa)
                    F_tilde[k] += 
                        -0.5*(1+(param.dt/param.dx/capa)*waveSpeed)
                        *limiting_factor*wave_state;
                else
                    F_tilde[k] += 
                        -0.5*(1+(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
                amdq[k] += wave_state;
            }
        }
        else if (waveSpeed > claw_zero)
        {
            real limiting_factor = limiting_shared_v_u<numStates, numWaves, blockSizeX, blockSizeY>(phi, fwaves, row, col, w);
            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                if (hasCapa)
                    F_tilde[k] += 
                        0.5*(1-(param.dt/param.dx/capa)*waveSpeed)*limiting_factor*wave_state;
                else
                    F_tilde[k] += 
                        0.5*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
                apdq[k] += wave_state;
            }
        }
        else
        {
#pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                real wave_state = getSharedWave(fwaves, row, col, w, k, numWaves, numStates, blockSizeY, blockSizeX);
                amdq[k] += 0.5*wave_state;
                apdq[k] += 0.5*wave_state;
            }
        }
    }

#pragma unroll
    for (int k = 0; k < numStates; k++) {
        amdq[k] += F_tilde[k];
        apdq[k] -= F_tilde[k];
    }
}

template <const unsigned int numStates, const unsigned int numWaves, const unsigned int numCoeff, const unsigned int blockSize, class Riemann, class Limiter>
__global__ void 
__launch_bounds__(HORIZONTAL_K_BLOCKSIZE,KERNEL_MIN_BLOCKS)
Riemann_horizontal_kernel(pdeParam param, Riemann Riemann_solver_h, Limiter phi)
{
    // Every threads gets a unique interface
    // a thread at (row, col) treats the interface between cells
    // (row,col)|(row,col+1)
    //          /\
    //		thread (row,col)
    //
    //  each block solves Riemann problems that is required for blockDim.x-3 columns of cells 
    
    int col = threadIdx.x + blockIdx.x*blockDim.x - 3*blockIdx.x;
    int row = threadIdx.y + blockIdx.y*blockDim.y;


    SharedMemory<real> smem;
    volatile real* wavesX = smem.getPointer();
    volatile real* waveSpeedsX = smem.getPointer(blockSize*numWaves*numStates);

    real apdq[numStates];
    real amdq[numStates];


#pragma unroll
    for (int i = 0; i < numWaves; i++)
        setSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, claw_zero);

    bool grid_valid = row < param.cellsY && col < param.cellsX-1;

    const real leftCell[numStates] = 
    {
        grid_valid ? param.getElement_qNew(row,col,0) : CUDART_NAN,
        grid_valid ? param.getElement_qNew(row,col,1) : CUDART_NAN,
        grid_valid ? param.getElement_qNew(row,col,2) : CUDART_NAN 
    };	

    const real rightCell[numStates] = 
    {
        grid_valid ? param.getElement_qNew(row,col+1,0) : CUDART_NAN,
        grid_valid ? param.getElement_qNew(row,col+1,1) : CUDART_NAN,
        grid_valid ? param.getElement_qNew(row,col+1,2) : CUDART_NAN 
    };	

    const real leftCoeff[numCoeff] = 
    {
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col,0) : CUDART_NAN,
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col,1) : CUDART_NAN,
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col,2) : CUDART_NAN
    };

    const real rightCoeff[numCoeff] = 
    {
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col+1,0) : CUDART_NAN,
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col+1,1) : CUDART_NAN,
        ((numCoeff > 0) && grid_valid) ? param.getElement_coeff(row,col+1,2) : CUDART_NAN
    };


    // Riemann solver
    // Let's say the cells in X direction are indexed from 0 to cellsX-1 (including ghost cells). 
    // The last thread is mapped to the second last cell, cell (row,cellsX-2), and will
    // solve the Riemann problem between cell (row,cellsX-2) and cell (row,cellsX-1) 
    if (grid_valid)	
    {
        Riemann_solver_h(leftCell, rightCell,
                leftCoeff, rightCoeff,
                threadIdx.y, threadIdx.x,
                numStates,
                numWaves,
                wavesX,
                waveSpeedsX,
                blockSize);
    }
    __syncthreads();

    // idle the first and the last thread since we don't need to limit the first and the last wave in this block
    // should also handle the edge case:  idle the last effective thread in the last block in x direction, in case
    // only part of the last block is mapped to internal of a grid patch
    grid_valid = grid_valid && 
        (0 < threadIdx.x && threadIdx.x < (HORIZONTAL_K_BLOCKSIZEX-1) && col < (param.cellsX-2));

    if (grid_valid)
    {
#pragma unroll
        for (int k = 0; k < numStates; k++)
        {
            amdq[k] = claw_zero;
            apdq[k] = claw_zero;
        }
        if (!param.second_order)	// simple first order scheme
        {
            first_order_update<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY>
                (wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x);
        }
        else
        {
#ifdef USE_FWAVES
#ifdef USE_CAPA
            first_second_order_update_fwave_horizontal<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY,true> 
                (param, wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x, phi,
                 0.5*(leftCoeff[param.mcapa]+rightCoeff[param.mcapa]));
#else
            first_second_order_update_fwave_horizontal<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY,false> 
                (param, wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x, phi);
#endif // USE_CAPA
#else
            first_second_order_update_horizontal<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY> 
                (param, wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x, phi);
#endif // USE_FWAVES
        }
    }
    __syncthreads(); // includes  __threadfence_block()

    // write the apdq to shared memory
    // overwrite some entries used by wavesX in shared memory
#pragma unroll
    for(int k = 0; k < numStates; k++)
    {
        setSharedWave(wavesX, threadIdx.y, threadIdx.x, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, apdq[k]);
    }

    // Some wave speeds should not be used in computing maximum CFL number,
    // e.g. wave speeds from the 1st and last Riemann problem in a row of cells
    // so we set them to 0.0
    // if (col == 0 || col == param.cellsX-2) {
    //     #pragma unroll
    //     for (int i = 0; i < numWaves; i++)
    //         setSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, claw_zero);
    // }

#ifdef DEBUG 
    // debug waveSpeeds
    __threadfence();
    if (grid_valid)
    {
        #pragma unroll
        for (int i = 0; i < numWaves; i++)
        {
            real ws =  getSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX );
            real cfl = fmax(-ws/leftCoeff[param.mcapa],ws/rightCoeff[param.mcapa])*param.dt/param.dx;
            if (cfl > 1.0)
            {
                printf("cfl larger than 1.0 in horizontal kernel. At row = %d, col = %d, wave number = %d. cfl = %E. ws = %E. local grid id = %d. \n", row, col, i, cfl, ws, param.id);
                printf("leftCell, h, hu, hv, capa in grid %d, row = %d, col = %d: %E, %E, %E, %E. \n", param.id, row, col, leftCell[0], leftCell[1], leftCell[2], leftCoeff[param.mcapa]);
                printf("rightCell, h, hu, hv, capa in grid %d, row = %d, col = %d: %E, %E, %E, %E. \n", param.id, row, col, rightCell[0], rightCell[1], rightCell[2], rightCoeff[param.mcapa]);
            }
        }
    }
#endif

    // Local Reduce over Wave Speeds
    // Stage 1
    // Bringing down the number of elements to compare to block size
    int tid = threadIdx.x + threadIdx.y*blockDim.x;
    real ws_max = 0.0;
#ifdef USE_CAPA
    // if s>0, use capa in right cell to compute CFL
    // if s<0, use capa in left cell to compute CFL
    if (grid_valid) {
#pragma unroll
        for (int i = 0; i < numWaves; i++)
        {
            real ws =  getSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, 
                    i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX );
            real ws_capa = fmax(-ws/leftCoeff[param.mcapa],ws/rightCoeff[param.mcapa]);
            ws_max = fmax(ws_max,ws_capa);
        }
    } 
#else
    for (int i = 0; i < numWaves; i++)
    {
        real ws =  getSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, 
                i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX );
        ws_max = fmax(ws_max,ws);
    }
#endif
    __syncthreads();
    // threads where (grid_valid == false) will have ws_max = 0.0
    waveSpeedsX[tid] = ws_max; 

    __syncthreads();	// unmovable syncthreads

    // Stage 2
    // Reducing over block size elements
    // use knowledge about your own block size:
    // I am assuming blocks to be of size no more than 512 will be used for the horizontal direction
    // There is a potential for a very subtle bug in the implementation below:
    // if the thread block has size non poer of 2, then there would be 2 threads reading/writing
    // from the same location, for example if the block size is 135, threads [0-67] will be active
    // and thread 0 will operate on (0,67) and thread 67 will operate on (67, 135). This can cause
    // an issue IF the SM somehow reads an unfinished write to the shared memory location. The
    // read data would be junk and could potentially hinder the simulation, either by a crash, or
    // behaving like a large number (behaving as a small number is no problem as this is a max reduce)
    // which could slow down the simulation. In any case block size should be multiples of 32 at least,
    // and never odd numbers. Also if the block size is between 32 and 64 warp reduce might access
    // off limit data. Rule of thumb keep block sizes to pwoers of 2.
    // At this stage there is no need to use fabs again, as all speeds were taken absolutely
    if (HORIZONTAL_K_BLOCKSIZE > 64 )
    {
        if (tid < (HORIZONTAL_K_BLOCKSIZE+1)/2)
            waveSpeedsX[tid] = fmax(waveSpeedsX[tid], waveSpeedsX[tid + HORIZONTAL_K_BLOCKSIZE/2]);
        __syncthreads();
    }
    if (HORIZONTAL_K_BLOCKSIZE/2 > 64 )
    {
        if (tid < (HORIZONTAL_K_BLOCKSIZE+3)/4)
            waveSpeedsX[tid] = fmax(waveSpeedsX[tid], waveSpeedsX[tid + HORIZONTAL_K_BLOCKSIZE/4]);
        __syncthreads();
    }
    if (HORIZONTAL_K_BLOCKSIZE/4 > 64 )
    {
        if (tid < (HORIZONTAL_K_BLOCKSIZE+7)/8)
            waveSpeedsX[tid] = fmax(waveSpeedsX[tid], waveSpeedsX[tid + HORIZONTAL_K_BLOCKSIZE/8]);
        __syncthreads();
    }
    if (HORIZONTAL_K_BLOCKSIZE/8 > 64 )
    {
        if (tid < (HORIZONTAL_K_BLOCKSIZE+15)/16)
            waveSpeedsX[tid] = fmax(waveSpeedsX[tid], waveSpeedsX[tid + HORIZONTAL_K_BLOCKSIZE/16]);
        __syncthreads();
    }
    if (tid < 32)
        warpReduce<blockSize>(waveSpeedsX, tid);

    if (grid_valid && threadIdx.x > 1)
    {
        #pragma unroll
        for (int k = 0; k < numStates; k++)
        {
#ifdef USE_CAPA
            param.setElement_q_tmp(row, col, k, 
                param.getElement_qNew(row, col, k) - param.dt/param.dx/leftCoeff[param.mcapa] * (amdq[k] + getSharedWave(wavesX, threadIdx.y, threadIdx.x-1, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX))
                    );
#else
            param.setElement_q_tmp(row, col, k, 
                param.getElement_qNew(row, col, k) - param.dt/param.dx * (amdq[k] + getSharedWave(wavesX, threadIdx.y, threadIdx.x-1, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX))
                    );
#endif
        }
    }
    else { // copy the remaining cells (ghost cells) to qNew, which is required by the vertical kernel
        if ( (row < param.cellsY) &&
             (col < 2 || (col > param.cellsX-3 && col < param.cellsX))
           ) {
            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                param.setElement_q_tmp(row, col, k, param.getElement_qNew(row, col, k));
            }
        }
    }

    __syncthreads();
    if (tid == 0) {
        // param.waveSpeedsX[blockIdx.x + blockIdx.y*gridDim.x] = waveSpeedsX[0];
        real cfl = waveSpeedsX[0]*param.dt/param.dx;
        atomicMax(param.cfl_grid, cfl);
#ifdef DEBUG
        if ( cfl > 1.0)
            printf("cfl larger than 1.0 in horizontal kernel. cfl = %E. \n blockIdx.x = %d, blockIdx.y = %d, local grid id = %d. \n", cfl, blockIdx.x, blockIdx.y, param.id);
//         int idx = blockIdx.x + blockIdx.y*gridDim.x;
//         if (idx == 0) {
//             printf("num of cells for grid %d: %d.\n", param.id, param.cellsX*param.cellsY);
//             printf("sizes of ws_x for grid %d: %d.\n", param.id, gridDim.x*gridDim.y);
//             printf("address of param.qNew for grid %d: %p.\n", param.id,(void*)param.qNew);
//             printf("address of param.q_tmp for grid %d: %p.\n", param.id,(void*)param.q_tmp);
//             printf("address of param.coefficients for grid %d: %p.\n", param.id,(void*)param.coefficients);
//             printf("address of param.waveSpeedsX for grid %d: %p.\n", param.id,(void*)param.waveSpeedsX);
//             printf("address of param.waveSpeedsY for grid %d: %p.\n", param.id,(void*)param.waveSpeedsY);
//         }
#endif
    }

}



template <const unsigned int numStates, const unsigned int numWaves, const unsigned int numCoeff,
         const unsigned int blockSize, class Riemann, class Limiter>
__global__ void 
__launch_bounds__(VERTICAL_K_BLOCKSIZE,KERNEL_MIN_BLOCKS)
Riemann_vertical_kernel(pdeParam param, Riemann Riemann_solver_v, Limiter phi)
{
    // Every threads gets a unique interface
    // a thread at (row, col) treats the interface between cells
    // (row+1,col)
    // -----------  << thread (row,col)
    // (row  ,col)
    int col = threadIdx.x + blockIdx.x*blockDim.x;
    int row = threadIdx.y + blockIdx.y*blockDim.y - 3*blockIdx.y;

    real apdq[numStates];
    real amdq[numStates];

    SharedMemory<real> smem;
    volatile real* wavesY = smem.getPointer();
    volatile real* waveSpeedsY = smem.getPointer(blockSize*numWaves*numStates);

    real upCell[numStates];		
    real downCell[numStates];	

//     staticArrayReal<numCoeff> upCoeff;
//     staticArrayReal<numCoeff> downCoeff;
    real upCoeff[numCoeff];		
    real downCoeff[numCoeff];	

#pragma unroll
    for (int i = 0; i < numWaves; i++)
        setSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX, claw_zero);

    bool grid_valid = row < param.cellsY-1 && col < param.cellsX;

    // Riemann solver
    if (grid_valid)
    {
#pragma unroll
        for (int i = 0; i < numStates; i++)
        {
            downCell[i] = param.getElement_q_tmp(row,col,i);	// use intermediate data
            upCell[i] = param.getElement_q_tmp(row+1,col,i);
        }
        if (numCoeff > 0) {
#pragma unroll
            for (int i = 0; i < numCoeff; i++)
            {
                downCoeff[i] = param.getElement_coeff(row,col,i);
                upCoeff[i] = param.getElement_coeff(row+1,col,i);
            }
        }
        Riemann_solver_v(	downCell,			// input comes from global memory
                upCell,				//
                downCoeff,			//
                upCoeff,			//
                threadIdx.y,		//
                threadIdx.x,		//
                numStates,			//
                numWaves,			//
                wavesY,				// output to shared memory
                waveSpeedsY,
                blockSize);		//

    }
    __syncthreads();

    // idle the first and the last thread since we don't need to limit the first and the last wave in this block
    // should also handle the edge case:  idle the last effective thread in the last block in y direction
    grid_valid = grid_valid && ( 0 < threadIdx.y && threadIdx.y < VERTICAL_K_BLOCKSIZEY-1 ) 
        && (row < (param.cellsY-2));

    if (grid_valid)
    {
#pragma unroll
        for (int k = 0; k < numStates; k++)
        {
            amdq[k] = claw_zero;
            apdq[k] = claw_zero;
        }
        if (!param.second_order)
        {
            first_order_update<numStates, numWaves, VERTICAL_K_BLOCKSIZEX, VERTICAL_K_BLOCKSIZEY>(wavesY, waveSpeedsY, amdq, apdq, threadIdx.y, threadIdx.x);
        }
        else	// simple first order scheme
        {
#ifdef USE_FWAVES
#ifdef USE_CAPA
            first_second_order_update_fwave_vertical<numStates, numWaves, VERTICAL_K_BLOCKSIZEX, VERTICAL_K_BLOCKSIZEY,true>
                (param, wavesY, waveSpeedsY, amdq, apdq, threadIdx.y, threadIdx.x, phi,
                 0.5*(upCoeff[param.mcapa]+downCoeff[param.mcapa]));
#else
            first_second_order_update_fwave_vertical<numStates, numWaves, VERTICAL_K_BLOCKSIZEX, VERTICAL_K_BLOCKSIZEY,false>
                (param, wavesY, waveSpeedsY, amdq, apdq, threadIdx.y, threadIdx.x, phi);
#endif // USE_CAPA
#else
            first_second_order_update_vertical<numStates, numWaves, VERTICAL_K_BLOCKSIZEX, VERTICAL_K_BLOCKSIZEY>
                (param, wavesY, waveSpeedsY, amdq, apdq, threadIdx.y, threadIdx.x, phi);
#endif // USE-FWAVES
        }
    }
    __syncthreads();

    // write the apdq to shared memory
#pragma unroll
    for(int k = 0; k < numStates; k++)
    {
        setSharedWave(wavesY, threadIdx.y, threadIdx.x, 0, k, numWaves, numStates, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX, apdq[k]);
    }

    // Some wave speeds should not be used in computing maximum CFL number,
    // e.g. wave speeds from the 1st and last Riemann problem in a col of cells
    // and from the first and last two columns
    // so we set them to 0.0
    // if (row == 0 || row == param.cellsY-2 || col < 2 || col > param.cellsX-3)
    //     #pragma unroll
    //     for (int i = 0; i < numWaves; i++)
    //         setSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX, claw_zero);

#ifdef DEBUG
    // debug waveSpeeds
    __threadfence();
    if (grid_valid)
    {
        #pragma unroll
        for (int i = 0; i < numWaves; i++)
        {
            real ws =  getSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX);
            real cfl = fmax(-ws/downCoeff[param.mcapa],ws/upCoeff[param.mcapa])*param.dt/param.dy;
            if (cfl > 1.0)
            {
                printf("cfl larger than 1.0 in vertical kernel. At row = %d, col = %d, wave number = %d. \ncfl = %E. local grid id = %d. \n", row, col, i, cfl, param.id);
            }
        }
    }
#endif

    // Local Reduce over Wave Speeds
    // Stage 1
    // Bringing down the number of elements to compare to block size
    int tid = threadIdx.x + threadIdx.y*blockDim.x;
    real ws_max = 0.0;
#ifdef USE_CAPA
    // if s>0, use capa in up cell to compute CFL
    // if s<0, use capa in down cell to compute CFL
    if (grid_valid) {
#pragma unroll
        for (int i = 0; i < numWaves; i++)
        {
            real ws =  getSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, 
                    i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX );
            real ws_capa = fmax(-ws/downCoeff[param.mcapa],ws/upCoeff[param.mcapa]);
            ws_max = fmax(ws_max,ws_capa);
        }
    }
#else
    for (int i = 0; i < numWaves; i++)
    {
        real ws =  getSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, 
                i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX );
        ws_max = fmax(ws_max,ws);
    }
#endif
    waveSpeedsY[tid] = ws_max;

    __syncthreads();

    if (VERTICAL_K_BLOCKSIZE > 64)
    {
        if (tid < (VERTICAL_K_BLOCKSIZE+1)/2)
            waveSpeedsY[tid] = fmax(waveSpeedsY[tid], waveSpeedsY[tid + VERTICAL_K_BLOCKSIZE/2]);
        __syncthreads();
    }
    if (VERTICAL_K_BLOCKSIZE/2 > 64)
    {
        if (tid < (VERTICAL_K_BLOCKSIZE+3)/4)
            waveSpeedsY[tid] = fmax(waveSpeedsY[tid], waveSpeedsY[tid + VERTICAL_K_BLOCKSIZE/4]);
        __syncthreads();
    }
    if (VERTICAL_K_BLOCKSIZE/4 > 64)
    {
        if (tid < (VERTICAL_K_BLOCKSIZE+7)/8)
            waveSpeedsY[tid] = fmax(waveSpeedsY[tid], waveSpeedsY[tid + VERTICAL_K_BLOCKSIZE/8]);
        __syncthreads();
    }
    if (VERTICAL_K_BLOCKSIZE/8 > 64)
    {
        if (tid < (VERTICAL_K_BLOCKSIZE+15)/16)
            waveSpeedsY[tid] = fmax(waveSpeedsY[tid], waveSpeedsY[tid + VERTICAL_K_BLOCKSIZE/16]);
        __syncthreads();
    }
    if (tid < 32)
        warpReduce<blockSize>(waveSpeedsY, tid);


    if (grid_valid && threadIdx.y > 1)
    {
        #pragma unroll
        for (int k = 0; k < numStates; k++)
        {
#ifdef USE_CAPA
            param.setElement_qNew(row, col, k, 
                param.getElement_q_tmp(row, col, k) - param.dt/param.dy/downCoeff[param.mcapa] * (amdq[k] + getSharedWave(wavesY, threadIdx.y-1, threadIdx.x, 0, k, numWaves, numStates, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX))
                    );
#else
            param.setElement_qNew(row, col, k, 
                param.getElement_q_tmp(row, col, k) - param.dt/param.dy * (amdq[k] + getSharedWave(wavesY, threadIdx.y-1, threadIdx.x, 0, k, numWaves, numStates, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX))
                    );
#endif
        }
    }

    __syncthreads();
    if (tid == 0)
    {
        // param.waveSpeedsY[blockIdx.x + blockIdx.y*gridDim.x] = waveSpeedsY[0];
        real cfl = waveSpeedsY[0]*param.dt/param.dy;
        atomicMax(param.cfl_grid, cfl);
#ifdef DEBUG
        if ( cfl > 1.0)
            printf("cfl larger than 1.0 in vertical kernel. cfl = %E. \n blockIdx.x = %d, blockIdx.y = %d, local grid id = %d. \n", cfl, blockIdx.x, blockIdx.y, param.id);
//         int idx = blockIdx.x + blockIdx.y*gridDim.x;
//         if (idx == 0) {
//             printf("sizes of ws_y for grid %d: %d.\n", param.id, gridDim.x*gridDim.y);
//         }
#endif
    }
}


#ifdef GEOCLAW
__global__ void 
b4xsweep_kernel(pdeParam param);
__global__ void 
b4ysweep_kernel(pdeParam param);

#endif

template <class Riemann_h, class Riemann_v, class Limiter>
int limited_Riemann_Update(pdeParam &param,						// Problem parameters
        Riemann_h Riemann_pointwise_solver_h,	//
        Riemann_v Riemann_pointwise_solver_v,	//
        Limiter limiter_phi,
        const cudaStream_t stream
        )
{
    {
        // before x-sweep
        const unsigned int blockDim_XR = HORIZONTAL_K_BLOCKSIZEX;
        const unsigned int blockDim_YR = HORIZONTAL_K_BLOCKSIZEY;
        unsigned int gridDim_XR = (param.cellsX + blockDim_XR - 1) / blockDim_XR;
        unsigned int gridDim_YR = (param.cellsY + blockDim_YR - 1) / blockDim_YR;
        dim3 dimGrid_hR(gridDim_XR, gridDim_YR);
        dim3 dimBlock_hR(blockDim_XR, blockDim_YR);
        b4xsweep_kernel<<<dimGrid_hR,dimBlock_hR,0,stream>>>(param);
    }
    {
        // RIEMANN, FLUCTUATIONS and UPDATES
        const unsigned int blockDim_XR = HORIZONTAL_K_BLOCKSIZEX;
        const unsigned int blockDim_YR = HORIZONTAL_K_BLOCKSIZEY;
        const unsigned int nCellsPerBlock = blockDim_XR-3; // each CUDA block will update this many cells in X direction
        unsigned int gridDim_XR = (param.cellsX - 2*param.numGhostCells + nCellsPerBlock-1) / nCellsPerBlock;
        unsigned int gridDim_YR = (param.cellsY + (blockDim_YR-1)) / (blockDim_YR);
        dim3 dimGrid_hR(gridDim_XR, gridDim_YR);
        dim3 dimBlock_hR(blockDim_XR, blockDim_YR);
        int shared_mem_size = HORIZONTAL_K_BLOCKSIZEX*HORIZONTAL_K_BLOCKSIZEY*NWAVES*(NEQNS+1)*sizeof(real);

        // allocate param.waveSpeedsX here
        // param.alloc_waveSpeedsX(dimGrid_hR.x * dimGrid_hR.y);

        Riemann_horizontal_kernel<NEQNS, NWAVES, NCOEFFS, HORIZONTAL_K_BLOCKSIZEX*HORIZONTAL_K_BLOCKSIZEY, Riemann_h, Limiter><<<dimGrid_hR, dimBlock_hR, shared_mem_size,stream>>>(param, Riemann_pointwise_solver_h, limiter_phi);

        // REDUCTION
        const unsigned int blockDim_X = 512;		// fine tune the best block size

        size_t SharedMemorySize = (blockDim_X)*sizeof(real);
        unsigned int gridDim_X1;

        gridDim_X1 = 1;

        dim3 dimGrid1(gridDim_X1);
        dim3 dimBlock1(blockDim_X);
//         check_cfl<blockDim_X><<< dimGrid1, dimBlock1, 0, stream>>>(param.waveSpeedsX, gridDim_XR*gridDim_YR, param, 0);
//         reduceMax_simplified<blockDim_X><<< dimGrid1, dimBlock1, SharedMemorySize, stream>>>(param.waveSpeedsX, gridDim_XR*gridDim_YR, param, 0,
//                 param.dx, param.dt, &(param.cfl_grid[0]) );

    }
    {
        // before y-sweep
        const unsigned int blockDim_XR = VERTICAL_K_BLOCKSIZEX;
        const unsigned int blockDim_YR = VERTICAL_K_BLOCKSIZEY;
        unsigned int gridDim_XR = (param.cellsX + blockDim_XR - 1) / blockDim_XR;
        unsigned int gridDim_YR = (param.cellsY + blockDim_YR - 1) / blockDim_YR;
        dim3 dimGrid_hR(gridDim_XR, gridDim_YR);
        dim3 dimBlock_hR(blockDim_XR,blockDim_YR);
        b4ysweep_kernel<<<dimGrid_hR,dimBlock_hR,0,stream>>>(param);
    }
    {
        // RIEMANN, FLUCTUATIONS and UPDATE
        const unsigned int blockDim_XR = VERTICAL_K_BLOCKSIZEX;
        const unsigned int blockDim_YR = VERTICAL_K_BLOCKSIZEY;
        const unsigned int nCellsPerBlock = blockDim_YR-3; // each CUDA block will update this many cells in Y direction
        unsigned int gridDim_XR = (param.cellsX + (blockDim_XR-1)) / (blockDim_XR);
        unsigned int gridDim_YR = (param.cellsY - 2*param.numGhostCells + nCellsPerBlock-1) / nCellsPerBlock;
        dim3 dimGrid_vR(gridDim_XR, gridDim_YR);
        dim3 dimBlock_vR(blockDim_XR, blockDim_YR);
        int shared_mem_size = VERTICAL_K_BLOCKSIZEX*VERTICAL_K_BLOCKSIZEY*NWAVES*(NEQNS+1)*sizeof(real);

        // allocate param.waveSpeedsY here
        // param.alloc_waveSpeedsY(dimGrid_vR.x * dimGrid_vR.y);
        Riemann_vertical_kernel<NEQNS, NWAVES, NCOEFFS, VERTICAL_K_BLOCKSIZEX*VERTICAL_K_BLOCKSIZEY, Riemann_v, Limiter><<<dimGrid_vR, dimBlock_vR,shared_mem_size,stream>>>(param, Riemann_pointwise_solver_v, limiter_phi);


        // REDUCTION
        const unsigned int blockDim_X = 512;		// fine tune the best block size

        size_t SharedMemorySize = (blockDim_X)*sizeof(real);
        unsigned int gridDim_X2 = 1;

        dim3 dimGrid2(gridDim_X2);
        dim3 dimBlock2(blockDim_X);

//         check_cfl<blockDim_X><<< dimGrid2, dimBlock2, 0, stream>>>(param.waveSpeedsY, gridDim_XR*gridDim_YR, param, 1);
//         reduceMax_simplified<blockDim_X><<< dimGrid2, dimBlock2, SharedMemorySize, stream>>>(param.waveSpeedsY, gridDim_XR*gridDim_YR, param, 1,
//                 param.dy, param.dt, &(param.cfl_grid[1]) );
    }
    // Let the memory manager free memory allocated by param later 
    // when the memory is not needed by GPU
    // param.register_allocated_memory();

    return 0;
}



#endif // __FUSED_SOLVERS_H__

