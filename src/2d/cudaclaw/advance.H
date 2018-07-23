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
#include <cstdio>


template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY>
inline __device__ void first_order_update(real* waves, real* waveSpeeds, real* amdq, real* apdq, const int row, const int col)
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



template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
inline __device__ void first_second_order_update_horizontal(pdeParam &param, real* waves, real* waveSpeeds, real* amdq, real* apdq, const int row, const int col, const Limiter phi)
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
                F_tilde[k] += -(real)0.5f*waveSpeed*(1+(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
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
                F_tilde[k] += (real)0.5f*waveSpeed*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
            }
        }
    }


#pragma unroll
    for (int k = 0; k < numStates; k++) {
        amdq[k] += F_tilde[k];
        apdq[k] -= F_tilde[k];
    }
}

template<const int numStates, const int numWaves, const int blockSizeX, const int blockSizeY, class Limiter>
inline __device__ void first_second_order_update_vertical(pdeParam &param, real* waves, real* waveSpeeds, real* amdq, real* apdq, const int row, const int col, const Limiter phi)
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
                F_tilde[k] += -(real)0.5f*waveSpeed*(1+(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
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
                F_tilde[k] += (real)0.5f*waveSpeed*(1-(param.dt/param.dx)*waveSpeed)*limiting_factor*wave_state;
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

template <const unsigned int numStates, const unsigned int numWaves, const unsigned int numCoeff, const unsigned int blockSize, class Riemann, class Limiter>
__global__ void Riemann_horizontal_kernel(pdeParam param, Riemann Riemann_solver_h, Limiter phi)
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
    real* wavesX = smem.getPointer();
    real* waveSpeedsX = smem.getPointer(blockSize*numWaves*numStates);

    real apdq[numStates];
    real amdq[numStates];

    real leftCell[numStates];	
    real rightCell[numStates];  

    staticArrayReal<numCoeff> leftCoeff;
    staticArrayReal<numCoeff> rightCoeff;


    // Parts of the CUDA block might be mapped to outside of grid patch. 
    // So we need to set waveSpeedsX that corresponds to those parts to 0.0
    // These 0.0's are automatically ignored when we do max reduction later
    if (blockIdx.x == gridDim.x-1 || blockIdx.y == gridDim.y-1) // last block in either direction
    {
#pragma unroll
        for (int i = 0; i < numWaves; i++)
            setSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, claw_zero);
    }

    bool grid_valid = row < param.cellsY && col < param.cellsX-1;

    // Riemann solver
    // Let's say the cells in X direction are indexed from 0 to cellsX-1 (including ghost cells). 
    // The last thread is mapped to the second last cell, cell (row,cellsX-2), and will
    // solve the Riemann problem between cell (row,cellsX-2) and cell (row,cellsX-1) 
    if (grid_valid)	
    {
#pragma unroll
        for (int i = 0; i < numStates; i++)
        {
            leftCell[i] = param.getElement_qNew(row,col,i);
            rightCell[i] = param.getElement_qNew(row,col+1,i);
        }
        if (numCoeff > 0) {
#pragma unroll
            for (int i = 0; i < numCoeff; i++)
            {
                leftCoeff[i] = param.getElement_coeff(row,col,i);
                rightCoeff[i] = param.getElement_coeff(row,col+1,i);
            }
        }
        Riemann_solver_h(leftCell, rightCell,
                leftCoeff, rightCoeff,
                threadIdx.y, threadIdx.x,
                numStates,
                numWaves,
                wavesX,
                waveSpeedsX);
    }
    __syncthreads();

    // idle the first and the last thread since we don't need to limit the first and the last wave in this block
    // should also handle the edge case:  idle the last effective thread in the last block in x direction
    grid_valid = grid_valid && 
        (0 < threadIdx.x && threadIdx.x < (HORIZONTAL_K_BLOCKSIZEX-1) && col < (param.cellsX-2));

    if (grid_valid)
    {
#pragma unroll
        for (int k = 0; k < numStates; k++)
        {
            amdq[k] = 0;
            apdq[k] = 0;
        }
        if (!param.second_order)	// simple first order scheme
        {
            first_order_update<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY>
                (wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x);
        }
        else
        {
            first_second_order_update_horizontal<numStates, numWaves, HORIZONTAL_K_BLOCKSIZEX, HORIZONTAL_K_BLOCKSIZEY> 
                (param, wavesX, waveSpeedsX, amdq, apdq, threadIdx.y, threadIdx.x, phi);
        }
    }
    __syncthreads();

    // write the apdq to shared memory
    // overwrite some entries used by wavesX in shared memory
#pragma unroll
    for(int k = 0; k < numStates; k++)
    {
        setSharedWave(wavesX, threadIdx.y, threadIdx.x, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, apdq[k]);
    }

    // Local Reduce over Wave Speeds
    // Stage 1
    // Bringing down the number of elements to compare to block size
    int tid = threadIdx.x + threadIdx.y*HORIZONTAL_K_BLOCKSIZEX;
    waveSpeedsX[tid] = fabs(waveSpeedsX[tid]);
#pragma unroll
    for (int i = 1; i < numWaves; i++)
    {
        // TODO: call getSharedSpeed to do this?
        // blockDim.x * blockDim.y is usally a multiple of 32, so there should be no bank conflicts.
        waveSpeedsX[tid] = fmax(waveSpeedsX[tid], fabs(waveSpeedsX[tid + i*blockSize]));
    }
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
            param.setElement_q(row, col, k, 
                param.getElement_qNew(row, col, k) - param.dt/param.dx * (amdq[k] + getSharedWave(wavesX, threadIdx.y, threadIdx.x-1, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX))
                    );
        }
    }
    else // copy the remaining cells (ghost cells) to qNew, which is required by the vertical kernel
    {
        if (col < 2 || (col > param.cellsX-3 && col < param.cellsX) )
        {
            #pragma unroll
            for (int k = 0; k < numStates; k++)
            {
                param.setElement_q(row, col, k, param.getElement_qNew(row, col, k));
            }
        }
    }

    if (tid == 0)
        param.waveSpeedsX[blockIdx.x + blockIdx.y*gridDim.x] =  waveSpeedsX[0];
}

// Riemann_horizontal_kernel using shuffle instruction within a warp
// Only 1st order
template <const int numStates, const int numWaves, const unsigned int blockSize, class Riemann, class Limiter>
__global__ void Riemann_horizontal_kernel_warp(pdeParam param, Riemann Riemann_solver_h, Limiter phi)
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
    real* wavesX = smem.getPointer();
    real* waveSpeedsX = smem.getPointer(blockSize*numWaves*numStates);

    real apdq[numStates];
    real amdq[numStates];

    real leftCell[numStates];	
    real rightCell[numStates];  

    real* leftCoeff = 0;
    real* rightCoeff = 0;


    if (blockIdx.x == gridDim.x-1 || blockIdx.y == gridDim.y-1) // last block in either direction
    {
#pragma unroll
        for (int i = 0; i < numWaves; i++)
            setSharedSpeed(waveSpeedsX, threadIdx.y, threadIdx.x, i, numWaves, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, claw_zero);
    }

    bool grid_valid = row < param.cellsY && col < param.cellsX-1;

    // Riemann solver
    // if there are 512 cells (including ghost cells) in X, the last thread is mapped to cell (row,X-1) and will
    // solve the Riemann problem between cell (row,X-1) and cell (row,X) 
    if (grid_valid)	
    {
#pragma unroll
        for (int i = 0; i < numStates; i++)
        {
            leftCell[i] = param.getElement_qNew(row,col,i);
            rightCell[i] = param.getElement_qNew(row,col+1,i);
        }
        Riemann_solver_h(leftCell, rightCell,
                leftCoeff, rightCoeff,
                threadIdx.y, threadIdx.x,
                numStates,
                numWaves,
                wavesX,
                waveSpeedsX,
                apdq,amdq);
    }

    // idle the first and the last thread since we don't need to limit the first and the last wave in this block
    grid_valid = grid_valid && ( 0 < threadIdx.x && threadIdx.x < HORIZONTAL_K_BLOCKSIZEX-1 && col < param.cellsX-2);
    __syncthreads();

    // write the apdq to shared memory
    // overwrite some entries used by wavesX in shared memory
#pragma unroll
    for(int k = 0; k < numStates; k++)
    {
        setSharedWave(wavesX, threadIdx.y, threadIdx.x, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX, apdq[k]);
    }

    // Local Reduce over Wave Speeds
    // Stage 1
    // Bringing down the number of elements to compare to block size
    int tid = threadIdx.x + threadIdx.y*HORIZONTAL_K_BLOCKSIZEX;
    waveSpeedsX[tid] = fabs(waveSpeedsX[tid]);
#pragma unroll
    for (int i = 1; i < numWaves; i++)
    {
        // blockDim.x * blockDim.y is usally a multiple of 32, so there should be no bank conflicts.
        waveSpeedsX[tid] = fmax(waveSpeedsX[tid], fabs(waveSpeedsX[tid + i*blockSize]));
    }
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


    int lane_id = threadIdx.x%warpSize;
    real apdq_k = claw_zero;
#pragma unroll
    for (int k = 0; k < numStates; k++)
    {
        apdq_k = __shfl_up_sync(FULL_MASK,apdq[k],1,warpSize);
        if (grid_valid && threadIdx.x > 1) {
            if (lane_id==0) 
                apdq_k = getSharedWave(wavesX, threadIdx.y, threadIdx.x-1, 0, k, numWaves, numStates, HORIZONTAL_K_BLOCKSIZEY, HORIZONTAL_K_BLOCKSIZEX);
            param.setElement_q(row, col, k, 
                    param.getElement_qNew(row, col, k) - param.dt/param.dx * (amdq[k] + apdq_k));
        }
    }

    if (tid == 0)
        param.waveSpeedsX[blockIdx.x + blockIdx.y*gridDim.x] =  waveSpeedsX[0];
}


template <const unsigned int numStates, const unsigned int numWaves, const unsigned int numCoeff, const unsigned int blockSize, class Riemann, class Limiter>
__global__ void Riemann_vertical_kernel(pdeParam param, Riemann Riemann_solver_v, Limiter phi)
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
    real* wavesY = smem.getPointer();
    real* waveSpeedsY = smem.getPointer(blockSize*numWaves*numStates);

    real upCell[numStates];		
    real downCell[numStates];	

    staticArrayReal<numCoeff> upCoeff;
    staticArrayReal<numCoeff> downCoeff;

    if (blockIdx.x == gridDim.x-1 || blockIdx.y == gridDim.y-1)
    {
#pragma unroll
        for (int i = 0; i < numWaves; i++)
            setSharedSpeed(waveSpeedsY, threadIdx.y, threadIdx.x, i, numWaves, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX, claw_zero);
    }

    bool grid_valid = row < param.cellsY-1 && col < param.cellsX;

    // Riemann solver
    if (grid_valid)
    {
#pragma unroll
        for (int i = 0; i < numStates; i++)
        {
            downCell[i] = param.getElement_q(row,col,i);	// use intermediate data
            upCell[i] = param.getElement_q(row+1,col,i);
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
                waveSpeedsY);		//

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
            first_second_order_update_vertical<numStates, numWaves, VERTICAL_K_BLOCKSIZEX, VERTICAL_K_BLOCKSIZEY>
                (param, wavesY, waveSpeedsY, amdq, apdq, threadIdx.y, threadIdx.x, phi);
        }
    }
    __syncthreads();

    // write the apdq to shared memory
#pragma unroll
    for(int k = 0; k < numStates; k++)
    {
        setSharedWave(wavesY, threadIdx.y, threadIdx.x, 0, k, numWaves, numStates, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX, apdq[k]);
    }

    // Local Reduce over Wave Speeds
    // Stage 1
    // Bringing down the number of elements to compare to block size
    int tid = threadIdx.x + threadIdx.y*VERTICAL_K_BLOCKSIZEX;
    // See horizontal version for comments
    waveSpeedsY[tid] = fabs(waveSpeedsY[tid]);
#pragma unroll
    for (int i = 1; i < numWaves; i++)
    {
        // blockDim.x * blockDim.y is usally a multiple of 32, so there should be no bank conflicts.
        waveSpeedsY[tid] = fmax(waveSpeedsY[tid], fabs(waveSpeedsY[tid + i*blockSize]));
    }
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
            param.setElement_qNew(row, col, k, 
                param.getElement_q(row, col, k) - param.dt/param.dy * (amdq[k] + getSharedWave(wavesY, threadIdx.y-1, threadIdx.x, 0, k, numWaves, numStates, VERTICAL_K_BLOCKSIZEY, VERTICAL_K_BLOCKSIZEX))
                    );
        }
    }

    if (tid == 0)
        param.waveSpeedsY[blockIdx.x + blockIdx.y*gridDim.x] =  waveSpeedsY[0];
}

template <class Riemann_h, class Riemann_v, class Limiter>
int limited_Riemann_Update(pdeParam &param,						// Problem parameters
        Riemann_h Riemann_pointwise_solver_h,	//
        Riemann_v Riemann_pointwise_solver_v,	//
        Limiter limiter_phi,
        const cudaStream_t stream
        )
{
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
        param.alloc_waveSpeedsX(dimGrid_hR.x * dimGrid_hR.y);

        Riemann_horizontal_kernel<NEQNS, NWAVES, NCOEFFS, HORIZONTAL_K_BLOCKSIZEX*HORIZONTAL_K_BLOCKSIZEY, Riemann_h, Limiter><<<dimGrid_hR, dimBlock_hR, shared_mem_size,stream>>>(param, Riemann_pointwise_solver_h, limiter_phi);
        // Riemann_horizontal_kernel_warp<NEQNS, NWAVES, HORIZONTAL_K_BLOCKSIZEX*HORIZONTAL_K_BLOCKSIZEY, Riemann_h, Limiter><<<dimGrid_hR, dimBlock_hR, shared_mem_size>>>(param, Riemann_pointwise_solver_h, limiter_phi);

        // REDUCTION
        const unsigned int blockDim_X = 512;		// fine tune the best block size

        size_t SharedMemorySize = (blockDim_X)*sizeof(real);
        unsigned int gridDim_X1;

        gridDim_X1 = 1;

        dim3 dimGrid1(gridDim_X1);
        dim3 dimBlock1(blockDim_X);
        reduceMax_simplified<blockDim_X><<< dimGrid1, dimBlock1, SharedMemorySize, stream>>>(param.waveSpeedsX, gridDim_XR*gridDim_YR,
                param.dx, param.dt, &(param.cfl_grid[0]) );
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
        param.alloc_waveSpeedsY(dimGrid_vR.x * dimGrid_vR.y);
        Riemann_vertical_kernel<NEQNS, NWAVES, NCOEFFS, VERTICAL_K_BLOCKSIZEX*VERTICAL_K_BLOCKSIZEY, Riemann_v, Limiter><<<dimGrid_vR, dimBlock_vR,shared_mem_size,stream>>>(param, Riemann_pointwise_solver_v, limiter_phi);


        // REDUCTION
        const unsigned int blockDim_X = 512;		// fine tune the best block size

        size_t SharedMemorySize = (blockDim_X)*sizeof(real);
        unsigned int gridDim_X2;

        gridDim_X2 = 1;

        dim3 dimGrid2(gridDim_X2);
        dim3 dimBlock2(blockDim_X);

        reduceMax_simplified<blockDim_X><<< dimGrid2, dimBlock2, SharedMemorySize, stream>>>(param.waveSpeedsY, gridDim_XR*gridDim_YR,
                param.dy, param.dt, &(param.cfl_grid[1]) );
    }
    // Let the memory manager free memory allocated by param later 
    // when the memory is not needed by GPU
    param.register_allocated_memory();

    return 0;
}



#endif // __FUSED_SOLVERS_H__

