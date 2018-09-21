#ifndef __REDUCE_MAX_H__
#define __REDUCE_MAX_H__

#include "real.H"
#include <assert.h>
#include <sm_20_atomic_functions.h>
#include "params.H"

template <unsigned int blockSize>
__global__ void 
check_cfl(real *globalWaves, int size, pdeParam param, int direction)
{
    int tid = threadIdx.x;
    if (tid < size) {
        real dx = (direction == 0) ? param.dx : param.dy; 
        real cfl = globalWaves[tid]*param.dt/dx;
        if (cfl > 1.0) {
            printf("cfl larger than 1.0 in check_cfl. cfl = %E. tid = %d. direction = %d. size = %d.\n", cfl, tid, direction, size);
            printf("local grid id: %d\n", param.id);
            printf("x-size: %d; y-size: %d\n", param.cellsX, param.cellsY);
            printf("address of globalWaves in check_cfl: %p.\n", (void*)globalWaves);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////   Reduction   //////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<unsigned int blockSize>
__device__ void 
warpReduce(volatile real *sdata, unsigned int tid) 
{
    if (blockSize >=  64) sdata[tid] = fmax(sdata[tid], sdata[tid + 32]);
    if (blockSize >=  32) sdata[tid] = fmax(sdata[tid], sdata[tid + 16]);
    if (blockSize >=  16) sdata[tid] = fmax(sdata[tid], sdata[tid +  8]);
    if (blockSize >=   8) sdata[tid] = fmax(sdata[tid], sdata[tid +  4]);
    if (blockSize >=   4) sdata[tid] = fmax(sdata[tid], sdata[tid +  2]);
    if (blockSize >=   2) sdata[tid] = fmax(sdata[tid], sdata[tid +  1]);
}


// Single block reduce Max. The assumption is that a single block reduces over the data.
// Technically gridDim = (1,0,0), blockDim = (blockSize,0,0)
// if there are more elements than threads, reduce elements into shared until all elements are checked
// else make sure the available data is put in shared and extra shared space is made 0
// Then reduce over all elements in shared memory
// return the first element in shared to the global input's first place
template <unsigned int blockSize>
__global__ void 
reduceMax_simplified(real *globalWaves, int size, pdeParam param, int direction,
        real dx = claw_zero, real dt = claw_zero, real* result = NULL)
{
    // size holds the exact number of wave speeds there are to be max reduced
    extern __shared__ real sharedWaveSpeeds[];	// This will be as big as the block size
    int tid = threadIdx.x;

    // test alternative: have a first phase to collect at least one layer in shared memory
    // then go block size by block size take the max of the elements
    if ( size > blockSize )
    {
        assert(globalWaves[tid] >= claw_zero);
        sharedWaveSpeeds[tid] = globalWaves[tid];

        int i = tid;

        while (i + blockSize < size)	// !!!! WE MUST MAKE SURE THAT WHEN ALLOCATING MEMORY FOR THE GLOBAL WAVES SPEEDS
            // !!!! WE ALLCOATE MORE MEMORY THAN NECESSARY (BY A MARGINAL AMMOUNT) AND WE
            // !!!! INITIALIZE THE MEMORY TO BE FILLED WITH ZEROS, AS i + blockSize COULD
            // !!!! READ OUTSIDE THE BOUNDS OF GLOBAL WAVE SPEEDS AND PRODUCE ERRONEOUS RESULTS
        {
            // all of these accesses will be cached,
            // in fact the whole reduced waves would fit in shared memory/L1 for medium sized problems
            // L2 has 768KB of memory, for a 1024*1024 grid we have about 16K elements to reduce
            // We do not need to take fabs as it is done in the local reduction stage
            assert(globalWaves[i+blockSize] >= claw_zero);
            sharedWaveSpeeds[tid] = fmax(sharedWaveSpeeds[tid], globalWaves[i+blockSize]);
            i += blockSize;
        }
    }
    else
    {
        if (tid < size) 
        {
            assert(globalWaves[tid] >= claw_zero);
            sharedWaveSpeeds[tid] = globalWaves[tid];

            real cfl = globalWaves[tid]*dt/dx;
            if (cfl > 1.0) {
                printf("cfl larger than 1.0 in early reduce_max. cfl = %E. tid = %d. direction = %d. size = %d.\n", cfl, tid, direction, size);
                printf("local grid id: %d\n", param.id);
                printf("x-size: %d; y-size: %d\n", param.cellsX, param.cellsY);
                printf("address of globalWaves at early reduce_max: %p.\n", (void*)globalWaves);
                printf("address of param.qNew at early reduce_max: %p.\n", (void*)param.qNew);
                printf("address of param.q_tmp at early reduce_max: %p.\n", (void*)param.q_tmp);
                printf("address of param.coefficients at early reduce_max: %p.\n", (void*)param.coefficients);
            }
        }
        else
            sharedWaveSpeeds[tid] = claw_zero;
    }
    __syncthreads();

    // ASSUMING BLOCK SIZES CAN BE ONLY 1024, 512, 256, 128, 64, OR 32, ... 1
    // At this stage we are sure that WHATEVER IS IN THE SHARED MEMORY WAS ALREADY TAKEN ABSOLUTELY
    // that is why we do not need to invoke the fabs function again.
    if (blockSize >= 1024)
    { 
        if (tid < 512)
            sharedWaveSpeeds[tid] = fmax(sharedWaveSpeeds[tid], sharedWaveSpeeds[tid + 512]);
        __syncthreads();
    }
    if (blockSize >= 512)
    { 
        if (tid < 256)
            sharedWaveSpeeds[tid] = fmax(sharedWaveSpeeds[tid], sharedWaveSpeeds[tid + 256]);
        __syncthreads();
    }
    if (blockSize >= 256)
    { 
        if (tid < 128)
            sharedWaveSpeeds[tid] = fmax(sharedWaveSpeeds[tid], sharedWaveSpeeds[tid + 128]);
        __syncthreads();
    }
    if (blockSize >= 128)
    { 
        if (tid < 64)
            sharedWaveSpeeds[tid] = fmax(sharedWaveSpeeds[tid], sharedWaveSpeeds[tid + 64]);
        __syncthreads();
    }
    if (tid < 32)
        warpReduce<blockSize>(sharedWaveSpeeds, tid);
    if (tid == 0) {
        real cfl = sharedWaveSpeeds[0]*dt/dx;
        if ( cfl > 1.0) {
            printf("cfl larger than 1.0 in reduce_max. cfl = %E. \n" , cfl);
            printf("array size: %d\n" , size);
            printf("local grid id: %d\n", param.id);
            printf("direction: %d\n", direction);
            printf("x-size: %d; y-size: %d\n", param.cellsX, param.cellsY);
            printf("address of globalWaves at later reduce_max: %p.\n", (void*)globalWaves);
            printf("address of param.qNew at later reduce_max: %p.\n", (void*)param.qNew);
            printf("address of param.q_tmp at later reduce_max: %p.\n", (void*)param.q_tmp);
            printf("address of param.coefficients at later reduce_max: %p.\n", (void*)param.coefficients);

            printf("globalWaves: \n");
            for (int i = 0; i < size; ++i) {
                printf("%d: %E\n", i, globalWaves[i]*dt/dx);
            }
            printf("sharedWaveSpeeds: \n");
            for (int i = 0; i < size; ++i) {
                printf("%d: %E\n", i, sharedWaveSpeeds[i]*dt/dx);
            }
        }

        globalWaves[blockIdx.x] = sharedWaveSpeeds[0];
        if (result)
            *result = sharedWaveSpeeds[0]*dt/dx;
    }
}

#endif

