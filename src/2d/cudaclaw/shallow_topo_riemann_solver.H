#ifndef __SHALLOW_WATER_RIEMANN_SOLVER_H__
#define __SHALLOW_WATER_RIEMANN_SOLVER_H__

#include "real.H"
#include "common.H"
#include "staticArray.H"
#include "fused_solvers_headers.H"
#include <cmath>
#include <cassert>

#ifdef CUDA
extern "C" 
{
    __device__
    void riemann_shallow_topo_x(
            const real* const q_left, const real* const q_right,
            const real* const aux_left, const real* const aux_right,
            volatile real* wave, volatile real* waveSpeeds, 
            const int* blocksize);
    __device__
    void riemann_shallow_topo_y(
            const real* const q_left, const real* const q_right,
            const real* const aux_left, const real* const aux_right,
            volatile real* wave, volatile real* waveSpeeds,
            const int* blocksize);
}
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
struct shallow_topo_horizontal
{
    // template<const size_t N>
    __device__ void 
    operator() (const real* const q_left, const real* const q_right, 
            // const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const real* const u_left, const real* const u_right,
            const int row, const int col, const int numStates, const int numWaves,
            volatile real* const wave, volatile real* const waveSpeeds,
            const int blocksize)	
    {
//         static_assert(N > 0, 
//                 "Number of coefficients must be larger than 0 for shallow_topo_horizontal");
        int i = 0;
        if (q_left[i] < 0.0)
            printf("q_left[%d] < 0.0 in x sweep: %E\n",i,q_left[i]);
        if (q_right[i] < 0.0)
            printf("q_right[%d] < 0.0 in x sweep: %E\n",i,q_right[i]);
        riemann_shallow_topo_x(
                // q_left, q_right, &(u_left[0]), &(u_right[0]), 
                q_left, q_right, u_left, u_right, 
                wave, waveSpeeds,
                &blocksize);

    }
};

struct shallow_topo_vertical
{
    // template<const size_t N>
    __device__ void 
    operator() (const real* const q_left, const real* const q_right, 
            // const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const real* const u_left, const real* const u_right,
            const int row, const int col, const int numStates, const int numWaves,
            volatile real* const wave, volatile real* const waveSpeeds,	
            const int blocksize)	
    {
//         static_assert(N > 0, 
//                 "Number of coefficients must be larger than 0 for shallow_topo_vertical");
        int i = 0;
        if (q_left[i] < 0.0)
            printf("q_left[%d] < 0.0 in y sweep: %E\n",i,q_left[i]);
        if (q_right[i] < 0.0)
            printf("q_right[%d] < 0.0 in y sweep: %E\n",i,q_right[i]);
        riemann_shallow_topo_y(
                // q_left, q_right, &(u_left[0]), &(u_right[0]), 
                q_left, q_right, u_left, u_right, 
                wave, waveSpeeds,
                &blocksize);
    }
};

#endif
#endif
