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
    void riemann_shallow_topo(const int ixy,
            const real* const q_left, const real* const q_right,
            const real* const aux_left, const real* const aux_right,
            real* wave, real* waveSpeeds);
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
    template<const size_t N>
    __device__ void 
    operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)	
    {
        static_assert(N == 1, 
                "Number of coefficients must be 1 for shallow_topo_horizontal");
        riemann_shallow_topo(1,
                q_left, q_right, &(u_left[0]), &(u_right[0]), 
                wave, waveSpeeds
                );

    }
};

struct shallow_topo_vertical
{
    template<const size_t N>
    __device__ void 
    operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)	
    {
        static_assert(N == 1, 
                "Number of coefficients must be 1 for shallow_topo_vertical");
        riemann_shallow_topo(2,
                q_left, q_right, &(u_left[0]), &(u_right[0]), 
                wave, waveSpeeds
                );
    }
};

#endif
#endif
