#ifndef __ACOUSTICS_RIEMANN_SOLVER_H__
#define __ACOUSTICS_RIEMANN_SOLVER_H__

#include "real.H"
#include "common.H"
#include "staticArray.H"
#include "fused_solvers_headers.H"
#include <cmath>
#include <cassert>

#ifdef CUDA

// This is for homogeneous acoustics.
// TODO: should let the user specify bulk and rho in case directory
constexpr real bulk = 4.0;
constexpr real rho = 1.0;

extern "C" 
{
    __device__
    void riemann_acoustics_homo_2d_fortran(const int ixy, const int num_eqn, 
            const int num_aux, const int num_waves,
            const real* const q_left, const real* const q_right,
            const real* const aux_left, const real* const aux_right,
            real* wave, real* waveSpeeds, 
            const real bulk, const real rho);
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
struct acoustics_homo_2d_horizontal
{
    // u_left and u_right are empty and not used here
    template <const size_t N>
    __device__ void operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)
    {
        static_assert(N == 0, 
                "Number of coefficients must be 0 for homogeneous acoustics");

        real wave_local[NWAVES*NEQNS];
        real waveSpeeds_local[NWAVES];
        riemann_acoustics_homo_2d_fortran(1, numStates, N, numWaves, 
                q_left, q_right, NULL, NULL, 
                wave_local, waveSpeeds_local, bulk, rho);

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, waveSpeeds_local[0]);
        setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x, waveSpeeds_local[1]);

        setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, wave_local[0]);
        setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, wave_local[1]);
        setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, wave_local[2]);

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, wave_local[3]);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, wave_local[4]);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, wave_local[5]);

    }
};

struct acoustics_homo_2d_horizontal_warp
{
    // u_left and u_right are empty and not used here
    __device__ void operator() (real* q_left, real* q_right, real* u_left, real* u_right,
            int row, int col, int numStates, int numWaves,
            real* wave, real* waveSpeeds, real* apdq, real* amdq)
    {
        assert(NCOEFFS == 0); // this is the case for homogeneous acoustics
        real c = sqrt(bulk/rho);
        real z = c*rho;

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c);
        setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c);

        real alpha1 = ( q_left[0] - q_right[0] + z*(q_right[1] - q_left[1])) / (2*z);
        real alpha2 = ( q_right[0] - q_left[0] + z*(q_right[1] - q_left[1])) / (2*z);

        amdq[0] = (-alpha1*z)*(-c);
        amdq[1] = alpha1*(-c);
        amdq[2] = claw_zero;

        apdq[0] = alpha2*z*c;
        apdq[1] = alpha2*c;
        apdq[2] = claw_zero;

    }
};

struct acoustics_homo_2d_vertical
{
    template<const size_t N>
    __device__ void operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)										// output
    {
        static_assert(N == 0, 
                "Number of coefficients must be 0 for homogeneous acoustics");

        real wave_local[NWAVES*NEQNS];
        real waveSpeeds_local[NWAVES];
        riemann_acoustics_homo_2d_fortran(2, numStates, N, numWaves, 
                q_left, q_right, NULL, NULL, 
                wave_local, waveSpeeds_local, bulk, rho);

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, waveSpeeds_local[0]);
        setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x, waveSpeeds_local[1]);

        setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, wave_local[0]);
        setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, wave_local[1]);
        setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, wave_local[2]);

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, wave_local[3]);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, wave_local[4]);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, wave_local[5]);

    }
};

struct acoustics_heterogeneous_2d_horizontal
{
    template<const size_t N>
    __device__ void operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)
    {
        static_assert(N>0, 
                "N must be larger than 0: heterogeneous acoustic riemann solver requires material property");
        real rho_l = u_left[0];
        real bulk_l = u_left[1];

        real rho_r = u_right[0];
        real bulk_r = u_right[1];

        real c_l = sqrt(bulk_l/rho_l); // sound speed
        real z_l = c_l*rho_l; // impedance

        real c_r = sqrt(bulk_r/rho_r);
        real z_r = c_r*rho_r;

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c_l);
        setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c_r);

        // ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
        real alpha1 = 
            ( q_left[0] - q_right[0] + z_r*(q_right[1] - q_left[1])) / (z_l+z_r); 
        // (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
        real alpha2 = 
            ( q_right[0] - q_left[0] + z_l*(q_right[1] - q_left[1])) / (z_l+z_r); 

        setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z_l);
        setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha1);
        setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z_r);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
    }
};

struct acoustics_heterogeneous_2d_vertical
{
    template<const size_t N>
    __device__ void operator() (const real* const q_left, const real* const q_right, 
            const staticArrayReal<N>& u_left, const staticArrayReal<N>& u_right,
            const int row, const int col, const int numStates, const int numWaves,
            real* const wave, real* const waveSpeeds)
    {
        static_assert(N>0, 
                "N must be larger than 0: heterogeneous acoustic riemann solver requires material property");
        real rho_l = u_left[0];
        real bulk_l = u_left[1];

        real rho_r = u_right[0];
        real bulk_r = u_right[1];

        real c_l = sqrt(bulk_l/rho_l);	// sound speed
        real z_l = c_l*rho_l;			// impedance

        real c_r = sqrt(bulk_r/rho_r);
        real z_r = c_r*rho_r;

        setSharedSpeed(waveSpeeds, row, col, 0, numWaves, blockDim.y, blockDim.x, -c_l);
        setSharedSpeed(waveSpeeds, row, col, 1, numWaves, blockDim.y, blockDim.x,  c_r);

        // ( -(pr-pl) + zr(vr-vl) )/ (zl+zr)
        real alpha1 = ( q_left[0] - q_right[0] + z_r*(q_right[2] - q_left[2])) / (z_l+z_r);
        // (  (pr-pl) + zl(vr-vl) )/ (zl+zr)
        real alpha2 = ( q_right[0] - q_left[0] + z_l*(q_right[2] - q_left[2])) / (z_l+z_r);

        setSharedWave(wave, row, col, 0, 0, numWaves, numStates, blockDim.y, blockDim.x, -alpha1*z_l);
        setSharedWave(wave, row, col, 0, 1, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
        setSharedWave(wave, row, col, 0, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha1);

        setSharedWave(wave, row, col, 1, 0, numWaves, numStates, blockDim.y, blockDim.x, alpha2*z_r);
        setSharedWave(wave, row, col, 1, 1, numWaves, numStates, blockDim.y, blockDim.x, claw_zero);
        setSharedWave(wave, row, col, 1, 2, numWaves, numStates, blockDim.y, blockDim.x, alpha2);
    }
};

#endif
#endif
