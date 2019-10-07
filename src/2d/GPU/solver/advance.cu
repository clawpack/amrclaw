#include "advance.H"
#include "real.H"
#include "common.H"

// some global constants that need to be set from Fortran


#ifdef GEOCLAW
// the dry_tolerance parameter in setrun.py 
__device__ __constant__  real dry_tolerance = claw_nan; // default value. Can be set from Fortran

void set_dry_tolerance (real value)
{
    cudaMemcpyToSymbol(dry_tolerance,&value,sizeof(real));
}
#endif




