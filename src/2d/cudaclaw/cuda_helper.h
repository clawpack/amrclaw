#ifndef __CUDA_HELPER_H__
#define __CUDA_HELPER_H__

#define CHKERR() if (cudaPeekAtLastError()) { printf("CUDA error in %s at %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(cudaPeekAtLastError())); return cudaPeekAtLastError(); }

#define CHKERRQ(n) if (n) { printf("CUDA error in %s at %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(n)); }

#define CHKERR_ABORT() if (cudaPeekAtLastError()) { printf("CUDA error in %s at %d: %s\n", __FILE__, __LINE__, cudaGetErrorString(cudaPeekAtLastError())); }

#define CHKERRQ_ABORT(n) if (n) { printf("CUDA error in %s at %d: %s\nIgnoring!\n", __FILE__, __LINE__, cudaGetErrorString(n)); }

#endif // __CUDA_HELPER_H__
