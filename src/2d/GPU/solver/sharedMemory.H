#ifndef __SHAREDMEM_H__
#define __SHAREDMEM_H__

#ifdef CUDA
template <typename T> 
class SharedMemory
{
    public:
        __device__ T *getPointer()
        {
            extern __shared__ T s_T[];
            return s_T;
        }

        __device__ T *getPointer(size_t offset)
        {
            extern __shared__ T s_T[];
            return &s_T[offset];
        }
};


#endif
#endif
