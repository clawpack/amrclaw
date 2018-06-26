#ifndef __SHAREDMEM_H__
#define __SHAREDMEM_H__

#ifdef CUDA
// Prevent instantiation of this class by only declaring it but not defining it
template <typename T> class SharedMemory;

// specializations for different types

template <>
class SharedMemory <int>
{
    public:
        __device__ int *getPointer()
        {
            extern __shared__ int s_int[];
            return s_int;
        }

        __device__ int *getPointer(size_t offset)
        {
            extern __shared__ int s_int[];
            return &s_int[offset];
        }
};

template <>
class SharedMemory <unsigned int>
{
    public:
        __device__ unsigned int *getPointer()
        {
            extern __shared__ unsigned int s_uint[];
            return s_uint;
        }

        __device__ unsigned int *getPointer(size_t offset)
        {
            extern __shared__ unsigned int s_uint[];
            return &s_uint[offset];
        }
};

template <>
class SharedMemory <float>
{
    public:
        __device__ float *getPointer()
        {
            extern __shared__ float s_float[];
            return s_float;
        }

        __device__ float *getPointer(size_t offset)
        {
            extern __shared__ float s_float[];
            return &s_float[offset];
        }
};

template <>
class SharedMemory <double>
{
    public:

        __device__ double *getPointer()
        {
            extern __shared__ double s_double[];
            return s_double;
        }

        __device__ double *getPointer(size_t offset)
        {
            extern __shared__ double s_double[];
            return &s_double[offset];
        }
};

#endif
#endif
