#ifndef CLAWPACK_DEVICE_H
#define CLAWPACK_DEVICE_H

#ifdef CUDA
extern "C" {
    void get_cuda_stream(const int id, const int device_id, cudaStream_t* stream); 
}
#endif

#endif
