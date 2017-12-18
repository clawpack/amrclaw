#ifndef CLAWPACK_MEMPOOL_H
#define CLAWPACK_MEMPOOL_H

extern "C" {
    void clawpack_mempool_init ();
    void clawpack_mempool_finalize();
    // void* clawpack_mempool_alloc (size_t n);
    // void  clawpack_mempool_free (void* p);
#ifdef CUDA
    void* clawpack_mempool_alloc_pinned (size_t n);
    void  clawpack_mempool_free_pinned (void* p);
    void* clawpack_mempool_alloc_gpu (size_t n, int dev_id);
    void  clawpack_mempool_free_gpu (void* p, int dev_id);
    // void* clawpack_mempool_alloc_gpu_hold (size_t n, int tag, int dev_id);
    int get_num_devices_used(); 
    int which_device_used(); 
#endif
    // void  clawpack_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  // min, max & tot in MB
    void  clawpack_real_array_init (double* p, size_t nelems); 
    void  clawpack_array_init_snan (double* p, size_t nelems);
}
#ifdef CUDA
    // void  clawpack_mempool_release_gpu (int tag, int dev_id);
    // void  CUDART_CB cudaCallback_release_gpu(cudaStream_t event, cudaError_t status, void *data);
#endif


#endif
