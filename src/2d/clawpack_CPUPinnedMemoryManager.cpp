#ifdef CUDA

#include <utility>
#include <cstring>
#include <cassert>
#include <iostream>
#include <mutex>

#include "clawpack_CPUPinnedMemoryManager.H"

#include "clawpack_CUDA_helper.H"

#include <cuda_runtime_api.h>
#include <cuda.h>

static bool local_verbose = false;
static std::mutex NL_mutex;

CPUPinnedMemoryManager::CPUPinnedMemoryManager (size_t hunk_size)
{
    //
    // Force alignment of hunksize.
    //
    m_hunk = CPUPinnedMemoryManager::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);
    m_used = 0;

    assert(m_hunk >= hunk_size);
    assert(m_hunk%CPUPinnedMemoryManager::align_size == 0);
}

CPUPinnedMemoryManager::~CPUPinnedMemoryManager ()
{
    CUresult cudaResult;
    CUcontext cuctx;
    cudaResult = cuCtxGetCurrent(&cuctx); 
    if (cudaResult != CUDA_ERROR_DEINITIALIZED) {
        for (unsigned int i = 0, N = m_alloc.size(); i < N; i++)
            checkCudaErrors(cudaFreeHost(m_alloc[i]));
    }
    // otherwise device has been reset and 
    // all memory has been free. So we should not call cudaFreeHost
    std::cout << "Summary from CPUPinnedMemoryManager: "<< std::endl;
    std::cout << "  Maximum CPU heap memory taken by the CPUPinnedMemoryManager during the execution: " << ((double) m_used)/(1024*1024*1024) << " GB." << std::endl;
    if (!m_busylist.empty()) 
        std::cout << "  Warning: memory leak might have occured in CPUPinnedMemoryManager!!!" << std::endl;
    else
        std::cout << "  All memory allocated through CPUPinnedMemoryManager has been freed." << std::endl;
}

void*
CPUPinnedMemoryManager::alloc_pinned (size_t nbytes)
{
    std::lock_guard<std::mutex> alloc_pinned_guard(NL_mutex);
    nbytes = CPUPinnedMemoryManager::align(nbytes == 0 ? 1 : nbytes);
    //
    // Find node in freelist at lowest memory address that'll satisfy request.
    //
    NL::iterator free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it)
        if ((*free_it).size() >= nbytes)
            break;

    void* vp = 0;

    if (free_it == m_freelist.end())
    {
        const size_t N = nbytes < m_hunk ? m_hunk : nbytes;

        if ( (m_used + N) > max_heap_size ) {
            std::cout << "Not enough GPU memory available in GPUMemoryManager." << std::endl;
            std::abort();
        }

        checkCudaErrors(cudaHostAlloc(&vp, N, cudaHostAllocDefault));

        m_used += N;

        m_alloc.push_back(vp);

        if (nbytes < m_hunk)
        {
            //
            // Add leftover chunk to free list.
            //
            // Insert with a hint -- should be largest block in the set.
            //
            void* block = static_cast<char*>(vp) + nbytes;

            m_freelist.insert(m_freelist.end(), Node(block, m_hunk-nbytes));
        }
    }
    else
    {
        assert((*free_it).size() >= nbytes);
        assert(m_busylist.find(*free_it) == m_busylist.end());

        vp = (*free_it).block();

        if ((*free_it).size() > nbytes)
        {
            //
            // Insert remainder of free block back into freelist.
            //
            // Insert with a hint -- right after the current block being split.
            //
            Node freeblock = *free_it;

            freeblock.size(freeblock.size() - nbytes);

            freeblock.block(static_cast<char*>(vp) + nbytes);

            m_freelist.insert(free_it, freeblock);
        }

        m_freelist.erase(free_it);
    }

    m_busylist.insert(Node(vp, nbytes));

    if (local_verbose) {
        std::cout << "Claim " << nbytes << " Bytes of memory at: " << (uintptr_t)vp << " from CPUPinnedMemoryManager. " << std::endl;
    }

    assert(!(vp == 0));

    return vp;
}

void
CPUPinnedMemoryManager::free_pinned (void* vp)
{
    std::lock_guard<std::mutex> free_pinned_guard(NL_mutex);
    if (vp == 0)
        //
        // Allow calls with NULL as allowed by C++ delete.
        //
        return;
    //
    // `vp' had better be in the busy list.
    //
    NL::iterator busy_it = m_busylist.find(Node(vp,0));

    assert(!(busy_it == m_busylist.end()));
    assert(m_freelist.find(*busy_it) == m_freelist.end());
    //
    // Put free'd block on free list and save iterator to insert()ed position.
    //
    std::pair<NL::iterator,bool> pair_it = m_freelist.insert(*busy_it);

    if (local_verbose) {
        std::cout << "Put pinned CPU memory block: " << (uintptr_t)vp << " of " << busy_it->size() << " back to CPUPinnedMemoryManager. " << std::endl;
    }

    assert(pair_it.second == true);

    NL::iterator free_it = pair_it.first;

    assert(free_it != m_freelist.end() && (*free_it).block() == (*busy_it).block());
    //
    // And remove from busy list.
    //
    m_busylist.erase(busy_it);
    //
    // Coalesce freeblock(s) on lo and hi side of this block.
    //
    if (!(free_it == m_freelist.begin()))
    {
        NL::iterator lo_it = free_it;

        --lo_it;

        void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

        if (addr == (*free_it).block())
        {
            bool merge = true;
            for (unsigned int i = 0, N = m_alloc.size(); i < N; i++) {
                // Don't merge two nodes if the merge will give a node 
                // whose memory block crosses the hunk boundary
                if (addr == m_alloc[i]) {
                    merge = false;
                }
            }
            //
            // This cast is needed as iterators to set return const values;
            // i.e. we can't legally change an element of a set.
            // In this case I want to change the size() of a block
            // in the freelist.  Since size() is not used in the ordering
            // relations in the set, this won't effect the order;
            // i.e. it won't muck up the ordering of elements in the set.
            // I don't want to have to remove the element from the set and
            // then reinsert it with a different size() as it'll just go
            // back into the same place in the set.
            //
            if (merge) {
                Node* node = const_cast<Node*>(&(*lo_it));
                assert(!(node == 0));
                node->size((*lo_it).size() + (*free_it).size());
                m_freelist.erase(free_it);
                free_it = lo_it;
            }
        }
    }

    NL::iterator hi_it = free_it;

    void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

    if (++hi_it != m_freelist.end() && addr == (*hi_it).block())
    {
        bool merge = true;
        for (unsigned int i = 0, N = m_alloc.size(); i < N; i++) {
            // Don't merge two nodes if the merge will give a node 
            // whose memory block crosses the hunk boundary
            if (addr == m_alloc[i]) {
                merge = false;
            }
        }
        if (merge) {
            //
            // Ditto the above comment.
            //
            Node* node = const_cast<Node*>(&(*free_it));
            assert(!(node == 0));
            node->size((*free_it).size() + (*hi_it).size());
            m_freelist.erase(hi_it);
        }
    }
}

std::size_t
CPUPinnedMemoryManager::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}


size_t
CPUPinnedMemoryManager::heap_space_used () const
{
    return m_used;
}


#endif
