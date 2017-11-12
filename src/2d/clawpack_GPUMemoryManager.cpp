#ifdef CUDA

#include <utility>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <mutex>

#include "clawpack_GPUMemoryManager.H"

#include "clawpack_CUDA_helper.H"

#include <cuda_runtime_api.h>
#include <cuda.h>

static bool local_verbose = false;

static std::mutex NL_mutex;
static std::mutex TagMemoryBlocksContainer_mutex;


GPUMemoryManager::GPUMemoryManager (int dev_id, size_t hunk_size)
    :device_id(dev_id)
{
    //
    // Force alignment of hunksize.
    //
    m_hunk = MemoryManager::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);
    m_used = 0;

    assert(m_hunk >= hunk_size);
    assert(m_hunk%MemoryManager::align_size == 0);
}

GPUMemoryManager::~GPUMemoryManager ()
{
    CUresult cudaResult;
    CUcontext cuctx;
    cudaResult = cuCtxGetCurrent(&cuctx); 
    if (local_verbose) 
        std::cout << "Maximum GPU heap memory used during the execution: " << m_used << std::endl;
    if (cudaResult != CUDA_ERROR_DEINITIALIZED) {
        checkCudaErrors(cudaSetDevice(device_id));
        for (unsigned int i = 0, N = m_alloc.size(); i < N; ++i)
            checkCudaErrors(cudaFree(m_alloc[i]));
        if (local_verbose) 
            std::cout << "Free all memory requested by GPUMemoryManager back to the system." << std::endl;
    } 
    // otherwise device has been reset and 
    // all device memory has been free. So we should not call cudaFree
}

void*
GPUMemoryManager::alloc_device (size_t nbytes, int dev_id)
{
    std::lock_guard<std::mutex> alloc_device_guard(NL_mutex);
    assert(device_id == dev_id);
    nbytes = GPUMemoryManager::align(nbytes == 0 ? 1 : nbytes);
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

        checkCudaErrors(cudaSetDevice(device_id));
        checkCudaErrors(cudaMalloc(&vp, N));
        m_used += N;
        if (local_verbose) {
            std::cout << "Allocate extra " << N << " Bytes of memory from the system. " << std::endl;
            std::cout << "Total memory allocated from the system: " << m_used << std::endl;
        }

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
        std::cout << "Claim " << nbytes << " Bytes of memory at: " << vp << " from GPUMemoryManager. " << std::endl;
    }

    assert(!(vp == 0));

    return vp;
}

void* 
GPUMemoryManager::alloc_device (size_t nbytes, int tag, int dev_id)
{
    void* vp = alloc_device(nbytes, dev_id);

    /*
     * Record this memory block and its tag here
     */
    std::lock_guard<std::mutex> alloc_device_tag_guard(TagMemoryBlocksContainer_mutex);
    TagMemoryBlocks idx_mem(tag);
    auto search = m_memory_blocks.find(idx_mem);
    if (search == m_memory_blocks.end()) { // new tag
        idx_mem.push_back(vp);
        m_memory_blocks.insert(idx_mem);
    }
    else { // already has memory blocks group associated with this tag in the container
        search->push_back(vp);
    }

    if (local_verbose) {
        std::cout << "Add pt: " << vp << " to m_memory_blocks. It's associated with tag: " << tag << std::endl;
        std::cout << "m_memory_blocks[tag] becomes: " << std::endl;
        auto search2 = m_memory_blocks.find(TagMemoryBlocks(tag));
        search2->print();
    }

    return vp;

}


void
GPUMemoryManager::free_device_tag (int tag, int dev_id)
{
    std::lock_guard<std::mutex> free_device_tag_guard(TagMemoryBlocksContainer_mutex);
    auto search = m_memory_blocks.find(TagMemoryBlocks(tag));
    assert( search != m_memory_blocks.end() );
    std::vector<void*>& all_blocks = *(search->data);
    for (auto it = all_blocks.begin(); it != all_blocks.end(); ++it) {
        if (local_verbose) {
            std::cout << "\tFree memory block at address: " << *it << std::endl;
        }
        free_device(*it, dev_id);
    }
    // remove idx from m_memory_blocks
    m_memory_blocks.erase(search);
}
void
GPUMemoryManager::free_device (void* vp, int dev_id)
{
    assert( device_id == dev_id );
    std::lock_guard<std::mutex> free_device_guard(NL_mutex);
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
        std::cout << "Put memory block: " << vp << " of " << busy_it->size() << " back to GPUMemoryManager. " << std::endl;
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
                if (addr == m_alloc[i]) merge = false; 
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
            if (addr == m_alloc[i])  merge = false;
        }
        //
        // Ditto the above comment.
        //
        if (merge) {
            Node* node = const_cast<Node*>(&(*free_it));
            assert(!(node == 0));
            node->size((*free_it).size() + (*hi_it).size());
            m_freelist.erase(hi_it);
        }
    }
}

std::size_t
GPUMemoryManager::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}


size_t
GPUMemoryManager::heap_space_used () const
{
    return m_used;
}


#endif //ifdef CUDA
