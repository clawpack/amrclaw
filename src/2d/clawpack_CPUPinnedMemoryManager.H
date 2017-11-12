#ifndef CLAWPACK_CPUPINNEDMEMORYMANAGER_H
#define CLAWPACK_CPUPINNEDMEMORYMANAGER_H

#ifdef CUDA

#include <cstddef>
#include <set>
#include <vector>

#include "clawpack_CPUPinnedMemoryManager.H"
#include "clawpack_MemoryManager.H"

/**
* \brief A Concrete Class for Dynamic unpageable (pinned) CPU Memory Management
* This is a coalescing memory manager.  It allocates (possibly) large
* chunks of heap space and apportions it out as requested.  It merges
* together neighboring chunks on each free().
*
* This class only manage 1 type of dynamic memory,
* CPU unpageable memory, from the memory pool, by using
* alloc_pinned(...) and free_pinned(...)
*
* This class should not be used directly. Instead, use clawpack_mempool_alloc_pinned
*/

class CPUPinnedMemoryManager
    :
    public MemoryManager
{
public:
    /**
    * \brief Construct a coalescing memory manager.  hunk_size is the
    * minimum size of hunks of memory to allocate from the heap.
    * If hunk_size == 0 we use DefaultHunkSize as specified below.
    */
    CPUPinnedMemoryManager (size_t hunk_size = 0);
    
    //! The destructor.
    virtual ~CPUPinnedMemoryManager () override;

    //! Allocate some memory.
    virtual void* alloc (size_t nbytes) override {
        void* pt = 0;
        std::cout << "CPUPinnedMemoryManager is for pinned memory. CPUPinnedMemoryManager::alloc() should not be called" << std::endl; 
        std::abort();
        return pt;
    }

    virtual void* alloc_pinned (std::size_t sz) override;

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free (void* ap) override {
        std::cout << "CPUPinnedMemoryManager is for pinned memory. CPUPinnedMemoryManager::free() should not be called" << std::endl; 
        std::abort();
    }

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free_pinned (void* ap) override;

    virtual void* alloc_device (size_t nbytes, int device_id) override {
        void* pt = 0;
        std::cout << "CPUPinnedMemoryManager is a CPU memory pool. CPUPinnedMemoryManager::alloc_device() should not be called" << std::endl;
        std::abort();
        return pt;
    }
#ifdef CUDA_ARRAY
    /*
     * do nothing for now
     */ 
    virtual void* alloc_device_2d (std::size_t& _pitch, std::size_t _isize, std::size_t _jsize, int device_id) override {
        void* pt = 0;
        std::cout << "CPUPinnedMemoryManager is a CPU memory pool. CPUPinnedMemoryManager::alloc_device() should not be called" << std::endl; 
        std::abort();
        return pt;
    }
#endif
    virtual void free_device (void* ap, int dev_id) override {
        std::cout << "CPUPinnedMemoryManager is a CPU memory pool. CPUPinnedMemoryManager::free_device() should not be called" << std::endl; 
        std::abort();
    }

    //! The current amount of heap space used by the CPUPinnedMemoryManager object.
    size_t heap_space_used () const;

    //! The default memory hunk size to grab from the heap.
    enum { DefaultHunkSize = 1024*1024*1024};

    /**
    * \brief Given a minimum required memory block of sz bytes, this returns
    * the next largest memory block of size that will align to align_size bytes
    * Override align in MemoryManager to use the align_size overridden in this class
    */
    static std::size_t align (std::size_t sz);

protected:
    //! The nodes in our free list and block list.
    class Node
    {
    public:
        //! The default constructor.
        Node ()
            :
            m_block(0), m_size(0) {}
        
        //! Another constructor.
        Node (void* block, size_t size)
            :
            m_block(block), m_size(size) {}

        //! The copy constructor.
        Node (const Node& rhs)
            :
            m_block(rhs.m_block), m_size(rhs.m_size) {}

        //! The copy assignment constructor.
        Node& operator= (const Node& rhs)
        {
            m_block = rhs.m_block;
            m_size  = rhs.m_size;
            return *this;
        }

        //! The "less-than" operator.
        bool operator< (const Node& rhs) const
        {
            return m_block < rhs.m_block;
        }

        //! The equality operator. 
        bool operator== (const Node& rhs) const
        {
            return m_block == rhs.m_block;
        }

        //! The block address.
        void* block () const { return m_block; }

        //! Set block address.
        void block (void* blk) { m_block = blk; }

        //! The size of the memory block.
        size_t size () const { return m_size; }

        //! Set size.
        void size (size_t sz) { m_size = sz; }


    private:
        //! The block of memory we reference.
        void* m_block;
        //! The size of the block we represent.
        size_t m_size;
    };

    /**
    * \brief The type of our freelist and blocklist.
    * We use a set sorted from lo to hi memory addresses.
    */
    typedef std::set<Node> NL;

    //! The list of blocks allocated via ::operator new().
    std::vector<void*> m_alloc;

    /**
    * \brief The free list of allocated but not currently used blocks.
    * Maintained in lo to hi memory sorted order.
    */
    NL m_freelist;

    /**
    * \brief The list of busy blocks.
    * A block is either on the freelist or on the blocklist, but not on both.
    */
    NL m_busylist;
    //! The minimal size of hunks to request via ::operator new().
    size_t m_hunk;
    //! The amount of heap space currently allocated.
    size_t m_used;

    //! overwrite align_size for GPU memory performance
    static const unsigned int align_size = 256;
    // How many memory in total we allow the CPUPinnedMemoryManager to allocate
    static const size_t max_heap_size = (size_t) 1024 * 1024 * 1024 * 32;

private:
    //! Disallowed.
    CPUPinnedMemoryManager (const CPUPinnedMemoryManager& rhs);
    CPUPinnedMemoryManager& operator= (const CPUPinnedMemoryManager& rhs);
};


#endif // CUDA
#endif /*CLAWPACK_CPUPINNEDMEMORYMANAGER_H*/
