#ifndef CLAWPACK_GPUMEMORYMANAGER_H
#define CLAWPACK_GPUMEMORYMANAGER_H

#ifdef CUDA

#include <cstddef>
#include <set>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "clawpack_MemoryManager.H"

/**
* \brief A Concrete Class for Dynamic **GPU*
* * Memory Management
* This is a coalescing memory manager.  It allocates (possibly) large
* chunks of heap space and apportions it out as requested.  It merges
* together neighboring chunks on each free().
*
* This class can only allocate 1 type of dynamic memory: GPU memory
* alloc_device(...) and free_device(...) can be called to get memory from 
* and put it back to the memory pool
*
* This class should not be used directly. Instead, use clawpack_mempool_alloc_gpu
*/

class GPUMemoryManager
    :
    public MemoryManager
{
public:
    /**
    * \brief Construct a coalescing memory manager.  hunk_size is the
    * minimum size of hunks of memory to allocate from the heap.
    * If hunk_size == 0 we use DefaultHunkSize as specified below.
    */
    GPUMemoryManager (int dev_id, size_t hunk_size = 0);
    
    //! The destructor.
    virtual ~GPUMemoryManager () override;

    //! Allocate some memory.
    virtual void* alloc (size_t nbytes) override {
        void* pt = 0;
        std::cout << "GPUMemoryManager is a GPU memory pool. GPUMemoryManager::alloc() should not be called" << std::endl;
        std::abort();
        return pt;
    }

    virtual void* alloc_pinned (std::size_t sz) override { 
        void* pt = 0;
        std::cout << "GPUMemoryManager is a GPU memory pool. GPUMemoryManager::alloc_pinned() should not be called" << std::endl;
        std::abort();
        return pt;
    }

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free (void* ap) override {
        std::cout << "GPUMemoryManager is a GPU memory pool. GPUMemoryManager::free() should not be called" << std::endl;
        std::abort();
    }

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free_pinned (void* ap) override { 
        std::cout << "GPUMemoryManager is a GPU memory pool. GPUMemoryManager::free_pinned() should not be called" << std::endl;
        std::abort();
    }

    virtual void* alloc_device (size_t nbytes, int dev_id) override; 

    void* alloc_device (size_t nbytes, int tag, int dev_id);
#ifdef CUDA_ARRAY
    /*
     * do nothing for now
     */ 
    virtual void* alloc_device_2d (std::size_t& _pitch, std::size_t _isize, std::size_t _jsize, int dev_id) override { void* pt = 0; return pt;};
#endif
    virtual void free_device (void* ap, int dev_id) override; 
    void free_device_tag (int tag, int dev_id);

    //! The current amount of heap space used by the GPUMemoryManager object.
    size_t heap_space_used () const;

    void useDevice(int dev_id) {device_id = dev_id;}

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
    // How many memory in total we allow the GPUMemoryManager to allocate
    static const size_t max_heap_size = (size_t) 1024 * 1024 * 1024 * 12;

    // which GPU's memory does this GPUMemoryManager manage
    int device_id;

    class TagMemoryBlocks {
    public:
        int tag;
        std::vector<void*>* data;

        TagMemoryBlocks() { data = new std::vector<void*>(); }
        TagMemoryBlocks(int tag_): tag(tag_) { data = new std::vector<void*>(); }
        ~TagMemoryBlocks() {
            delete data;
        }

        // copy constructor
        TagMemoryBlocks( const TagMemoryBlocks& src) {
            tag = src.tag;
            data = new std::vector<void*>(*(src.data));
        }

        // assignment operator
        TagMemoryBlocks& operator= (const TagMemoryBlocks& other) = delete;

        void push_back(void* memory_block) const {
            data->push_back(memory_block);
        }

        // TODO: remove this 
        // for debuggging
        void print() const {
            std::cout << "tag: " << tag << std::endl;
            std::cout << "data size: " << data->size() << std::endl;
            std::cout << "data: " << std::endl;
            for (int j = 0; j < data->size(); ++j) {
                std::cout << "\t" << (*data)[j];
            }
            std::cout << std::endl;
        }
    };
    class TagMemoryBlocksComparator {
        public:
        bool operator()( const TagMemoryBlocks& lhs, const TagMemoryBlocks& rhs ) const {
            return lhs.tag < rhs.tag;
        }
    };

    typedef std::set<TagMemoryBlocks, TagMemoryBlocksComparator> TagMemoryBlocksContainer;

    TagMemoryBlocksContainer m_memory_blocks;



private:
    //! Disallowed.
    GPUMemoryManager (const GPUMemoryManager& rhs);
    GPUMemoryManager& operator= (const GPUMemoryManager& rhs);
};




#endif // CUDA
#endif /*CLAWPACK_GPUMEMORYMANAGER_H*/
