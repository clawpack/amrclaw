#ifndef CLAWPACK_MEMORYMANAGER_H
#define CLAWPACK_MEMORYMANAGER_H

#include <cstddef>

// class Arena;

// Arena* The_Arena ();

/**
* \brief 
* A virtual base class for objects that manage their own dynamic
* memory allocation.  
*/

class MemoryManager
{
public:

    virtual ~MemoryManager();
    /**
    * Allocate a dynamic memory block
    * \param sz size of the memory request
    * \return a pointer to the allocated memory
    */
    virtual void* alloc (std::size_t sz) = 0;
    /**
    * Allocate a dynamic pinned memory 
    * \param sz size of the memory request
    * \return a pointer to the allocated memory
    */
    virtual void* alloc_pinned (std::size_t sz) = 0;
    /**
    * \brief A pure virtual function for deleting the memory block pointed to by pt
    */
    virtual void free (void* pt) = 0;
    /**
    * \brief A pure virtual function for deleting the memory block pointed to by pt
    */
    virtual void free_pinned (void* pt) = 0;
    /**
    * Allocate a dynamic device memory block
    * \param sz size of the memory request
    * \return a pointer to the allocated memory
    */
    virtual void* alloc_device (std::size_t sz, int device_id) = 0;
    /**
    * \brief A pure virtual function for deleting the device memory block pointed to by pt
    */
    virtual void free_device (void* pt, int device_id) = 0;
    /**
    * \brief Given a minimum required memory size of sz bytes, this returns
    * the next largest memory size that will align to align_size bytes
    */
    static std::size_t align (std::size_t sz);

protected:
    static const unsigned int align_size = 16;
};


#endif /*CLAWPACK_MEMORYMANAGER_H*/
