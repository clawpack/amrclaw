#include "clawpack_MemoryManager.H"

const unsigned int MemoryManager::align_size;

MemoryManager::~MemoryManager () {}

std::size_t
MemoryManager::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}
