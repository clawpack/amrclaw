#ifndef __ARRAY_H__
#define __ARRAY_H__

template <typename T, const size_t N>
class staticArray { 
    public:
        T x[N]; 

        __host__
        __device__
        T& operator[](const size_t i)
        {
            return x[i];
        }

        __host__
        __device__
        const T& operator[] (const size_t i) const
        {
            return x[i];
        }
};

// specialization 
// staticArray of size 0 for type T
// the [] and = operators are defined in an
// unconventional way such that size 0 
// staticArray can be used the same way as 
// non-zero sized staticArray
template <typename T>
class staticArray<T,0> { 
    public:
        __host__
        __device__
        staticArray<T,0>& operator[](const size_t i)
        {
            return *this;
        }

        __host__
        __device__
        const staticArray<T,0>& operator[](const size_t i) const
        {
            return *this;
        }

        __host__
        __device__
        staticArray<T,0>& operator=(const T& right) {
            // overide the default assignment operator 
            // for 0-sized array. 
            // The = operator does nothing in this case
            return *this;
        }


};

// define type staticArrayReal for storing 
// real type
template <const size_t N>
using staticArrayReal = staticArray<real,N>;



#endif
