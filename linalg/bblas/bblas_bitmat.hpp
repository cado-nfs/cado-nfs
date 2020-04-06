#ifndef BBLAS_BITMAT_HPP_
#define BBLAS_BITMAT_HPP_

/* very generic code on bit matrices. We don't keep much here. More
 * concrete types such as mat64 provide most of the actual code.
 */
#include <cstdint>
#include <cstring>
#include <vector>
#include "macros.h"
#include "memory.h"     // malloc_aligned in utils

namespace bblas_bitmat_details {

    template<typename T> struct bblas_bitmat_type_supported {
        static constexpr const bool value = false;
    };
    template<typename D> struct bitmat_ops;
}


template<typename T>
class bitmat : public bblas_bitmat_details::bitmat_ops<bitmat<T>>
{
    typedef bblas_bitmat_details::bblas_bitmat_type_supported<T> S;
    static_assert(S::value, "bblas bitmap must be built on uintX_t");

    public:
    static constexpr const int width = S::width;
    typedef T datatype;
    typedef std::vector<bitmat, aligned_allocator<bitmat, S::alignment>> vector_type;
    // typedef std::vector<bitmat> vector_type;

    private:
    T x[width]; // ATTRIBUTE((aligned(64))) ;

    public:
    static inline bitmat * alloc(size_t n) {
        return (bitmat *) malloc_aligned(n * sizeof(bitmat), S::alignment);
    }
    static inline void free(bitmat * p) {
        free_aligned(p);
    }

    inline T* data() { return x; }
    inline const T* data() const { return x; }
    T& operator[](int i) { return x[i]; }
    T operator[](int i) const { return x[i]; }
    inline bool operator==(bitmat const& a) const
    {
        return memcmp(x, a.x, sizeof(x)) == 0;
    }
    inline bool operator!=(bitmat const& a) const { return !operator==(a); }
    bitmat() {}
    inline bitmat(bitmat const& a) { memcpy(x, a.x, sizeof(x)); }
    inline bitmat& operator=(bitmat const& a)
    {
        memcpy(x, a.x, sizeof(x));
        return *this;
    }
    inline bitmat(int a) { *this = a; }
    inline bitmat& operator=(int a)
    {
        if (a & 1) {
            T mask = 1;
            for (int j = 0; j < width; j++, mask <<= 1)
                x[j] = mask;
        } else {
            memset(x, 0, sizeof(x));
        }
        return *this;
    }
};

#endif	/* BBLAS_BITMAT_HPP_ */
