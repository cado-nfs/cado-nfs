#ifndef CADO_BBLAS_BITMAT_HPP
#define CADO_BBLAS_BITMAT_HPP

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

/* very generic code on bit matrices. We don't keep much here. More
 * concrete types such as mat64 provide most of the actual code.
 */
#include <cstdint>
#include <cstring>
#include <vector>
#include "macros.h"
#include "memory.h"     // malloc_aligned in utils

template<typename T> class bitmat;

namespace bblas_bitmat_details {

    template<typename T> struct bblas_bitmat_type_supported {
        static constexpr const bool value = false;
    };
    template<typename T> struct bitmat_ops {
        typedef bitmat<T> matrix;
        static void fill_random(matrix & w, gmp_randstate_t rstate);
        static void add(matrix & C, matrix const & A, matrix const & B);
        static void transpose(matrix & C, matrix const & A);
        static void mul(matrix & C, matrix const & A, matrix const & B);
        static void mul_lt_ge(matrix & C, matrix const & A, matrix const & B) {
            mul(C, A, B);
        }
        /* do C[0] = A[0]*B, C[block_stride]=A[block_stride]*B, etc */
        static void mul_blocks(matrix * C, matrix const * A, matrix const& B, size_t nblocks, size_t Cstride, size_t Astride);
        static void addmul_blocks(matrix * C, matrix const * A, matrix const& B, size_t nblocks, size_t Cstride, size_t Astride);
        static void addmul(matrix & C, matrix const & A, matrix const & B);
        static void addmul(matrix & C,
                   matrix const & A,
                   matrix const & B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1);
        static void trsm(matrix const & L,
                matrix & U,
                unsigned int yi0,
                unsigned int yi1);
        static void trsm(matrix const & L, matrix & U);
        static void extract_uppertriangular(matrix & a, matrix const & b);
        static void extract_lowertriangular(matrix & a, matrix const & b);
        /* Keeps only the upper triangular part in U, and copy the lower
         * triangular, together with a unit block, to L */
        static void extract_LU(matrix & L, matrix & U);
        protected:
        /* these are accessed as _member functions_ in the matrix type */
        static bool is_lowertriangular(matrix const & a);
        static bool is_uppertriangular(matrix const & a);
        static bool triangular_is_unit(matrix const & a);
        static void make_uppertriangular(matrix & a);
        static void make_lowertriangular(matrix & a);
        static void make_unit_uppertriangular(matrix & a);
        static void make_unit_lowertriangular(matrix & a);
        static void triangular_make_unit(matrix & a);
    };
} /* namespace bblas_bitmat_details */

template<typename T>
class
    alignas(bblas_bitmat_details::bblas_bitmat_type_supported<T>::alignment)
    bitmat
    : public bblas_bitmat_details::bitmat_ops<T>
{
    typedef bblas_bitmat_details::bitmat_ops<T> ops;
    typedef bblas_bitmat_details::bblas_bitmat_type_supported<T> S;
    static_assert(S::value, "bblas bitmap must be built on uintX_t");

    typedef aligned_allocator<bitmat, S::alignment> allocator_type;
    public:
    static constexpr const int alignment = S::alignment;
    static constexpr const int width = S::width;
    typedef T datatype;
    typedef std::vector<bitmat, allocator_type> vector_type;
    // typedef std::vector<bitmat> vector_type;

    private:
    T x[width]; // ATTRIBUTE((aligned(64))) ;

    public:
    static bitmat * alloc(size_t n) {
        bitmat * p = allocator_type().allocate(n);
        return new(p) bitmat[n];
    }
    static void free(bitmat * p, size_t n) {
        ::operator delete (p, p);
        allocator_type().deallocate(p, n);
    }

    T* data() { return x; }
    const T* data() const { return x; }
    T& operator[](int i) { return x[i]; }
    T operator[](int i) const { return x[i]; }
    T& operator[](unsigned int i) { return x[i]; }
    T operator[](unsigned int i) const { return x[i]; }
    bool operator==(bitmat const& a) const
    {
        /* anyway we're not going to do it any smarter in instantiations,
         * so let's rather keep this as a simple and stupid inline */
        return memcmp(x, a.x, sizeof(x)) == 0;
    }
    bitmat() { memset(x, 0, sizeof(x)); }
    bitmat(bitmat const& a) { memcpy(x, a.x, sizeof(x)); }
    bitmat& operator=(bitmat const& a)
    {
        memcpy(x, a.x, sizeof(x));
        return *this;
    }
    bitmat(int a) { *this = a; }
    bitmat& operator=(int a)
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
    bool operator==(int a) const
    {
        if (a&1) {
            T mask = a&1;
            for (int j = 0; j < width; j++, mask <<= 1)
                if (x[j]&~mask) return false;
        } else {
            for (int j = 0; j < width; j++)
                if (x[j]) return false;
        }
        return true;
    }

    bool is_lowertriangular() const { return ops::is_lowertriangular(*this); }
    bool is_uppertriangular() const { return ops::is_uppertriangular(*this); }
    bool triangular_is_unit() const { return ops::triangular_is_unit(*this); }
    void make_uppertriangular() { ops::make_uppertriangular(*this); }
    void make_lowertriangular() { ops::make_lowertriangular(*this); }
    void make_unit_uppertriangular() { ops::make_unit_uppertriangular(*this); }
    void make_unit_lowertriangular() { ops::make_unit_lowertriangular(*this); }
    void triangular_make_unit() { ops::triangular_make_unit(*this); }
};

#endif	/* BBLAS_BITMAT_HPP_ */
