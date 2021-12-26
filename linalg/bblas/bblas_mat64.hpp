#ifndef BBLAS_MAT64_HPP_
#define BBLAS_MAT64_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"
// IWYU pragma: friend "bpack.hpp"

#include <stddef.h>          // for size_t
#include <stdint.h>          // for uint64_t
#include "macros.h"          // for ATTRIBUTE
#include "bblas_bitmat.hpp"

/* The attribute seems to be really necessary, and actually abided by by
 * clang at least for stack-defined variables. On freebsd, this matters,
 * too.
 *
 * 20211226: added explicit alignas in the definition of bitmat. I
 * _think_ that this is enough to make the extra alignment constraint
 * here unnecessary.
 */
typedef bitmat<uint64_t> mat64 ATTRIBUTE((aligned(64)));

namespace bblas_bitmat_details {

    template<> struct bblas_bitmat_type_supported<uint64_t> {
        static constexpr const bool value = true;
        static constexpr const int width = 64;
        static constexpr const int alignment = 64;  // in bytes
    };

    /* warn the compiler that we have some specializations */
    template<> void bitmat_ops<uint64_t>::add(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<uint64_t>::transpose(mat64 & C, mat64 const & A);
    template<> void bitmat_ops<uint64_t>::mul(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<uint64_t>::addmul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t Cstride, size_t Astride);
    template<> void bitmat_ops<uint64_t>::mul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t Cstride, size_t Astride);
    template<> void bitmat_ops<uint64_t>::mul_lt_ge(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<uint64_t>::addmul(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<uint64_t>::addmul(mat64 & C,
            mat64 const & A,
            mat64 const & B,
            unsigned int i0,
            unsigned int i1,
            unsigned int yi0,
            unsigned int yi1);
    template<> void bitmat_ops<uint64_t>::trsm(mat64 const & L,
            mat64 & U,
            unsigned int yi0,
            unsigned int yi1);

    extern template struct bitmat_ops<uint64_t>;
}

#endif	/* BBLAS_MAT64_HPP_ */
