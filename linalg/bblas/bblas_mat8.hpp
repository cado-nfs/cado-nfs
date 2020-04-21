#ifndef BBLAS_MAT8_HPP_
#define BBLAS_MAT8_HPP_

#include "bblas_bitmat.hpp"

/* we provide a matrix type for 8*8 matrices. This should not be
 * understood as something pretending any sort of efficiency, but rather
 * as a means to better test some corner cases that might be hard to see
 * with the 64*64 case
 *
 * This is also an illustration of the fact that our code does not do
 * wild assumptions on the bit size.
 */
typedef bitmat<uint8_t> mat8;

namespace bblas_bitmat_details {

    template<> struct bblas_bitmat_type_supported<uint8_t> {
        static constexpr const bool value = true;
        static constexpr const int width = 8;
        static constexpr const int alignment = 64;  // in bytes
    };

    /* warn the compiler that we have some specializations */
    template<> void bitmat_ops<mat8>::add(mat8 & C, mat8 const & A, mat8 const & B);
    template<> void bitmat_ops<mat8>::transpose(mat8 & C, mat8 const & A);
    template<> void bitmat_ops<mat8>::mul(mat8 & C, mat8 const & A, mat8 const & B);
    template<> void bitmat_ops<mat8>::addmul_blocks(mat8 * C, mat8 const * A, mat8 const& B, size_t nblocks, size_t block_stride);
    template<> void bitmat_ops<mat8>::mul_blocks(mat8 * C, mat8 const * A, mat8 const& B, size_t nblocks, size_t block_stride);
    template<> void bitmat_ops<mat8>::addmul(mat8 & C,
            mat8 const & A,
            mat8 const & B,
            unsigned int i0,
            unsigned int i1,
            unsigned int yi0,
            unsigned int yi1);
    template<> void bitmat_ops<mat8>::trsm(mat8 const & L,
            mat8 & U,
            unsigned int yi0,
            unsigned int yi1);

    extern template struct bitmat_ops<mat8>;
}

#endif	/* BBLAS_MAT8_HPP_ */
