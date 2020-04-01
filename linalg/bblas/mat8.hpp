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

#include "level3a.hpp"
#include "level3b.hpp"

namespace bblas_bitmat_details {

    template<> struct bblas_bitmat_type_supported<uint8_t> {
        static constexpr const bool value = true;
        static constexpr const int width = 8;
        static constexpr const int alignment = 64;  // in bytes
    };

    template<> struct bitmat_ops<mat8> {
        static void fill_random(mat8 & w, gmp_randstate_t rstate);
        static void add(mat8 & C, mat8 const & A, mat8 const & B);
        static void transpose(mat8 & C, mat8 const & A);
        static void mul(mat8 & C, mat8 const & A, mat8 const & B);
        static void addmul(mat8 & C, mat8 const & A, mat8 const & B);
        static void addmul(mat8 & C,
                   mat8 const & A,
                   mat8 const & B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1);
        static void trsm(mat8 const & L,
                mat8 & U,
                unsigned int yi0,
                unsigned int yi1);
        static void trsm(mat8 const & L, mat8 & U);
    };
}

#endif	/* BBLAS_MAT8_HPP_ */
