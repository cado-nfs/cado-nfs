#ifndef BBLAS_MAT64_HPP_
#define BBLAS_MAT64_HPP_

#include "bblas_bitmat.hpp"

typedef bitmat<uint64_t> mat64; // ATTRIBUTE((aligned(64)));

#include "level3a.hpp"
#include "level3b.hpp"

namespace bblas_bitmat_details {

    template<> struct bblas_bitmat_type_supported<uint64_t> {
        static constexpr const bool value = true;
        static constexpr const int width = 64;
        static constexpr const int alignment = 64;  // in bytes
    };

    template<> struct bitmat_ops<mat64> {
        static void fill_random(mat64 &, gmp_randstate_t rstate);
        static void add(mat64 & C, mat64 const & A, mat64 const & B);
        static void transpose(mat64 & C, mat64 const & A);
        static void mul(mat64 & C, mat64 const & A, mat64 const & B);
        static void addmul(mat64 & C, mat64 const & A, mat64 const & B);
        static void addmul(mat64 & C,
                   mat64 const & A,
                   mat64 const & B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1);
        static void trsm(mat64 const & L,
                mat64 & U,
                unsigned int yi0,
                unsigned int yi1);
        static void trsm(mat64 const & L,
                mat64 & U);
    };
}

#endif	/* BBLAS_MAT64_HPP_ */
