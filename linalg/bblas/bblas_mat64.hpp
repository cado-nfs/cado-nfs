#ifndef BBLAS_MAT64_HPP_
#define BBLAS_MAT64_HPP_

#include "bblas_bitmat.hpp"

typedef bitmat<uint64_t> mat64; // ATTRIBUTE((aligned(64)));

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
        static bool is_lowertriangular(mat64 const & A);
        static bool is_uppertriangular(mat64 const & A);
        static bool triangular_is_unit(mat64 const & A);
        static void extract_uppertriangular(mat64 & a, mat64 const & b);
        static void extract_lowertriangular(mat64 & a, mat64 const & b);
        static void make_uppertriangular(mat64 & u);
        static void make_lowertriangular(mat64 & u);
        static void make_unit_uppertriangular(mat64 & u);
        static void make_unit_lowertriangular(mat64 & u);
        static void triangular_make_unit(mat64 & u);
    };
}

#include "bblas_level3a.hpp"
#include "bblas_level3b.hpp"

#endif	/* BBLAS_MAT64_HPP_ */
