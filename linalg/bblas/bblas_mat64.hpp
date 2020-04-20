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

    /* warn the compiler that we have some specializations */
    template<> void bitmat_ops<mat64>::add(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<mat64>::transpose(mat64 & C, mat64 const & A);
    template<> void bitmat_ops<mat64>::mul(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<mat64>::mul_lt_ge(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<mat64>::addmul(mat64 & C, mat64 const & A, mat64 const & B);
    template<> void bitmat_ops<mat64>::addmul(mat64 & C,
            mat64 const & A,
            mat64 const & B,
            unsigned int i0,
            unsigned int i1,
            unsigned int yi0,
            unsigned int yi1);
    template<> void bitmat_ops<mat64>::trsm(mat64 const & L,
            mat64 & U,
            unsigned int yi0,
            unsigned int yi1);

    extern template struct bitmat_ops<mat64>;
}

#endif	/* BBLAS_MAT64_HPP_ */
