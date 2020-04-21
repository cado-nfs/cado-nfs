#include "cado.h"

#include "bblas.hpp"
#include "bblas_bitmat_inl.hpp"
#include "bblas_level3b.hpp"
#include "bblas_level3d.hpp"
#include "bblas_mat64.hpp"

using namespace bblas_bitmat_details;

template<>
void
bitmat_ops<mat64>::add(mat64& C, mat64 const& A, mat64 const& B)
{
    mat64_add(C, A, B);
}

template<>
void
bitmat_ops<mat64>::transpose(mat64& C, mat64 const& A)
{
    mat64_transpose(C, A);
}

template<>
void
bitmat_ops<mat64>::trsm(mat64 const& L,
                        mat64& U,
                        unsigned int yi0,
                        unsigned int yi1)
{
    trsm64_general(L, U, yi0, yi1);
}

template<>
void
bitmat_ops<mat64>::addmul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t block_stride)
{
    addmul_6464_blocks(C, A, B, nblocks, block_stride);
}

template<>
void
bitmat_ops<mat64>::mul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t block_stride)
{
    mul_6464_blocks(C, A, B, nblocks, block_stride);
}

template<>
void
bitmat_ops<mat64>::mul(mat64& C, mat64 const& A, mat64 const& B)
{
    mul_6464_6464(C, A, B);
}

template<> void bitmat_ops<mat64>::mul_lt_ge(mat64 & C, mat64 const & A, mat64 const & B)
{
    mul_6464lt_6464(C, A, B);
}

template<>
void
bitmat_ops<mat64>::addmul(mat64& C, mat64 const& A, mat64 const& B)
{
    addmul_6464_6464(C, A, B);
}

template<>
void
bitmat_ops<mat64>::addmul(mat64& C,
                          mat64 const& A,
                          mat64 const& B,
                          unsigned int i0,
                          unsigned int i1,
                          unsigned int yi0,
                          unsigned int yi1)
{
    addmul_6464_6464_fragment_lookup4(C, A, B, i0, i1, yi0, yi1);
}

template struct bblas_bitmat_details::bitmat_ops<mat64>;
