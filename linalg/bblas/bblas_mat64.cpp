#include "cado.h" // IWYU pragma: keep
#include "bblas_bitmat.hpp"   // for bitmat_ops, bblas_bitmat_d...
#include "bblas_bitmat_inl.hpp" // many generic inlines of the bitmat type // IWYU pragma: keep
#include "bblas_level3a.hpp"  // for mat64_add, mat64_transpose
#include "bblas_level3b.hpp"
#include "bblas_level3c.hpp"  // for addmul_6464_blocks, mul_64...
#include "bblas_level3d.hpp"
#include "bblas_mat64.hpp"

using namespace bblas_bitmat_details;

template<>
void
bitmat_ops<uint64_t>::add(mat64& C, mat64 const& A, mat64 const& B)
{
    mat64_add(C, A, B);
}

template<>
void
bitmat_ops<uint64_t>::transpose(mat64& C, mat64 const& A)
{
    mat64_transpose(C, A);
}

template<>
void
bitmat_ops<uint64_t>::trsm(mat64 const& L,
                        mat64& U,
                        unsigned int yi0,
                        unsigned int yi1)
{
    trsm64_general(L, U, yi0, yi1);
}

template<>
void
bitmat_ops<uint64_t>::addmul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t Cstride, size_t Astride)
{
    addmul_6464_blocks(C, A, B, nblocks, Cstride, Astride);
}

template<>
void
bitmat_ops<uint64_t>::mul_blocks(mat64 * C, mat64 const * A, mat64 const& B, size_t nblocks, size_t Cstride, size_t Astride)
{
    mul_6464_blocks(C, A, B, nblocks, Cstride, Astride);
}

template<>
void
bitmat_ops<uint64_t>::mul(mat64& C, mat64 const& A, mat64 const& B)
{
    mul_6464_6464(C, A, B);
}

template<> void bitmat_ops<uint64_t>::mul_lt_ge(mat64 & C, mat64 const & A, mat64 const & B)
{
    mul_6464lt_6464(C, A, B);
}

template<>
void
bitmat_ops<uint64_t>::addmul(mat64& C, mat64 const& A, mat64 const& B)
{
    addmul_6464_6464(C, A, B);
}

template<>
void
bitmat_ops<uint64_t>::addmul(mat64& C,
                          mat64 const& A,
                          mat64 const& B,
                          unsigned int i0,
                          unsigned int i1,
                          unsigned int yi0,
                          unsigned int yi1)
{
    addmul_6464_6464_fragment_lookup4(C, A, B, i0, i1, yi0, yi1);
}

template struct bblas_bitmat_details::bitmat_ops<uint64_t>;
