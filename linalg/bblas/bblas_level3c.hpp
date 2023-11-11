#ifndef BBLAS_LEVEL3C_HPP_
#define BBLAS_LEVEL3C_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "cado_config.h"    // for HAVE_SSE2, ULONG_BITS, HAVE_AVX2
#include <cstddef>         // for size_t
#include <cstdint>         // for uint64_t, uint32_t
#include "bblas_mat64.hpp"

/**********************************************************************/
/* level 3c: operations on matrices with one arbitrary length (sometimes 2)
 *      copy_N64
 *      cmp_N64
 *      mul_N64_6464   (several variants)
 *      mul_N64_T6464   (several variants)
 *      addmul_N64_6464   (several variants)
 *      mul_TN64_N64
 *      addmul_TN64_N64
 * with varying left length:
 *      mul_TN32_N64
 *      mul_TN64K_N64
 *      addmul_TN64K_N64
 */

/* implemented here:
 *    - mul_N64_6464: Given an N*64 matrix A (N uint64_t's) and a 64*64
 *      matrix B, compute the product (in place).  With respect to
 *      endianness, we match the column of (A[i]&1)'s with B[0].
 *    - mul_N64_T6464: same, but multiply by the transpose of the second
 *      operand.
 *    - mul_TN64_N64: rank-N update. More wordy: this takes, in row major
 *      order, an Nx64 matrix A (transpose of a 64xN matrix), together
 *      with another Nx64 matrix B, and xors the output matrix with the
 *      product transpose(A)*B -- this may as well be seen as the block
 *      dot product of A and B.
 *    - mul_TN32_N64: rank-N update, but creates a matrix of size 32*N
 *    - mul_TN64K_N64: rank-N update, but creates a matrix of size 64K*N
 */

void copy_N64(uint64_t * dst, uint64_t const * src, size_t m);

int cmp_N64(uint64_t const * dst, uint64_t const * src, size_t m);

/* implementation details, variants */
//
///////////////////////////////////////////////////////////////////////
void mul_N64_6464_lookup4(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);
void mul_6464lt_6464_lookup4(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B);
void mul_N64_6464_lookup8(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);
void mul_N64_6464_vec(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);
void mul_N64_6464_transB(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);

#if defined(HAVE_AVX2)
void mul_N64_6464_avx2(uint64_t *C,
		 uint64_t const *A,
		 mat64 const & B, size_t m);
#endif

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void mul_N64_6464_sse(uint64_t *C,
		 uint64_t const *A,
		 mat64 const & B, size_t m);
#endif

///////////////////////////////////////////////////////////////////////

void mul_N64_T6464_vec(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);
void mul_N64_T6464_transB(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);

///////////////////////////////////////////////////////////////////////

void addmul_N64_6464_lookup4(uint64_t *C, 
                   uint64_t const *A,
                   mat64 const & B, size_t m);
#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_N64_6464_sse(uint64_t *C,
		 uint64_t const *A,
		 mat64 const & B, size_t m);
#endif

///////////////////////////////////////////////////////////////////////

void mul_TN64_N64_addmul(mat64 & r, uint64_t const * a, uint64_t const * w, size_t n);
void mul_TN64_N64_C(mat64 &  b, uint64_t const * A, uint64_t const * x, unsigned int ncol);

///////////////////////////////////////////////////////////////////////

void addmul_TN64_N64_C(mat64 &  b, uint64_t const * A, uint64_t const * x, unsigned int ncol);

///////////////////////////////////////////////////////////////////////

void mul_TN32_N64_C(uint64_t * b, uint32_t const * A, uint64_t const * x, unsigned int ncol);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_TN64K_N64_sse2(mat64 * w, uint64_t const * u, uint64_t const * v, unsigned int n, unsigned int K);
#endif

void addmul_TN64K_N64_C(mat64 * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K);

///////////////////////////////////////////////////////////////////////

/* final exported choices */
void mul_N64_6464(uint64_t *C,
		 uint64_t const *A,
		 mat64 const & B, size_t m);
void addmul_6464_blocks(mat64 *C,
                   mat64 const *A,
                   mat64 const & B, size_t nblocks, size_t Cstride, size_t Astride);
void mul_6464_blocks(mat64 *C,
                   mat64 const *A,
                   mat64 const & B, size_t nblocks, size_t Cstride, size_t Astride);
void addmul_N64_6464(uint64_t *C,
		 uint64_t const *A,
		 mat64 const & B, size_t m);
void mul_N64_T6464(uint64_t *C,
                   uint64_t const *A,
                   mat64 const & B, size_t m);
void mul_TN64_N64(mat64 &  b, uint64_t const * A, uint64_t const * x, unsigned int ncol);
void addmul_TN64_N64(mat64 &  b, uint64_t const * A, uint64_t const * x, unsigned int ncol);
void mul_TN32_N64(uint64_t * b, uint32_t const * A, uint64_t const * x, unsigned int ncol);
void mul_TN64K_N64(uint64_t * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K);
void addmul_TN64K_N64(uint64_t * b, uint64_t const * A, uint64_t const * x, unsigned int ncol, unsigned int K);

#endif	/* BBLAS_LEVEL3C_HPP_ */
