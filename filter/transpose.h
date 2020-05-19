#ifndef TRANSPOSE_H_
#define TRANSPOSE_H_

#include "typedefs.h"

// #define TRANSPOSE_EASY_WAY

#ifdef __cplusplus
extern "C" {
#endif


/* copy 64 bytes from src to dst using non-temporal store instructions
   if available (this bypasses the cache). */
static inline void store_nontemp_64B(void * dst, void * src);


/* converts a sparse matrix in COOrdinate format to the CSR format.
   INPUT:  COO sparse matrix in Ai, Aj (both of size nnz), with n rows

   OUTPUT: CSR sparse matrix in Rp, Ri.
   
   Rp and Ri MUST be preallocated (of sizes n+1 and nnz, respectively). 
   The "row pointers" Rp MUST be already computed. 

   Ai, Aj, Ri MUST be aligned on a 64-byte boundary (for good cache behavior).
   The input arrays are expendable (i.e. they might be destroyed). 
   The current code only reads them though. */
void transpose(uint64_t nnz, index_t *Ai, index_t *Aj, index_t n, index_t *Rp, index_t *Ri);

#ifdef __cplusplus
}
#endif

/* L1 cache line has size 64 on most CPUs */
#define CACHELINE_SIZE ((int) (64 / sizeof(index_t)))

#if __AVX__
#include <immintrin.h>
static inline void store_nontemp_64B(void * dst, void * src)
{
  register __m256i * d1 = (__m256i*) dst;
  register __m256i s1 = *((__m256i*) src);
  register __m256i * d2 = d1+1;
  register __m256i s2 = *(((__m256i*) src)+1);
  _mm256_stream_si256(d1, s1);
  _mm256_stream_si256(d2, s2);
  /* note : it can also be done using SSE for non-AVX machines */
}
#else
static inline void store_nontemp_64B(void * dst, void * src)
{
  index_t *in = src;
  index_t *out = dst;
  for (int i = 0; i < CACHELINE_SIZE; i++)
    out[i] = in[i];
}
#endif

#endif	/* TRANSPOSE_H_ */
