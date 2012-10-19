#ifndef LAS_COORDINATES_H_
#define LAS_COORDINATES_H_

#include <stdint.h>
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Forward declarations of conversion functions */
void xToIJ(int *i, unsigned int *j, const unsigned int X, sieve_info_srcptr si);
void NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, sieve_info_srcptr si);

void IJTox(unsigned int * x, int i, unsigned int j, sieve_info_srcptr si);
void IJToNx(unsigned int *N, unsigned int * x, int i, unsigned int j, sieve_info_srcptr si);
void IJToAB(int64_t *a, uint64_t *b, const int i, const unsigned int j, sieve_info_srcptr si);
static inline void xToAB(int64_t *a, uint64_t *b, const unsigned int x, sieve_info_srcptr si);
static inline void NxToAB(int64_t *a, uint64_t *b, const unsigned int N, const unsigned int x, sieve_info_srcptr si);
int ABToIJ(int *i, unsigned int *j, const int64_t a, const uint64_t b, sieve_info_srcptr si);
int ABTox(unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si);
int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si);
/*  */

/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void xToAB(int64_t *a, uint64_t *b, const unsigned int x, sieve_info_srcptr si)
{
    int i, j;
    int64_t c;
    uint32_t I = si->I;

    i = (x & (I - 1)) - (I >> 1);
    j = x >> si->logI;
    *a = (int64_t) i * si->a0 + (int64_t) j * si->a1;
    c =  (int64_t) i * si->b0 + (int64_t) j * si->b1;
    if (c >= 0)
      *b = c;
    else
      {
        *a = -*a;
        *b = -c;
      }
}
static inline void NxToAB(int64_t *a, uint64_t *b, const unsigned int N, const unsigned int x, sieve_info_srcptr si)
{
    xToAB(a, b, (N << LOG_BUCKET_REGION) + x, si);
}

#ifdef __cplusplus
}
#endif

#endif	/* LAS_COORDINATES_H_ */
