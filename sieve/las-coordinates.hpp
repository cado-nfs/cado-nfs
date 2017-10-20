#ifndef LAS_COORDINATES_HPP_
#define LAS_COORDINATES_HPP_

#include <stdint.h>
#include "las-types.hpp"

/*  Forward declarations of conversion functions */
void NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, sieve_info const & si);
void adjustIJsublat(int *i, unsigned int *j, sieve_info const & si);

void IJTox(uint64_t * x, int i, unsigned int j, sieve_info const & si);
void IJToNx(unsigned int *N, unsigned int * x, int i, unsigned int j, sieve_info const & si);
void IJToAB(int64_t *a, uint64_t *b, const int i, const unsigned int j, sieve_info const & si);
static inline void xToAB(int64_t *a, uint64_t *b, const uint64_t x, sieve_info const & si);
static inline void NxToAB(int64_t *a, uint64_t *b, const unsigned int N, const unsigned int x, sieve_info const & si);
int ABToIJ(int *i, unsigned int *j, const int64_t a, const uint64_t b, sieve_info const & si);
int ABTox(uint64_t *x, const int64_t a, const uint64_t b, sieve_info const & si);
int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, sieve_info const & si);
/*  */

/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void xToAB(int64_t *a, uint64_t *b, const uint64_t x, sieve_info const & si)
{
    int i, j;
    int64_t c;
    uint32_t I = si.I;

    i = (x & (I - 1)) - (I >> 1);
    j = x >> si.conf.logI_adjusted;
    if (si.conf.sublat.m != 0) {
        i = i*si.conf.sublat.m + si.conf.sublat.i0;
        j = j*si.conf.sublat.m + si.conf.sublat.j0;
    }
    *a = (int64_t) i * si.qbasis.a0 + (int64_t) j * si.qbasis.a1;
    c =  (int64_t) i * si.qbasis.b0 + (int64_t) j * si.qbasis.b1;
    if (c >= 0)
      *b = c;
    else
      {
        *a = -*a;
        *b = -c;
      }
}
static inline void NxToAB(int64_t *a, uint64_t *b, const unsigned int N, const unsigned int x, sieve_info const & si)
{
    xToAB(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, si);
}

#ifdef SUPPORT_LARGE_Q
/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void xToABmpz(mpz_ptr a, mpz_ptr b,
        const uint64_t x, sieve_info const & si)
{
    int i, j;
    uint32_t I = si.I;

    i = (x & (I - 1)) - (I >> 1);
    j = x >> si.conf.logI_adjusted;
    if (si.conf.sublat.m != 0) {
        i = i*si.conf.sublat.m + si.conf.sublat.i0;
        j = j*si.conf.sublat.m + si.conf.sublat.j0;
    }

    mpz_t aux_i, aux_j;
    mpz_t aux;
    mpz_init(aux);
    mpz_init_set_si(aux_i, i);
    mpz_init_set_ui(aux_j, j);
    
    mpz_mul_si(a, aux_i, si.qbasis.a0);
    mpz_mul_si(aux, aux_j, si.qbasis.a1);
    mpz_add(a, a, aux);

    mpz_mul_si(b, aux_i, si.qbasis.b0);
    mpz_mul_si(aux, aux_j, si.qbasis.b1);
    mpz_add(b, b, aux);

    if (mpz_sgn(b) < 0) {
        mpz_neg(a, a);
        mpz_neg(b, b);
    }
    mpz_clear(aux);
    mpz_clear(aux_i);
    mpz_clear(aux_j);
}

static inline void NxToABmpz(mpz_ptr a, mpz_ptr b,
        const unsigned int N, const unsigned int x, sieve_info const & si)
{
    xToABmpz(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, si);
}
#endif

#endif	/* LAS_COORDINATES_HPP_ */
