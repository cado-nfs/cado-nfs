#ifndef LAS_COORDINATES_HPP_
#define LAS_COORDINATES_HPP_

#include <stdint.h>
#include "las-types.hpp"

/*  Forward declarations of conversion functions */
void xToIJ(int *i, unsigned int *j, const uint64_t X, sieve_info const & si);
void NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, sieve_info const & si);

void IJTox(uint64_t * x, int i, unsigned int j, sieve_info const & si);
void IJToNx(unsigned int *N, unsigned int * x, int i, unsigned int j, sieve_info const & si);
void IJToAB(int64_t *a, int64_t *b, const int i, const unsigned int j, sieve_info const & si);
static inline void xToAB(int64_t *a, int64_t *b, const uint64_t x, sieve_info const & si);
static inline void NxToAB(int64_t *a, int64_t *b, const unsigned int N, const unsigned int x, sieve_info const & si);
int ABToIJ(int *i, unsigned int *j, const int64_t a, const int64_t b, sieve_info const & si);
int ABTox(uint64_t *x, const int64_t a, const int64_t b, sieve_info const & si);
int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const int64_t b, sieve_info const & si);
/*  */

/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void xToAB(int64_t *a, int64_t *b, const uint64_t x, sieve_info const & si)
{
    uint32_t I = si.I;
    int i = (x & (I - 1)) - (I >> 1);
    int j = x >> si.conf.logI_adjusted;
    int64_t xa = (int64_t) i * si.qbasis.a0 + (int64_t) j * si.qbasis.a1;
    int64_t xb = (int64_t) i * si.qbasis.b0 + (int64_t) j * si.qbasis.b1;

    /* normalization choice: b>=0 */
    if (xb >= 0) {
        *a = xa;
        *b = xb;
    } else {
        *a = -xa;
        *b = -xb;
    }
}
static inline void NxToAB(int64_t *a, int64_t *b, const unsigned int N, const unsigned int x, sieve_info const & si)
{
    xToAB(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, si);
}

#ifdef SUPPORT_LARGE_Q
/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void xToABmpz(mpz_t a, mpz_t b,
        const uint64_t x, sieve_info const & si)
{
    uint32_t I = si.I;

    int i = (x & (I - 1)) - (I >> 1);
    int j = x >> si.conf.logI_adjusted;
    
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

    /* normalization choice: b>=0 */
    if (mpz_sgn(b) < 0) {
        mpz_neg(a, a);
        mpz_neg(b, b);
    }
    mpz_clear(aux);
    mpz_clear(aux_i);
    mpz_clear(aux_j);
}

static inline void NxToABmpz(mpz_t a, mpz_t b,
        const unsigned int N, const unsigned int x, sieve_info const & si)
{
    xToABmpz(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, si);
}
#endif

#endif	/* LAS_COORDINATES_HPP_ */
