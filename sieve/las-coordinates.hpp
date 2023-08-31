#ifndef LAS_COORDINATES_HPP_
#define LAS_COORDINATES_HPP_

#include <gmp.h>             // for mpz_mul_si, mpz_ptr, mpz_clear, mpz_t
#include <cstdint>       // for uint64_t, int64_t
#include "fb-types.h"    // for sublat_t
#include "las-config.h"  // for LOG_BUCKET_REGION
#include "las-qlattice.hpp"

/* See las-coordinates.cpp for documentation on the various coordinate
 * systems */

/*  Forward declarations of conversion functions */
void convert_X_to_ij(int & i, unsigned int & j, const uint64_t x, int logI);
void convert_Nx_to_ij(int & i, unsigned int & j, const unsigned int N, const unsigned int x, int logI);
void adjustIJsublat(int & i, unsigned int & j, sublat_t const & S);

void convert_ij_to_X(uint64_t & x, int i, unsigned int j, int logI);
void convert_ij_to_Nx(unsigned int & N, unsigned int &  x, int i, unsigned int j, int logI);
void convert_ij_to_ab(int64_t & a, uint64_t & b, const int i, const unsigned int j, qlattice_basis const & Q);
static inline void convert_X_to_ab(int64_t & a, uint64_t & b, const uint64_t x, int logI, qlattice_basis const & Q);
static inline void convert_Nx_to_ab(int64_t & a, uint64_t & b, const unsigned int N, const unsigned int x, int logI, qlattice_basis const & Q);
int convert_ab_to_ij(int & i, unsigned int & j, const int64_t a, const uint64_t b, qlattice_basis const & Q);
// code exists in las-coordinates.cpp, but is unused, so untested.
// int convert_ab_to_X(uint64_t *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q);
// int convert_ab_to_Nx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q);

/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void convert_X_to_ab(int64_t & a, uint64_t & b, const uint64_t x, int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    convert_X_to_ij(i, j, x, logI);
    convert_ij_to_ab(a, b, i, j, Q);
}
static inline void convert_Nx_to_ab(int64_t & a, uint64_t & b, const unsigned int N, const unsigned int x, int logI, qlattice_basis const & Q)
{
    convert_X_to_ab(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, logI, Q);
}

// this is only used with SUPPORT_LARGE_Q, obviously, but having them
// doesn't hurt
// #ifdef SUPPORT_LARGE_Q
/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void convert_X_to_abmpz(mpz_ptr a, mpz_ptr b,
        const uint64_t x,
        int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    convert_X_to_ij(i, j, x, logI);
    adjustIJsublat(i, j, Q.sublat);

    mpz_t aux_i, aux_j;
    mpz_t aux;
    mpz_init(aux);
    mpz_init_set_si(aux_i, i);
    mpz_init_set_ui(aux_j, j);
    
    mpz_mul_si(a, aux_i, Q.a0);
    mpz_mul_si(aux, aux_j, Q.a1);
    mpz_add(a, a, aux);

    mpz_mul_si(b, aux_i, Q.b0);
    mpz_mul_si(aux, aux_j, Q.b1);
    mpz_add(b, b, aux);

    if (mpz_sgn(b) < 0) {
        mpz_neg(a, a);
        mpz_neg(b, b);
    }
    mpz_clear(aux);
    mpz_clear(aux_i);
    mpz_clear(aux_j);
}

static inline void convert_Nx_to_abmpz(mpz_ptr a, mpz_ptr b,
        const unsigned int N, const unsigned int x,
        int logI, qlattice_basis const & Q)
{
    convert_X_to_abmpz(a, b, (((uint64_t)N) << LOG_BUCKET_REGION) + (uint64_t)x, logI, Q);
}
// #endif

#endif	/* LAS_COORDINATES_HPP_ */
