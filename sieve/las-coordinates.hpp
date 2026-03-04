#ifndef CADO_LAS_COORDINATES_HPP
#define CADO_LAS_COORDINATES_HPP

#include <cstdint>

#include <gmp.h>

#include "fb-types.hpp"
#include "las-config.hpp"
#include "las-qlattice.hpp"

/* See las-coordinates.cpp for documentation on the various coordinate
 * systems */

/* The functions convert_ab_to_ij and convert_ij_to_ab are special-q dependant
 * and are implemented as methods of the special-q classes in
 * las-qlattice.[ch]pp.
 */

/*  Forward declarations of conversion functions */
void convert_X_to_ij(int & i, unsigned int & j, uint64_t x, int logI);
void convert_Nx_to_ij(int & i, unsigned int & j, unsigned int N, unsigned int x, int logI);

void convert_ij_to_X(uint64_t & x, int i, unsigned int j, int logI);
void convert_ij_to_Nx(unsigned int & N, unsigned int &  x, int i, unsigned int j, int logI);
// code exists in las-coordinates.cpp, but is unused, so untested.
// int convert_ab_to_X(uint64_t *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q);
// int convert_ab_to_Nx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q);

static inline uint64_t convert_Nx_to_X(
        const unsigned int N,
        const unsigned int x)
{
    return (((uint64_t) N) << LOG_BUCKET_REGION) + (uint64_t) x;
}

/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void convert_X_to_ab(
        int64_t & a,
        uint64_t & b,
        const uint64_t x,
        int logI,
        special_q_data_base const & Q)
{
    int i;
    unsigned int j;
    convert_X_to_ij(i, j, x, logI);
    Q.convert_ij_to_ab(a, b, i, j);
}

static inline void convert_Nx_to_ab(
        int64_t & a,
        uint64_t & b,
        const unsigned int N,
        const unsigned int x,
        int logI,
        special_q_data_base const & Q)
{
    convert_X_to_ab(a, b, convert_Nx_to_X(N, x), logI, Q);
}

// this is only used with SUPPORT_LARGE_Q, obviously, but having them
// doesn't hurt
/* Warning: b might be negative, in which case we return (-a,-b) */
static inline void convert_X_to_abmpz(
        cxx_mpz & a,
        cxx_mpz & b,
        const uint64_t x,
        int logI,
        special_q_data_base const & Q)
{
    int i;
    unsigned int j;
    convert_X_to_ij(i, j, x, logI);
    Q.convert_ij_to_ab(a, b, i, j);
}

static inline void convert_Nx_to_abmpz(
        cxx_mpz & a,
        cxx_mpz & b,
        const unsigned int N,
        const unsigned int x,
        int logI,
        special_q_data_base const & Q)
{
    convert_X_to_abmpz(a, b, convert_Nx_to_X(N, x), logI, Q);
}

#endif	/* CADO_LAS_COORDINATES_HPP */
