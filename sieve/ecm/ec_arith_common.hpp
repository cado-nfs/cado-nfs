#ifndef CADO_SIEVE_ECM_EC_ARITH_COMMON_HPP
#define CADO_SIEVE_ECM_EC_ARITH_COMMON_HPP

#include <cstdio>
#include <cstdlib>

#include <utility>

#include "cxx_mpz.hpp"
#include "fmt/base.h"

/* Types of coordinates */
enum ec_point_coord_type_t {
    SHORT_WEIERSTRASS_aff,
    SHORT_WEIERSTRASS_proj,
    MONTGOMERY_xz,
    TWISTED_EDWARDS_proj,
    TWISTED_EDWARDS_ext,
};

/* A point on an elliptic curve, over an arithxx layer.
 *   The significant coordinates are:
 *      - x,y     for SHORT_WEIERSTRASS_aff
 *      - x,y,z   for SHORT_WEIERSTRASS_proj
 *      - x,  z   for MONTGOMERY_xz
 *      - x,y,z   for TWISTED_EDWARDS_proj
 *      - x,y,z,t for TWISTED_EDWARDS_ext
 */
template<typename layer>
struct ec_point {
    typename layer::Residue x, y, z, t;
    explicit ec_point(typename layer::Modulus const & m)
        : x(m), y(m), z(m), t(m)
    {}
};

template<typename layer>
void ec_point_set(ec_point<layer> & Q, ec_point<layer> const & P,
        typename layer::Modulus const & m,
        ec_point_coord_type_t const coord)
{
    m.set(Q.x, P.x);
    if (coord != MONTGOMERY_xz)
        m.set(Q.y, P.y);
    if (coord != SHORT_WEIERSTRASS_aff)
        m.set(Q.z, P.z);
    if (coord == TWISTED_EDWARDS_ext)
        m.set(Q.t, P.t);
}

template<typename layer>
void ec_point_swap(ec_point<layer> & Q, ec_point<layer> & P,
        ec_point_coord_type_t const coord)
{
    std::swap(Q.x, P.x);
    if (coord != MONTGOMERY_xz)
        std::swap(Q.y, P.y);
    if (coord != SHORT_WEIERSTRASS_aff)
        std::swap(Q.z, P.z);
    if (coord == TWISTED_EDWARDS_ext)
        std::swap(Q.t, P.t);
}

/* The value of the residue r, as an integer (for printing) */
template<typename layer>
cxx_mpz ec_residue_get_mpz(typename layer::Residue const & r,
        typename layer::Modulus const & m)
{
    return cxx_mpz(m.get(r));
}

template<typename layer>
void ec_point_fprintf(FILE * out, ec_point<layer> const & P,
        ec_point_coord_type_t const coord,
        typename layer::Modulus const & m)
{
    cxx_mpz const x = ec_residue_get_mpz<layer>(P.x, m);
    cxx_mpz y, z, t;
    if (coord != MONTGOMERY_xz)
        y = ec_residue_get_mpz<layer>(P.y, m);
    if (coord != SHORT_WEIERSTRASS_aff)
        z = ec_residue_get_mpz<layer>(P.z, m);
    if (coord == TWISTED_EDWARDS_ext)
        t = ec_residue_get_mpz<layer>(P.t, m);

    switch (coord) {
        case SHORT_WEIERSTRASS_aff: /* (x, y) */
            fmt::print(out, "({:#x}, {:#x})", x, y);
            break;
        case TWISTED_EDWARDS_proj: /* (X : Y : Z) */
        case SHORT_WEIERSTRASS_proj: /* (X : Y : Z) */
            fmt::print(out, "({:#x} : {:#x} : {:#x})", x, y, z);
            break;
        case MONTGOMERY_xz: /* (X :: Z) */
            fmt::print(out, "({:#x} :: {:#x})", x, z);
            break;
        case TWISTED_EDWARDS_ext: /* (X : Y : Z : T) */
            fmt::print(out, "({:#x} : {:#x} : {:#x} : {:#x})", x, y, z, t);
            break;
        default:
            fmt::print(stderr, "{}: unknown coordinates system\n", __func__);
            abort();
    }
}

#endif /* CADO_SIEVE_ECM_EC_ARITH_COMMON_HPP */
