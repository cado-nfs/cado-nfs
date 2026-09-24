#ifndef CADO_SIEVE_ECM_ECM_HPP
#define CADO_SIEVE_ECM_ECM_HPP

#include <cstdint>

#include "facul_ecm.h"

/* ECM with the given plan, modulo m, with an arithxx layer. The factor
 * found, or 1, is written to f. Returns 1 if backtracking was used, 0
 * otherwise. ecm.cpp instantiates this for modredc64, modredc96,
 * modredc126 and mod_mpz_new.
 */
template<typename layer>
int ecm(typename layer::Integer & f, typename layer::Modulus const & m,
        ecm_plan_t const & plan);

/* The order of the starting point of the curve given by parameterization
 * and parameter, modulo the prime p. If the curve order is known to be ==
 * known_r (mod known_m), this can be supplied. Returns 0 if an inversion
 * failed (which indicates that the curve is not defined modulo p).
 */
uint64_t ec_parameterization_point_order(ec_parameterization_t parameterization,
        unsigned long parameter, uint64_t known_m, uint64_t known_r,
        uint64_t p, int verbose);

/* The order of the curve given by parameterization and parameter, modulo
 * the prime p. Returns 0 if an inversion failed. This takes O(p) time. */
uint64_t ec_parameterization_curve_order(ec_parameterization_t parameterization,
        unsigned long parameter, uint64_t p);

#endif /* CADO_SIEVE_ECM_ECM_HPP */
