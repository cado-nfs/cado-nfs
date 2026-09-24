#ifndef CADO_SIEVE_ECM_ECM_HPP
#define CADO_SIEVE_ECM_ECM_HPP

#include <cstdint>

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"
#include "facul_ecm.h"

/* ECM with the given plan. The arithmetic is done with the arithxx layer
 * that suits the size of the modulus. The factor found, or 1, is written
 * to f. Returns 1 if backtracking was used, 0 otherwise.
 */
int ecm(modintredcul_t f, const modulusredcul_t m, const ecm_plan_t * plan);
int ecm(modintredc15ul_t f, const modulusredc15ul_t m, const ecm_plan_t * plan);
int ecm(modintredc2ul2_t f, const modulusredc2ul2_t m, const ecm_plan_t * plan);
int ecm(modintmpz_t f, const modulusmpz_t m, const ecm_plan_t * plan);

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
