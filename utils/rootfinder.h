#ifndef CADO_UTILS_ROOTFINDER_H_
#define CADO_UTILS_ROOTFINDER_H_

#include <gmp.h>
#include <stdint.h>
#include "mpz_poly.h"
#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif



/* This is the entry point for the root finding routines.
 *
 * It relies on either the mod_ul inline assembly layer
 */


unsigned long mpz_poly_roots_gen(mpz_t **r, mpz_poly_srcptr F, mpz_srcptr p, gmp_randstate_ptr rstate);
int mpz_poly_roots(mpz_t * r, mpz_poly_srcptr F, mpz_srcptr p, gmp_randstate_ptr rstate);
int mpz_poly_roots_ulong(unsigned long * r, mpz_poly_srcptr F, unsigned long p, gmp_randstate_ptr rstate);
int mpz_poly_roots_uint64(uint64_t * r, mpz_poly_srcptr F, uint64_t p, gmp_randstate_ptr rstate);
int mpz_poly_roots_mpz (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_ptr rstate);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* Some of the prototypes are available only from C++ */
#include <vector>

/* instatiations are defined in rootfinder.cpp */
template<typename T>
std::vector<T> mpz_poly_roots(cxx_mpz_poly const & f, T const & q, gmp_randstate_ptr rstate);

/* When the factorization is known, compute the result via crt. The
 * factors need not be of the same type as q itself.
 */
template<typename T, typename F = T>
std::vector<T> mpz_poly_roots(cxx_mpz_poly const & f, T const & q, std::vector<F> const & qfac, gmp_randstate_ptr rstate);

#endif  /* __cplusplus */

#endif	/* CADO_UTILS_ROOTFINDER_H_ */
