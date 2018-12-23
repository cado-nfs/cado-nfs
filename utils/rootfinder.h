#ifndef CADO_UTILS_ROOTFINDER_H_
#define CADO_UTILS_ROOTFINDER_H_

#include <gmp.h>
#include <stdint.h>
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif



/* This is the entry point for the root finding routines.
 *
 * It relies on either the mod_ul inline assembly layer
 */

/* Structure to hold factorization of an unsigned long, and to iterate through
   its proper divisors. An integer < 2^64 has at most 15 prime factors. The 
   smallest integer with 16 prime factors is 53# =~ 2^64.8. */
typedef struct {
  unsigned long p[15];
  unsigned char e[15];
  unsigned char c[15];
  unsigned int ndiv;
} enumeratediv_t;


unsigned long mpz_poly_roots_gen(mpz_t **r, mpz_poly_srcptr F, mpz_srcptr p);
int mpz_poly_roots(mpz_t * r, mpz_poly_srcptr F, mpz_srcptr p);
int mpz_poly_roots_ulong(unsigned long * r, mpz_poly_srcptr F, unsigned long p);
int mpz_poly_roots_uint64(uint64_t * r, mpz_poly_srcptr F, uint64_t p);
int mpz_poly_roots_mpz (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p);

extern int roots_mod_uint64 (uint64_t * r, uint64_t a, int d, uint64_t p);

unsigned char factor_ul (unsigned long *, unsigned char *, unsigned long n);
void enumeratediv_init (enumeratediv_t *, unsigned long n);
unsigned long enumeratediv (enumeratediv_t *);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* Some of the prototypes are available only from C++ */
#include <vector>
#include "cxx_mpz.hpp"

/* instatiations are defined in rootfinder.cpp */
template<typename T>
std::vector<T> mpz_poly_roots(cxx_mpz_poly const & f, T const & q);

/* When the factorization is known, compute the result via crt. The
 * factors need not be of the same type as q itself.
 */
template<typename T, typename F = T>
std::vector<T> mpz_poly_roots(cxx_mpz_poly const & f, T const & q, std::vector<F> const & qfac);

#endif  /* __cplusplus */

#endif	/* CADO_UTILS_ROOTFINDER_H_ */
