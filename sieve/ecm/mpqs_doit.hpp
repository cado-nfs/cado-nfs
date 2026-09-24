#ifndef CADO_SIEVE_ECM_MPQS_DOIT_HPP
#define CADO_SIEVE_ECM_MPQS_DOIT_HPP

#include <gmp.h>

/* Put in f a factor of N0 using MPQS. N0 must be odd. */
void mpqs_doit (mpz_ptr f, mpz_srcptr N0, int verbose);

/* Statistics on the norms tested so far; a negative argument prints them. */
void smooth_stat (int);

#endif	/* CADO_SIEVE_ECM_MPQS_DOIT_HPP */
