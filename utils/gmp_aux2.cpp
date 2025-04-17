#include "cado.h" // IWYU pragma: keep

#include "gmp_aux.h"
#include "cado_math_aux.hpp"

/* This compilation unit is a complement to gmp_aux.c ; some of the
 * utilities there are actually implemented by c++ code
 */
long double
mpz_get_ld (mpz_srcptr z)
{
#if 0
    long double ld;
    double d;
    mpz_t t;

    d = mpz_get_d (z);
    mpz_init (t);
    mpz_set_d (t, d);
    mpz_sub (t, z, t);
    ld = (long double) d + (long double) mpz_get_d (t);
    mpz_clear (t);
    return ld;
#endif
    return cado_math_aux::mpz_get<long double>(z);
}

