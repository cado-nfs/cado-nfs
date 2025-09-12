#include "cado.h" // IWYU pragma: keep

#include <gmp.h>

#include "gmp_aux.h"
#include "cado_mp_conversions.hpp"

/* This compilation unit is a complement to gmp_aux.c ; some of the
 * utilities there are actually implemented by c++ code
 */
long double
mpz_get_ld (mpz_srcptr z)
{
    return cado_math_aux::mpz_get<long double>(z);
}

