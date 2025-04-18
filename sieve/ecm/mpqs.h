#ifndef CADO_MPQS_H
#define CADO_MPQS_H

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"

#ifdef __cplusplus
extern "C" {
#endif

int mpqs_ul (modint_t, const modulus_t);
int mpqs_15ul (modint_t, const modulus_t);
int mpqs_2ul2 (modint_t, const modulus_t);
int mpqs_mpz (modint_t, const modulus_t);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MPQS_H */

