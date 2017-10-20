#ifndef MPQS_H_
#define MPQS_H_

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"

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

#endif	/* MPQS_H_ */

