#ifndef MODUL_POLY_H_
#define MODUL_POLY_H_

/* This file provides polynomial arithmetic using the mod_ul layer.
 * mod_ul is rather tightly bound to providing arithmetic on unsigned
 * longs (hence the name), so on first approximation, we're also
 * providing unsigned longs here.
 */

#include <gmp.h>
#include <stdint.h>
#include "mod_ul.h"
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  residueul_t *coeff; /* coefficient list */
} __modul_poly_struct;
typedef __modul_poly_struct modul_poly_t[1];

/* The exported interface contains only the used functions, while in fact a
 * complete set of functions could be made available */

/* The type at the end of the name merely indicates the return type for
 * the stored values */

/* [LI] Do we want these functions here or in rootfinder? */
int modul_poly_roots (residueul_t*, mpz_poly_srcptr, modulusul_t);
int modul_poly_roots_ulong  (unsigned long*, mpz_poly_srcptr, modulusul_t);
int modul_poly_cantor_zassenhaus (residueul_t*, modul_poly_t, modulusul_t);


void modul_poly_init (modul_poly_t, int);
int modul_poly_set_mod (modul_poly_t, mpz_poly_srcptr, modulusul_t);
int modul_poly_set_mod_raw (modul_poly_t, mpz_poly_srcptr, modulusul_t);
void modul_poly_set_immediate (modul_poly_t, int, modulusul_t, ...);
void modul_poly_eval (residueul_t, modul_poly_t, residueul_t, modulusul_t);
void modul_poly_clear (modul_poly_t);
int modul_poly_is_irreducible(modul_poly_t, modulusul_t);
int modul_poly_is_squarefree (modul_poly_t, modulusul_t);


#ifdef __cplusplus
}
#endif

#endif	/* MODUL_POLY_H_ */
