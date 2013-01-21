#ifndef CADO_UTILS_PLAIN_POLY_H_
#define CADO_UTILS_PLAIN_POLY_H_

#include <gmp.h>
#include <stdint.h>

/* These functions provide plain C polynomial arithmetic using a C
 * integral type.
 *
 * The constraint on the modulus (and degree) is that dp^2 must fit
 * within the type.
 */


/* The type used for plain_poly coefficients. Must be a POD (plain old
 * datatype). Note that when asserting whether a given modulus is
 * eligible for use with a plain_poly layer, we consider the modulus as
 * an uint64_t when calling plain_poly_fits. Afterwards, the prime
 * requirement is that the modulus itself fits in a plain_poly_coeff_t
 * (plain_poly_fits is actually more precise than that). */

#include <stdint.h>
typedef int64_t plain_poly_coeff_t;
#define PLAIN_POLY_COEFF_MAX    INT64_MAX

/* For INT64_MAX: it is a pain, but for C++ usage, this macro had better
 * be defined on top of the implementation unit -- and not here, since we
 * can't be sure that stdint.h hasn't been included already (ok,
 * strictly speaking, we don't need PLAIN_POLY_COEFF_MAX really). */
#if defined(__cplusplus) && !defined(__STDC_LIMIT_MACROS)
#error "Please define __STDC_LIMIT_MACROS in the implementation unit"
#endif



typedef struct {
  int alloc;    /* number of allocated coefficients */
  int degree;   /* degree < alloc */
  plain_poly_coeff_t *coeff;  /* coefficient list */
} plain_poly_struct;
typedef plain_poly_struct plain_poly_t[1];
typedef int (*sortfunc_t) (const void *, const void *);

#ifdef __cplusplus
extern "C" {
#endif

extern void plain_poly_init (plain_poly_t f, int d);
extern void plain_poly_clear (plain_poly_t f);
extern int plain_poly_set_mod (plain_poly_t fp, mpz_t *f, int d, plain_poly_coeff_t p);
extern int plain_poly_fits (unsigned int, uint64_t);
extern int plain_poly_is_irreducible(plain_poly_t fp, const plain_poly_coeff_t p);

/* These do not require the caller to create the plain_poly_t structure
 * beforehand.
 *
 * XXX: Use the rootfinding wrappers poly_roots defined in rootfinder.h
 */
extern int plain_poly_roots (plain_poly_coeff_t*, mpz_t*, int, const plain_poly_coeff_t);
extern int plain_poly_roots_long (long*, mpz_t*, int, const plain_poly_coeff_t);
extern int plain_poly_roots_int64 (int64_t*, mpz_t*, int, const plain_poly_coeff_t);


static inline int roots_mod_long (plain_poly_coeff_t* r, mpz_t* c, int d, const plain_poly_coeff_t p)
    __attribute__((deprecated));
static inline int roots_mod_long (plain_poly_coeff_t* r, mpz_t* c, int d, const plain_poly_coeff_t p)
{ return plain_poly_roots(r,c,d,p); }
static inline int isirreducible_mod_long(plain_poly_t fp, const plain_poly_coeff_t p)
    __attribute__((deprecated));
static inline int isirreducible_mod_long(plain_poly_t fp, const plain_poly_coeff_t p)
{ return plain_poly_is_irreducible(fp,p); }

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_PLAIN_POLY_H_ */
