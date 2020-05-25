#ifndef INT64_POLY_H
#define INT64_POLY_H

#include <stdint.h>
#include "sieving_bound.h"
#include "mpz_poly.h"
#include "macros.h"

typedef struct {
  int alloc;
  int deg;
  int64_t * coeff;
} s_int64_poly_t;

typedef s_int64_poly_t int64_poly_t[1];
typedef s_int64_poly_t * int64_poly_ptr;
typedef const s_int64_poly_t * int64_poly_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Initialize a polynomial.
 *
 * f: the polynomial.
 * d: degree of the polyniomial.
*/
void int64_poly_init(int64_poly_ptr f, int d);

/*
 * Realloc f to (at least) nc coefficients.
 *
 * f: the polynomial.
 * nc: number of coefficients.
 */
void int64_poly_realloc(int64_poly_ptr f, int nc);

/*
 * Copy f to g.
 *
 * g: the receiver polynomial.
 * f: the source polynomial
 */
void int64_poly_copy(int64_poly_ptr g, int64_poly_srcptr f);

/*
 * Swap f and g. g must just be initialized.
 *
 * g: a polynomial.
 * f: a polynomial
 */
void int64_poly_swap(int64_poly_ptr f, int64_poly_ptr g);

/*
  Delete the polynomial.

  f: the polynomial.
*/
void int64_poly_clear(int64_poly_ptr f);

/*
 * Find polynomial degree, deg is an upper bound for the degree.
 *
 * f: the polynomial.
 * deg: upper bound for the degree of f.
 */
void int64_poly_cleandeg(int64_poly_ptr f, int deg);

/*
 * Set f to a zero polynomial.
 *
 * f: the polynomial.
 */
void int64_poly_set_zero(int64_poly_ptr f);


/*
 * Set the ith coefficient of a polynomial.
 *
 * f: the polynomial.
 * i: index of the coefficient.
 * coeff: new value of the ith coefficient.
 */
static inline void int64_poly_setcoeff(int64_poly_ptr f, int i, int64_t coeff);
/*
 * Return the ith coefficient of a polynomial.
 *
 * i: index of the coefficient.
 * f: the polynomial.
 */
static inline int64_t int64_poly_getcoeff(int i, int64_poly_srcptr f);
/*
 * Set f to x^i.
 *
 * f: the polynomial.
 * i: degree of the polynomial.
 */
void int64_poly_set_xi(int64_poly_ptr f, int i);

/*
 * Set f to b*x^i.
 *
 * f: the polynomial.
 * i: degree of the polynomial.
 * b: multiplicative coefficient.
 */
void int64_poly_set_bxi(int64_poly_ptr f, int i, int64_t b);

/*
 * Return 1 if f and g are equal, 0 otherwise. Assumes f and g are normalized.
 *
 * a: a polynomial.
 * b: a polynowial.
 */
int int64_poly_equal(int64_poly_srcptr a, int64_poly_srcptr b);

/*
 * Return 1 if f is normalized, i.e. f[deg] != 0, or the null polynomial.
 *
 * f: the tested polynomial.
 */
int int64_poly_normalized(int64_poly_srcptr f);

/*
 * To write a polynomial in a file.
 *
 * file: the file.
 * f: the polynomial.
 */
void int64_poly_fprintf(FILE * file, int64_poly_srcptr f);

/*
 * Return the maximum of the coefficients of the polynomial.
 *
 * f: the polynomial.
 */
int64_t int64_poly_max(int64_poly_srcptr f);

/*
 * Find the infinite norm of f and set in in.
 *
 * f: the polynomial.
 */
uint64_t int64_poly_infinity_norm(int64_poly_srcptr f);

/*
 * Make a mpz_poly from the int64_poly.
 *
 * a: the new mpz_poly.
 * b: the old int64_poly.
 */
void int64_poly_to_mpz_poly(mpz_poly_ptr a, int64_poly_srcptr b);

#ifdef __cplusplus
}
#endif

static inline void int64_poly_setcoeff(int64_poly_ptr f, int i, int64_t coeff)
{
  int64_poly_realloc (f, i + 1);
  f->coeff[i] = coeff;
  if (i >= f->deg) {
    int64_poly_cleandeg (f, i);
  }
}

static inline int64_t int64_poly_getcoeff(int i, int64_poly_srcptr f)
{
  ASSERT_ALWAYS(f->deg == -1 ||  f->deg >= i);
  if (i > f->deg) {
    return 0;
  } else {
    return f->coeff[i];
  }
}

#endif /* INT64_POLY_H */
