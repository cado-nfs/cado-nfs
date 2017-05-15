/*
   These files (mpz_poly.*) implement arithmetics of polynomials whose
   coefficients are in multiprecision integers (using mpz_t from GNU MP).

   We use them in sqrt/algsqrt.c to represent rings of integers.

   Please see the file README_POLYNOMIALS for more details and
   comparisons to other poly formats.
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <sstream>
#include "portability.h"
#include "mpz_poly.h"
#include "lll.h"
#include "rootfinder.h"
#include "gmp_aux.h"
#include "misc.h"
/* and just because we expose a proxy to usp.c's root finding... */
#include "usp.h"
#include "double_poly.h"

#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif

static inline mpz_ptr mpz_poly_lc(mpz_poly_ptr f)
{
    assert(f->deg >= 0);
    return f->coeff[f->deg];
}

/* --------------------------------------------------------------------------
   Static functions
   -------------------------------------------------------------------------- */

/* Return f=g*h, where g has degree r, and h has degree s. */
static int
mpz_poly_mul_basecase (mpz_t *f, mpz_t *g, int r, mpz_t *h, int s) {
  int i, j;
  assert(f != g && f != h);
  for (i = 0; i <= r + s; i++)
    mpz_set_ui (f[i], 0);
  for (i = 0; i <= r; ++i)
    for (j = 0; j <= s; ++j)
      mpz_addmul(f[i+j], g[i], h[j]);
  return r + s;
}

/* f <- g^2 where g has degree r */
static int
mpz_poly_sqr_basecase (mpz_t *f, mpz_t *g, int r) {
  int i, j;
  assert(f != g);
  for (i = 0; i <= 2 * r; i++)
    mpz_set_ui (f[i], 0);
  for (i = 0; i <= r; ++i)
    for (j = 0; j < i; ++j)
      mpz_addmul(f[i+j], g[i], g[j]);
  for (i = 0; i < 2*r ; i++)
      mpz_mul_2exp(f[i], f[i], 1);
  for (i = 0; i < r ; i++) {
      mpz_mul(f[2*r], g[i], g[i]);
      mpz_add(f[2*i], f[2*i], f[2*r]);
  }
  mpz_mul(f[2*r], g[r], g[r]);
  return 2 * r;
}

/* To interpolate a polynomial of degree <= MAX_TC_DEGREE, we use the following
   MAX_TC_DEGREE evaluation points, plus the point at infinity.
   The first point should always be 0. */
static const long tc_points[MAX_TC_DEGREE] = {0, 1, -1, 2, -2, 3, -3, 4, -4,
					    5, -5, 6, -6, 7, -7, 8, -8, 9, -9};

#if GNUC_VERSION_ATLEAST(7,0,0)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow" 
#endif
/* Given f[0]...f[t] that contain respectively f(0), ..., f(t),
   put in f[0]...f[t] the coefficients of f. Assumes t <= MAX_TC_DEGREE.

   In the square root, with an algebraic polynomial of degree d,
   we have to multiply polynomials of degree d-1, thus we have t=2(d-1).
*/
static void
mpz_poly_mul_tc_interpolate (mpz_t *f, int t) {
  int64_t M[MAX_TC_DEGREE+1][MAX_TC_DEGREE+1], g, h;
  int i, j, k, l;
  /* G[] is the list of gcd's that appear in the forward Gauss loop, in the
     order they appear (they don't depend on t, since we start with the low
     triangular submatrix for t-1) */
  static const int64_t G[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000};

  ASSERT (t <= MAX_TC_DEGREE); /* Ensures that all M[i][j] fit in uint64_t,
				  and similarly for all intermediate
				  computations on M[i][j]. This avoids the
				  use of mpz_t to store the M[i][j]. */

  /* initialize M[i][j] = tc_points[i]^j */
  for (i = 0; i <= t; i++)
    for (j = 0; j <= t; j++)
      {
	if (i < t)
	  M[i][j] = (j == 0) ? 1 : tc_points[i] * M[i][j-1];
	else /* evaluation point is infinity */
	  M[i][j] = (j == t);
      }

  /* Forward Gauss: zero the under-diagonal coefficients while going down.
     Since the last point of evaluation is infinity, the last row is already
     reduced, thus we go up to i=t-1. */
  for (i = 1; i < t; i++)
    {
      for (j = 0, l = 0; j < i; j++)
	{
	  g = G[l++]; /* same as gcd_int64 (M[i][j], M[j][j]) */
	  h = M[i][j] / g;
	  g = M[j][j] / g;
	  /* f[i] <- g*f[i] - h*f[j] */
	  ASSERT(g > 0);
	  if (g != 1)
	    mpz_mul_uint64 (f[i], f[i], g);
	  mpz_submul_int64 (f[i], f[j], h);
	  for (k = j; k <= t; k++)
	    M[i][k] = g * M[i][k] - h * M[j][k];
	}
    }

  /* now zero upper-diagonal coefficients while going up */
  for (i = t - 1; i >= 0; i--)
  {
    for (j = i + 1; j <= t; j++)
      /* f[i] = f[i] - M[i][j] * f[j] */
      mpz_submul_int64 (f[i], f[j], M[i][j]);
    ASSERT (mpz_divisible_uint64_p (f[i], M[i][i]));
    mpz_divexact_uint64 (f[i], f[i], M[i][i]);
  }
}
#if GNUC_VERSION_ATLEAST(7,0,0)
#pragma GCC diagnostic pop
#endif

/* v <- g(i) */
static void
mpz_poly_mul_eval_si (mpz_t v, mpz_t *g, int r, long i)
{
  mpz_set (v, g[r]);
  for (int j = r - 1; j >= 0; j--)
    {
      mpz_mul_si (v, v, i);
      mpz_add (v, v, g[j]);
    }
}

/* Generic Toom-Cook implementation: stores in f[0..r+s] the coefficients
   of g*h, where g has degree r and h has degree s, and their coefficients
   are in g[0..r] and h[0..s].
   Assumes f differs from g and h, and f[0..r+s] are already allocated.
   Returns the degree of f.
*/
static int
mpz_poly_mul_tc (mpz_t *f, mpz_t *g, int r, mpz_t *h, int s)
{
  int t, i;

  if ((r == -1) || (s == -1)) /* g or h is 0 */
    return -1;

  if (r < s)
    return mpz_poly_mul_tc (f, h, s, g, r);

  /* now r >= s */

  if (s == 0)
    {
      for (i = 0; i <= r; i++)
	mpz_mul (f[i], g[i], h[0]);
      return r;
    }

  t = r + s; /* degree of f, which has t+1 coefficients */

  if (t == 2) /* necessary r = s = 1, use vanilla Karatsuba with 3 MUL,
		 this is always optimal */
    {
      mpz_add (f[0], g[0], g[1]);
      mpz_add (f[2], h[0], h[1]);
      mpz_mul (f[1], f[0], f[2]);
      mpz_mul (f[0], g[0], h[0]);
      mpz_mul (f[2], g[1], h[1]);
      mpz_sub (f[1], f[1], f[0]);
      mpz_sub (f[1], f[1], f[2]);
      return 2;
    }

  if (t == 3) /* r = 2, s = 1, the following code in 4 MUL is optimal */
    {
      mpz_add (f[1], g[2], g[0]);
      mpz_add (f[2], f[1], g[1]); /* g(1) */
      mpz_sub (f[1], f[1], g[1]); /* g(-1) */
      mpz_add (f[0], h[0], h[1]); /* h(1) */
      mpz_mul (f[2], f[2], f[0]); /* f(1) = g(1)*h(1) */
      mpz_sub (f[0], h[0], h[1]); /* h(-1) */
      mpz_mul (f[0], f[1], f[0]); /* f(-1) = g(-1)*h(-1) */
      mpz_sub (f[1], f[2], f[0]); /* f(1) - f(-1) = 2*g2*h1+2*g1*h0+2*g0*h1 */
      mpz_add (f[2], f[2], f[0]); /* f(1) + f(-1) = 2*g2*h0+2*g1*h1+2*g0*h0 */
      mpz_div_2exp (f[1], f[1], 1); /* g2*h1 + g1*h0 + g0*h1 */
      mpz_div_2exp (f[2], f[2], 1); /* g2*h0 + g1*h1 + g0*h0 */
      mpz_mul (f[0], g[0], h[0]);
      mpz_sub (f[2], f[2], f[0]); /* g2*h0 + g1*h1 */
      mpz_mul (f[3], g[2], h[1]);
      mpz_sub (f[1], f[1], f[3]); /* g1*h0 + g0*h1 */
      return 3;
    }

  if (r == 2) /* Necessarily s = 2 since r >= s and the cases s = 0 and
		 (r,s) = (2,1) have already been treated. */
    {
      /* we use here the code from Appendix Z from "Towards Optimal Toom-Cook
	 Multiplication for Univariate and Multivariate Polynomials in
	 Characteristic 2 and 0" from Marco Bodrato, WAIFI 2007, LNCS 4547,
	 with the change: W3 -> f[1], W2 -> f[3], W1 -> f[2]. */
      mpz_add (f[0], g[2], g[0]); mpz_add (f[4], h[2], h[0]);
      mpz_sub (f[3], f[0], g[1]); mpz_sub (f[2], f[4], h[1]);
      mpz_add (f[0], f[0], g[1]); mpz_add (f[4], f[4], h[1]);
      mpz_mul (f[1], f[3], f[2]); mpz_mul (f[2], f[0], f[4]);
      mpz_add (f[0], f[0], g[2]);
      mpz_mul_2exp (f[0], f[0], 1);
      mpz_sub (f[0], f[0], g[0]);
      mpz_add (f[4], f[4], h[2]);
      mpz_mul_2exp (f[4], f[4], 1);
      mpz_sub (f[4], f[4], h[0]);
      mpz_mul (f[3], f[0], f[4]);
      mpz_mul (f[0], g[0], h[0]);
      mpz_mul (f[4], g[2], h[2]);
      /* interpolation */
      mpz_sub (f[3], f[3], f[1]);
      ASSERT (mpz_divisible_ui_p (f[3], 3));
      mpz_divexact_ui (f[3], f[3], 3);
      mpz_sub (f[1], f[2], f[1]);
      ASSERT (mpz_divisible_ui_p (f[1], 2));
      mpz_div_2exp (f[1], f[1], 1); /* exact */
      mpz_sub (f[2], f[2], f[0]);
      mpz_sub (f[3], f[3], f[2]);
      ASSERT (mpz_divisible_ui_p (f[3], 2));
      mpz_div_2exp (f[3], f[3], 1); /* exact */
      mpz_submul_ui (f[3], f[4], 2);
      mpz_sub (f[2], f[2], f[1]);
      mpz_sub (f[2], f[2], f[4]);
      mpz_sub (f[1], f[1], f[3]);
      return 4;
    }

  if (t > MAX_TC_DEGREE) {
    /* naive product */
    /* currently we have to resort to this for larger degree, because
     * the generic Toom implementation is bounded in degree.
     */
    return mpz_poly_mul_basecase (f, g, r, h, s);
  }

  /* now t <= MAX_TC_DEGREE */

  /* store g(tc_points[i])*h(tc_points[i]) in f[i] for 0 <= i <= t */

  /* the first evaluation point is 0 */
  ASSERT(tc_points[0] == 0);
  mpz_mul (f[0], g[0], h[0]);

  for (i = 1; i < t; i++)
    {
      /* f[i] <- g(i) */
      mpz_poly_mul_eval_si (f[i], g, r, tc_points[i]);
      /* f[t] <- h(i) */
      mpz_poly_mul_eval_si (f[t], h, s, tc_points[i]);
      /* f[i] <- g(i)*h(i) */
      mpz_mul (f[i], f[i], f[t]);
    }

  /* last evaluation point is infinity */
  mpz_mul (f[t], g[r], h[s]);

  mpz_poly_mul_tc_interpolate (f, t);

  return t;
}

/* Same as mpz_poly_mul_tc for the squaring: store in f[0..2r] the coefficients
   of g^2, where g has degree r, and its coefficients are in g[0..r].
   Assumes f differs from g, and f[0..2r] are already allocated.
   Returns the degree of f.
*/
static int
mpz_poly_sqr_tc (mpz_t *f, mpz_t *g, int r)
{
  int t, i;
  size_t nbits;

  if (r == -1) /* g is 0 */
    return -1;

  if (r == 0)
    {
      mpz_mul (f[0], g[0], g[0]);
      return 0;
    }

  if (r == 1) /* 3 SQR: always optimal */
    {
      mpz_mul (f[0], g[0], g[0]);
      mpz_mul (f[2], g[1], g[1]);
      mpz_add (f[1], g[0], g[1]);
      mpz_mul (f[1], f[1], f[1]);
      mpz_sub (f[1], f[1], f[0]);
      mpz_sub (f[1], f[1], f[2]);
      return 2;
    }

  /* we assume the number of bits of f[0] is representative of the size of
     the coefficients of f */
  nbits = mpz_sizeinbase (f[0], 2);

  if (r == 2 && nbits < 4096)
    /* (g2*x^2+g1*x+g0)^2 = g2^2*x^4 + (2*g2*g1)*x^3 +
       (g1^2+2*g2*g0)*x^2 + (2*g1*g0)*x + g0^2 in 3 SQR + 2 MUL.
       Experimentally this is faster for less than 4096 bits than the
       algorithm in 5 SQR from mpz_poly_mul_tc_interpolate. */
    {
      mpz_mul (f[4], g[2], g[2]); /* g2^2 */
      mpz_mul (f[3], g[2], g[1]);
      mpz_mul_2exp (f[3], f[3], 1); /* 2*g2*g1 */
      mpz_mul (f[1], g[1], g[0]);
      mpz_mul_2exp (f[1], f[1], 1); /* 2*g1*g0 */
      mpz_mul (f[0], g[0], g[0]); /* g0^2 */
      mpz_add (f[2], g[2], g[1]);
      mpz_add (f[2], f[2], g[0]);
      mpz_mul (f[2], f[2], f[2]); /* (g2+g1+g0)^2 */
      mpz_sub (f[2], f[2], f[4]);
      mpz_sub (f[2], f[2], f[0]); /* g1^2 + 2*g2*g0 + 2*g2*g1 + 2*g1*g0 */
      mpz_sub (f[2], f[2], f[3]); /* g1^2 + 2*g2*g0 + 2*g1*g0 */
      mpz_sub (f[2], f[2], f[1]); /* g1^2 + 2*g2*g0 */
      return 4;
    }

  if (r == 3 && nbits < 4096)
    /* (g3*x^3+g2*x^2+g1*x+g0)^2 = g3^2*x^6 + (2*g3*g2)*x^5
       + (g2^2+2*g3*g1)*x^4 + (2*g3*g0+2*g2*g1)*x^3
       + (g1^2+2*g2*g0)*x^2 + (2*g1*g0)*x + g0^2 in 4 SQR + 6 MUL.
       Experimentally this is faster for less than 4096 bits than the
       algorithm in 7 SQR from mpz_poly_mul_tc_interpolate. */
    {
      mpz_mul (f[6], g[3], g[3]); /* g3^2 */
      mpz_mul (f[5], g[3], g[2]);
      mpz_mul_2exp (f[5], f[5], 1); /* 2*g3*g2 */
      mpz_mul (f[1], g[1], g[0]);
      mpz_mul_2exp (f[1], f[1], 1); /* 2*g1*g0 */
      mpz_mul (f[2], g[1], g[1]); /* g1^2 */
      mpz_mul (f[0], g[2], g[0]);
      mpz_addmul_ui (f[2], f[0], 2); /* g1^2+2*g2*g0 */
      mpz_mul (f[4], g[2], g[2]); /* g2^2 */
      mpz_mul (f[0], g[3], g[1]);
      mpz_addmul_ui (f[4], f[0], 2); /* g2^2+2*g3*g1 */
      mpz_mul (f[3], g[3], g[0]);
      mpz_mul (f[0], g[2], g[1]);
      mpz_add (f[3], f[3], f[0]);
      mpz_mul_2exp (f[3], f[3], 1); /* 2*g3*g0+2*g2*g1 */
      mpz_mul (f[0], g[0], g[0]); /* g0^2 */
      return 6;
    }

  t = 2 * r; /* product has degree t thus t+1 coefficients */
  if (t > MAX_TC_DEGREE) {
    /* naive product */
    /* currently we have to resort to this for larger degree, because
     * the generic toom implementation is bounded in degree.
     */
    return mpz_poly_sqr_basecase (f, g, r);
  }

  ASSERT (t <= MAX_TC_DEGREE);

  /* store g(tc_points[i])^2 in f[i] for 0 <= i <= t */

  /* first evaluation point is 0 */
  ASSERT(tc_points[0] == 0);
  mpz_mul (f[0], g[0], g[0]);

  for (i = 1; i < t; i++)
    {
      /* f[i] <- g(i) */
      mpz_poly_mul_eval_si (f[i], g, r, tc_points[i]);
      mpz_mul (f[i], f[i], f[i]);
    }

  /* last evaluation point is infinity */
  mpz_mul (f[t], g[r], g[r]);

  mpz_poly_mul_tc_interpolate (f, t);

  return t;
}

/* Pseudo-reduce a plain polynomial p modulo a non-monic polynomial F.
   The result is of type polymodF_t P, and satisfies:
   P->p = lc(F)^P->v * p mod F.
   WARNING: this function destroys its input p !!! */
static void
mpz_poly_reducemodF(polymodF_t P, mpz_poly_ptr p, mpz_poly_srcptr F) {
  int v = 0;

  if (p->deg < F->deg) {
    mpz_poly_set(P->p, p);
    P->v = 0;
    return;
  }

  const int d = F->deg;

  while (p->deg >= d) {
    const int k = p->deg;
    int i;

    /* We compute F[d]*p - p[k]*F. In case F[d] divides p[k], we can simply
       compute p - p[k]/F[d]*F. However this will happen rarely with
       Kleinjung's polynomial selection, since lc(F) is large. */

    /* FIXME: in msieve, Jason Papadopoulos reduces by F[d]^d*F(x/F[d])
       instead of F(x). This might avoid one of the for-loops below. */

    // temporary hack: account for the possibility that we're indeed
    // using f_hat instead of f.
    if (mpz_cmp_ui(F->coeff[d], 1) != 0) {
      v++; /* we consider p/F[d]^v */
      for (i = 0; i < k; ++i)
        mpz_mul (p->coeff[i], p->coeff[i], F->coeff[d]);
    }

    for (i = 0; i < d; ++i)
      mpz_submul (p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

    mpz_poly_cleandeg (p, k-1);
  }

  mpz_poly_set(P->p, p);
  P->v = v;
}



/* --------------------------------------------------------------------------
   Public functions
   -------------------------------------------------------------------------- */


/* Management of the structure, set and print coefficients. */


/* Allocate a polynomial that can contain 'd+1' coefficients and set to zero.
   We allow d < 0, which is equivalent to d = -1.
 */
void mpz_poly_init(mpz_poly_ptr f, int d)
{
  f->deg = -1;
  if (d < 0)
  {
    f->alloc = 0;
    f->coeff = (mpz_t *) NULL;
  }
  else
  {
    int i;
    f->alloc = d+1;
    f->coeff = (mpz_t *) malloc ((d+1)*sizeof(mpz_t));
    FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
    for (i = 0; i <= d; ++i)
      mpz_init (f->coeff[i]);
  }
}


/* realloc f to (at least) nc coefficients */
void mpz_poly_realloc (mpz_poly_ptr f, int nc)
{
  int i;
  if (f->alloc < nc)
    {
      f->coeff = (mpz_t*) realloc (f->coeff, nc * sizeof (mpz_t));
      FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
      for (i = f->alloc; i < nc; i++)
        mpz_init (f->coeff[i]);
      f->alloc = nc;
    }
}

/* Copy f to g, where g must be initialized (but there is no need it has
   enough allocated coefficients, since mpz_poly_setcoeff reallocates if
   needed). */
void mpz_poly_set(mpz_poly_ptr g, mpz_poly_srcptr f) 
{
  int i;

  if (g == f)
    return; /* nothing to do */

  g->deg = f->deg;
  for (i = f->deg; i >= 0; --i)
    mpz_poly_setcoeff (g, i, f->coeff[i]);
}

void mpz_poly_set_double_poly(mpz_poly_ptr g, double_poly_srcptr f) 
{
  mpz_poly_realloc (g, f->deg + 1);
  g->deg = f->deg;
  for (int i = f->deg; i >= 0; --i)
    mpz_set_d (g->coeff[i], f->coeff[i]);
}

/* Init polynomial rel and set it to a - b*x */
void
mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b)
{
    mpz_poly_init(rel, 1);
    mpz_poly_setcoeff_int64(rel, 0, a);
    mpz_poly_setcoeff_int64(rel, 1, -b);
}

/* rel <- a-b*x */
void
mpz_poly_init_set_mpz_ab (mpz_poly_ptr rel, mpz_t a, mpz_t b)
{
    mpz_poly_init(rel, 1);
    mpz_poly_setcoeff(rel, 0, a);
    mpz_t mb;
    mpz_init(mb);
    mpz_neg(mb, b);
    mpz_poly_setcoeff(rel, 1, mb);
    mpz_clear(mb);
}

/* swap f and g */
void
mpz_poly_swap (mpz_poly_ptr f, mpz_poly_ptr g)
{
  int i;
  mpz_t *t;

  i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->deg;
  f->deg = g->deg;
  g->deg = i;
  t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

/* Free polynomial f in mpz_poly. */
void mpz_poly_clear(mpz_poly_ptr f) 
{
  int i;
  for (i = 0; i < f->alloc; ++i)
    mpz_clear(f->coeff[i]);
  if (f->coeff != NULL)
    free(f->coeff);
  f->coeff = NULL; /* to avoid a double-free */
  memset(f, 0, sizeof(mpz_poly));
  f->deg = -1;
  f->alloc = 0; /* to avoid a double-free */
}

/* Return 0 if f[i] is zero, -1 is f[i] is negative and +1 if f[i] is positive,
   like mpz_sgn function. */
static inline int mpz_poly_coeff_sgn (mpz_poly_srcptr f, int i)
{
  if (i >= f->alloc)
    return 0;
  else
    return mpz_sgn (f->coeff[i]);
}

/* removed mpz_poly_set_deg, as for all purposes there is no reason to
 * not use the more robust mpz_poly_cleandeg */

/* Find polynomial degree. */
void mpz_poly_cleandeg(mpz_poly_ptr f, int deg)
{
  ASSERT(deg >= -1);
  while ((deg >= 0) && (mpz_poly_coeff_sgn (f, deg)==0))
    deg--;
  f->deg = deg;
}

/* Sets f to the polynomial of degree d, of coefficients
   given by coeffs. */
void mpz_poly_setcoeffs (mpz_poly_ptr f, mpz_t * coeffs, int d)
{
  int i;
  for (i=d; i>=0; --i)
    mpz_poly_setcoeff(f, i, coeffs[i]);
  mpz_poly_cleandeg(f, d);
}

/* Set a zero polynomial. */
void mpz_poly_set_zero(mpz_poly_ptr f) 
{
  f->deg = -1;
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff (mpz_poly_ptr f, int i, mpz_srcptr z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* Set signed int coefficient for the i-th term. */
void mpz_poly_setcoeff_si(mpz_poly_ptr f, int i, long z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set_si (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* Set unsigned int coefficient for the i-th term. */
void mpz_poly_setcoeff_ui(mpz_poly_ptr f, int i, unsigned long z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set_ui (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* Set int64 coefficient for the i-th term. */
void mpz_poly_setcoeff_int64(mpz_poly_ptr f, int i, int64_t z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set_int64 (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* Set uint64 coefficient for the i-th term. */
void mpz_poly_setcoeff_uint64(mpz_poly_ptr f, int i, uint64_t z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set_uint64 (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* f[i] <- z */
void mpz_poly_setcoeff_double(mpz_poly_ptr f, int i, double z)
{
  mpz_poly_realloc (f, i + 1);
  mpz_set_d (f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg (f, i);
}

/* Get coefficient for the i-th term. */
void mpz_poly_getcoeff(mpz_t res, int i, mpz_poly_srcptr f)
{
    if (i > f->deg)
        mpz_set_ui (res, 0);
    else
        mpz_set (res, f->coeff[i]);
}

/* x^i is often useful */
void mpz_poly_set_xi(mpz_poly_ptr f, int i)
{
    mpz_poly_realloc (f, i + 1);
    for(int j = 0 ; j <= i ; j++) {
        mpz_set_ui(f->coeff[j], j == i);
    }
    f->deg = i;
}

/* g <- quo (f, x^i) */
void mpz_poly_div_xi(mpz_poly_ptr g, mpz_poly_srcptr f, int i)
{
    if (f->deg < i) {
        mpz_poly_set_zero(g);
        return;
    }
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }
    if (g == f) {
        /* rotate the coefficients, don't do any freeing */
        mpz_t * temp = (mpz_t*) malloc(i * sizeof(mpz_t));
        memcpy(temp, g->coeff, i * sizeof(mpz_t));
        memmove(g->coeff, g->coeff + i, (g->deg + 1 - i) * sizeof(mpz_t));
        memcpy(g->coeff + (g->deg + 1 - i), temp, i * sizeof(mpz_t));
        g->deg -= i;
        free(temp);
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg - i);
    for(int j = i ; j <= f->deg ; j++) {
        mpz_set(g->coeff[j - i], f->coeff[j]);
    }
    g->deg = f->deg - i;
}

/* g <- f * x^i */
void
mpz_poly_mul_xi (mpz_poly_ptr g, mpz_poly_srcptr f, int i)
{
    if (i == 0) {
        mpz_poly_set(g, f);
        return;
    }

    mpz_poly_realloc(g, 1 + f->deg + i);

    if (g == f) {
        for(int j = g->deg ; j >= 0 ; j--) {
            mpz_swap(g->coeff[j + i], g->coeff[j]);
        }
    } else {
        for(int j = g->deg ; j >= 0 ; j--) {
            mpz_set(g->coeff[j + i], g->coeff[j]);
        }
    }
    for(int j = 0 ; j < i ; j++) {
        mpz_set_ui(g->coeff[j], 0);
    }
    g->deg = f->deg + i;
}

/* g <- f * (x + a) */
void mpz_poly_mul_xplusa(mpz_poly_ptr g, mpz_poly_srcptr f, mpz_srcptr a)
{
    mpz_t aux;
    mpz_init_set_ui(aux, 0);
    mpz_poly_realloc(g, f->deg + 2);
    for(int i = 0 ; i <= f->deg ; i++) {
        /* aux is is coeff of degree [i-1]. We want
         * (coeff_i, aux) <- (coeff_{i-1} + a * coeff_i, coeff_i)
         *                   (aux + a * coeff_i, coeff_i)
         */
        if (a) mpz_addmul(aux, f->coeff[i], a);
        mpz_swap(g->coeff[i], aux);
        if (f != g) {
            mpz_set(aux, f->coeff[i]);
        }
    }
    /* last coefficient */
    mpz_swap(g->coeff[f->deg + 1], aux);
    /* This is just as valid as it was for f */
    g->deg = f->deg + 1;
    mpz_clear(aux);
}

/* return the valuation of f */
int mpz_poly_valuation(mpz_poly_srcptr f)
{
    int n = 0;
    assert(f->deg >= 0);
    for( ; n < f->deg  && mpz_cmp_ui(f->coeff[n], 0) == 0 ; n++) ;
    return n;
}


int mpz_poly_asprintf(char ** res, mpz_poly_srcptr f)
{
    size_t alloc = 0;
    size_t size = 0;
    int rc;
    static const size_t batch = 4;
    *res = NULL;

    alloc += batch;
    *res = (char*) realloc(*res, alloc);
    if (*res == NULL) goto oom;

#define SNPRINTF_FRAGMENT(fmt, arg) do {				\
    rc = gmp_snprintf(*res + size, alloc - size, fmt, arg);     	\
    if (size + (size_t) rc + 1 >= alloc) {				\
        alloc = MAX(alloc + batch, size + (size_t) rc + 1);             \
        *res = (char*) realloc(*res, alloc);				\
        if (*res == NULL) goto oom;					\
        rc = gmp_snprintf(*res + size, alloc - size, fmt, arg); 	\
    }									\
    size += rc;								\
    ASSERT_ALWAYS(size < alloc);					\
} while (0)

#define PUTS_FRAGMENT(arg) do {						\
    rc = strlen(arg);                                                   \
    if (size + (size_t) rc + 1 >= alloc) {      			\
        alloc = MAX(alloc + batch, size + (size_t) rc + 1);		\
        *res = (char*) realloc(*res, alloc);				\
        if (*res == NULL) goto oom;					\
    }									\
    size += strlcpy(*res + size, arg, alloc-size);			\
    ASSERT_ALWAYS(size < alloc);					\
} while (0)

    if (f->deg == -1) {
        SNPRINTF_FRAGMENT("%d", 0);
        return size;
    }
    for (int i = 0, printed = 0; i <= f->deg; ++i) {
        if (mpz_cmp_ui(f->coeff[i], 0) == 0) continue;
        if (printed++ && mpz_cmp_ui(f->coeff[i], 0) > 0)
            PUTS_FRAGMENT ("+");
        if (i && mpz_cmp_ui(f->coeff[i], 1) == 0) {
            PUTS_FRAGMENT ("x");
        } else if (i && mpz_cmp_si(f->coeff[i], -1) == 0) {
            PUTS_FRAGMENT ("-x");
        } else {
            SNPRINTF_FRAGMENT ("%Zd", f->coeff[i]);
            if (i) {
                PUTS_FRAGMENT ("*x");
            }
        }

        if (i > 1) {
            SNPRINTF_FRAGMENT ("^%d", i);
        }
    }

    return size;
#undef SNPRINTF_FRAGMENT
#undef PUTS_FRAGMENT
oom:
    free(res);
    return -1;
}

/* Print coefficients of f.
 * endl = 1 if "\n" at the end of fprintf. */
void mpz_poly_fprintf_endl (FILE *fp, mpz_poly_srcptr f, int endl)
{
    char * res;
    int rc = mpz_poly_asprintf(&res, f);
    ASSERT_ALWAYS(rc >= 0);
    fprintf(fp, "%s", res);
    if (endl) {
      fprintf(fp, "\n");
    }
    free(res);
}

/* Print coefficients of f. */
void mpz_poly_fprintf (FILE *fp, mpz_poly_srcptr f)
{
  mpz_poly_fprintf_endl (fp, f, 1);
}

/* Print f of degree d with the following format
    f0<sep>f1<sep>...<sep>fd\n
   Print only '\n' if f = 0 (ie deg(f) = -1)
*/
void mpz_poly_fprintf_coeffs (FILE *fp, mpz_poly_srcptr f, const char sep)
{
  if (f->deg >= 0)
  {
    gmp_fprintf (fp, "%Zd", f->coeff[0]);
    for (int i = 1; i <= f->deg; i++)
      gmp_fprintf (fp, "%c%Zd", sep, f->coeff[i]);
  }
  fprintf (fp, "\n");
}

/* Print f of degree d with the following format
    <pre><letter>0: f0\n
    <pre><letter>1: f1\n
    ...
    <pre><letter>d: fd\n
   Print nothing if f = 0 (ie deg(f) = -1)
*/
void
mpz_poly_fprintf_cado_format (FILE *fp, mpz_poly_srcptr f, const char letter,
                              const char *prefix)
{
  for (int i = 0; i <= f->deg; i++)
  {
    if (prefix)
      fputs (prefix, fp);
    gmp_fprintf (fp, "%c%d: %Zd\n", letter, i, f->coeff[i]);
  }
}

void mpz_poly_print_raw(mpz_poly_srcptr f){
    cxx_mpz_poly F;
    mpz_poly_set(F, f);
    std::string s = F.print_poly("x");
    printf("%s\n", s.c_str());
}

/* -------------------------------------------------------------------------- */

/* Tests and comparison functions */


/* return 0 if f and g are equal,
 * -1 if f is "smaller" and 1 if f is "bigger", for some arbitrary
 * ordering (lowest degree first, then lex order).
 *
 * Assumes f and g are normalized */
int mpz_poly_cmp (mpz_poly_srcptr a, mpz_poly_srcptr b)
{
    int r = (a->deg > b->deg) - (b->deg > a->deg);
    if (r) return r;
    for(int d = a->deg; d >= 0 ; d--) {
        r = mpz_cmp(a->coeff[d], b->coeff[d]);
        if (r) return r;
    }
    return 0;
}

/* return 1 if f is normalized, i.e. f[deg] != 0, or the null polynomial.  */
int mpz_poly_normalized_p (mpz_poly_srcptr f)
{
  return (f->deg == -1) || mpz_cmp_ui (f->coeff[f->deg], 0) != 0;
}

/* return 1 if f is nonmonic, i.e. f[deg] != 1, return 0 otherwise (null
 * polynomial is considered monic).
 */
int mpz_poly_is_nonmonic (mpz_poly_srcptr f)
{
  if (f->deg == -1)
    return 0;
  else if (mpz_cmp_ui (f->coeff[f->deg], 1) == 0)
    return 0;
  else
    return 1;
}

/* -------------------------------------------------------------------------- */

/* Set f=-g.
   Note: f can be the same as g. */
void mpz_poly_neg (mpz_poly_ptr f, mpz_poly_srcptr g) {
  mpz_poly_realloc(f, g->deg + 1);
  for (int i = 0 ; i <= g->deg ; i++) {
    mpz_neg (f->coeff[i], g->coeff[i]);
  }
}

/* Set f=g+h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_add(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h) {
  int i, maxdeg;
  mpz_t z;
  mpz_init(z);
  maxdeg = max(g->deg, h->deg);
  mpz_poly_realloc(f, maxdeg + 1);
  for (i = 0 ; i <= maxdeg ; i++) {
    if (i <= g->deg)
      mpz_set(z, g->coeff[i]);
    else
      mpz_set_ui(z, 0);
    if (i <= h->deg)
      mpz_add(z, z, h->coeff[i]);
    mpz_set(f->coeff[i], z);
  }
  f->deg = maxdeg;
  mpz_clear(z);
  mpz_poly_cleandeg(f, maxdeg);
}

/* Set f=g-h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_sub(mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h) {
  int i, maxdeg;
  mpz_t z;
  mpz_init(z);
  maxdeg = max(g->deg, h->deg);
  mpz_poly_realloc(f, maxdeg + 1);
  for (i = 0 ; i <= maxdeg ; i++) {
    if (i <= g->deg)
      mpz_set(z, g->coeff[i]);
    else
      mpz_set_ui(z, 0);
    if (i <= h->deg)
      mpz_sub(z, z, h->coeff[i]);
    mpz_set(f->coeff[i], z);
  }
  f->deg = maxdeg;
  mpz_clear(z);
  mpz_poly_cleandeg(f, maxdeg);
}

/* Set g=f+a where a is unsigned long. */
void
mpz_poly_add_ui (mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_add_ui(g->coeff[0], f->coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->coeff[0], a);
        g->deg = 0;
    }
}

/* Set g=f-a where a is unsigned long. */
void
mpz_poly_sub_ui (mpz_poly_ptr g, mpz_poly_srcptr f, unsigned long a)
{
    mpz_poly_set(g, f);
    mpz_poly_realloc(g, 1);
    if (f->deg >= 0) {
        mpz_sub_ui(g->coeff[0], f->coeff[0], a);
        g->deg = f->deg;
        if (g->deg == 0)
            mpz_poly_cleandeg(g, g->deg);
    } else {
        mpz_set_ui(g->coeff[0], a);
        mpz_neg(g->coeff[0], g->coeff[0]);
        g->deg = 0;
    }
}

/* Set f=g-h (mod m). */
void
mpz_poly_sub_mod_mpz (mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h, mpz_srcptr m)
{
    mpz_poly_sub(f, g, h);
    mpz_poly_mod_mpz(f, f, m, NULL);
}

/* Set f=g*h. Note: f might equal g or h.
   Assumes the g[g->deg] and h[h->deg[] are not zero. */
void
mpz_poly_mul (mpz_poly_ptr f, mpz_poly_srcptr g, mpz_poly_srcptr h) {
  if (f == h || f == g)
    {
      mpz_poly aux;
      mpz_poly_init (aux, -1);
      mpz_poly_mul (aux, g, h);
      mpz_poly_set (f, aux);
      mpz_poly_clear (aux);
      return;
    }

  /* now f differs from g and h */

  if ((g->deg == -1) || (h->deg == -1)) {
    f->deg = -1;
    return;
  }

  mpz_poly_realloc (f, g->deg + h->deg + 1);

  if (g == h) /* this is a square */
    {
      f->deg = mpz_poly_sqr_tc (f->coeff, g->coeff, g->deg);
      return;
    }

  if (g->deg == 0) {
    mpz_poly_mul_mpz (f, h, g->coeff[0]);
    return;
  }

  if (h->deg == 0) {
    mpz_poly_mul_mpz (f, g, h->coeff[0]);
    return;
  }

  ASSERT(mpz_cmp_ui (g->coeff[g->deg], 0) != 0);
  ASSERT(mpz_cmp_ui (h->coeff[h->deg], 0) != 0);
  ASSERT(f != g);
  ASSERT(f != h);

#if 1
  f->deg = mpz_poly_mul_tc (f->coeff, g->coeff, g->deg, h->coeff, h->deg);
#else /* segmentation, this code has problem with huge runs, for example
         degree 5 with lifting to 631516975 bits */
  {
    mpz_t G, H;
    size_t sg, sh, s;
    int i;

    mpz_init (G);
    mpz_init (H);
    sg = mpz_poly_sizeinbase (g, 2);
    sh = mpz_poly_sizeinbase (h, 2);
    /* the +1 accounts for a possible sign */
    for (s = sg + sh + 1, i = (g->deg >= h->deg) ? h->deg + 1 : g->deg + 1;
	 i > 1; i = (i + 1) / 2, s++);
    mpz_set (G, g->coeff[g->deg]);
    for (i = g->deg - 1; i >= 0; i--)
    {
      mpz_mul_2exp (G, G, s);
      mpz_add (G, G, g->coeff[i]);
    }
    /* sanity check: G should have sizeinbase(lc(g))+d*s bits (or -1) */
    size_t size_g = mpz_sizeinbase (g->coeff[g->deg], 2) + g->deg * s;
    ASSERT(mpz_sizeinbase (G, 2) == size_g ||
	   mpz_sizeinbase (G, 2) == size_g - 1);
    mpz_set (H, h->coeff[h->deg]);
    for (i = h->deg - 1; i >= 0; i--)
    {
      mpz_mul_2exp (H, H, s);
      mpz_add (H, H, h->coeff[i]);
    }
    /* sanity check: H should have sizeinbase(lc(h))+d*s bits (or -1) */
    size_t size_h = mpz_sizeinbase (h->coeff[h->deg], 2) + h->deg * s;
    ASSERT(mpz_sizeinbase (H, 2) == size_h ||
	   mpz_sizeinbase (H, 2) == size_h - 1);
    size_g = mpz_sizeinbase (G, 2);
    size_h = mpz_sizeinbase (H, 2);
    /* sanity check: we verify that the product agrees both mod B and B-1 */
    mp_limb_t g0 = mpz_getlimbn (G, 0);
    mpz_mul (G, G, H);
    ASSERT(mpz_getlimbn (G, 0) == g0 * mpz_getlimbn (H, 0));
    ASSERT(mpz_sizeinbase (G, 2) == size_g + size_h ||
	   mpz_sizeinbase (G, 2) == size_g + size_h - 1);
    for (i = 0; i < g->deg + h->deg; i++)
    {
      mpz_fdiv_r_2exp (f->coeff[i], G, s);
      if (mpz_sizeinbase (f->coeff[i], 2) == s)
      {
        mpz_cdiv_r_2exp (f->coeff[i], G, s);
        mpz_cdiv_q_2exp (G, G, s);
      }
      else
        mpz_fdiv_q_2exp (G, G, s);
      ASSERT(mpz_sizeinbase (f->coeff[i], 2) < s);
    }
    mpz_set (f->coeff[i], G);
    mpz_clear (G);
    mpz_clear (H);
    f->deg = g->deg + h->deg;
  }
#endif
  /* there is no need to run mpz_poly_cleandeg since g[g->deg] <> 0
     and h[h->deg] <> 0 */
  ASSERT(mpz_cmp_ui (f->coeff[f->deg], 0) != 0);
}

/* Set Q=a*P, where a is an mpz_t */
void
mpz_poly_mul_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a)
{
  int i;
  mpz_t aux;

  mpz_init (aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i)
    {
      mpz_mul (aux, P->coeff[i], a);
      mpz_poly_setcoeff (Q, i, aux);
    }
  mpz_clear (aux);
}

/* Set Q=P/a, where a is an mpz_t. Assume a divides the content of P (the
   division is done with mpz_divexact). Otherwise the result is not correct. */
void
mpz_poly_divexact_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr a)
{
  mpz_t aux;
  mpz_init (aux);
  Q->deg = P->deg;
  for (int i = 0; i <= P->deg; ++i)
    {
      mpz_divexact (aux, P->coeff[i], a);
      mpz_poly_setcoeff (Q, i, aux);
    }
  mpz_clear (aux);
}

/* Test whether a divides P, where a is an mpz_t. */
int
mpz_poly_divisible_mpz (mpz_poly_srcptr P, mpz_srcptr a)
{
  for (int i = 0; i <= P->deg; ++i)
      if (!mpz_divisible_p(P->coeff[i], a)) return 0;
  return 1;
}

/* Set ft(x) = f(x+k) where k is an mpz_t, ft and f can be the same poly. */
void
mpz_poly_translation (mpz_poly_ptr ft, mpz_poly_srcptr f, const mpz_t k)
{
  int i, j;
  int d = f->deg;

  mpz_poly_set (ft, f);
  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul (ft->coeff[j], ft->coeff[j+1], k);
}

/* Set fr = f + k * x^t * g such that t+deg(g) <= deg(f) and t >= 0 (those two
 * assumptions are not checked). fr and f can be the same poly */
void mpz_poly_rotation (mpz_poly_ptr fr, mpz_poly_srcptr f, mpz_poly_srcptr g,
                        const mpz_t k, int t)
{
  mpz_poly_set (fr, f);
  for (int i = 0; i <= g->deg; i++)
    mpz_addmul (fr->coeff[i+t], g->coeff[i], k);
}

/* Set f = f + k * g such that deg(g) <= deg(f) (this assumption is not
   checked). */
void
mpz_poly_addmul_si (mpz_poly_ptr f, mpz_poly_srcptr g, long k)
{
  for (int i = 0; i <= g->deg; i++)
    mpz_addmul_si (f->coeff[i], g->coeff[i], k);
}

/* Set f = k * g such that deg(g) <= deg(f) (this assumption is not
   checked). */
void
mpz_poly_mul_si (mpz_poly_ptr f, mpz_poly_srcptr g, long k)
{
  for (int i = 0; i <= g->deg; i++)
    mpz_mul_si (f->coeff[i], g->coeff[i], k);
}

/* Set h = fr + k * x^t * g such that t+deg(g) <= deg(f) and t >= 0 (those two
 * assumptions are not checked). fr and f can be the same poly */
void mpz_poly_rotation_int64 (mpz_poly_ptr fr, mpz_poly_srcptr f,
                              mpz_poly_srcptr g, const int64_t k, int t)
{
  mpz_poly_set (fr, f);
  for (int i = 0; i <= g->deg; i++)
    mpz_addmul_int64 (fr->coeff[i+t], g->coeff[i], k);
}

/* h=rem(h, f) mod N, f not necessarily monic, N not necessarily prime */
/* Coefficients of f must be reduced mod N
 * Coefficients of h need not be reduced mod N on input, but are reduced
 * on output.
 */
static int
mpz_poly_pseudodiv_r (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_srcptr N, mpz_t factor)
{
  int i, d = f->deg, dh = h->deg;
  mpz_t tmp, inv;

  mpz_init_set_ui (inv, 1);
  mpz_init_set_ui (tmp, 1);

  mpz_set (tmp, f->coeff[d]);
  if (mpz_cmp_ui(tmp, 1) != 0) {
    /* inv is 1/f[d] mod N */
    if (!mpz_invert (inv, tmp, N)) {
      if (factor != NULL)
        mpz_gcd(factor, tmp, N);
      mpz_clear(inv);
      mpz_clear(tmp);
      return 0;
    }
  }

  while (dh >= d)
  {
    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(inv, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], inv);
      mpz_mod (h->coeff[dh], h->coeff[dh], N);
    }

    for (i = 0; i < d; i++) {
      mpz_mul (tmp, h->coeff[dh], f->coeff[i]);
      mpz_mod (tmp, tmp, N);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], tmp);
      mpz_mod (h->coeff[dh - d + i], h->coeff[dh - d + i], N);
    }

    do {
      dh --;
    }
    while (dh >= 0 && mpz_divisible_p (h->coeff[dh], N));

    h->deg = dh;
  }

  mpz_clear (inv);
  mpz_clear (tmp);
  return 1;
}

/* h=rem(h, f) mod p, f not necessarily monic. */
/* returns 0 if an inverse of the leading coeff could not be found. */
int
mpz_poly_div_r (mpz_poly_ptr h, mpz_poly_srcptr f, mpz_srcptr p)
{
  return mpz_poly_pseudodiv_r (h, f, p, NULL);
}

/*
   computes q, r such that f = q*g + r mod p, with deg(r) < deg(g)
   and p in mpz_t
   q and r must be allocated!
   the only aliasing allowed is f==r.
*/
int mpz_poly_div_qr (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_srcptr p)
{
  int k, j, df = f->deg, dg = g->deg, dq = df - dg;
  mpz_t lg, invlg;

  if (df < dg) /* f is already reduced mod g */
    {
      mpz_poly_set_zero (q);
      mpz_poly_set (r, f);
      return 1;
    }

  /* now df >= dg */
  mpz_poly_realloc(q, dq + 1);

  mpz_init(lg);
  mpz_init_set_ui(invlg, 1);

  mpz_poly_set(r, f);
  q->deg = dq;

  mpz_set (lg, g->coeff[dg]);
  mpz_mod (lg, lg, p);
  /* invlg = 1/g[dg] mod p */
  if (mpz_cmp_ui(lg, 1) != 0)
    if (!mpz_invert(invlg, lg, p)) {
        mpz_clear(lg);
        mpz_clear(invlg);
        return 0;
    }


  for (k = df-dg ; k >=0 ; k--) {
    mpz_mul(q->coeff[k], r->coeff[k+dg], invlg);
    mpz_mod(q->coeff[k], q->coeff[k], p);
    for (j = dg+k ; j >= k ; j--) {
      mpz_submul(r->coeff[j], q->coeff[k], g->coeff[j-k]);
      mpz_mod(r->coeff[j], r->coeff[j], p);
    }
  }
  mpz_poly_cleandeg(r, r->deg);

  mpz_clear(invlg);
  mpz_clear(lg);
  return 1;
}

/* This also computes q and r such that f = q * g + r, but over Z, not
 * modulo a prime. Also, we do not assume that g is monic. Of course, if
 * it is not, then most often the result will be undefined (over Z). We
 * return something well-defined if q and r happen to be integer
 * polynomials.
 * We return 0 if this is not the case (in which case q and r are
 * undefined).
 * r==f is allowed.
 */
int mpz_poly_div_qr_z (mpz_poly_ptr q, mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g)
{
  int k, j, df = f->deg, dg = g->deg, dq = df - dg;

  if (df < dg) /* f is already reduced mod g */
    {
      mpz_poly_set_zero (q);
      mpz_poly_set (r, f);
      return 1;
    }

  /* now df >= dg */
  mpz_poly_realloc(q, dq + 1);

  mpz_poly_set(r, f);
  q->deg = dq;

  mpz_srcptr lg = g->coeff[dg];

  for (k = df-dg ; k >=0 ; k--) {
    if (!mpz_divisible_p(r->coeff[k+dg], lg))
        return 0;
    mpz_divexact(q->coeff[k], r->coeff[k+dg], lg);
    for (j = dg+k ; j >= k ; j--) {
      mpz_submul(r->coeff[j], q->coeff[k], g->coeff[j-k]);
    }
  }
  mpz_poly_cleandeg(r, r->deg);
  return 1;
}
int mpz_poly_div_r_z (mpz_poly_ptr r, mpz_poly_srcptr f, mpz_poly_srcptr g)
{
    mpz_poly quo;
    mpz_poly_init(quo,-1);
    int ret = mpz_poly_div_qr_z(quo,r,f,g);
    mpz_poly_clear(quo);
    return ret;
}


/* q=divexact(h, f) mod p, f not necessarily monic.
   Assumes lc(h) <> 0 mod p.
   Clobbers h. */
/* Coefficients of f must be reduced mod p
 * Coefficients of h need not be reduced mod p
 * Coefficients of q are reduced mod p
 */
static int
mpz_poly_divexact_clobber (mpz_poly_ptr q, mpz_poly_ptr h, mpz_poly_srcptr f,
                   mpz_srcptr p)
{
  int i, d = f->deg, dh = h->deg;
  mpz_t t, aux;

  mpz_init (t);
  mpz_init (aux);
  ASSERT (d >= 0);
  ASSERT (dh >= 0);
  ASSERT (dh >= d);
  ASSERT (mpz_divisible_p (h->coeff[dh], p) == 0);

  mpz_poly_realloc (q, dh + 1 - d);
  q->deg = dh - d;
  /* t is 1/f[d] mod p */
  mpz_set (aux, f->coeff[d]);
  mpz_mod (aux, aux, p);
  if (!mpz_invert (t, aux, p)) {
      mpz_clear(t);
      mpz_clear(aux);
      return 0;
  }


  while (dh >= d) {

    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(t, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], t);
      mpz_mod (h->coeff[dh], h->coeff[dh], p);
    }
    mpz_set (q->coeff[dh-d], h->coeff[dh]);
    mpz_mod (q->coeff[dh-d], q->coeff[dh - d], p);

    /* we only need to update the coefficients of degree >= d of h,
       i.e., we want i >= 2d - dh */
    for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
      mpz_mul (aux, h->coeff[dh], f->coeff[i]);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], aux);
      mpz_mod (h->coeff[dh - d + i], h->coeff[dh - d + i], p);
    }
    dh --;
  }
  /* since lc(h) <> 0 mod p, q is normalized */

  mpz_clear (t);
  mpz_clear (aux);
  return 1;
}

/* q <- divexact(h, f) mod p, f not necessarily monic. */
/* Coefficients of f must be reduced mod p
 * Coefficients of h need not be reduced mod p
 * Coefficients of q are reduced mod p
 */
int
mpz_poly_divexact (mpz_poly_ptr q, mpz_poly_srcptr h, mpz_poly_srcptr f, mpz_srcptr p)
{
    mpz_poly hh;
    mpz_poly_init(hh, h->deg);
    mpz_poly_set(hh, h);
    int r = mpz_poly_divexact_clobber(q, hh, f, p);
    mpz_poly_clear(hh);
    return r;
}

/* Set f=g/2 (mod m), where f might equal g.
   Assumes m is odd. */
/* If coefficients of g are reduced mod m, then coefficients of f are
 * reduced.
 */
void mpz_poly_div_2_mod_mpz (mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr m)
{
  int i;

  ASSERT_ALWAYS(mpz_scan1 (m, 0) == 0);

  mpz_poly_realloc (f, g->deg + 1);

  for (i = g->deg; i >= 0; --i)
    {
      if (mpz_scan1 (g->coeff[i], 0) == 0) /* g[i] is odd */
        {
          mpz_add (f->coeff[i], g->coeff[i], m);
          mpz_div_2exp (f->coeff[i], f->coeff[i], 1);
        }
      else
        mpz_div_2exp (f->coeff[i], g->coeff[i], 1);
    }
}

/* Set res=f(x). Assumes res and x are different variables. */
void mpz_poly_eval(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x) {
  int i, d;
  d = f->deg;
  if (d == -1) {
    mpz_set_ui(res, 0);
    return;
  }
  mpz_set(res, f->coeff[d]);
  for (i = d-1; i>=0; --i) {
    mpz_mul(res, res, x);
    mpz_add(res, res, f->coeff[i]);
  }
}

/* Set res=f(x) where x is an unsigned long. */
void mpz_poly_eval_ui (mpz_t res, mpz_poly_srcptr f, unsigned long x)
{
  int d = f->deg;

  mpz_set (res, f->coeff[d]);
  for (int i = d - 1; i >= 0; i--)
  {
    mpz_mul_ui (res, res, x);
    mpz_add (res, res, f->coeff[i]);
  }
}

/* Set res=f'(x), where x is an unsigned long */
void
mpz_poly_eval_diff_ui (mpz_t res, mpz_poly_srcptr f, unsigned long x)
{
  int d = f->deg;

  mpz_mul_ui (res, f->coeff[d], d);
  for (int i = d - 1; i >= 1; i--)
    {
      mpz_mul_ui (res, res, x);
      mpz_addmul_ui (res, f->coeff[i], i); /* res <- res + i*f[i] */
    }
}


/* Set res=f(x) (mod m).  Assume res and x are different variables. */
/* Coefficients of f(x) need not be reduced mod m.
 * The computed value res is reduced mod m
 */
void mpz_poly_eval_mod_mpz(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x,
                       mpz_srcptr m)
{
    mpz_poly_eval_mod_mpz_barrett(res, f, x, m, NULL);
}

/* Return 1 if poly(root) % modulus == 0, return 0 otherwise */
/* Coefficients of f(x) need not be reduced mod m */
int mpz_poly_is_root(mpz_poly_srcptr poly, mpz_t root, mpz_t modulus)
{
    mpz_t x;
    mpz_init(x);
    mpz_poly_eval_mod_mpz(x, poly, root, modulus);
    int is_root = (mpz_cmp_ui(x, 0) == 0);
    mpz_clear(x);
    return is_root;
}

/* Set res=f(x) (mod m).  Assume res and x are different variables. */
/* Coefficients of f(x) need not be reduced mod m.
 * The computed value res is reduced mod m
 */
void mpz_poly_eval_mod_mpz_barrett(mpz_t res, mpz_poly_srcptr f, mpz_srcptr x,
                       mpz_srcptr m, mpz_srcptr mx) {
  int i, d;

  d = f->deg;
  if (d == -1) {
    mpz_set_ui(res, 0);
    return;
  }
  barrett_mod(res, f->coeff[d], m, mx);
  for (i = d-1; i>=0; --i) {
    mpz_mul(res, res, x);
    mpz_add(res, res, f->coeff[i]);
    barrett_mod(res, res, m, mx);
  }
}

/* This evaluates several polynomials at the same point w. It is possible
 * to do fewer modular reductions in this case.
 *
 * When k==1, we use mpz_poly_eval_mod_mpz instead, since it's faster
 * then.
 */
/* Coefficients of f[i](x) need not be reduced mod m.
 * The computed values r[i] are reduced mod m
 */
void mpz_poly_eval_several_mod_mpz_barrett(mpz_ptr * r, mpz_poly_srcptr * f, int k, mpz_srcptr x,
                       mpz_srcptr m, mpz_srcptr mx)
{
    int i;

    if (k == 1) {
        mpz_poly_eval_mod_mpz_barrett(r[0], f[0], x, m, mx);
        return;
    }

    mpz_t w;
    mpz_init(w);
    int maxdeg = -1;
    for(int j = 0 ; j < k ; j++) {
        if (f[j]->deg >= 0)
            mpz_set(r[j], f[j]->coeff[0]);
        else
            mpz_set_ui(r[j], 0);
        if (f[j]->deg > maxdeg)
            maxdeg = f[j]->deg;
    }

    mpz_set(w, x);
    for(int j = 0 ; j < k ; j++) {
        if (f[j]->deg >= 1)
            mpz_addmul(r[j], w, f[j]->coeff[1]);
    }
    for(i = 2 ; i <= maxdeg ; i++) {
        mpz_mul(w, w, x);
        barrett_mod(w, w, m, mx);
        for(int j = 0 ; j < k ; j++)
            if (f[j]->deg >= i)
                mpz_addmul(r[j], w, f[j]->coeff[i]);
    }
    for(int j = 0 ; j < k ; j++) {
        barrett_mod(r[j], r[j], m, mx);
    }
    mpz_clear(w);
}

/* Set Q=P1*P2 (mod F). Warning: Q might equal P1 (or P2). */
void polymodF_mul (polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
                   mpz_poly_srcptr F) {
  mpz_poly prd;
  int v;

  /* beware: if P1 and P2 are zero, P1->p->deg + P2->p->deg = -2 */
  mpz_poly_init (prd, (P1->p->deg == -1) ? -1 : P1->p->deg + P2->p->deg);

  ASSERT_ALWAYS(mpz_poly_normalized_p (P1->p));
  ASSERT_ALWAYS(mpz_poly_normalized_p (P2->p));

  mpz_poly_mul (prd, P1->p, P2->p);
  v = P1->v + P2->v;

  mpz_poly_reducemodF (Q, prd, F);
  Q->v += v;
  mpz_poly_clear (prd);
}

/* Set Q = P/lc(P) (mod m). Q and P might be identical. */
/* Coefficients of P need not be reduced mod m
 * Coefficients of Q are reduced mod m
 */
void
mpz_poly_makemonic_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_srcptr m)
{
  int i;
  mpz_t aux;
  mpz_init(aux);
  for(i = P->deg; i>=0; i--) {
      mpz_mod(aux, P->coeff[i], m);
      if (mpz_cmp_ui(aux, 0) != 0) break;
  }
  /* i is the degree of the leading monomial */
  Q->deg = i;
  if (i < 0) {
      /* if i == -1, then Q is the zero polynomial, there's nothing to do */
      mpz_clear(aux);
      return;
  }

  mpz_t aux2;
  mpz_init(aux2);
  mpz_invert(aux2, aux, m);
  for (i = 0; i < Q->deg; ++i) {
      mpz_mul(aux, aux2, P->coeff[i]);
      mpz_mod(aux, aux, m);
      mpz_poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux2);

  /* we can directly set the leading coefficient to 1 */
  mpz_poly_setcoeff_si (Q, Q->deg, 1);
  mpz_clear(aux);
}

/* Coefficients of A need not be reduced mod m
 * Coefficients of R are reduced mod m
 * invm may be NULL (or computed by barrett_init)
 */
int
mpz_poly_mod_mpz (mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m, mpz_srcptr invm)
{
  /* reduce lower coefficients */
  mpz_poly_realloc(R, A->deg + 1);
  for (int i = 0; i <= A->deg; ++i)
    barrett_mod (R->coeff[i], A->coeff[i], m, invm);

  mpz_poly_cleandeg(R, A->deg);
  return R->deg;
}

/* reduce non-negative coefficients in [0, p-1], negative ones in [1-p, -1] */
int
mpz_poly_mod_mpz_lazy (mpz_poly_ptr R, mpz_poly_srcptr A, mpz_srcptr m,
		       mpz_srcptr invm)
{
  mpz_poly_realloc(R, A->deg + 1);
  for (int i = 0; i <= A->deg; ++i)
    {
      if (mpz_sgn (A->coeff[i]) >= 0)
	barrett_mod (R->coeff[i], A->coeff[i], m, invm);
      else
	{
	  mpz_neg (R->coeff[i], A->coeff[i]);
	  barrett_mod (R->coeff[i], R->coeff[i], m, invm);
	  mpz_neg (R->coeff[i], R->coeff[i]);
	}
    }

  mpz_poly_cleandeg(R, A->deg);
  return R->deg;
}

/* Reduce R[d]*x^d + ... + R[0] mod f[df]*x^df + ... + f[0] modulo m.
   Return the degree of the remainder.
   Assume invm = floor(B^(2k)/m), m having k limbs, and B is the limb base.
   Coefficients of f must be reduced mod m on input.
   Coefficients of R need not be reduced mod m on input, but are reduced
   on output.
   invm may be NULL (or computed by barrett_init).
   If invf is not NULL, it should be 1/m mod lc(f). */
int
mpz_poly_mod_f_mod_mpz (mpz_poly_ptr R, mpz_poly_srcptr f, mpz_srcptr m,
			mpz_srcptr invm, mpz_srcptr invf)
{
  mpz_t aux, c;

  if (f == NULL)
    goto reduce_R;

  if (invf == NULL)
    {
      mpz_init (aux);
      /* aux = 1/m mod lc(f) */
      mpz_invert (aux, m, f->coeff[f->deg]);
    }

  mpz_init (c);
  // FIXME: write a subquadratic variant
  while (R->deg >= f->deg)
    {
      /* Here m is large (thousand to million bits) and lc(f) is small
       * (typically one word). We first add to R lambda * m * x^(dR-df) ---
       * which is zero mod m --- such that the new coefficient of degree dR
       * is divisible by lc(f), i.e., lambda = -lc(R)/m mod lc(f). Then if
       * c = (lc(R) + lambda * m) / lc(f), we subtract c * x^(dR-df) * f. */
      mpz_mod (c, R->coeff[R->deg], f->coeff[f->deg]); /* lc(R) mod lc(f) */
      mpz_mul (c, c, (invf == NULL) ? aux : invf);
      mpz_mod (c, c, f->coeff[f->deg]);    /* lc(R)/m mod lc(f) */
      mpz_submul (R->coeff[R->deg], m, c);  /* lc(R) - m * (lc(R) / m mod lc(f)) */
      ASSERT (mpz_divisible_p (R->coeff[R->deg], f->coeff[f->deg]));
      mpz_divexact (c, R->coeff[R->deg], f->coeff[f->deg]);
      /* If R[deg] has initially size 2n, and f[deg] = O(1), then c has size
	 2n here. */
      for (int i = R->deg - 1; i >= R->deg - f->deg; --i)
	mpz_submul (R->coeff[i], c, f->coeff[f->deg-R->deg+i]);
      R->deg--;
    }
  
  mpz_clear (c);
  if (invf == NULL)
    mpz_clear (aux);

 reduce_R:
  mpz_poly_mod_mpz (R, R, m, invm);

  return R->deg;
}

/* Stores in g the polynomial linked to f by:
 *  - alpha is a root of f if and only if lc(f)*alpha is a root of g
 *  - g is monic
 */
void mpz_poly_to_monic(mpz_poly_ptr g, mpz_poly_srcptr f)
{
    mpz_t fd,temp;
    mpz_init(fd);
    mpz_init(temp);
    mpz_poly_getcoeff(fd,f->deg,f);

    mpz_poly_set(g,f);
    for (int k = 0 ; k < g->deg ; k++) {
        mpz_poly_getcoeff(temp,k,g);
        for(int j = 1 ; j <= g->deg-1-k ; j++){
            mpz_mul(temp,temp,fd);
        }
        mpz_poly_setcoeff(g,k,temp);
    }
    mpz_poly_setcoeff_ui(g,g->deg,1);

    mpz_clear(temp);
    mpz_clear(fd);
}

/*  Reduce frac (= num / denom) mod F mod m ,
    i.e. compute num * denom^-1 mod F mod m .
    The return value is in num, denom is set to constant polynomial 1
    */
/* TODO: We must state the input / output requirements with respect to
 * reduction mod m */
void
mpz_poly_reduce_frac_mod_f_mod_mpz (mpz_poly_ptr num, mpz_poly_ptr denom,
                                    mpz_poly_srcptr F, mpz_srcptr m, mpz_srcptr invm)
{
  if (denom->deg == 0)
  {
    mpz_t inv;
    mpz_init (inv);
    mpz_poly_getcoeff (inv, 0, denom); /* inv <- denom[0] */
    mpz_invert (inv, inv, m);          /* inv <- denom[0]^-1 */
    mpz_poly_mul_mpz (num, num, inv);  /* num <- num * inv */
    mpz_poly_mod_mpz (num, num, m, invm); /* num <- num * inv mod m */
    mpz_clear (inv);
  } else {
    mpz_poly g, U, V;
    mpz_poly_init (g, 0);
    mpz_poly_init (U, 0);
    mpz_poly_init (V, 0);
    mpz_poly_xgcd_mpz (g, F, denom, U, V, m);
    mpz_poly_mul (num, num, V);
    mpz_poly_mod_f_mod_mpz (num, F, m, invm, NULL);
    mpz_poly_clear (g);
    mpz_poly_clear (U);
    mpz_poly_clear (V);
  }
  mpz_poly_set_zero (denom);
  mpz_poly_setcoeff_si (denom, 0, 1);
}


/* Q = P1*P2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Assume invm = floor(B^(2k)/m), m having k limbs, and B is the limb base.
   Coefficients of P1 and P2 need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If invf is not NULL, it is 1/m mod lc(f). */
void
mpz_poly_mul_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2,
			    mpz_poly_srcptr f, mpz_srcptr m, mpz_srcptr invm,
			    mpz_srcptr invf)
{
  int d1 = P1->deg;
  int d2 = P2->deg;
  int d = d1+d2;
  mpz_poly R;

  mpz_poly_init(R, d);

  d = mpz_poly_mul_tc (R->coeff, P1->coeff, d1, P2->coeff, d2);
  mpz_poly_cleandeg(R, d);

  // reduce mod f
  mpz_poly_mod_f_mod_mpz (R, f, m, invm, invf);

  mpz_poly_set(Q, R);
  mpz_poly_clear(R);
}

/* Q = P1*P2 mod f, assuming f is monic */
void
mpz_poly_mul_mod_f (mpz_poly_ptr Q, mpz_poly_srcptr P1, mpz_poly_srcptr P2,
                        mpz_poly_srcptr f)
{
    mpz_poly_mul(Q,P1,P2);
    mpz_poly_div_r_z(Q,Q,f);
}

/* Q = P^2 mod f, mod m
   f is the original algebraic polynomial (non monic but small coefficients)
   Assume invm = floor(B^(2k)/m), m having k limbs, and B is the limb base.
   Coefficients of P need not be reduced mod m.
   Coefficients of Q are reduced mod m.
   If not NULL, invf = 1/m mod lc(f). */
void
mpz_poly_sqr_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f,
			    mpz_srcptr m, mpz_srcptr invm, mpz_srcptr invf)
{
  int d1 = P->deg;
  int d = d1 + d1;
  mpz_poly R;

  mpz_poly_init(R, d);

  /* Fast squaring in 2d1+1 squares, i.e., 2d-1 squares.
     For d=5, this gives 9 squares. */
  d = mpz_poly_sqr_tc (R->coeff, P->coeff, d1);
  mpz_poly_cleandeg(R, d);

  // reduce mod f
  mpz_poly_mod_f_mod_mpz (R, f, m, invm, invf);

  mpz_poly_set(Q, R);
  mpz_poly_clear(R);
}

/* Affects the derivative of f to df. Assumes df different from f.
   Assumes df has been initialized with degree at least f->deg-1. */
void mpz_poly_derivative(mpz_poly_ptr df, mpz_poly_srcptr f) {
  int n;

  df->deg = (f->deg <= 0) ? -1 : f->deg - 1;
  for (n = 0; n <= f->deg - 1; n++)
    mpz_mul_si (df->coeff[n], f->coeff[n + 1], n + 1);
}

/* B = A^n mod f, assuming f is monic */
void mpz_poly_pow_ui_mod_f(mpz_poly_ptr B, mpz_poly_srcptr A, unsigned long n, mpz_poly_srcptr f)/*{{{*/
{
    if (n == 0) {
        mpz_poly_set_xi(B, 0);
        return;
    }
    if (B == A || B  == f) {
        mpz_poly C;
        mpz_poly_init(C, f->deg-1);
        mpz_poly_pow_ui_mod_f(C, A, n, f);
        mpz_poly_swap(B, C);
        mpz_poly_clear(C);
        return;
    }
    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > n ; k >>= 1);
    mpz_poly_set(B, A);
    mpz_poly Q;
    mpz_poly_init(Q, 3*f->deg-3);
    for( ; k >>= 1 ; ) {
        mpz_poly_mul(Q, B, B);
        mpz_poly_swap(Q, B);
        if (n & k) {
            mpz_poly_mul(Q, B, A);
            mpz_poly_swap(Q, B);
        }
        mpz_poly_div_r_z(B, B, f);
    }
    mpz_poly_clear(Q);
}/*}}}*/

/* Q = P^a mod f, mod p (f is the algebraic polynomial, non monic) */
/* Coefficients of f must be reduced mod m on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void
mpz_poly_pow_mod_f_mod_ui (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f,
                         mpz_srcptr a, unsigned long p)
{
    mpz_t m;
    mpz_init_set_ui(m, p);
    mpz_poly_pow_mod_f_mod_mpz(Q, P, f, a, m);
    mpz_clear(m);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t */
/* f may be NULL, in case there is only reduction mod p */
/* Coefficients of f must be reduced mod m on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void
mpz_poly_pow_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f,
                          mpz_srcptr a, mpz_srcptr p)
{
  mpz_t invp;

  mpz_init (invp);
  barrett_init (invp, p);
  mpz_poly_pow_mod_f_mod_mpz_barrett(Q, P, f, a, p, invp);
  mpz_clear (invp);
}

/* Coefficients of f must be reduced mod m on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p.
 */
void
mpz_poly_pow_ui_mod_f_mod_mpz (mpz_poly_ptr Q, mpz_poly_srcptr P, mpz_poly_srcptr f,
                          unsigned long a, mpz_srcptr p)
{
    mpz_t az;
    mpz_init_set_ui(az, a);
    mpz_poly_pow_mod_f_mod_mpz(Q, P, f, az, p);
    mpz_clear(az);
}

/* store in invm the value of floor(B^(2k)/m), where m has k limbs,
   and B is the limb base */
void barrett_init (mpz_ptr invm, mpz_srcptr m)
{
  size_t k = mpz_size (m);

  mpz_set_ui (invm, 1);
  mpz_mul_2exp (invm, invm, 2 * k * GMP_NUMB_BITS);
  mpz_tdiv_q (invm, invm, m);
}

/* a <- b mod m
   Assume invm = floor(B^(2k)/m), m having k limbs, and B is the limb base */
void barrett_mod (mpz_ptr a, mpz_srcptr b, mpz_srcptr m, mpz_srcptr invm)
{
  size_t k = mpz_size (m), sizeb, l;
  mpz_t c;

  sizeb = mpz_size (b);
  if (invm == NULL || sizeb <= k || k == 1)
  {
    mpz_mod (a, b, m);
    return;
  }

  /* now sizeb > k */
  l = sizeb - k;
  if (l > k)
    l = k;
  /* l > 0 is the number of limbs to reduce */
  mpz_init (c);
  mpz_div_2exp (c, b, (sizeb - l) * GMP_NUMB_BITS); /* c has l limbs */
  mpz_mul (c, c, invm);                             /* c has k+l limbs */
  mpz_div_2exp (c, c, k * GMP_NUMB_BITS);           /* c has l limbs */
  mpz_mul (c, c, m);                                /* c has k+l limbs */
  mpz_mul_2exp (c, c, (sizeb - k - l) * GMP_NUMB_BITS);
  mpz_sub (a, b, c);
  while (mpz_size (a) > k)
  {
    sizeb = mpz_size (a);
    l = sizeb - k;
    if (l > k)
      l = k;
    mpz_div_2exp (c, a, (sizeb - l) * GMP_NUMB_BITS);
    mpz_mul (c, c, invm);
    mpz_div_2exp (c, c, k * GMP_NUMB_BITS);
    mpz_mul (c, c, m);
    mpz_mul_2exp (c, c, (sizeb - k - l) * GMP_NUMB_BITS);
    mpz_sub (a, a, c);
  }
  mpz_clear (c);
  mpz_mod (a, a, m);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t */
/* Same as mpz_poly_pow_mod_f_mod_mpz but use barrett for reduction mod p */
/* For f == NULL, there is no reduction mod f */
/* Coefficients of f must be reduced mod p on input
 * Coefficients of P need not be reduced mod p.
 * Coefficients of Q are reduced mod p
 * invp may be NULL, of must have been computed by barrett_init
 */
void
mpz_poly_pow_mod_f_mod_mpz_barrett (mpz_poly_ptr Q, mpz_poly_srcptr P,
				    mpz_poly_srcptr f, mpz_srcptr a,
				    mpz_srcptr p, mpz_srcptr invp)
{
  int k = mpz_sizeinbase(a, 2), l, L = 0, j;
  mpz_poly R, *T = NULL;
  mpz_t invf;

  if (mpz_cmp_ui(a, 0) == 0) {
    mpz_poly_set_xi (Q, 0);
    return;
  }

  mpz_poly_init(R, f ? 2*f->deg : -1);

  // Initialize R to P
  mpz_poly_set(R, P);

  /* compute invf = 1/p mod lc(f) */
  if (f != NULL)
    {
      mpz_init (invf);
      mpz_invert (invf, p, f->coeff[f->deg]);
    }

  /* We use base-2^l exponentiation with sliding window,
     thus we need to precompute P^2, P^3, ..., P^(2^l-1).
     The expected average cost k squarings plus k/(l+1) + 2^(l-1) multiplies.
  */

  for (l = 1; k / (l + 1) + (1 << (l - 1)) > k / (l + 2) + (1 << l); l++);
  /* this gives (roughly) l=1 for k < 8, l=2 for k < 27, l=3 for k < 84,
     l=4 for k < 245 */
  if (l > 1)
    {
      L = 1 << (l - 1);
      T = (mpz_poly*) malloc (L * sizeof (mpz_poly));
      /* we store P^2 in T[0], P^3 in T[1], ..., P^(2^l-1) in T[L-1] */
      for (j = 0; j < L; j++)
	mpz_poly_init (T[j], f ? 2*f->deg : -1);
      mpz_poly_sqr_mod_f_mod_mpz (T[0], R, f, p, invp, invf);       /* P^2 */
      mpz_poly_mul_mod_f_mod_mpz (T[1], T[0], R, f, p, invp, invf); /* P^3 */
      for (j = 2; j < L; j++)
	mpz_poly_mul_mod_f_mod_mpz (T[j], T[j-1], T[0], f, p, invp, invf);
    }

  // Horner
  for (k -= 2; k >= 0;)
  {
    while (k >= 0 && mpz_tstbit (a, k) == 0)
      {
	mpz_poly_sqr_mod_f_mod_mpz (R, R, f, p, invp, invf);
	k --;
      }
    if (k < 0)
      break;
    j = mpz_scan1 (a, (k >= l) ? k - (l - 1) : 0);
    /* new window starts at bit k, and ends at bit j <= k */
    int e = 0;
    while (k >= j)
      {
	mpz_poly_sqr_mod_f_mod_mpz (R, R, f, p, invp, invf);
	e = 2 * e + mpz_tstbit (a, k);
	k --;
      }
    mpz_poly_mul_mod_f_mod_mpz (R, R, (e == 1) ? P : T[e/2], f, p, invp, invf);
  }

  mpz_poly_swap(Q, R);
  mpz_poly_clear(R);
  if (f != NULL)
    mpz_clear (invf);
  if (l > 1)
    {
      for (k = 0; k < L; k++)
	mpz_poly_clear (T[k]);
      free (T);
    }
}

/* Return a list of polynomials P[0], P[1], ..., P[l] such that
   P0 = Q[l-1] + p^K[1]*P[l]
   Q[l-1] = Q[l-2] + p^K[2]*P[l-1]
   ...
   Q[2] = Q[1] + p^K[l-1]*P[2]
   Q[1] = Q[0] + p*P[1]         where K[l] = 1
   Q[0] = P[0]
   ...
   With all coefficients of P[i] smaller than p^(K[l-i]-K[l-(i-1)]).
   Assume K[l]=1, K[l-1]=2, ..., K[i] = 2*K[i+1] or 2*K[i+1]-1.
   P0 = P[0] + p*P[1] + p^2*P[2] + ... + p^K[1]*P[l] < p^K[0]
   The end of the list is P[l+1]=0.
   Assume l > 0.
*/
mpz_poly*
mpz_poly_base_modp_init (mpz_poly_srcptr P0, int p, int *K, int l)
{
  mpz_poly *P;
  int k, i, j;
  mpz_t *pk;

  ASSERT_ALWAYS (l > 0);
  ASSERT_ALWAYS (K[l] == 1);

  /* initialize pk[i] = p^K[l-i] for 0 <= i < l */
  pk = (mpz_t*) malloc (l * sizeof (mpz_t));
  FATAL_ERROR_CHECK (pk == NULL, "not enough memory");
  mpz_init_set_ui (pk[0], p);
  for (i = 1; i < l; i++)
  {
    mpz_init (pk[i]);
    mpz_mul (pk[i], pk[i-1], pk[i-1]);
    if (K[l-i] & 1)
    {
      ASSERT_ALWAYS (K[l-i] == 2 * K[l-i+1] - 1);
      mpz_div_ui (pk[i], pk[i], p);
    }
    else
      ASSERT_ALWAYS (K[l-i] == 2 * K[l-i+1]);
  }

  /* now decompose P0: we need P[0], P[1] for factor p, P[2] for p^2,
     ..., P[l] for p^K[1], and one for the end of list,
     thus l+2 polynomials */
  P = (mpz_poly*) malloc ((l + 2) * sizeof(mpz_poly));
  FATAL_ERROR_CHECK (P == NULL, "not enough memory");
  for (i = 0; i < l + 2; i++)
    mpz_poly_init (P[i], P0->deg);
  /* P[l+1] is initialized to 0 by mpz_poly_init */

  /* initialize P[k], and put remainder in P[k-1] */
  for (k = l, i = 0; i <= P0->deg; i++)
    mpz_tdiv_qr (P[k]->coeff[i], P[k-1]->coeff[i], P0->coeff[i], pk[k-1]);
  mpz_poly_cleandeg (P[k], P0->deg);

  /* now go down */
  for (j = l-1; j >= 1; j--)
  {
    /* reduce P[j] into P[j] and P[j-1] */
    for (i = 0; i <= P0->deg; i++)
      mpz_tdiv_qr (P[j]->coeff[i], P[j-1]->coeff[i], P[j]->coeff[i],
                   pk[j-1]);
    mpz_poly_cleandeg (P[j], P0->deg);
  }
  mpz_poly_cleandeg (P[0], P0->deg);

  for (i = 0; i < l; i++)
    mpz_clear (pk[i]);
  free (pk);

  return P;
}

/* a <- a + pk*P[k] */
void
mpz_poly_base_modp_lift (mpz_poly_ptr a, mpz_poly *P, int k, mpz_srcptr pk)
{
  int i;

  /* first check P[k] exists and is not zero */
  if (P[k]->deg == -1)
    return;

  mpz_poly_realloc (a, P[k]->deg + 1);

  for (i = 0; i <= P[k]->deg && i <= a->deg; i++)
    mpz_addmul (a->coeff[i], P[k]->coeff[i], pk);

  for (; i <= P[k]->deg; i++)
    mpz_mul (a->coeff[i], P[k]->coeff[i], pk);

  mpz_poly_cleandeg (a, (a->deg >= P[k]->deg) ? a->deg : P[k]->deg);
}

void
mpz_poly_base_modp_clear (mpz_poly *P, int l)
{
  for (int i = 0; i < l + 2; i++)
    mpz_poly_clear (P[i]);
  free (P);
}

/* return the maximal size of the coefficients of f in base b */
size_t
mpz_poly_sizeinbase (mpz_poly_srcptr f, int b)
{
  size_t S = 0, s;
  int i;
  int d = f->deg;

  for (i = 0; i <= d; i++)
  {
    s = mpz_sizeinbase (f->coeff[i], b);
    if (s > S)
      S = s;
  }
  return S;
}

void
mpz_poly_infinity_norm (mpz_ptr in, mpz_poly_srcptr f)
{
  if (f->deg == -1) {
    mpz_set_ui(in, 0);
  } else {
  mpz_abs (in, f->coeff[0]);
  for (int i = 1; i <= f->deg; i++)
    {
      if (mpz_cmpabs (f->coeff[i], in) > 0)
	mpz_abs (in, f->coeff[i]);
    }
  }
}

/* return the total size (in bytes) to store the polynomial f */
size_t
mpz_poly_totalsize (mpz_poly_srcptr f)
{
  int i;
  size_t s = 0;

  for (i = 0; i <= f->deg; i++)
    s += mpz_size (f->coeff[i]);
  return s * sizeof (mp_limb_t);
}

/* f=gcd(f, g) mod p, with p in mpz_t */
/* clobbers g */
/* Coefficients of f and g need not be reduced mod p on input.
 * Coefficients of f are reduced mod p on output */
static void
mpz_poly_gcd_mpz_clobber (mpz_poly_ptr f, mpz_poly_ptr g, mpz_srcptr p)
{
    /* First reduce mod p */
    mpz_poly_mod_mpz(f, f, p, NULL);
    mpz_poly_mod_mpz(g, g, p, NULL);
  while (g->deg >= 0)
    {
      mpz_poly_div_r (f, g, p);
      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (f, g);
    }
}

/* f <- gcd(a, b) mod p. */
/* Coefficients of a and b need not be reduced mod p
 * Coefficients of f are reduced mod p */
void mpz_poly_gcd_mpz(mpz_poly_ptr f, mpz_poly_srcptr a, mpz_poly_srcptr b, mpz_srcptr p)
{
    mpz_poly hh;
    if (f == b) {
        mpz_poly_init(hh, a->deg);
        mpz_poly_set(hh, a);
    } else {
        mpz_poly_init(hh, b->deg);
        mpz_poly_set(hh, b);
        mpz_poly_set(f, a); /* will do nothing if f = a */
    }
    mpz_poly_gcd_mpz_clobber(f, hh, p);
    mpz_poly_clear(hh);
}


/* Attempt to compute the f=gcd(f, g) mod N, where N is not necessarily a
 * prime. If at some point a division fails, this gives a proper factor
 * of N that is put in the corresponding argument.
 * The return value tells whether the process was successful (1 means
 * that no inversion failed, 0 means that a factor was found).
 * WARNING: this function destroys its input.
 */
/* Coefficients of f and g need not be reduced mod p on input.
 * Coefficients of f are reduced mod p on output */
int
mpz_poly_pseudogcd_mpz(mpz_poly_ptr f, mpz_poly_ptr g, mpz_srcptr N, mpz_t factor)
{
    mpz_poly_mod_mpz(f, f, N, NULL);
    mpz_poly_mod_mpz(g, g, N, NULL);
  while (g->deg >= 0)
  {
    int ret = mpz_poly_pseudodiv_r(f, g, N, factor);
    if (!ret)
        return ret;
    /* now deg(f) < deg(g): swap f and g */
    mpz_poly_swap (f, g);
  }
  // success: all inversions mod N worked.
  return 1;
}


/* computes d = gcd(f, g) = u*f + v*g mod p, with p in mpz_t */
/* Coefficients of f and g need not be reduced mod p.
 * Coefficients of d, u, v are reduced mod p */
void
mpz_poly_xgcd_mpz (mpz_poly_ptr d, mpz_poly_srcptr f, mpz_poly_srcptr g, mpz_poly_ptr u, mpz_poly_ptr v, mpz_srcptr p)
{
  mpz_poly q, tmp;
  mpz_poly gg;

  if (f->deg < g->deg) {
      mpz_poly_xgcd_mpz(d, g, f, v, u, p);
      return;
  }
  mpz_poly_init(gg, g->alloc);
  mpz_poly_set(d, f);
  mpz_poly_set(gg, g);
  mpz_poly_mod_mpz(d, d, p, NULL);
  mpz_poly_mod_mpz(gg, gg, p, NULL);

  mpz_poly uu, vv;
  mpz_poly_init (uu, 0);
  mpz_poly_init (vv, 0);

  mpz_poly_set_xi(u, 0);
  mpz_poly_set_zero(uu);

  mpz_poly_set_xi(vv, 0);
  mpz_poly_set_zero(v);

  mpz_poly_init(q, d->deg);
  mpz_poly_init(tmp, d->deg + gg->deg);

  while (gg->deg >= 0)
    {

      /* q, r := f div g mod p */
      mpz_poly_div_qr (q, d, d, gg, p);

      /* u := u - q * uu mod p */
      mpz_poly_mul(tmp, q, uu);
      mpz_poly_sub_mod_mpz(u, u, tmp, p);
      mpz_poly_swap (u, uu);

      /* v := v - q * vv mod p */
      mpz_poly_mul(tmp, q, vv);
      mpz_poly_sub_mod_mpz(v, v, tmp, p);
      mpz_poly_swap (v, vv);

      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (d, gg);
    }

  /* make monic */
  mpz_t inv;
  if (mpz_cmp_ui(d->coeff[d->deg], 1) != 0)
    {
      mpz_init(inv);
      mpz_invert(inv, d->coeff[0], p);
      mpz_poly_mul_mpz(d, d, inv);
      mpz_poly_mod_mpz(d, d, p, NULL);
      mpz_poly_mul_mpz(u, u, inv);
      mpz_poly_mod_mpz(u, u, p, NULL);
      mpz_poly_mul_mpz(v, v, inv);
      mpz_poly_mod_mpz(v, v, p, NULL);
      mpz_clear(inv);
    }

  mpz_poly_clear(gg);
  mpz_poly_clear(uu);
  mpz_poly_clear(vv);
  mpz_poly_clear(q);
  mpz_poly_clear(tmp);
}

/*  Homographic transform on polynomials */
/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
*/
void
mpz_poly_homography (mpz_poly_ptr Fij, mpz_poly_srcptr F, int64_t H[4])
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;
  mpz_t *f = F->coeff;
  int d = F->deg;

  mpz_poly_realloc (Fij, d + 1);

  mpz_t *fij = Fij->coeff;
  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  Fij->deg = d;

  g = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  FATAL_ERROR_CHECK (g == NULL, "not enough memory");
  for (k = 0; k <= d; k++)
    mpz_init (g[k]);
  mpz_init (f0);

  /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
     F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
     We use that formula recursively. */

  mpz_set_ui (g[0], 1); /* g = 1 */

  for (k = d - 1; k >= 0; k--)
    {
      /* invariant: we have already translated coefficients of degree > k,
         in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
         g[0..d - (k+1)]:
         f[k] <- a1*f[k+1]
         ...
         f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
         ...
         f[d] <- a0*f[d] */
      mpz_swap (f0, fij[k]); /* save the new constant coefficient */
      mpz_mul_si (fij[k], fij[k + 1], H[2]);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], H[0]);
          mpz_addmul_si (fij[l], fij[l + 1], H[2]);
        }
      mpz_mul_si (fij[d], fij[d], H[0]);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], H[1]);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], H[3]);
          mpz_addmul_si (g[l], g[l-1], H[1]);
        }
      mpz_mul_si (g[0], g[0], H[3]);

      /* now g has degree d-k, and we add f0*g */
      for (l = 0; l <= d-k; l++)
          mpz_addmul(fij[l+k], g[l], f0);

    }

  mpz_clear (f0);
  for (k = 0; k <= d; k++)
    mpz_clear (g[k]);
  free (g);
}

/* v <- |f(i,j)|, where f is homogeneous of degree d */
void mpz_poly_homogeneous_eval_siui (mpz_t v, mpz_poly_srcptr f, const int64_t i, const uint64_t j)
{
  unsigned int k = f->deg;
  mpz_t jpow;

  ASSERT(k > 0);
  mpz_set (v, f->coeff[f->deg]);
  mpz_mul_si (v, f->coeff[k], i);
  mpz_init (jpow);
  mpz_set_uint64 (jpow, j);
  mpz_addmul (v, f->coeff[--k], jpow); /* v = i*f[d] + j*f[d-1] */
  for (; k-- > 0;)
      {
        /* this test will be resolved at compile time by most compilers */
        if ((uint64_t) ULONG_MAX >= UINT64_MAX)
          { /* hardcode since this function is critical in las */
            mpz_mul_si (v, v, i);
            mpz_mul_ui (jpow, jpow, j);
          }
        else
          {
            mpz_mul_int64 (v, v, i);
            mpz_mul_uint64 (jpow, jpow, j);
          }
        mpz_addmul (v, f->coeff[k], jpow);
      }
  mpz_abs (v, v); /* avoids problems with negative norms */
  mpz_clear (jpow);
}

/* put in c the content of f */
void
mpz_poly_content (mpz_t c, mpz_poly_srcptr F)
{
  int i;
  mpz_t *f = F->coeff;
  int d = F->deg;

  mpz_set (c, f[0]);
  for (i = 1; i <= d; i++)
    mpz_gcd (c, c, f[i]);
  mpz_abs (c, c);
}

/*
 * Compute the pseudo division of a and b such that
 *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that deg(a) >= deg(b) and b is not the zero polynomial.
 */
static void mpz_poly_pseudo_division(mpz_poly_ptr q, mpz_poly_ptr r,
    mpz_poly_srcptr a, mpz_poly_srcptr b)
{
  ASSERT(a->deg >= b->deg);
  ASSERT(b->deg != -1);

  int m = a->deg;
  int n = b->deg;
  mpz_t d;
  int e;
  mpz_poly s;

#ifndef NDEBUG
  MAYBE_UNUSED mpz_poly q_tmp;
#endif // NDEBUG

  mpz_init(d);
  mpz_set(d, mpz_poly_lc_const(b));

  if (q != NULL) {
    mpz_poly_set_zero(q);
  }

#ifndef NDEBUG
  else {
    mpz_poly_init(q_tmp, 0);
  }
#endif // NDEBUG

  mpz_poly_set(r, a);

  e = m - n + 1;

  while (r->deg >= n) {
    mpz_poly_init(s, r->deg - n);
    mpz_poly_setcoeff(s, r->deg - n, mpz_poly_lc(r));
    if (q != NULL) {
      mpz_poly_mul_mpz(q, q, d);
      mpz_poly_add(q, s, q);
    }

#ifndef NDEBUG
    else {
      mpz_poly_mul_mpz(q_tmp, q_tmp, d);
      mpz_poly_add(q_tmp, s, q_tmp);
    }
#endif // NDEBUG

    mpz_poly_mul_mpz(r, r, d);
    mpz_poly_mul(s, b, s);
    mpz_poly_sub(r, r, s);
    mpz_poly_clear(s);
    e--;
  }

  ASSERT(e >= 0);

  mpz_pow_ui(d, d, (unsigned long int) e);
  if (q != NULL) {
    mpz_poly_mul_mpz(q, q, d);
  }

#ifndef NDEBUG
  else {
    mpz_poly_mul_mpz(q_tmp, q_tmp, d);
  }
#endif // NDEBUG

  mpz_poly_mul_mpz(r, r, d);

#ifndef NDEBUG
  mpz_poly f, g;
  mpz_poly_init(f, a->deg);
  mpz_poly_init(g, b->deg);
  mpz_poly_set(f, a);
  mpz_poly_set(g, b);

  mpz_set(d, mpz_poly_lc_const(g));

  ASSERT(m - n + 1 >= 0);

  mpz_pow_ui(d, d, (unsigned long int) (m - n + 1));
  mpz_poly_mul_mpz(f, f, d);

  if (q != NULL) {
    mpz_poly_mul(g, g, q);
  } else {
    mpz_poly_mul(g, g, q_tmp);
  }
  mpz_poly_add(g, g, r);

  ASSERT(mpz_poly_cmp(f, g) == 0);

  mpz_poly_clear(f);
  mpz_poly_clear(g);
  if (q == NULL) {
    mpz_poly_clear(q_tmp);
  }
#endif // NDEBUG

  mpz_clear(d);
}

/*
 * Like mpz_poly_pseudo_division, but give only the remainder.
 */
static void mpz_poly_pseudo_remainder(mpz_poly_ptr r, mpz_poly_srcptr a,
    mpz_poly_srcptr b)
{
  mpz_poly_pseudo_division(NULL, r, a, b);
}

/*
 * Compute the resultant of p and q and set the resultat in res.
 *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
 *  for more information.
 *
 * Assume that the polynomials are normalized.
 */
void mpz_poly_resultant(mpz_ptr res, mpz_poly_srcptr p, mpz_poly_srcptr q)
{
  if (p->deg == -1 || q->deg == -1) {
    mpz_set_ui(res, 0);
    return;
  }

  ASSERT(mpz_cmp_ui(p->coeff[p->deg], 0) != 0);
  ASSERT(mpz_cmp_ui(q->coeff[q->deg], 0) != 0);

  long int s = 1;
  mpz_t g;
  mpz_t h;
  mpz_t t;
  mpz_t tmp;
  int d;
  mpz_poly r;
  mpz_poly a;
  mpz_poly b;

  mpz_init(g);
  mpz_init(h);
  mpz_init(t);
  mpz_init(tmp);
  mpz_poly_init(r, -1);
  mpz_poly_init(a, p->deg);
  mpz_poly_init(b, q->deg);

  mpz_poly_set(a, p);
  mpz_poly_set(b, q);

  mpz_poly_content(g, a);
  mpz_poly_content(h, b);

  mpz_poly_divexact_mpz(a, a, g);
  mpz_poly_divexact_mpz(b, b, h);

#ifndef NDEBUG
  mpz_poly_content(tmp, a);
  ASSERT(mpz_cmp_ui(tmp, 1) == 0);
  mpz_poly_content(tmp, b);
  ASSERT(mpz_cmp_ui(tmp, 1) == 0);
#endif // NDEBUG

  mpz_pow_ui(t, g, (unsigned long int) b->deg);
  mpz_pow_ui(tmp, h, (unsigned long int) a->deg);
  mpz_mul(t, t, tmp);

  mpz_set_ui(g, 1);
  mpz_set_ui(h, 1);

  if (a->deg < b->deg) {
    mpz_poly_swap(a, b);

    if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
      s = -1;
    }
  }

  while (b->deg > 0) {
    d = a->deg - b->deg;

    if ((a->deg % 2) == 1 && (b->deg % 2) == 1) {
      s = -s;
    }

    mpz_poly_pseudo_remainder(r, a, b);
    mpz_poly_set(a, b);

    ASSERT(d >= 0);

    mpz_pow_ui(tmp, h, (unsigned long int) d);
    mpz_mul(tmp, g, tmp);

    mpz_poly_divexact_mpz(b, r, tmp);

    mpz_set(g, mpz_poly_lc(a));

#ifdef NDEBUG
    if (d == 0) {
      ASSERT(mpz_cmp_ui(h, 1) == 0);
    }
#endif // NDEBUG
    mpz_pow_ui(h, h, (unsigned long int) (d - 1));
    mpz_pow_ui(tmp, g, (unsigned long int) d);
    mpz_divexact(h, tmp, h);
  }

  //Prevent an error if b = 0.
  if (b->deg == -1) {
    mpz_set_ui(res, 0);
  } else {
    ASSERT(a->deg > 0);
    ASSERT(b->deg == 0);

    mpz_pow_ui(h, h, (unsigned long int) (a->deg - 1));

    ASSERT(a->deg >= 0);

    mpz_pow_ui(tmp, b->coeff[0], (unsigned long int) a->deg);
    mpz_divexact(h, tmp, h);

    mpz_mul_si(t, t, s);
    mpz_mul(h, h, t);
    mpz_set(res, h);
  }

  mpz_clear(g);
  mpz_clear(h);
  mpz_clear(t);
  mpz_clear(tmp);
  mpz_poly_clear(a);
  mpz_poly_clear(b);
  mpz_poly_clear(r);
}

void mpz_poly_discriminant(mpz_ptr res, mpz_poly_srcptr f)
{
    mpz_poly df;
    mpz_poly_init(df, f->deg);
    mpz_poly_derivative(df, f);
    mpz_poly_resultant(res, f, df);
    ASSERT(mpz_divisible_p(res, mpz_poly_lc_const(f)));
    mpz_divexact(res, res, mpz_poly_lc_const(f));
    mpz_poly_clear(df);
}

/* returns non-zero iff f is square-free in Z[x] */
int
mpz_poly_squarefree_p (mpz_poly_srcptr f)
{
  mpz_poly df;
  mpz_t res;
  int ret;

  mpz_poly_init (df, f->deg);
  mpz_poly_derivative (df, f);
  mpz_init (res);
  mpz_poly_resultant (res, f, df);
  ret = mpz_cmp_ui (res, 0);
  mpz_clear (res);
  mpz_poly_clear (df);
  return ret;
}

/**
 * Quick-and-dirty test if the polynomial f is irreducible over Z.
 * This function is extracted from dlpolyselect.c written by PZ
 *
 * Let p be a prime such that f has a root r modulo p.
 * Search by LLL a small linear combination between 1, r, ..., r^(d-1).
 * If f is not irreducible, then p will be root of a factor of degree <= d-1,
 * which will yield a linear dependency of the same order as the coefficients
 * of f, otherwise if f is irreducible the linear dependency will be of order
 * p^(1/d). 
 *
 * \return 0 by default, 1 if the poly is irreducible.
 * 
 * this function is without any warranty and at your own risk
 */
int mpz_poly_is_irreducible_z(mpz_poly_srcptr f)
{
  mpz_t p;
  int d = f->deg;
  int dg = d - 1;
  int i, j, nr;
  mpz_t a, b, det, r, *roots;
  size_t normf;
  int ret = 0; // init
  mat_Z g;

  mpz_init (p);
  mpz_poly_infinity_norm (p, f);
  normf = mpz_sizeinbase (p, 2);
  /* The following table might be useful to optimize the value of MARGIN:
     degree d    infinity_norm(f)    optimal MARGIN
        3               8                 5
        4               4                 5
  */
#define MARGIN 16
  /* add some margin bits */
  mpz_mul_2exp (p, p, MARGIN);
  mpz_pow_ui (p, p, d);

  roots = (mpz_t*) malloc (d * sizeof (mpz_t));
  for (i = 0; i < d; i++)
    mpz_init (roots[i]);

  do {
    mpz_nextprime (p, p);
    nr = mpz_poly_roots_mpz (roots, f, p);
    /* If f has no root mod p and degree <= 3, it is irreducible,
       since a degree 2 polynomial can only factor into 1+1 or 2,
       and a degree 3 polynomial can only factor into 1+2 or 3. */
    if (nr == 0 && d <= 3)
      {
        ret = 1;
        goto clear_and_exit;
      }
  } while (nr == 0);

  g.coeff = (mpz_t **) malloc ((dg + 2)*sizeof(mpz_t*));
  g.NumRows = g.NumCols = dg + 1;
  for (i = 0; i <= dg + 1; i ++) {
    g.coeff[i] = (mpz_t *) malloc ((dg + 2)*sizeof(mpz_t));
    for (j = 0; j <= dg + 1; j ++) {
      mpz_init (g.coeff[i][j]);
    }
  }

  mpz_init (det);
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 1);
  mpz_init_set (r, roots[0]);
  for (i = 0; i <= dg + 1; i ++) {
    for (j = 0; j <= dg + 1; j ++) {
      mpz_set_ui (g.coeff[i][j], 0);
    }
  }

  for (i = 1;  i <= dg + 1; i++) {
    for (j = 1; j <= dg + 1; j++) {
      if (i == 1) {
        if (j == 1) {
          mpz_set (g.coeff[j][i], p);
        }
        else {
          mpz_neg (g.coeff[j][i], r);
          mpz_mul (r, r, roots[0]);
        }
      }
      else
        mpz_set_ui (g.coeff[j][i], i==j);
    }
  }

  LLL (det, g, NULL, a, b);

  for (j = 1; j <= dg + 1; j ++)
    {
      /* the coefficients of vector j are in g.coeff[j][i], 1 <= i <= dg + 1 */
      mpz_abs (a, g.coeff[j][1]);
      for (i = 2; i <= dg + 1; i++)
        if (mpz_cmpabs (g.coeff[j][i], a) > 0)
          mpz_abs (a, g.coeff[j][i]);
      /* now a = max (|g.coeff[j][i]|, 1 <= i <= dg+1) */
      if (j == 1 || mpz_cmpabs (a, b) < 0)
        mpz_set (b, a);
    }
  /* now b is the smallest infinity norm */
  if (mpz_sizeinbase (b, 2) < normf + MARGIN / 2)
    ret = 0;
  else
    ret = 1;

  for (i = 0; i <= dg + 1; i ++) {
    for (j = 0; j <= dg + 1; j ++)
      mpz_clear (g.coeff[i][j]);
    free(g.coeff[i]);
  }
  free (g.coeff);

  mpz_clear (det);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (r);
 clear_and_exit:
  for (i = 0; i < d; i++)
    mpz_clear (roots[i]);
  free (roots);
  mpz_clear (p);

  return ret;
}

  
/* factoring polynomials */

void mpz_poly_factor_list_init(mpz_poly_factor_list_ptr l)
{
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_clear(mpz_poly_factor_list_ptr l)
{
    for(int i = 0 ; i < l->alloc ; i++)
        mpz_poly_clear(l->factors[i]->f);
    free(l->factors);
    memset(l, 0, sizeof(*l));
}

void mpz_poly_factor_list_flush(mpz_poly_factor_list_ptr l)
{
    /* There's a design choice here. We may elect to free everything.
     * Instead, we'll simply mark everything as zero, but keep all
     * allocated space.
     */
    for(int i = 0 ; i < l->alloc ; i++)
        l->factors[i]->f->deg = -1;
    l->size = 0;
}


void mpz_poly_factor_list_prepare_write(mpz_poly_factor_list_ptr l, int index)
{
    if (index >= l->alloc) {
        l->alloc = index + 4 + l->alloc / 4;
        l->factors = (mpz_poly_with_m*) realloc(l->factors, l->alloc * sizeof(mpz_poly_with_m));
        /* We need to set something. A zero polynomial has NULL storage
         * area, so that will do (realloc behaves as needed).  */
        for(int i = l->size ; i < l->alloc ; i++) {
            mpz_poly_init(l->factors[i]->f, -1);
            l->factors[i]->m = 0;
        }
    }
    if (l->size <= index)
        l->size = index + 1;
}

void mpz_poly_factor_list_push(mpz_poly_factor_list_ptr l, mpz_poly_srcptr f, int m)
{
    mpz_poly_factor_list_prepare_write(l, l->size);
    mpz_poly_set(l->factors[l->size - 1]->f, f);
    l->factors[l->size - 1]->m = m;
}

void mpz_poly_factor_list_fprintf(FILE* fp, mpz_poly_factor_list_srcptr l)
{
    for (int i = 0 ; i < l->size ; i++){
        char * res;
        int rc = mpz_poly_asprintf(&res,l->factors[i]->f);
        ASSERT_ALWAYS(rc >= 0);
        if (i) fprintf(fp, "*");
        fprintf(fp, "(%s)^%d", res, l->factors[i]->m);
        free(res);
    }
    fprintf(fp, "\n");
}
/* Squarefree factorization */

/* This auxiliary function almost does the sqf. It fills
 * lf->factors[stride*i] (i from 1 to deg(f)) with the factors with
 * multiplicity i in f.  lf->factors[0] is filled with the product whose
 * multiplicity is a multiple of the field characteristic.  returns max
 * multiplicity stored (multiplied by the stride value).
 *
 * (stride has an importance for the recursive call, for p small. E.g. on
 * f=g^p, we get called with g=f^(1/p) and stride=p).
 *
 * We make no effort to check that lf is clean on input, which is to be
 * guaranteed by the caller (e.g. sufficiently many polynomials, all
 * equal to 1 -- or 0, which can easily be identified as something
 * unset).
 *
 * Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_sqf_inner(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int stride, mpz_srcptr p)
{
    int r = 0;

    mpz_poly g, mi, mi1;
    mpz_poly t0,t1, T, tmp;
    mpz_poly_init(g, f->deg);
    mpz_poly_init(mi, f->deg);
    mpz_poly_init(mi1, f->deg);
    mpz_poly_init(t0, f->deg);
    mpz_poly_init(t1, f->deg);
    mpz_poly_init(T, f->deg);
    mpz_poly_init(tmp, f->deg);

    mpz_poly_derivative(t0, f);
    mpz_poly_gcd_mpz(g, f, t0, p);
    mpz_poly_divexact(mi, f, g, p);
    /* mi is f/gcd(f,f') == all repeated prime factors of f whose
     * multiplicity isn't a multiple of the field characteristic.
     */

    mpz_poly_set_xi(T, 0);
    for(int i = 1 ; mi->deg > 0 ; i++) {
        /* note the weird argument ordering */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, mi, NULL, i, p);
        mpz_poly_divexact(t1, f, t0, p);
        /* t1 = all polynomials in mi taken out from f with multiplicity i */
        mpz_poly_gcd_mpz(mi1, t1, mi, p);
        /* mi1 = almost like mi, but since factors with multiplicity i
         * are no longer in t1, there's absent from mi1 too. Whence
         * mi/mi1 is exactly the product of factors of multiplicity 1.
         */
        mpz_poly_factor_list_prepare_write(lf, i * stride);
        /* Use tmp so that we don't absurdly keep storage within
         * lf->factors */
        mpz_poly_divexact(tmp, mi, mi1, p);
        mpz_poly_set(lf->factors[i * stride]->f, tmp);
        /* multiplicity field still unused at this point */
        mpz_poly_pow_ui_mod_f_mod_mpz(t0, lf->factors[i * stride]->f, NULL, i, p);
        mpz_poly_mul(T, T, t0);
        mpz_poly_mod_mpz(T, T, p, NULL);
        mpz_poly_swap(mi, mi1);
        r = i * stride;
    }

    mpz_poly_factor_list_prepare_write(lf, 0);
    mpz_poly_divexact(lf->factors[0]->f, f, T, p);

    mpz_poly_clear(g);
    mpz_poly_clear(tmp);
    mpz_poly_clear(mi);
    mpz_poly_clear(mi1);
    mpz_poly_clear(t0);
    mpz_poly_clear(t1);
    mpz_poly_clear(T);
    return r;
}

/* Fills lf->factors[i] (i from 1 to deg(f)) with product of the factors
 * with multiplicity i in f, and to 1 if there are none. lf->factors[0]
 * is set to 1.  return the largest multiplicity stored. */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_sqf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0, mpz_srcptr p)
{
    /* factoring 0 doesn't make sense, really */
    assert(f0->deg >= 0);

    /* We'll call mpz_poly_factor_sqf_inner, possibly several times if
     * we are in small characteristic.
     */
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    mpz_poly_makemonic_mod_mpz(f, f0, p);
    assert(mpz_cmp_ui(mpz_poly_lc_const(f), 1) == 0);

    int m = 0;
    int pu = mpz_get_ui(p);  // see below
    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    for(int stride = 1 ; ; stride *= pu) {
        int r = mpz_poly_factor_sqf_inner(lf, f, stride, p);
        if (r > m) m = r;
        if (lf->factors[0]->f->deg == 0) {
            // if p is LAAAARGE, then of course we'll never have a linear
            // polynomial out of sqf_inner, thus we'll break early here.
            break;
        }
        /* divide coefficients */
        for(int i = 0 ; i <= lf->factors[0]->f->deg ; i++) {
            if (i % pu == 0) {
                mpz_set(f->coeff[i / pu], lf->factors[0]->f->coeff[i]);
            } else {
                assert (mpz_cmp_ui(lf->factors[0]->f->coeff[i], 0) == 0);
            }
        }
        f->deg = lf->factors[0]->f->deg / pu;
        mpz_poly_set_xi(lf->factors[0]->f, 0);
    }
    /* Now make sure that all factors in the factor list are non-zero */
    for(int i = 0 ; i < lf->size ; i++) {
        if (lf->factors[i]->f->deg < 0) {
            mpz_poly_set_xi(lf->factors[i]->f, 0);
        }
    }
    mpz_poly_clear(f);
    return m;
}




/* This performs distinct degree factorization */
/* Input polynomial must be squarefree -- otherwise repeated factors
 * probably won't show up in the factor list, or maybe at the wrong place
 * as parasites. */
/* Returns max degree of factors found (i.e. if the largest factors we
 * have are two factors of degree 7, we return 7, not 14). */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_ddf_inner(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0, mpz_srcptr p, int only_check_irreducible)
{
    mpz_poly g, gmx, x, tmp;
    mpz_poly f;
    mpz_poly_init(f, f0->deg);
    int i;

    /* factoring 0 doesn't make sense, really */
    assert(f0->deg >= 0);

    mpz_poly_makemonic_mod_mpz(f, f0, p);
    assert(mpz_cmp_ui(mpz_poly_lc_const(f), 1) == 0);

    /* reset the factor list completely */
    mpz_poly_factor_list_flush(lf);

    mpz_poly_init (g, 2 * f->deg - 1);
    mpz_poly_init (gmx, 2 * f->deg - 1);
    mpz_poly_init (x, 1);
    mpz_poly_init (tmp, f->deg);

    mpz_poly_set_xi(x, 1);
    mpz_poly_set_xi(g, 1);

    for (i = 1; i <= f->deg ; ++i) {
        if (2 * i > f->deg) {
            /* Then we know that the remaining f is irreducible.  */
            mpz_poly_factor_list_prepare_write(lf, f->deg);
            for( ; i < f->deg ; i++) {
                /* multiplicity field still unused at this point */
                mpz_poly_set_xi(lf->factors[i]->f, 0);
            }
            /* multiplicity field still unused at this point */
            mpz_poly_swap(lf->factors[f->deg]->f, f);
            break;
        }

        /* g <- g^p mod fp */
        mpz_poly_pow_mod_f_mod_mpz (g, g, f, p, p);

        /* subtract x */
        mpz_poly_sub (gmx, g, x);
        mpz_poly_mod_mpz(gmx, gmx, p, NULL);

        /* lf[i] <- gcd (f, x^(p^i)-x) */
        mpz_poly_factor_list_prepare_write(lf, i);
        /* multiplicity field still unused at this point */

        /* see remark in _sqf regarding the relevance of tmp for storage */
        mpz_poly_gcd_mpz(tmp, f, gmx, p);
        mpz_poly_divexact(f, f, tmp, p);
        mpz_poly_set(lf->factors[i]->f, tmp);

        /* Note for a mere irreducibility test: the length of the loop in
         * the irreducible case would still be deg(f)/2, and the penalty
         * caused by storing factors can be neglected.
         */
        if (only_check_irreducible && lf->factors[i]->f->deg > 0)
            break;

        if (f->deg == 0)
            break;
    }

    mpz_poly_clear (g);
    mpz_poly_clear (x);
    mpz_poly_clear (gmx);
    mpz_poly_clear (tmp);
    mpz_poly_clear (f);

    return i;
}

int mpz_poly_factor_ddf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f0, mpz_srcptr p)
{
    return mpz_poly_factor_ddf_inner(lf, f0, p, 0);
}

/* Note that this also works for non squarefree polynomials -- the factor
 * list returned by mpz_poly_factor_ddf will be rubbish, but the m ==
 * f->deg test will tell the truth. */
int mpz_poly_is_irreducible(mpz_poly_srcptr f, mpz_srcptr p)
{
    mpz_poly_factor_list lf;
    mpz_poly_factor_list_init(lf);
    int m = mpz_poly_factor_ddf_inner(lf, f, p, 1);
    mpz_poly_factor_list_clear(lf);
    return m == f->deg;
}

/* Naive factorization of polynomials over GF(2). We assume we're always
 * dealing with polynomials of degree at most 10, so we're beter off
 * simply enumerating potential factors...
 */

/*
 * Add 1 to f. If the constant term is equal to 1, set this term to 0 and
 *  propagate the addition of 1 to the next coefficient, and so on.
 *
 * f: the polynomial on which the addition is computed, the modifications are
 *  made on f.
 */
static void mpz_poly_add_one_in_F2(mpz_poly_ptr f)
{
  ASSERT(f->deg >= 1);

  int i = 0;
  while (mpz_cmp_ui(f->coeff[i], 1) == 0) {
    mpz_poly_setcoeff_si(f, i, 0);
    i++;
    if (i > f->deg) {
      break;
    }
  }
  mpz_poly_setcoeff_si(f, i, 1);
}

/*
 * Factorize naively a mpz_poly mod 2.
 *
 * Return the number of factors found.
 */
static int mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly_srcptr f,
    mpz_srcptr p)
{
  ASSERT(mpz_cmp_ui(p, 2) == 0);

  //make a copy of f.
  mpz_poly fcopy;
  mpz_poly_init(fcopy, f->deg);
  mpz_poly_set(fcopy, f);

  //reduce all the coefficient mod p.
  mpz_t coeff;
  mpz_init(coeff);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  //Purge list.
  mpz_poly_factor_list_flush(list);

  //If deg(f) in F2 is less than 1, we have the factor.
  if (fcopy->deg < 1) {
    mpz_clear(coeff);
    mpz_poly_clear(fcopy);
    ASSERT(list->size == 0);
    return list->size;
  } else if (fcopy->deg == 1) {
    mpz_clear(coeff);
    mpz_poly_factor_list_push(list, fcopy, 1);
    mpz_poly_clear(fcopy);
    ASSERT(list->size == 1);
    return list->size;
  }

  if (mpz_poly_is_irreducible(f, p)) {
    //If f is irreducible mod 2, fcopy is the factorisation of f mod 2.
    mpz_poly_factor_list_push(list, fcopy, 1);
  } else {
    //Create the first possible factor.
    mpz_poly tmp;
    mpz_poly_init(tmp, 1);
    mpz_poly_setcoeff_int64(tmp, 1, 1);

    //Enumerate all the possible factor of f mod 2.
    while (tmp->deg <= fcopy->deg) {
      //tmp is a possible factor.
      if (mpz_poly_is_irreducible(tmp, p)) {
        mpz_poly q;
        mpz_poly_init(q, 0);
        mpz_poly r;
        mpz_poly_init(r, 0);
        //Euclidean division of fcopy
        mpz_poly_div_qr(q, r, fcopy, tmp, p);
        //Power of the possible factor.
        unsigned int m = 0;
        //While fcopy is divisible by tmp.
        while (r->deg == -1) {
          //Increase the power of tmp.
          m++;
          mpz_poly_set(fcopy, q);
          if (fcopy->deg == 0 || fcopy->deg == -1) {
            //No other possible factor.
            break;
          }
          mpz_poly_div_qr(q, r, fcopy, tmp, p);
        }
        if (m != 0) {
          //Push tmp^m as a factor of f mod 2.
          mpz_poly_factor_list_push(list, tmp, m);
        }
        mpz_poly_clear(q);
        mpz_poly_clear(r);
      }
      //Go to the next possible polynomial in F2.
      mpz_poly_add_one_in_F2(tmp);
    }
    mpz_poly_clear(tmp);
  }

#ifndef NDBEBUG
  //Verify if the factorisation is good.
  mpz_poly_cleandeg(fcopy, -1);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  mpz_poly fmul;
  mpz_poly_init(fmul, -1);
  mpz_poly_set(fmul, list->factors[0]->f);
  for (int j = 1; j < list->factors[0]->m; j++) {
    mpz_poly_mul(fmul, fmul, list->factors[0]->f);
  }
  for (int i = 1; i < list->size ; i++) {
    for (int j = 0; j < list->factors[i]->m; j++) {
      mpz_poly_mul(fmul, fmul, list->factors[i]->f);
    }
  }
  for (int i = 0; i <= fmul->deg; i++) {
    mpz_poly_getcoeff(coeff, i, fmul);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fmul, i, coeff);
  }

  ASSERT(mpz_poly_cmp(fcopy, fmul) == 0);

  mpz_poly_clear(fmul);
#endif // NDBEBUG

  mpz_clear(coeff);
  mpz_poly_clear(fcopy);

  return list->size;
}

/*
 * this tries to split between squares and non-squares -- it's the most
 * basic split of course, and the other splits merely randomize on top of
 * this. Note that this building block must be changed for characteristic
 * two
 *
 * returns non-zero if a non-trivial or split is obtained.
 *
 * k is the degree of factors we are looking for.
 *
 * We split in two parts:
 * 
 *  - factors whose roots are squares in GF(p^k).
 *  - factors whose roots are non squares in GF(p^k).
 *
 * Note that we do not find the factor X this way; this is to be done by
 * the caller.
 *
 * Obviously, we include some shift, and hope that eventually there is a
 * shift that works. Large characteristic is generally happy with some
 * translation shift. Small characteristic may need more general shifts.
 *
 * Coefficients of f0 need not be reduced mod p.
 * Coefficients of g[0] and g[1] are reduced mod p.
 */
static void mpz_poly_factor_edf_pre(mpz_poly g[2], mpz_poly_srcptr f, int k, mpz_srcptr p)
{
    int nontrivial = 0;
    mpz_poly_set_xi(g[0], 0);
    mpz_poly_set_xi(g[1], 0);

    ASSERT_ALWAYS (f->deg > k);

    mpz_poly xplusa;
    mpz_poly_init(xplusa, 1);

    mpz_t half_pk;
    mpz_init(half_pk);
    mpz_pow_ui(half_pk, p, k);
    mpz_fdiv_q_ui(half_pk, half_pk, 2); /* (p^k-1)/2 */

    for(unsigned long a = 0 ; ! nontrivial ; a++) {
        /* we want to iterate on monic polynomials of degree <= k-1. */
        /* In order to bear in mind what happens in large enough
         * characteristic, we'll name these polynomials xplusa, although
         * it does not have to be x+a (and it can't be restricted to only
         * x+a if the characteristic is small -- that does not give
         * enough legroom).
         */
        if (mpz_fits_ulong_p(p)) {
            unsigned long pz = mpz_get_ui(p);
            if (a == 0) {
                /* special case, really */
                mpz_poly_set_xi(xplusa, 1);
            } else {
                /* write a in base p, and add 1 */
                int i = 0;
                for(unsigned long tmp = a ; tmp ; i++, tmp /= pz) {
                    mpz_poly_setcoeff_ui(xplusa, i, tmp % pz);
                }
                mpz_poly_setcoeff_ui(xplusa, i, 1);
            }
        } else {
            /* take the polynomial x+a */
            mpz_poly_set_xi(xplusa, 1);
            mpz_poly_setcoeff_ui(xplusa, 0, a);
        }

        mpz_poly_pow_mod_f_mod_mpz(g[0], xplusa, f, half_pk, p);

        mpz_poly_add_ui(g[1], g[0], 1);     /* (x+a)^((p^k-1)/2) + 1 */
        mpz_poly_sub_ui(g[0], g[0], 1);     /* (x+a)^((p^k-1)/2) - 1 */

        mpz_poly_mod_mpz(g[0], g[0], p, NULL);
        mpz_poly_mod_mpz(g[1], g[1], p, NULL);

        mpz_poly_gcd_mpz(g[0], g[0], f, p);
        mpz_poly_gcd_mpz(g[1], g[1], f, p);

        if (g[0]->deg + g[1]->deg < f->deg) {
            /* oh, we're lucky. x+a is a factor ! */
            int s = g[0]->deg > g[1]->deg;
            /* multiply g[s] by (x+a) */
            mpz_poly_mul_mod_f_mod_mpz(g[s], g[s], xplusa, f, p, NULL, NULL);
        }
        assert(g[0]->deg + g[1]->deg == f->deg);
        assert(g[0]->deg % k == 0);
        assert(g[1]->deg % k == 0);

        nontrivial += g[0]->deg != 0 && g[0]->deg != f->deg;
        nontrivial += g[1]->deg != 0 && g[1]->deg != f->deg;
    }
    mpz_clear(half_pk);
    mpz_poly_clear(xplusa);
}

/* This factors f, and for each factor q found, store q in lf.
 * Return the number of distinct factors found. */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
static int mpz_poly_factor_edf_inner(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int k, mpz_srcptr p, gmp_randstate_t rstate)
{
    if (f->deg == k) {
        mpz_poly_factor_list_push(lf, f, 1);
        return 1;
    }
    if (f->deg == 0) {
        return 0;
    }

    mpz_poly h[2];

    mpz_poly_init(h[0], f->deg);
    mpz_poly_init(h[1], f->deg);

    mpz_poly_factor_edf_pre(h, f, k, p);

    int n = 0;

    n += mpz_poly_factor_edf_inner(lf, h[0], k, p, rstate);
    mpz_poly_clear(h[0]);

    n += mpz_poly_factor_edf_inner(lf, h[1], k, p, rstate);
    mpz_poly_clear(h[1]);

    return n;
}

// returns f0->deg / d
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor_edf(mpz_poly_factor_list_ptr lf, mpz_poly_srcptr f, int k, mpz_srcptr p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list_flush(lf);

    if (mpz_cmp_ui(p, 2) == 0) {
        /* we need some other code for edf. Currently we have very naive
         * code, but that's good enough for small degree. */
        return mpz_poly_factor2(lf, f, p);
    }

    int v = mpz_poly_valuation(f);
    if (v) {
        /* Since our input is square-free, then we expect v==1.
         * Furthermore, k prescribes the extension field where the
         * expected roots live, thus for 0 it should really be 1. */
        ASSERT_ALWAYS(v == 1 && k == 1);
        mpz_poly_factor_list_prepare_write(lf, lf->size);
        mpz_poly_set_xi(lf->factors[lf->size-1]->f, 1);

        mpz_poly f1;
        mpz_poly_init(f1, f->deg - 1);
        mpz_poly_div_xi(f1, f, v);
        int n = 1 + mpz_poly_factor_edf_inner(lf, f1, k, p, rstate);
        mpz_poly_clear(f1);
        return n;
    }

    return mpz_poly_factor_edf_inner(lf, f, k, p, rstate);
}

typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_with_m_cmp(
        const mpz_poly_with_m * a,
        const mpz_poly_with_m * b)
{
    int r = mpz_poly_cmp((*a)->f, (*b)->f);
    if (r) return r;
    return ((*a)->m > (*b)->m) - ((*b)->m > (*a)->m);
}


/* putting it all together */
/* Coefficients of f need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
int mpz_poly_factor(mpz_poly_factor_list lf, mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_t rstate)
{
    mpz_poly_factor_list sqfs, ddfs, edfs;
    mpz_poly_factor_list_init(sqfs);
    mpz_poly_factor_list_init(ddfs);
    mpz_poly_factor_list_init(edfs);

    mpz_poly_factor_list_flush(lf);

    int maxmult = mpz_poly_factor_sqf(sqfs, f, p);
    for(int m = 1 ; m <= maxmult ; m++) {
        ASSERT_ALWAYS(sqfs->factors[m]->f->deg >= 0);
        if (sqfs->factors[m]->f->deg == 0) continue;
        int maxdeg = mpz_poly_factor_ddf(ddfs, sqfs->factors[m]->f, p);
        for(int k = 1 ; k <= maxdeg ; k++) {
            ASSERT_ALWAYS(ddfs->factors[k]->f->deg >= 0);
            if(ddfs->factors[k]->f->deg == 0) continue;
            mpz_poly_factor_edf(edfs, ddfs->factors[k]->f, k, p, rstate);
            for(int j = 0 ; j < edfs->size ; j++) {
                /* register this factor, with multiplicity m */
                /* cheat a bit... */
                mpz_poly_factor_list_prepare_write(lf, lf->size);
                mpz_poly_swap(lf->factors[lf->size-1]->f, edfs->factors[j]->f);
                mpz_poly_makemonic_mod_mpz(lf->factors[lf->size-1]->f, lf->factors[lf->size-1]->f, p);
                lf->factors[lf->size-1]->m = m;
            }
        }
    }
    mpz_poly_factor_list_clear(edfs);
    mpz_poly_factor_list_clear(ddfs);
    mpz_poly_factor_list_clear(sqfs);

    /* sort factors by degree and lexicographically */
    qsort(lf->factors, lf->size, sizeof(mpz_poly_with_m), (sortfunc_t) &mpz_poly_with_m_cmp);

    return lf->size;
}

int mpz_poly_number_of_real_roots(mpz_poly_srcptr f)
{
    /* This is coded in usp.c, with an interface pretty different from
     * what we have here.
     *
     * 0.0 in usp.c means: find a bound by yourself.
     */
    return numberOfRealRoots(f->coeff, f->deg, 0.0, 0, NULL);
}

int mpz_poly_factor_list_lift(mpz_poly_factor_list_ptr fac, mpz_poly_srcptr f, mpz_srcptr ell, mpz_srcptr ell2)
{
    mpz_poly f1; /* f - the product of its factors mod ell */

    mpz_poly_init(f1, -1);
    mpz_poly_set(f1, f);

    /* Lift all factors. Take g0 h0 unitary, g1 h1 of lesser
     * degree.
     * want (g0+ell*g1) * (h0 + ell*h1) = f
     *    g1 * h0 + h1 * g0 = f1 = (f - g0 * h0) / ell
     *
     * say a*g0 + b*h0 = 1
     *
     *    a*f1 = a*g1*h0+a*h1*g0 = h1 mod h0
     *    b*f1 = b*g1*h0+b*h1*g0 = g1 mod g0
     *
     */
    mpz_poly_set_xi(f1, 0);
    for(int i = 0 ; i < fac->size ; i++) {
        mpz_poly_srcptr g0 = fac->factors[i]->f;
        mpz_poly_mul_mod_f_mod_mpz(f1, f1, g0, 0, ell2, NULL, NULL);
    }
    mpz_poly_sub_mod_mpz(f1, f, f1, ell2);
    mpz_poly_divexact_mpz(f1, f1, ell);

    for(int i = 0 ; i < fac->size ; i++) {
        mpz_poly g1, d, a, b, h0;
        mpz_poly_ptr g = fac->factors[i]->f;
        mpz_poly_srcptr g0 = g;   /* alias */
        mpz_poly_init(g1, g0->deg - 1);
        mpz_poly_init(h0, f->deg - g0->deg);
        mpz_poly_init(d, 0);
        mpz_poly_init(a, f->deg - g0->deg - 1);
        mpz_poly_init(b, g0->deg - 1);
        if (fac->factors[i]->m != 1) {
            fprintf(stderr, "Ramified ell not supported\n");
            exit(EXIT_FAILURE);
        }
        ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_lc_const(g), 1) == 0);
        /* compute h0 = product of other factors */
        mpz_poly_divexact(h0, f, g0, ell);

        /* say a*g0 + b*h0 = 1 */
        mpz_poly_xgcd_mpz(d, g0, h0, a, b, ell);
        ASSERT_ALWAYS(d->deg == 0);
        ASSERT_ALWAYS(mpz_cmp_ui(d->coeff[0], 1) == 0);

        /* b*f1 = b*g1*h0+b*h1*g0 = g1 mod g0 */
        mpz_poly_mul_mod_f_mod_mpz(g1, f1, b, g0, ell, NULL, NULL);

        /* now update g */
        mpz_poly_mul_mpz(g1, g1, ell);
        mpz_poly_add(g, g, g1);

        mpz_poly_clear(g1);
        mpz_poly_clear(a);
        mpz_poly_clear(b);
        mpz_poly_clear(d);
        mpz_poly_clear(h0);
    }

    mpz_poly_clear(f1);

    return 1;
}

int mpz_poly_factor_and_lift_padically(mpz_poly_factor_list_ptr fac, mpz_poly_srcptr f, mpz_srcptr ell, int prec, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_lc_const(f), 1) == 0);

    mpz_poly_factor(fac, f, ell, rstate);
    mpz_t ellx;
    mpz_init_set(ellx, ell);

    /* keep for a rainy day: compute the list of prime powers we want to
     * pass by. (the mpz_poly_factor_list_lift function is happy to take
     * ell^m and ell^n, with n being 2m or 2m-1). */
    ASSERT_ALWAYS(prec == 2);
    mpz_mul(ellx, ell, ell);
    mpz_poly_factor_list_lift(fac, f, ell, ellx);
    mpz_clear(ellx);

    return 1;
}

std::string cxx_mpz_poly::print_poly(std::string const& var) const
{
    std::ostringstream os;
    for(int i = 0 ; i <= x->deg ; i++) {
        int r = mpz_cmp_ui(x->coeff[i], 0);
        if (r == 0) continue;
        if (r > 0 && os.str().size())
            os << "+";
        if (i == 0) {
            os << x->coeff[i];
        } else {
            if (mpz_cmp_ui(x->coeff[i], -1) == 0) {
                os << "-";
            } else if (mpz_cmp_ui(x->coeff[i], 1) != 0) {
                os << x->coeff[i] << "*";
            }
            os << var;
            if (i > 1) os << "^" << i;
        }

    }
    return os.str();
}

/* functions for Joux--Lercier and Generalized Joux--Lercier */
/**
 * \brief set the (mpz_t) coefficients of f from ::counter
 * \param[out] f               polynomial
 * \param[out] max_abs_coeffs  largest absolute value of coefficients of f
 * \param[out] next_counter    the next counter to try after the present one, larger than counter+1 sometimes
 * \param[in]  counter         counter storing the coefficients of f
 * \param[in]  bound           max absolute value of coefficients of f
 *
 * \return 1 if the poly is valid (content = 1, not a duplicate, and +/-1 is not a root) 
 *           (in fact the value returned is != 0) and in this case,
 *           max_abs_coeffs is set to the largest absolute value of the coefficients of f
 * \return 0 if the poly corresponding to the counter is a duplicate, and in that case,
 *           max_abs_coeffs is set to an error_code (poly duplicate, +/- 1 root, content != 1...)
 *
 * Given ::counter, compute the coefficients of ::f.
 * Assume that the leading coefficient of f is > 0, the next one is >= 0, and the last one is not 0
 * 1 <= f[deg] <= bound
 * 0 <= f[deg-1] <= bound
 * -bound <= f[i] <= bound for 0 < i < deg-1
 * no other test is made: you need to check after if the poly is square-free, irreducible, etc.
 */
int mpz_poly_setcoeffs_counter(mpz_poly_ptr f, int* max_abs_coeffs, unsigned long *next_counter, int deg, unsigned long counter, unsigned int bound){
    unsigned int i;
    unsigned long idx = counter, mod_j=0, mod_k=0;
    int j;
    int *fint, ok = 1;
    mpz_t content;
    //unsigned long idx_next_poly_ok = idx;
    int error_code = 0; // because I want to know what is the pb
#define POLY_EQUIV_INV_X -1 // the poly is the same as (+/-)f(1/x)
#define POLY_EQUIV_MINUS_X -2 // the poly is the same as (+/-)f(-x)
#define POLY_EQUIV_MINUS_COEFFS -3
#define POLY_ROOT_ONE -4
#define POLY_ROOT_MINUS_ONE -5
#define POLY_CONTENT -6

    fint = (int *) malloc ((deg + 1)*sizeof(int));
    // by default the next counter is counter+1
    *next_counter = counter+1;

    /* compute polynomial f and fint of index idx */
    /* assume that f was already initialized with mpz_poly_init(f, deg) */
    f->deg = deg;

    /* we take 1 <= f[deg] <= bound */
    fint[deg] = 1 + (idx % bound);
    idx = idx / bound;
    /* we take 0 <= f[deg-1] <= bound */
    fint[deg-1] = idx % (bound + 1);
    idx = idx / (bound + 1);
    for (i = deg-2; i > 0; i --)
      {
	/* we take -bound <= f[i] <= bound */
	fint[i] = (idx % (2 * bound + 1)) - bound;
        idx = idx / (2 * bound + 1);
    }
    
    /* we take -bound <= f[0] <= bound, f[0] <> 0,
       which makes 2*bound possible values */
    ASSERT(idx < 2 * bound);
    fint[0] = (idx < bound) ? idx - bound : idx - (bound - 1);
    
    
    /* since f and the reversed polynomial are equivalent, we can assume
       |f[deg]| < |f[0]| or (|f[deg]| = |f[0]| and |f[deg-1]| < |f[1]|) or ... */

    /* other point of view
    ok = 1;
    j = 0;
    while ((abs(fint[deg-j]) == abs(fint[j])) && (2*j < deg)){
      j++;
    }
    ok = ((j*2 >= deg) || (abs(fint[deg-j]) < abs(fint[j])));
    */
    ok = 1;
    for (int i = 0; 2 * i < deg; i++)
      {
        if (abs(fint[deg-i]) != abs(fint[i]))
          {
            /* we want |f[deg-i]| < |f[i]| */
            ok = abs(fint[deg-i]) < abs(fint[i]);
            break;
          }
      }
    if(!ok){
      error_code = POLY_EQUIV_INV_X; // the poly is an equivalent to (+/-)x^d*f(1/x)
    }
    if(abs(fint[deg]) >= abs(fint[0]) && (bound > 1)){
      // next valid poly: the next one is an increment of the leading coefficient.
      // but either leading coeff  = |constant coeff|, and incrementing by one will lead to
      // leading coeff > |constant coeff|, and this is not a valid poly,
      // or we already are in the case:
      // leading coeff > |constant coeff|. In both cases,
      // since the leading coefficient is > 0, it means:
      // set the leading coeff to 1 and increment the next coefficient.
      *next_counter = counter - (counter % bound) + bound;
      // then, what if f[0] = +/- 1?
      if(abs(fint[0]) == 1){
	// means that setting ld=1 might not be enough: what if f[d-1] > |f[1]| now ?
	// in that case, counter++ is not enough because f[d]=1 is the only possibility,
	// so the next counter we are looking for is with f[d-1]++
	if(abs(fint[deg-1]+1) > abs(fint[1])){
	  // setting ld=1 then incrementing the second high deg coeff is not enough
	  idx = *next_counter;
	  // set the second high deg coeff to 0
	  *next_counter = idx - (idx % (bound*(bound+1))) + (idx % bound);
	  // and increment the third high deg coeff
	  *next_counter += bound*(bound+1);
	  j = 2;
	  mod_j = bound*(bound+1)*(2*bound+1);
	  mod_k = mod_j*(2*bound+1);
	  while ((deg-j > j) && (fint[j-1] == 0)) {
	    //setting the deg-j+1 coeff to 0 and incrementing the next deg-j coeff migh not be enough
	    if (abs(fint[deg-j]+1) > abs(fint[j])){
	      if ((fint[deg-j]+1) < fint[j]){// negative values
		//set fint[deg-j] to fint[j]
		idx = *next_counter;
		*next_counter = idx - (idx % mod_j) + (-abs(fint[j])+bound)*mod_j;
	      }else{//(fint[deg-j]+1) is too large anyway: set it to its minimal value which is -|fint[j]| and do fint[d-j-1]++
		// set fint[deg-j] to -|fint[j]| and
		// fint[deg-j-1]++
		idx = *next_counter;
		*next_counter = idx - (idx % mod_j) + (-abs(fint[j])+bound)*mod_j + mod_k;
	      }
	    }
	    j++;
	    mod_j = mod_k;
	    mod_k *= (2*bound+1);
	  }// end while
	}
      }
    }// end of the computation of the next valid non-duplicate counter

    /* since f(x) is equivalent to (-1)^deg*f(-x), if f[deg-1] = 0, then the
       largest i = deg-3, deg-5, ..., such that f[i] <> 0 should have f[i] > 0 */
    if (ok && fint[deg-1] == 0)
      {
        for (int i = deg - 3; i >= 0; i -= 2)
          {
            if (fint[i])
              {
                ok = fint[i] > 0;
                break;
              }
          }
      }
    if((!ok) && (error_code == 0)){
      error_code = POLY_EQUIV_MINUS_X; // the poly is an equivalent to (-1)^deg*f(-x)
    }

    /* if |f[i]| = |f[deg-i]| for all i, then [f[deg], f[deg-1], ..., f[1], f[0]]
       is equivalent to [s*f[0], s*t*f[1], s*t^2*f[2], ...], where
       s = sign(f[0]), and s*t^i is the sign of f[i] where i is the smallest
       odd index > 0 such that f[i] <> 0.  */
    if (ok && fint[deg] == abs(fint[0])) /* f[deg] > 0 */
      {
        int s = (fint[0] > 0) ? 1 : -1, t = 0, i;
        for (i = 1; i <= deg; i++)
          {
            if (2 * i < deg && abs(fint[deg-i]) != abs(fint[i]))
              break;
            if (t == 0 && fint[i] != 0 && (i & 1))
              t = (fint[i] > 0) ? s : -s;
          }
        if (2 * i >= deg) /* |f[i]| = |f[deg-i]| for all i */
          {
            /* if t=0, then all odd coefficients are zero, but since
               |f[deg-i]| = |f[i]| this can only occur for deg even, but then
               we can set t to 1, since t only affects the odd coefficients */
            if (t == 0)
              t = 1;
            for (i = 0; i <= deg; i++)
              {
                if (fint[deg-i] != s * fint[i])
                  {
                    ok = fint[deg-i] > s * fint[i];
                    break;
                  }
                s = s * t;
              }
          }
      }
    if((!ok) && (error_code == 0)){
      error_code = POLY_EQUIV_MINUS_COEFFS;
    }

    if (ok){
      /* check if +-1 is a root of f (very common): compute the sum of the coefficients, and the alterned sum: it should be != 0 */
      // evaluating at 1 means sum the coefficients
      int sum_even_coeffs = 0;
      int sum_odd_coeffs = 0;
      for(j=0; j<=deg; j += 2){
	sum_even_coeffs += fint[j];
      }
      for(j=1; j<=deg; j += 2){
	sum_odd_coeffs += fint[j];
      }
      ok = (sum_even_coeffs + sum_odd_coeffs) != 0;
      if (!ok){
	error_code = POLY_ROOT_ONE;
      }else{// eval at -1: alterning sum of coefficients
	ok = (sum_even_coeffs - sum_odd_coeffs) != 0;
	if (!ok){
	  error_code = POLY_ROOT_MINUS_ONE;
	}
      }
    }

    /* set mpz_poly */
    f->deg = deg;
    for (j = 0; j <= deg; j ++)
      mpz_set_si (f->coeff[j], fint[j]);

    if (ok){
      /* content test */
      mpz_init (content);
      mpz_poly_content (content, f);
      ok = mpz_cmp_ui (content, 1) == 0; /* duplicate with f/t */
      if (ok == 0){
	error_code = POLY_CONTENT;
      }
      mpz_clear (content);
    }

    // compute the max absolute value of coefficients,
    // usefull for computing the lognorm after
    if(ok){
      *max_abs_coeffs = fint[0]; // this is positive anyway
      for(j=0; j <= deg; j++){
	if ((fint[j] < 0) && (*max_abs_coeffs < -fint[j])){
	  *max_abs_coeffs = -fint[j];
	}else{// fint[i] >= 0
	  if (*max_abs_coeffs < fint[j]){
	    *max_abs_coeffs = fint[j];
	  }
	}
      }
    }else{
      *max_abs_coeffs = error_code;
    }

    free(fint);
    return ok;

#undef POLY_EQUIV_INV_X
#undef POLY_EQUIV_MINUS_X
#undef POLY_EQUIV_MINUS_COEFFS
#undef POLY_ROOT_ONE
#undef POLY_ROOT_MINUS_ONE
#undef POLY_CONTENT
}

/**
 return the total number of polys f that satisfy:
 deg(f) = deg exactly (leading coefficient != 0)
 |f_i| <= bound: the coefficients are bounded by bound
 f_deg > 0: leading coeff strictly positive
 f_{deg-1} > 0
 f_0 != 0: constant coeff non-zero

This corresponds to the number of polynomials that can be enumerated with the function 
mpz_poly_setcoeffs_counter
*/
unsigned long mpz_poly_cardinality(int deg, unsigned int bound){
  /* we take 1 <= f[d] <= bound */
  /* we take 0 <= f[d-1] <= bound */
  /* we take -bound <= f[i] <= bound for d-2 >= i > 0 */
  /* we take -bound <= f[0] < bound, f[0] <> 0,
     which makes 2*bound possible values */
  unsigned long number_polys=0;
  int i;
  if(deg >= 3){
    number_polys = bound*(bound+1);
    for (i = deg-2; i > 0; i --){
      number_polys *= 2*bound + 1;
    }
    number_polys *= 2*bound;
  }
  return number_polys;
}

/**
 * return the counter in basis ::bound corresponding to f
 *
 */
unsigned long mpz_poly_getcounter(mpz_poly_ptr f, unsigned int bound){
  unsigned int counter;
  int i;
  // leading coeff: 1 <= f->coeff[deg] <= bound
  counter = mpz_get_ui(f->coeff[f->deg]);
  counter *= bound;
  // next leading coeff: 0 <= f->coeff[deg-1] <= bound
  counter += mpz_get_ui(f->coeff[f->deg -1]);
  counter *= (bound+1);
  // next coeffs: -bound <= f->coeff[i] <= bound
  for(i=f->deg-2; i > 0; i--){
    counter += mpz_get_si(f->coeff[i]) + bound;
    counter *= 2*bound + 1;
  }
  if(mpz_sgn(f->coeff[0]) < 0){
    counter += mpz_get_si(f->coeff[0]) + bound;
  }else{// f->coeff[0] > 0
    counter += mpz_get_ui(f->coeff[0]) + bound - 1;
  }
  return counter;
}

void  mpz_poly_setcoeffs_counter_print_error_code(int error_code){
  switch(error_code){
  case -1: printf("-1: f <-> +/- x^d f(1/x) \n"); break;
  case -2: printf("-2: f <-> +/- f(-x) \n"); break;
  case -3: printf("-3: f <-> another one with permuted/sign coeffs\n"); break;
  case -4: printf("-4: f(1) = 0\n"); break;
  case -5: printf("-5: f(-1) = 0\n"); break;
  case -6: printf("-6: content(f) > 1\n"); break;
  default: printf(" 0: undetermined error.\n");
  }
}
