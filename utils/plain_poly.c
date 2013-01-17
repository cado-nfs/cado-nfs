/* Arithmetic on polynomials over Z/pZ, with coefficients represented by
   the 'plain_poly_coeff_t' type, defined in plain_poly.h

   Since products are performed with the '*' operator, and several products
   are accumulated before a reduction by '% p', p should be at most half the
   size of the 'long' type.

   Assumes that (p-1)+deg(f)*(p-1)^2 fits in a long [for plain_poly_div_r],
   where f is the algebraic polynomial. This gives the following bounds:

   deg(f)      plain_poly_coeff_t=int32_t        plain_poly_coeff_t=int64_t
     2         p <= 32749          p <= 2147483647
     3         p <= 26737          p <= 1753413037
     4         p <= 23167          p <= 1518500213
     5         p <= 20719          p <= 1358187913
     6         p <= 18917          p <= 1239850223
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "macros.h"
#include "plain_poly.h"
#include "gmp_aux.h"

/* return non-zero if (p-1)+d*(p-1)^2 fits in type plain_poly_coeff_t (defined in utils.h) */
int
plain_poly_fits (unsigned int d, uint64_t p)
{
  uint64_t max = 0;

#if ! ( PLAIN_POLY_COEFF_MAX <= UINT64_MAX)
#error "We have a (small) problem here."
/* It's really a matter of deciding something ``optimal'' or not */
#endif

  max = (uint64_t) PLAIN_POLY_COEFF_MAX - (p - 1); /* > 0 since p <= LONG_MAX */

  /* computes floor(max/d) */
  max = max / d;
  
  if (max < p - 1)
    return 0;
  
  max = max / (p - 1);
  return (max >= p - 1);
}

/* allocate an array of d coefficients, and initialize it */
static plain_poly_coeff_t*
alloc_long_array (int d)
{
  plain_poly_coeff_t *f;

  f = (plain_poly_coeff_t*) malloc (d * sizeof (plain_poly_coeff_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* reallocate an array to d coefficients */
static plain_poly_coeff_t*
realloc_long_array (plain_poly_coeff_t *f, int d)
{
  f = (plain_poly_coeff_t*) realloc (f, d * sizeof (plain_poly_coeff_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* initialize a polynomial with maximal degree d */
void
plain_poly_init (plain_poly_t f, int d)
{
  ASSERT (d >= 0);
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_long_array (d + 1);
  f->alloc = d + 1;
}

/* clear a polynomial */
void
plain_poly_clear (plain_poly_t f)
{
  free (f->coeff);
  f->alloc = 0;
}

/* realloc f to (at least) n coefficients */
static void
plain_poly_realloc (plain_poly_t f, int n)
{
  if (f->alloc < n)
    {
      f->coeff = realloc_long_array (f->coeff, n);
      f->alloc = n;
    }
}

/* Return 1/s mod t.
   Assumes t > 0, |s| < t, and gcd(s,t) = 1.
 */
static plain_poly_coeff_t
invert_si (plain_poly_coeff_t s, plain_poly_coeff_t t)
{
  plain_poly_coeff_t u1, v1, q;
  plain_poly_coeff_t u2, v2;

  ASSERT (t > 0);
 
  u1 = 1;
  v1 = 0;
  u2 = (s > 0) ? s : -s;
  v2 = t;

  /* invariant: u2 = u1*s (mod t)
                v2 = v1*s (mod t) */

  while (v2 != 0)
    {
      /* (u2,v2) are elements of the remainder sequence starting from
         |s|, thus remain nonnegative and in the interval 0..t */

      /* unroll twice and swap u/v */
      q = u2 / v2;
      u1 = u1 - q * v1;
      u2 = u2 - q * v2;
      
      if (u2 == 0)
	{
	  u1 = v1;
	  break;
	}

      q = v2 / u2;
      v1 = v1 - q * u1;
      v2 = v2 - q * u2;
    }

  /* now u2 = 1 = u1*s + w1*t with |u1| < t */

  if (u1 < 0)
    u1 = u1 + t;
  
  return (s > 0) ? u1 : t - u1; /* since u1 <> 0, 0 < t - u1 < t */
}

/* f <- g */
static void
plain_poly_set (plain_poly_t f, const plain_poly_t g)
{
  int i, d = g->degree;

  plain_poly_realloc (f, d + 1);
  f->degree = d;
  for (i = 0; i <= d; i++)
    f->coeff[i] = g->coeff[i];
}

/* f <- f/lc(f) mod p.
   Assumes the coefficients of f are in [-p+1, p-1].
   Requires that (p-1)^2 fits in a long.
 */
static void
plain_poly_make_monic (plain_poly_t f, plain_poly_coeff_t p)
{
  int d = f->degree, i;
  plain_poly_coeff_t ilc;

  if (f->coeff[d] == 1)
    return;

  ilc = invert_si (f->coeff[d], p);
  for (i = 0; i < d; i++)
    f->coeff[i] = (f->coeff[i] * ilc) % p;
  f->coeff[d] = 1;
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
plain_poly_set_mod (plain_poly_t fp, mpz_t *f, int d, plain_poly_coeff_t p)
{
  int i;

  while (d >= 0 && mpz_divisible_ui_p (f[d], p))
    d --;
  /* Let's allow N that have a small factor and remove this assertion:
  ASSERT (d >= 0); // f is 0 mod p: should not happen in the CADO-NFS context
                   // since otherwise p would divide N, indeed f(m)=N
  */
  if (d<0)
      return d;
  plain_poly_realloc (fp, d + 1);
  fp->degree = d;
  for (i = 0; i <= d; i++)
    fp->coeff[i] = mpz_fdiv_ui (f[i], p); /* 0 <= c[i] < p */
  plain_poly_make_monic (fp, p);

  return d;
}

/* f <- a*x+b, a <> 0 */
static void
plain_poly_set_linear (plain_poly_t f, plain_poly_coeff_t a, plain_poly_coeff_t b)
{
  plain_poly_realloc (f, 2);
  f->degree = 1;
  f->coeff[1] = a;
  f->coeff[0] = b;
}

/* swap f and g */
static void
plain_poly_swap (plain_poly_t f, plain_poly_t g)
{
  int i;
  plain_poly_coeff_t *t;

  i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->degree;
  f->degree = g->degree;
  g->degree = i;
  t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

/* h <- f*g mod p
   Assumes the coefficients of f and g are in [-p+1, p-1].
   Work properly if min(deg(f), deg(g))*(p-1)^2 fits in a long.
 */
static void
plain_poly_mul (plain_poly_t h, const plain_poly_t f, const plain_poly_t g, plain_poly_coeff_t p)
{
  int df = f->degree;
  int dg = g->degree;
  int i, j;
  plain_poly_t res;

  plain_poly_init (res, df + dg);
  for (i = 0; i <= df + dg; ++i)
    res->coeff[i] = 0;
  for (i = 0; i <= df; ++i)
    for (j = 0; j <= dg; ++j)
      res->coeff[i+j] += f->coeff[i] * g->coeff[j];
  /* reduce mod p */
  for (i = 0; i <= df + dg; i++)
    res->coeff[i] %= p;
  res->degree = df + dg;
  plain_poly_set (h, res);
  plain_poly_clear (res);
}


/* h <- g^2 mod p, g and h must differ.
   Assumes the coefficients of g are in [-p+1, p-1].
   Assumes deg(g)*(p-1)^2 fits in a long.
 */
static void
plain_poly_sqr (plain_poly_t h, const plain_poly_t g, plain_poly_coeff_t p)
{
  int i, j, dg = g->degree;
  plain_poly_coeff_t *gc, *hc;

  ASSERT (dg >= -1);
  if (dg == -1) /* g is zero */
    {
      h->degree = -1;
      return;
    }
  plain_poly_realloc (h, 2 * dg + 1);
  gc = g->coeff;
  hc = h->coeff;
  for (i = 0; i <= dg; i++)
    for (j = i + 1; j <= dg; j++)
      if (i == 0 || j == dg)
	hc[i + j] = gc[i] * gc[j];
      else
	hc[i + j] += gc[i] * gc[j];
  for (i = 1; i < 2 * dg; i++)
    hc[i] <<= 1;
  hc[0] = gc[0] * gc[0];
  hc[2 * dg] = gc[dg] * gc[dg];
  for (i = 1; i < dg; i++)
    hc[2 * i] += gc[i] * gc[i];
  /* reduce mod p */
  for (i = 0; i <= 2 * dg; i++)
    hc[i] %= p;
  h->degree = 2 * dg;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
static void
plain_poly_normalize (plain_poly_t h)
{
  int dh = h->degree;

  while (dh >= 0 && h->coeff[dh] == 0)
    dh --;
  h->degree = dh;
}

/* g <- h mod p */
static void
plain_poly_mod_ui (plain_poly_t g, const plain_poly_t h, plain_poly_coeff_t p)
{
  int i;
  int dh = h->degree;

  plain_poly_realloc (g, dh + 1);
  for (i = 0; i <= dh; i++)
    g->coeff[i] = h->coeff[i] % p;
  g->degree = dh;
  plain_poly_normalize (g);
}

/* h <- (x+a)*h mod p with 0 <= a < p.
   Requires that p*(p-1) fits in a long.
 */
static void
plain_poly_mul_x (plain_poly_t h, plain_poly_coeff_t a, plain_poly_coeff_t p)
{
  int i, d = h->degree;
  plain_poly_coeff_t *hc;

  if (d == -1) /* h = 0 */
    return;

  plain_poly_realloc (h, d + 2); /* (x+a)*h has degree d+1, thus d+2 coeffs */
  hc = h->coeff;
  hc[d + 1] = hc[d];
  for (i = d - 1; i >= 0; i--)
    hc[i + 1] = (hc[i] + a * hc[i + 1]) % p;
  hc[0] = (a * hc[0]) % p;
  h->degree = d + 1;
}

/* h <- g - x mod p */
static void
plain_poly_sub_x (plain_poly_t h, const plain_poly_t g, plain_poly_coeff_t p)
{
  int i, d = g->degree;

  /* if g has degree d >= 2, then g-x has degree d too;
     otherwise g-x has degree <= 1 */
  plain_poly_realloc (h, ((d < 2) ? 1 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 1; i++)
    h->coeff[i] = 0;
  h->coeff[1] = (h->coeff[1] - 1) % p;
  h->degree = (d < 1) ? 1 : d;
  plain_poly_normalize (h);
}

/* h <- g - 1 mod p */
static void
plain_poly_sub_1 (plain_poly_t h, const plain_poly_t g, plain_poly_coeff_t p)
{
  int i, d = g->degree;

  /* g-1 has degree d if d >= 1, and degree 0 otherwise */
  plain_poly_realloc (h, ((d < 1) ? 0 : d) + 1);
  for (i = 0; i <= d; i++)
    h->coeff[i] = g->coeff[i];
  for (i = d + 1; i <= 0; i++)
    h->coeff[i] = 0;
  h->coeff[0] = (h->coeff[0] - 1) % p;
  h->degree = (d < 0) ? 0 : d;
  plain_poly_normalize (h);
}

/* h <- rem(h, f) mod p, f not necessarily monic.
   Requires that (p-1)+deg(f)*(p-1)^2 fits in a long.
*/
static void
plain_poly_div_r (plain_poly_t h, const plain_poly_t f, plain_poly_coeff_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  plain_poly_coeff_t *hc = h->coeff, t = 1;

  monic = f->coeff[d] == 1;
  if (!monic)
    t = invert_si (f->coeff[d], p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
	hc[dh] = (hc[dh] * t) % p;
      for (i = 0; i < d; i++)
	hc[dh - d + i] -= hc[dh] * f->coeff[i];
      /* it is not necessary to reduce h[j] mod p here for j < dh */
      do
	{
	  dh --;
	  if (dh >= 0)
	    hc[dh] = hc[dh] % p;
	}
      while (dh >= 0 && hc[dh] == 0);
    }
  h->degree = dh;
  /* reduce mod p */
  plain_poly_mod_ui (h, h, p);
}

/* q <- divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
static void
plain_poly_divexact (plain_poly_t q, plain_poly_t h, const plain_poly_t f, plain_poly_coeff_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  plain_poly_coeff_t *hc = h->coeff, t = 1;

  ASSERT (d >= 0);
  ASSERT (dh >= 0);
  ASSERT (dh >= d);
  plain_poly_realloc (q, dh + 1 - d);
  q->degree = dh - d;
  monic = f->coeff[d] == 1;
  if (!monic)
    t = invert_si (f->coeff[d], p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
	hc[dh] = (hc[dh] * t) % p;
      q->coeff[dh - d] = hc[dh];
      /* we only need to update the coefficients of degree >= d of h,
	 i.e., we want i >= 2d - dh */
      for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++)
	hc[dh - d + i] -= hc[dh] * f->coeff[i];
      h->degree = dh - 1;
      /* it is not necessary to reduce h[j] mod p here for j < dh */
      dh --;
      if (dh >= 0)
	hc[dh] = hc[dh] % p;
    }
  /* no need to normalize q, since deg(q) = deg(h) - deg(f) */
}

/* fp <- gcd (fp, g), clobbers g */
static void
plain_poly_gcd (plain_poly_t fp, plain_poly_t g, plain_poly_coeff_t p)
{
  while (g->degree >= 0)
    {
      plain_poly_div_r (fp, g, p);
      /* now deg(fp) < deg(g): swap f and g */
      plain_poly_swap (fp, g);
    }
}

#if 0
static void
plain_poly_out (FILE *fp, plain_poly_t f)
{
  if (f->degree < 0)
    fprintf (fp, "0\n");
  else
    {
      int i;
      for (i = 0; i <= f->degree; i++)
	if (f->coeff[i] > 0)
	  fprintf (fp, "+%lld*x^%d", (long long int) f->coeff[i], i);
	else if (f->coeff[i] < 0)
	  fprintf (fp, "%lld*x^%d", (long long int) f->coeff[i], i);
      fprintf (fp, ";\n");
    }
}
#endif

/* returns f(x) mod p */
static plain_poly_coeff_t
plain_poly_eval (plain_poly_t f, plain_poly_coeff_t x, plain_poly_coeff_t p)
{
  int i, d = f->degree;
  plain_poly_coeff_t v;
  
  ASSERT (d >= 0);
  v = f->coeff[d];
  for (i = d - 1; i >= 0; i--)
    v = (v * x + f->coeff[i]) % p;
  return v;
}

/* Return the number n of roots of f mod p using a naive algorithm.
   If r is not NULL, put the roots in r[0], ..., r[n-1]. */
static int
plain_poly_roots_naive (plain_poly_coeff_t *r, plain_poly_t f, const plain_poly_coeff_t p)
{
  int n = 0;
  plain_poly_coeff_t x;
  
  for (x = 0; x < p; x++)
    if (plain_poly_eval (f, x, p) == 0)
      {
	if (r != NULL)
	  r[n] = x;
	n ++;
      }
  return n;
}

/* g <- (x+a)^e mod (fp, p), using auxiliary polynomial h.
   Assumes 0 <= a < p.
 */
static void
plain_poly_powmod_ui (plain_poly_t g, plain_poly_t fp, plain_poly_t h,
                      plain_poly_coeff_t a, plain_poly_coeff_t e,
                      plain_poly_coeff_t p)
{
  int k = nbits (e);

  plain_poly_make_monic (fp, p);

  /* initialize g to x */
  plain_poly_set_linear (g, 1, a);

  ASSERT (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      plain_poly_sqr (h, g, p);             /* h <- g^2 */
      if (e & (1L << k))
        plain_poly_mul_x (h, a, p);            /* h <- x*h */

      plain_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      plain_poly_swap (g, h);            /* h is already reduced mod p */
    }
}

/* g <- g^e mod (fp, p), using auxiliary polynomial h
   Assumes deg(g) < deg(fp) on input.
   Assumes (p-1)+deg(fp)*(p-1)^2 fits in a long [for plain_poly_div_r].
 */
static void
plain_poly_general_powmod_ui (plain_poly_t g, plain_poly_t fp, plain_poly_t h,
		     plain_poly_coeff_t e, plain_poly_coeff_t p)
{
  int k = nbits (e);
  plain_poly_t g_sav;

  plain_poly_init (g_sav, g->degree);
  plain_poly_set (g_sav, g);

  plain_poly_make_monic (fp, p);

  ASSERT (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      plain_poly_sqr (h, g, p);             /* h <- g^2 */
      plain_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      if (e & (1L << k))
        plain_poly_mul (h, h, g_sav, p);            /* h <- g_sav*h */
      plain_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      plain_poly_swap (g, h);            /* h is already reduced mod p */
    }
  plain_poly_clear(g_sav);
}
/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. 
   Return number of roots found (should be degree of f).
   Assumes p is odd, and deg(f) >= 1.
*/
static int
plain_poly_cantor_zassenhaus (plain_poly_coeff_t *r, plain_poly_t f,
                              plain_poly_coeff_t p, int depth)
{
  plain_poly_coeff_t a;
  plain_poly_t q, h, ff;
  int d = f->degree, dq, n, m;

  ASSERT (p & 1);
  ASSERT (d >= 1);

  if (d == 1) /* easy case: linear factor x + a ==> root -a mod p */
    {
      a = invert_si (-f->coeff[1], p); /* -1/f[1] mod p */
      r[0] = (a * f->coeff[0]) % p;
      return 1;
    }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the powering
     algorithm */
  plain_poly_init (q, 2 * d - 1);
  plain_poly_init (h, 2 * d - 1);
  plain_poly_init (ff, d);
  a = lrand48 () % p;
  for (;;)
    {
      plain_poly_powmod_ui (q, f, h, a, (p - 1) / 2, p);
      plain_poly_sub_1 (q, q, p);
      plain_poly_set (h, f);
      plain_poly_gcd (q, h, p);
      dq = q->degree;
      ASSERT (dq >= 0);
      if (0 < dq && dq < d)
	{
	  n = plain_poly_cantor_zassenhaus (r, q, p, depth + 1);
	  ASSERT (n == dq);
	  plain_poly_set (ff, f); /* plain_poly_divexact clobbers its 2nd arg */
	  plain_poly_divexact (h, ff, q, p);
	  m = plain_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
	  ASSERT (m == h->degree);
	  n += m;
          break;
	}
      if (++a >= p)
        a -= p;
    }
  plain_poly_clear (q);
  plain_poly_clear (h);
  plain_poly_clear (ff);
  return n;
}

static int plain_poly_coeff_cmp(
        const plain_poly_coeff_t * a,
        const plain_poly_coeff_t * b)
{
    if (*a < *b) return -1;
    if (*b < *a) return 1;
    return 0;
}


#define ROOTS_MOD_THRESHOLD  43 /* if only the number of roots is needed */
#define ROOTS_MOD_THRESHOLD2 67 /* if the roots are needed too */

/* The following function returns the number of roots of f(x) mod p,
   where f has degree d.
   If r is not NULL, put the n roots of f(x) mod p in r[0]...r[n-1].

   For p <= ROOTS_MOD_THRESHOLD, determine the number of roots of f mod p
   by evaluating f on 0, 1, ..., p-1. In theory, this threshold should 
   depend on the degree of f. The above thresholds seem to be optimal
   on a Pentium M. We must have ROOTS_MOD_THRESHOLD2 >= 2.

   The following ideas are from Emmanuel Thome:
   FIXME1: instead of first computing f1 = gcd(x^p-x, f) which is the
   product of linear factors, and then splitting f1, we can directly
   compute f1 = gcd(x^((p-1)/2)-1, f) and f2 = gcd(x^((p-1)/2)+1, f),
   assuming the roots 0 has been taken out.
   FIXME2: in the EDF, instead of dividing out f by gcd((x+a)^((p-1)/2)-1, f),
   we can simply compute gcd((x+a)^((p-1)/2)+1, f).
   FIXME3: instead of computing (x+a)^((p-1)/2), we can translate f by -a,
   and keep x^((p-1)/2) unchanged. [But then why not translate x^((p-1)/2),
   which has lower degree? Warning: if f has root -a, we might miss it.]
*/
int
plain_poly_roots (plain_poly_coeff_t *r, mpz_t *f, int d, const plain_poly_coeff_t p)
{
  plain_poly_t fp, g, h;
  int df;

  if (plain_poly_fits (d, p) == 0)
    {
      fprintf (stderr, "Error, (p-1)+deg(f)(p-1)^2 does not fits in a long\n");
      exit (EXIT_FAILURE);
    }

  /* the number of roots is the degree of gcd(x^p-x,f) */

  plain_poly_init (fp, d);
  d = plain_poly_set_mod (fp, f, d, p);
  /* d is the degree of fp (-1 if fp=0) */

  if (d <= 0)
    {
      df = 0;
      goto clear_fp;
    }

  ASSERT (d > 0);

  if ((r == NULL && p <= ROOTS_MOD_THRESHOLD) ||
      (r != NULL && p <= ROOTS_MOD_THRESHOLD2))
    {
      df = plain_poly_roots_naive (r, fp, p);
      goto clear_fp;
    }

  plain_poly_init (g, 2 * d - 1);
  plain_poly_init (h, 2 * d - 1);

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d: a square give at most degree
     2d-2, and multiplication by x gives 2d-1. */

  /* g <- x^p mod fp */
  plain_poly_powmod_ui (g, fp, h, 0, p, p);

  /* subtract x */
  plain_poly_sub_x (g, g, p);

  /* fp <- gcd (fp, g) */
  plain_poly_gcd (fp, g, p);

  /* now fp contains gcd(x^p-x, f) */

  df = fp->degree;
  ASSERT (df >= 0);
  
  /* If r is NULL, then we're not interested in finding the exact roots.
   * Only their number is important, and it's easily read as the degree
   * of the gcd.
   */
  if (r != NULL && df > 0)
    {
      int n MAYBE_UNUSED = plain_poly_cantor_zassenhaus (r, fp, p, 0);
      ASSERT (n == df);
    }

  plain_poly_clear (g);

  plain_poly_clear (h);

 clear_fp:
  plain_poly_clear (fp);

  /* Sort the roots */
  if (r && df)
  qsort(r, df, sizeof(plain_poly_coeff_t), (sortfunc_t) &plain_poly_coeff_cmp);

  return df;
}

int
plain_poly_roots_long (long *r, mpz_t *f, int d, const plain_poly_coeff_t p)
{
    plain_poly_coeff_t * pr;
    int i, n;

    pr = malloc(d * sizeof(plain_poly_coeff_t));
    n = plain_poly_roots(pr,f,d,p);
    for(i = 0 ; i < n ; i++) {
        r[i] = pr[i];
    }
    free(pr);
    return n;
}

int
plain_poly_roots_int64 (int64_t *r, mpz_t *f, int d, const plain_poly_coeff_t p)
{
    plain_poly_coeff_t * pr;
    int i, n;

    pr = malloc(d * sizeof(plain_poly_coeff_t));
    n = plain_poly_roots(pr,f,d,p);
    for(i = 0 ; i < n ; i++) {
        r[i] = pr[i];
    }
    free(pr);
    return n;
}

    

int plain_poly_is_irreducible(plain_poly_t fp, const plain_poly_coeff_t p) {
  plain_poly_t g, gmx, h;
  int d, i;

  plain_poly_make_monic (fp, p);
  d = fp->degree;

  if (d == 0)
    return 1;

  ASSERT (d > 0);

  plain_poly_init (g, 2 * d - 1);
  plain_poly_init (gmx, 2 * d - 1);
  plain_poly_init (h, 2 * d - 1);

  for (i = 1; 2*i <= d; ++i) {
    /* we first compute x^(p^i) mod fp; since fp has degree d, all operations can
       be done with polynomials of degree < 2d: a square give at most degree
       2d-2, and multiplication by x gives 2d-1. */

    if (i == 1) {
      /* g <- x^p mod fp */
      plain_poly_powmod_ui (g, fp, h, 0, p, p);
    } else {
      /* g <- g^p mod fp */
      plain_poly_general_powmod_ui (g, fp, h, p, p);
    }

    /* subtract x */
    plain_poly_sub_x (gmx, g, p);

    /* h <- gcd (fp, x^(p^i)-x) */
    plain_poly_set (h, fp);
    plain_poly_gcd (h, gmx, p);

    if (h->degree > 0)
      return 0;
  }

  plain_poly_clear (g);
  plain_poly_clear (gmx);
  plain_poly_clear (h);
  return 1;
}
