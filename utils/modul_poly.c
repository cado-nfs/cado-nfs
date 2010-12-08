/* arithmetic on polynomials over Z/pZ, with coefficients represented by
   the same type as in mod_ul.[ch] ; most presumably an unsigned long */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "mod_ul.h"
#include "modul_poly.h"
#include "gmp_aux.h"

/* allocate an array of d coefficients, and initialize it */
static residueul_t*
alloc_long_array (int d)
{
  residueul_t *f;

  f = (residueul_t*) malloc (d * sizeof (residueul_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* reallocate an array to d coefficients */
static residueul_t*
realloc_long_array (residueul_t *f, int d)
{
  f = (residueul_t*) realloc (f, d * sizeof (residueul_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
  return f;
}

/* initialize a polynomial with maximal degree d */
void
modul_poly_init (modul_poly_t f, int d)
{
  ASSERT (d >= 0);
  f->degree = -1; /* initialize to 0 */
  f->coeff = alloc_long_array (d + 1);
  f->alloc = d + 1;
}

/* clear a polynomial */
void
modul_poly_clear (modul_poly_t f)
{
  free (f->coeff);
  f->alloc = 0;
}

/* realloc f to (at least) n coefficients */
void
modul_poly_realloc (modul_poly_t f, int n)
{
  if (f->alloc < n)
    {
      f->coeff = realloc_long_array (f->coeff, n);
      f->alloc = n;
    }
}

/* f <- g */
void
modul_poly_set (modul_poly_t f, const modul_poly_t g, modulusul_t p)
{
  int i, d = g->degree;

  modul_poly_realloc (f, d + 1);
  f->degree = d;
  for (i = 0; i <= d; i++)
    modul_set(f->coeff[i], g->coeff[i], p);
}

/* f <- f/lc(f) mod p */
void
modul_poly_make_monic (modul_poly_t f, modulusul_t p)
{
  int d = f->degree, i;
  residueul_t ilc;

  if (modul_is1(f->coeff[d],p))
    return;

  modul_inv(ilc, f->coeff[d], p);
  for (i = 0; i < d; i++)
      modul_mul(f->coeff[i], f->coeff[i], ilc, p);
  modul_set1(f->coeff[d], p);
}

/* fp <- f/lc(f) mod p. Return degree of fp (-1 if fp=0). */
int
modul_poly_set_mod (modul_poly_t fp, mpz_t *f, int d, modulusul_t p)
{
  d = modul_poly_set_mod_raw (fp, f, d, p);
  modul_poly_make_monic (fp, p);

  return d;
}

/* fp <- f mod p. Return degree of fp (-1 if fp=0). */
int
modul_poly_set_mod_raw (modul_poly_t fp, mpz_t *f, int d, modulusul_t p)
{
  int i;

  while (d >= 0 && mpz_divisible_ui_p (f[d], modul_getmod_ul (p)))
    d --;
  ASSERT (d >= 0); /* f is 0 mod p: should not happen in the CADO-NFS context
                      since otherwise p would divide N, indeed f(m)=N */
  modul_poly_realloc (fp, d + 1);
  fp->degree = d;
  for (i = 0; i <= d; i++)
    modul_set_ul (fp->coeff[i], mpz_fdiv_ui (f[i], modul_getmod_ul (p)), p);

  return d;
}

/* f <- a*x+b, a <> 0 */
void
modul_poly_set_linear (modul_poly_t f, residueul_t a, residueul_t b, modulusul_t p)
{
  modul_poly_realloc (f, 2);
  f->degree = 1;
  modul_set(f->coeff[1], a, p);
  modul_set(f->coeff[0], b, p);
}

/* swap f and g */
void
modul_poly_swap (modul_poly_t f, modul_poly_t g)
{
  int i;
  residueul_t *t;

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

/* h <- f*g mod p */
void
modul_poly_mul(modul_poly_t h, const modul_poly_t f, const modul_poly_t g, modulusul_t p) {
  int df = f->degree;
  int dg = g->degree;
  int i, j;
  modul_poly_t res;
  residueul_t aux;

  modul_poly_init(res, df+dg);
  for (i = 0; i <= df+dg; ++i)
    modul_set0(res->coeff[i], p);
  for (i = 0; i <= df; ++i)
    for (j = 0; j <= dg; ++j) {
        modul_mul(aux, f->coeff[i], g->coeff[j], p);
        modul_add(res->coeff[i+j], res->coeff[i+j], aux, p);
    }
  res->degree=df+dg;
  modul_poly_set(h, res, p);
  modul_poly_clear(res);
}


/* h <- g^2 mod p, g and h must differ */
void
modul_poly_sqr (modul_poly_t h, const modul_poly_t g, modulusul_t p)
{
  int i, j, dg = g->degree;
  residueul_t *gc, *hc;
  residueul_t aux;

  ASSERT (dg >= -1);
  if (dg == -1) /* g is zero */
    {
      h->degree = -1;
      return;
    }
  modul_poly_realloc (h, 2 * dg + 1);
  gc = g->coeff;
  hc = h->coeff;
  for (i = 0; i <= dg; i++)
    for (j = i + 1; j <= dg; j++)
      if (i == 0 || j == dg) 
          modul_mul(hc[i + j], gc[i], gc[j], p);
      else {
          modul_mul(aux, gc[i], gc[j], p);
          modul_add(hc[i + j], hc[i + j], aux, p);
      }
  for (i = 1; i < 2 * dg; i++)
      modul_add(hc[i], hc[i], hc[i], p);
  modul_mul(hc[0], gc[0], gc[0], p);
  modul_mul(hc[2*dg], gc[dg], gc[dg], p);
  for (i = 1; i < dg; i++) {
      modul_mul(aux, gc[i], gc[i], p);
      modul_add(hc[2*i], hc[2*i], aux, p);
  }
  h->degree = 2 * dg;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
static void
modul_poly_normalize (modul_poly_t h, modulusul_t p)
{
  int dh = h->degree;

  while (dh >= 0 && modul_is0(h->coeff[dh],p))
    dh --;
  h->degree = dh;
}

/* h <- (x+a)*h mod p */
void
modul_poly_mul_x (modul_poly_t h, residueul_t a, modulusul_t p)
{
  int i, d = h->degree;
  residueul_t *hc;
  residueul_t aux;

  modul_poly_realloc (h, d + 2); /* (x+a)*h has degree d+1, thus d+2 coeffs */
  hc = h->coeff;
  modul_set(hc[d + 1], hc[d], p);
  for (i = d - 1; i >= 0; i--) {
      modul_mul(aux, a, hc[i+1], p);
      modul_add(hc[i+1], aux, hc[i], p);
  }
  modul_mul(hc[0], hc[0], a, p);
  h->degree = d + 1;
}

/* h <- g - x mod p */
void
modul_poly_sub_x (modul_poly_t h, const modul_poly_t g, modulusul_t p)
{
  int i, d = g->degree;

  /* if g has degree d >= 2, then g-x has degree d too;
     otherwise g-x has degree <= 1 */
  modul_poly_realloc (h, ((d < 2) ? 1 : d) + 1);
  for (i = 0; i <= d; i++)
    modul_set(h->coeff[i], g->coeff[i], p);
  for (i = d + 1; i <= 1; i++)
    modul_set0(h->coeff[i], p);
  modul_sub_ul(h->coeff[1], h->coeff[1], 1, p);
  h->degree = (d < 1) ? 1 : d;
  modul_poly_normalize (h, p);
}

/* h <- g - 1 mod p */
void
modul_poly_sub_1 (modul_poly_t h, const modul_poly_t g, modulusul_t p)
{
  int i, d = g->degree;

  /* g-1 has degree d if d >= 1, and degree 0 otherwise */
  modul_poly_realloc (h, ((d < 1) ? 0 : d) + 1);
  for (i = 0; i <= d; i++)
    modul_set(h->coeff[i], g->coeff[i], p);
  for (i = d + 1; i <= 0; i++)
    modul_set0(h->coeff[i], p);
  modul_sub_ul(h->coeff[0], h->coeff[0], 1, p);
  h->degree = (d < 0) ? 0 : d;
  modul_poly_normalize (h, p);
}

/* h <- rem(h, f) mod p, f not necessarily monic */
void
modul_poly_div_r (modul_poly_t h, const modul_poly_t f, modulusul_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  residueul_t *hc = h->coeff, t;
  residueul_t aux;
  modul_set1(t, p);

  monic = modul_is1(f->coeff[d], p);
  if (!monic)
      modul_inv(t, f->coeff[d], p); /* 1/f[d] mod p */
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
          modul_mul(hc[dh], hc[dh], t, p);
      for (i = 0; i < d; i++) {
          modul_mul(aux, hc[dh], f->coeff[i], p);
          modul_sub(hc[dh - d + i], hc[dh - d + i], aux, p);
      }
      do {
	  dh --;
      } while (dh >= 0 && modul_is0(hc[dh],p));
      h->degree = dh;
    }
}

/* q <- divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
void
modul_poly_divexact (modul_poly_t q, modul_poly_t h, const modul_poly_t f, modulusul_t p)
{
  int i, d = f->degree, dh = h->degree, monic;
  residueul_t *hc = h->coeff, t;
  residueul_t aux;
  modul_set1(t, p);

  ASSERT (d >= 0);
  ASSERT (dh >= 0);
  ASSERT (dh >= d);
  modul_poly_realloc (q, dh + 1 - d);
  q->degree = dh - d;
  monic = modul_is1(f->coeff[d], p);
  if (!monic)
      modul_inv(t, f->coeff[d], p);
  while (dh >= d)
    {
      /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
      if (!monic)
          modul_mul(hc[dh], hc[dh], t, p);
      modul_set(q->coeff[dh - d], hc[dh], p);
      /* we only need to update the coefficients of degree >= d of h,
	 i.e., we want i >= 2d - dh */
      for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
          modul_mul(aux, hc[dh], f->coeff[i], p);
          modul_sub(hc[dh - d + i], hc[dh - d + i], aux, p);
      }
      dh --;
    }
  modul_poly_normalize (q, p);
}

/* fp <- gcd (fp, g), clobbers g */
void
modul_poly_gcd (modul_poly_t fp, modul_poly_t g, modulusul_t p)
{
  while (g->degree >= 0)
    {
      modul_poly_div_r (fp, g, p);
//      modul_poly_mod_ui (fp, fp, p);
      /* now deg(fp) < deg(g): swap f and g */
      modul_poly_swap (fp, g);
    }
}

void
modul_poly_out (FILE *fp, modul_poly_t f, modulusul_t p)
{
  if (f->degree < 0)
    fprintf (fp, "0\n");
  else
    {
      int i;
      for (i = 0; i <= f->degree; i++)
	if (!modul_is0(f->coeff[i],p))
	  gmp_fprintf (fp, "+%Nu*x^%d", f->coeff[i], MODUL_SIZE, i);
      fprintf (fp, ";\n");
    }
}

/* returns f(x) mod p */
void modul_poly_eval (residueul_t r, modul_poly_t f, residueul_t x, modulusul_t p)
{
  int i, d = f->degree;
  residueul_t aux;
  
  ASSERT (d >= 0);
  modul_set(r, f->coeff[d], p);
  for (i = d - 1; i >= 0; i--) {
      modul_mul(aux, r, x, p);
      modul_add(r, aux, f->coeff[i], p);
  }
}

/* Return the number n of roots of f mod p using a naive algorithm.
   If r is not NULL, put the roots in r[0], ..., r[n-1]. */
int
modul_poly_roots_naive (residueul_t *r, modul_poly_t f, modulusul_t p)
{
  int n = 0;
  residueul_t x;
  modul_set0(x, p);
  for ( ; !modul_finished(x, p) ; modul_next(x, p) ) {
      residueul_t v;
    modul_poly_eval (v, f, x, p);
    if (modul_is0(v, p))
      {
	if (r != NULL)
	  modul_set(r[n], x, p);
	n ++;
      }
  }
  return n;
}

/* g <- (x+a)^e mod (fp, p), using auxiliary polynomial h */
void
modul_poly_powmod_ui (modul_poly_t g, modul_poly_t fp, modul_poly_t h, residueul_t a,
		     unsigned long e, modulusul_t p)
{
  int k = nbits (e);

  modul_poly_make_monic (fp, p);

  residueul_t one;
  modul_set1(one, p);
  /* initialize g to x */
  modul_poly_set_linear (g, one, a, p);

  ASSERT (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      modul_poly_sqr (h, g, p);             /* h <- g^2 */
      if (e & (1UL << k))
        modul_poly_mul_x (h, a, p);            /* h <- x*h */

      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h, p);       /* g <- h  */
    }
}

/* g <- g^e mod (fp, p), using auxiliary polynomial h */
void
modul_poly_general_powmod_ui (modul_poly_t g, modul_poly_t fp, modul_poly_t h,
		     unsigned long e, modulusul_t p)
{
  int k = nbits (e);
  modul_poly_t g_sav;

  modul_poly_init(g_sav, g->degree);
  modul_poly_set(g_sav, g, p);

  modul_poly_make_monic (fp, p);

  ASSERT (e > 0);
  for (k -= 2; k >= 0; k--)
    {
      modul_poly_sqr (h, g, p);             /* h <- g^2 */
      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h, p);       /* g <- h */
      if (e & (1UL << k))
        modul_poly_mul (h, h, g_sav, p);            /* h <- g_sav*h */
      modul_poly_div_r (h, fp, p);       /* h -> rem(h, fp) */
      modul_poly_set (g, h, p);       /* g <- h */
    }
  modul_poly_clear(g_sav);
}
/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. 
   Return number of roots found (should be degree of f).
   Assumes p is odd, and deg(f) >= 1.
*/
int
modul_poly_cantor_zassenhaus (residueul_t *r, modul_poly_t f, modulusul_t p, int depth)
{
  residueul_t a;
  modul_poly_t q, h, ff;
  int d = f->degree, dq, n, m;

  ASSERT (p[0] & 1);
  ASSERT (d >= 1);

  if (d == 1) /* easy case: linear factor x + a ==> root -a mod p */
    {
        residueul_t aux;
        modul_neg(aux, f->coeff[1], p);
        modul_inv(a, aux, p);
        modul_mul(r[0], a, f->coeff[0], p);
        return 1;
    }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the powering
     algorithm */
  modul_poly_init (q, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);
  modul_poly_init (ff, d);
  modul_set_ul(a, lrand48(), p);
  for (;;)
    {
      modul_poly_powmod_ui (q, f, h, a, (modul_getmod_ul(p) - 1) / 2, p);
      modul_poly_sub_1 (q, q, p);
      modul_poly_set (h, f, p);
      modul_poly_gcd (q, h, p);
      dq = q->degree;
      ASSERT (dq >= 0);
      if (0 < dq && dq < d)
	{
	  n = modul_poly_cantor_zassenhaus (r, q, p, depth + 1);
	  ASSERT (n == dq);
	  modul_poly_set (ff, f, p); /* modul_poly_divexact clobbers its 2nd arg */
	  modul_poly_divexact (h, ff, q, p);
	  m = modul_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
	  ASSERT (m == h->degree);
	  n += m;
          break;
	}
      modul_next(a, p);
      if (modul_finished(a, p))
          modul_set0(a, p);
    }
  modul_poly_clear (q);
  modul_poly_clear (h);
  modul_poly_clear (ff);
  return n;
}

typedef int (*sortfunc_t) (const void *, const void *);

int coeff_cmp(
        const unsigned long * a,
        const unsigned long * b)
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
modul_poly_roots(residueul_t *r, mpz_t *f, int d, modulusul_t p)
{
  modul_poly_t fp, g, h;
  int df, n;

  /* the number of roots is the degree of gcd(x^p-x,f) */

  modul_poly_init (fp, d);
  d = modul_poly_set_mod (fp, f, d, p);

   /* d is the degree of fp (-1 if fp=0) */

  if (d == 0)
    {
      df = 0;
      goto clear_fp;
    }

  ASSERT (d > 0);

  if ((r == NULL && p[0] <= ROOTS_MOD_THRESHOLD) ||
      (r != NULL && p[0] <= ROOTS_MOD_THRESHOLD2))
    {
      df = modul_poly_roots_naive (r, fp, p);
      goto clear_fp;
    }

  modul_poly_init (g, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);

  /* we first compute x^p mod fp; since fp has degree d, all operations can
     be done with polynomials of degree < 2d: a square give at most degree
     2d-2, and multiplication by x gives 2d-1. */

  /* g <- x^p mod fp */
  residueul_t zero;
  modul_init(zero,p);
  modul_poly_powmod_ui (g, fp, h, zero, p[0], p);
  modul_clear(zero,p);


  /* subtract x */
  modul_poly_sub_x (g, g, p);

  /* fp <- gcd (fp, g) */
  modul_poly_gcd (fp, g, p);

  /* now fp contains gcd(x^p-x, f) */

  df = fp->degree;
  ASSERT (df >= 0);
  
  /* If r is NULL, then we're not interested in finding the exact roots.
   * Only their number is important, and it's easily read as the degree
   * of the gcd.
   */
  if (r != NULL && df > 0)
    {
      n = modul_poly_cantor_zassenhaus (r, fp, p, 0);
      ASSERT (n == df);
    }

  modul_poly_clear (g);

  modul_poly_clear (h);

 clear_fp:
  modul_poly_clear (fp);

  if (r && df)
  qsort(r, df, sizeof(unsigned long), (sortfunc_t) &coeff_cmp);

  return df;
}

int
modul_poly_roots_ulong (unsigned long *r, mpz_t *f, int d, modulusul_t p)
{
    residueul_t * pr;
    int i, n;

    pr = malloc(d * sizeof(residueul_t));
    n = modul_poly_roots(pr,f,d,p);
    for(i = 0 ; i < n ; i++) {
        r[i] = modul_getmod_ul (pr[i]);
    }
    free(pr);
    return n;
}

int
modul_poly_roots_int64 (int64_t *r, mpz_t *f, int d, modulusul_t p)
{
    residueul_t *pr;
    int i, n;

    pr = malloc(d * sizeof(residueul_t));
    n = modul_poly_roots(pr,f,d,p);
    for(i = 0 ; i < n ; i++) {
      r[i] = modul_getmod_ul (pr[i]);
    }
    free(pr);
    return n;
}


int modul_poly_is_irreducible(modul_poly_t fp, modulusul_t p)
{
  modul_poly_t g, gmx, h;
  int d, i;

  modul_poly_make_monic (fp, p);
  d = fp->degree;

  if (d == 0)
    return 1;

  ASSERT (d > 0);

  modul_poly_init (g, 2 * d - 1);
  modul_poly_init (gmx, 2 * d - 1);
  modul_poly_init (h, 2 * d - 1);

  for (i = 1; 2*i <= d; ++i) {
    /* we first compute x^(p^i) mod fp; since fp has degree d, all operations can
       be done with polynomials of degree < 2d: a square give at most degree
       2d-2, and multiplication by x gives 2d-1. */

    if (i == 1) {
      /* g <- x^p mod fp */
      modul_poly_powmod_ui (g, fp, h, 0, p[0], p);
    } else {
      /* g <- g^p mod fp */
      modul_poly_general_powmod_ui (g, fp, h, p[0], p);
    }

    /* subtract x */
    modul_poly_sub_x (gmx, g, p);

    /* h <- gcd (fp, x^(p^i)-x) */
    modul_poly_set(h, fp, p);
    modul_poly_gcd (h, gmx, p);

    if (h->degree > 0)
      return 0;
  }

  modul_poly_clear (g);
  modul_poly_clear (gmx);
  modul_poly_clear (h);
  return 1;
}
