/* This implements Algorithm 1 from "Root Optimization of Polynomials in the
   Number Field Sieve" by Shi Bai, Richard P. Brent and Emmanuel Thomé,
   Mathematics of Computation, 2015.
   More precisely:
   * when the macro ORIGINAL is defined, it implements the original version
     described in the above reference;
   * when ORIGINAL is not defined (which is the default), it implements a
     variant which gives better results. Namely, when p^k is the largest
     power of p < B, the contribution is multiplied by p/(p-1) to take into
     account lifted roots mod p^(k+1), p^(k+2), ...
     On the RSA-768 polynomial, with B=V=W=200 and ORIGINAL defined, we get a
     maximal difference of 0.55 between the affine alpha-value and the
     computed estimation, and an average difference of 0.067, for all
     polynomials of the [-200,200]^2 grid.
     With ORIGINAL undefined, we get a maximal difference of 0.42 only,
     and an average of 0.0045 only.
*/

#include "cado.h" // IWYU pragma: keep

#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gmp.h>

#include "area.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "cado_poly.h"
#include "gcd.h"
#include "gmp_aux.h"    // ulong_isprime
#include "macros.h"
#include "mpz_poly.h"
#include "murphyE.h"
#include "omp_proxy.h"
#include "polyselect_alpha.h"
#include "polyselect_norms.h"
#include "size_optimization.h"
#include "timing.h"     // seconds

/* define ORIGINAL if you want the original algorithm from the paper */
// #define ORIGINAL

// #define TRACE_V 7
// #define TRACE_W 3

/* The algorithm is very sensitive to GUARD_ALPHA: with GUARD_ALPHA=1.0,
   almost 87% of the time is spent checking potential records.
   With GUARD_ALPHA=0.5, only about 20% of the time is spent for that. */
#define GUARD_ALPHA 0.5

/* global variables */
int verbose = 0;                /* verbosity level */
long *Primes, nprimes;          /* primes less than B */
long *Q;                        /* largest p^k < B */
long bestu = 0, bestv = 0;      /* current best rotation */
mpz_t bestw;                    /* current best rotation in w */
double best_alpha = DBL_MAX;    /* alpha of best rotation */
double best_E = 0;              /* E of best rotation (with -E) */
double tot_pols = 0;            /* number of sieved polynomial */
long u0 = 0, v0 = 0, w0 = 0;    /* initial translation */
int optimizeE = 0;              /* if not zero, optimize E instead of alpha */
double guard_alpha = 0.0;       /* guard when -E */
long mod = 0;                   /* consider class of u,v,w % mod, 0 = undef */
double tot_alpha = 0;           /* sum of alpha's */
int keep = 10;                  /* number of best classes kept */
double effort = DBL_MAX;        /* total effort */

typedef struct sieve_data {
  uint16_t q;
  uint16_t s;
  float nu;
} sieve_data;

static void
usage_and_die (const char *argv0)
{
  fprintf (stderr, "usage: %s [-area a] [-I n] [-Bf b] [-Bg c] [-margin x] [-v] [-sopt] [-V vmax] [-W wmax] [-B bbb] poly\n", argv0);
  fprintf (stderr, "  poly: filename of polynomial\n");
  fprintf (stderr, "  -area a      area parameter for Murphy-E computation (default %.2e)\n", AREA);
  fprintf (stderr, "  -I nnn       I-value for Murphy-E computation (overrides area)\n");
  fprintf (stderr, "  -Bf b        Bf bound for Murphy-E computation (default %.2e)\n", BOUND_F);
  fprintf (stderr, "  -Bg c        Bg bound for Murphy-E computation (default %.2e)\n", BOUND_G);
  fprintf (stderr, "  -margin x    allows a lognorm increase of x (default %.2f)\n", NORM_MARGIN);
  fprintf (stderr, "  -v           verbose toggle\n");
  fprintf (stderr, "  -sopt        first size-optimize the given polynomial\n");
  fprintf (stderr, "  -B nnn       parameter for alpha computation (default %d)\n", ALPHA_BOUND);
  fprintf (stderr, "  -E           optimize E instead of alpha\n");
  exit (1);
}

/* return x mod m, with 0 <= x < m */
static long
get_mod (long x, long m)
{
  x = x % m;
  return (x >= 0) ? x : x + m;
}

static unsigned long
initPrimes (unsigned long B)
{
  unsigned long nprimes = 0, p, q, l;

  Primes = malloc (B * sizeof (long));
  ASSERT_ALWAYS(Primes != NULL);
  for (p = 2; p < B; p += 1 + (p > 2))
    if (ulong_isprime (p))
      Primes[nprimes++] = p;
  CHECKED_REALLOC(Primes, nprimes, long);

  /* compute prime powers */
  Q = malloc (nprimes * sizeof (long));
  ASSERT_ALWAYS(Q != NULL);
  for (l = 0; l < nprimes; l++)
    {
      p = Primes[l];
      for (q = p; q * p < B; q *= p);
      Q[l] = q;
    }

  return nprimes;
}

/* Put in roots[0], roots[1], ... the roots of f + g * w = 0 mod q,
   and return the number of roots.
   Assume 0 <= f, g < q. */
static unsigned long
get_roots (unsigned long *roots, unsigned long f, unsigned long g,
           unsigned long q)
{
  unsigned long nroots = 0;
  unsigned long h = gcd_ul (g, q);
  if (h == 1) /* only one root, namely -f/g mod q */
    {
      unsigned long invg = invert_ul (g, q);
      roots[0] = get_mod (-f * invg, q);
      nroots = 1;
    }
  else if ((f % h) != 0)
    nroots = 0;
  else
    {
      f /= h;
      g /= h;
      q /= h;
      unsigned long invg = invert_ul (g, q);
      roots[0] = get_mod (-f * invg, q);
      for (unsigned long j = 1; j < h; j++)
        roots[j] = roots[j-1] + q;
      nroots = h;
    }
  return nroots;
}

#define TRIES 10

/* Return the average value of alpha in the class (v,w) = (modv,modw) % mod.
   Assume q is a prime power. */
static double
average_alpha (cado_poly_srcptr poly0, long modv, long modw, long q, gmp_randstate_ptr rstate)
{
  cado_poly poly;
  double s = 0.0, alpha;
  long v0 = 0, w0 = 0, p, t;

  /* check if q is a prime power */
  ASSERT_ALWAYS (q >= 2);
  for (p = 2; p * p <= q && q % p != 0; p += 1 + (p > 2));
  if (q % p != 0)
    p = q; /* q is prime */
  /* now p is the smallest prime factor of q */
  for (t = q; t % p == 0; t = t / p);
  ASSERT_ALWAYS (t == 1);

  /* first make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  /* first rotate by modv*x+modw */
  rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, modv, 1);
  v0 = modv;
  rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, modw, 0);
  w0 = modw;

  for (long j = 0; j < TRIES; j++)
    {
      rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], v0, q * j + modv, 1);
      v0 = q * j + modv;
      for (long k = 0; k < TRIES; k++)
        {
          rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], w0, q*k + modw, 0);
          w0 = q * k + modw;
          alpha = get_alpha_affine_p (poly->pols[ALG_SIDE], p, rstate);
          s += alpha;
        }
    }

  cado_poly_clear (poly);
  return s / pow ((double) TRIES, 2.0);
}

/* return 1/p mod q */
static long
invert (long p, long q)
{
  for (long t = 1; t < q; t++)
    if ((t * p) % q == 1)
      return t;
  ASSERT_ALWAYS(0);
}

/* Return c such that c = a mod p and c = b mod q.
   Assume invp = 1/p mod q. */
static long
crt (long a, long b, long p, long q, long invp)
{
  /* assume c = a + t*p, then t = (b-a)/p mod q */
  long t = (b - a) % q;
  if (t < 0)
    t += q;
  t = (t * invp) % q;
  return a + t * p;
}

typedef struct
{
  long vmod, wmod;
  double alpha;
} class;

/* Insert alpha into c, where c has already n entries (maximum is keep).
   Return the new value of n. */
static int
insert_class (class *c, int n, int keep, double alpha, long v, long w,
	      long vmin, long vmax, long mod)
{
  int i;

  /* check if this class has at least one representative in [vmin,vmax] */
  long t = get_mod (v - vmin, mod);
  if (vmin + t > vmax)
    return n; /* no representative in [vmin,vmax] */

  /* if alpha exceeds the best alpha value + guard_alpha, then it cannot
     yield A[j] < best_alpha + guard_alpha in rotate_v(). Note that this
     remains true after crt: if mod=mod1*mod2, and alpha1 > best_alpha1
     + guard_alpha, then since alpha = alpha1 + alpha2, then
     alpha > best_alpha1 + best_alpha2 + guard_alpha */
  if (n > 0 && alpha > c[0].alpha + guard_alpha)
    return n;

  for (i = n; i > 0 && alpha < c[i-1].alpha; i--)
    c[i] = c[i-1];
  /* now i = 0 or alpha >= c[i-1].alpha */
  if (i < keep)
    {
      c[i].vmod = v;
      c[i].wmod = w;
      c[i].alpha = alpha;
    }
  n += (n < keep);
  return n;
}

/* Return the (at most keep) best classes (v,w) mod 'mod'.
   Put in *nc the number of returned classes. */
static class*
best_classes (cado_poly poly0, long mod, int keep, long vmin, long vmax,
              int *nc, long u, gmp_randstate_ptr rstate)
{
  int nfactors = 0;
  unsigned long *factors = NULL, p;
  class *c, *d, *e;
  int nd, ne;
  int i;
  cado_poly poly;
  long q, Q = 1;

  if (mod == 1)
    {
      c = malloc (sizeof (class));
      c[0].vmod = c[0].wmod = 0;
      c[0].alpha = 0; /* value does not matter */
      *nc = 1;
      return c;
    }

  *nc = 0;
  /* first determine the prime factors of mod */
  for (long t = mod, p = 2; t != 1; p += 1 + (p & 1))
    {
      if ((t % p) == 0)
        {
          nfactors ++;
          CHECKED_REALLOC(factors, nfactors, unsigned long);
          q = 1;
          while ((t % p) == 0)
            {
              t /= p;
              q *= p;
            }
          factors[nfactors - 1] = q;
        }
    }

  /* make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  c = malloc ((keep + 1) * sizeof (class));
  ASSERT_ALWAYS(c != NULL);
  d = malloc ((keep + 1) * sizeof (class));
  ASSERT_ALWAYS(d != NULL);
  e = malloc ((keep + 1) * sizeof (class));
  ASSERT_ALWAYS(e != NULL);

  for (i = 0; i < nfactors; i++)
    {
      nd = 0; /* number of elements in 'd' */
      q = factors[i];
      /* determine p such that q=p^k */
      for (p = 2; q % p; p++);
      for (long v = 0; v < q; v++)
        {
          for (long w = 0; w < q; w++)
            {
              double alpha;
              alpha = average_alpha (poly, v, w, q, rstate);
              nd = insert_class (d, nd, keep, alpha, v, w, vmin, vmax, q);
            }
        }
      if (i == 0) /* copy d into c */
        {
          memcpy (c, d, nd * sizeof (class));
          *nc = nd;
        }
      else /* merge c and d into e */
        {
          long inv = invert (Q, q);
          ne = 0;
          for (int ic = 0; ic < *nc; ic++)
            for (int id = 0; id < nd; id++)
              {
                double alpha = c[ic].alpha + d[id].alpha;
		/* if alpha is larger (i.e., worse) than the last element,
		   since d[] is sorted by increasing values of alpha, we
		   assume all further values will be worse */
		if (ne == keep && e[ne-1].alpha < alpha)
		  break;
                long v = crt (c[ic].vmod, d[id].vmod, Q, q, inv);
                long w = crt (c[ic].wmod, d[id].wmod, Q, q, inv);
                ne = insert_class (e, ne, keep, alpha, v, w, vmin, vmax, Q * q);
              }
          /* copy back e to c */
          memcpy (c, e, ne * sizeof (class));
          *nc = ne;
        }
      Q *= q;
    }

  /* if u = -u0, check the class of the initial polynomial (-v0,-w0) */
  if (u == -u0)
    {
      int included = -1;
      for (i = 0; i < *nc; i++)
        if (get_mod (-v0, mod) == c[i].vmod && get_mod (-w0, mod) == c[i].wmod)
          included = i;
      if (included >= 0)
        printf ("class of initial polynomial has rank %d (%.2f)\n",
                included, c[included].alpha);
      else
        {
          double alpha = 0;
          for (int j = 0; j < nfactors; j++)
            alpha += average_alpha (poly, get_mod (-v0, factors[j]),
                                    get_mod (-w0, factors[j]), factors[j], rstate);
          printf ("class of initial polynomial is not included");
          if (*nc > 0)
            printf (" (last %.2f wrt %.2f)\n", c[*nc - 1].alpha, alpha);
          else
            printf (" (%.2f)\n", alpha);
        }
    }

  cado_poly_clear (poly);
  free (factors);
  free (d);
  free (e);
  return c;
}

/* rotation for a fixed value of v */
static void
rotate_v (cado_poly_srcptr poly0, long v, long B,
          double maxlognorm, double Bf, double Bg, double area, long u,
          long modw)
{
  long w, wmin, wmax;
  cado_poly poly;
  long l;
  mpz_t wminz, wmaxz;
  double tot_pols_local = 0;
  double tot_alpha_local = 0;

  mpz_init (wminz);
  mpz_init (wmaxz);

  /* first make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  /* compute f + (v*x)*g */
  rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, v, 1);

  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0,
                   maxlognorm, poly->skew);
  mpz_set_d (wminz, r.kmin);
  mpz_set_d (wmaxz, r.kmax);

  if (verbose > 1)
#pragma omp critical
    {
      gmp_printf ("v=%ld: wmin=%Zd wmax=%Zd\n", v, wminz, wmaxz);
      fflush (stdout);
    }

  /* Ensure wminz % mod = modw. Since mod <= MAX_LONG, we have
     s := wminz % mod < MAX_LONG thus modw - s fits in a long and
     there is no underflow below. */
  long t = get_mod (modw - mpz_fdiv_ui (wminz, mod), mod);
  ASSERT_ALWAYS(0 <= t && t < mod);
  mpz_add_ui (wminz, wminz, t);

  if (mpz_cmp (wminz, wmaxz) >= 0)
    goto end;

  /* if mod != 1, we have f + (k*mod+modw)*g = (f+modw*g) + k*(mod*g) */
  if (mod > 1)
    {
      rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, modw, 0); /* f <- f+modw*g */
      mpz_poly_mul_si (poly->pols[RAT_SIDE], poly->pols[RAT_SIDE], mod);
      /* wmin -> (wmin - modw) / mod */
      mpz_sub_ui (wminz, wminz, modw);
      ASSERT_ALWAYS(mpz_divisible_ui_p (wminz, mod));
      mpz_divexact_ui (wminz, wminz, mod);
      /* wmax -> (wmax - modw) / mod */
      mpz_sub_ui (wmaxz, wmaxz, modw);
      mpz_cdiv_q_ui (wmaxz, wmaxz, mod);
    }

  /* compute the expected value sum(log(p)/(p-1), p < B) and the
     sum of largest prime powers sum(p^floor(log(B-1)/log(p)), p < B) */
  double expected = 0.0;
  unsigned long sum_of_prime_powers = 0;
  for (l = 0; l < nprimes; l++)
    {
      long p = Primes[l];
      expected += log ((double) p) / (double) (p - 1);
      sum_of_prime_powers += Q[l];
    }

  ASSERT_ALWAYS (mpz_fits_slong_p (wmaxz));
  ASSERT_ALWAYS (mpz_fits_slong_p (wminz));

  wmax = mpz_get_si (wmaxz);
  wmin = mpz_get_si (wminz);

  mpz_t ump;
  mpz_init (ump);
  unsigned long *roots = malloc (B * sizeof (long));
  ASSERT_ALWAYS(roots != NULL);
  float *L;
  double nu;
  L = malloc (B * sizeof (float));

  /* sieve data: sieve_q contains the largest p^k < B for each prime p,
     it thus fits in an uint16_t if B does;
     sieve_s contains the first index i multiple of q where the contribution
     sieve_nu should be added, it is thus smaller than q and thus fits too. */
  sieve_data *sieve_d = malloc (sum_of_prime_powers * sizeof (sieve_data));
  unsigned long sieve_n = 0;
  for (l = 0; l < nprimes; l++)
    {
      long p = Primes[l], s, t, q;
      double logp = log ((double) p);
      memset (L, 0, B * sizeof (float));
      for (q = p; q <= Q[l]; q *= p)
        {
          /* the contribution is log(p)/p^(k-1)/(p+1) when the exponent k
             is not the largest one, and log(p)/p^(k-1)/(p+1)*p/(p-1) for the
             largest exponent k */
          nu = logp / (double) q * (double) p / (double) (p + 1);
#ifndef ORIGINAL
          if (q == Q[l])
            nu *= (double) p / (double) (p - 1);
#endif
          for (long x = 0; x < q; x++)
            {
              /* compute f(x) and g(x) mod p^k, where q = p^k */
              unsigned long fx, gx;
              mpz_poly_eval_ui (ump, poly->pols[ALG_SIDE], x);
              fx = mpz_fdiv_ui (ump, q);
              mpz_poly_eval_ui (ump, poly->pols[RAT_SIDE], x);
              gx = mpz_fdiv_ui (ump, q);
              /* search roots w of fx + w*gx = 0 mod q */
              unsigned long nroots = get_roots (roots, fx, gx, q);
              for (unsigned long i = 0; i < nroots; i++)
                {
                  long w = roots[i];
                  /* update for w+t*q */
                  for (t = 0; t < Q[l] / q; t++)
                    L[w + t * q] += nu;
                }
            }
        }

      /* prepare data for the sieve */
      q = Q[l];
      for (w = 0; w < q; w++)
        {
          nu = L[w];
          if (nu == 0.0)
            continue;
          /* compute s = w+t*q-wmin such that s - q < 0 <= s, i.e.,
             (t-1)*q < wmin-w <= t*q: t = ceil((wmin-w)/q) */
          if (wmin - w < 0)
            t = (wmin - w) / q;
          else
            t = (wmin - w + q - 1) / q;
          s = w + t * q - wmin;
          sieve_d[sieve_n].q = q;
          sieve_d[sieve_n].s = s;
          sieve_d[sieve_n].nu = nu;
          sieve_n ++;
        }
    }

  ASSERT_ALWAYS(sieve_n <= sum_of_prime_powers);

#define LEN (1<<14) /* length of the sieve array */

  float *A = malloc (LEN * sizeof (float));
  ASSERT_ALWAYS(A != NULL);

  /* we sieve by chunks of LEN cells at a time */
  long wcur = wmin;
  while (wcur < wmax)
    {
      /* A[j] corresponds to w = wcur + j */
      for (long j = 0; j < LEN; j++)
        A[j] = expected;

#if defined(TRACE_V) && defined(TRACE_W)
      if (v == TRACE_V && (mod * wcur + modw <= TRACE_W &&
                           TRACE_W < mod * (wcur + LEN) + modw))
        printf ("initialized A[%d] to %f\n", TRACE_W, A[(TRACE_W - modw) / mod - wcur]);
#endif

      /* now perform the sieve */
      for (unsigned long i = 0; i < sieve_n; i++)
        {
          long q = sieve_d[i].q;
          long s = sieve_d[i].s;
          float nu = sieve_d[i].nu;
          /* if mod=1, a given value of s corresponds to w = wcur + s
             if mod>1, then s corresponds to w = mod*(wcur + s) + modw */
          while (s < LEN)
            {
#if defined(TRACE_V) && defined(TRACE_W)
              if (v == TRACE_V && TRACE_W == mod * (wcur + s) + modw)
                printf ("q=%ld: update A[%d] from %f to %f\n",
                        q, TRACE_W, A[s], A[s] - nu);
#endif
              A[s] -= nu;
              s += q;
            }
          sieve_d[i].s = s - LEN;
        }

      /* check for the smallest A[s] */
      /* if wcur + LEN > wmax, we check only wmax - wcur entries */
      long maxj = (wcur + LEN <= wmax) ? LEN : wmax - wcur;
      for (long j = 0; j < maxj; j++)
        {
          tot_alpha_local += A[j];
          tot_pols_local += 1;
          /* print alpha and E of original polynomial */
          if (u == -u0 && v == -v0 && mod * (wcur + j) + modw == -w0)
            {
              w = wcur + j; /* local value of w, the global one is mod * w + modw */
              rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, w, 0);
              double skew = poly->skew; /* save skewness */
              poly->skew = L2_skewness (poly->pols[ALG_SIDE]);
              double lognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew);
              /* to compute E, we need to divide g by mod */
              mpz_poly_divexact_ui (poly->pols[RAT_SIDE], poly->pols[RAT_SIDE], mod);
              double E = MurphyE (poly, Bf, Bg, area, MURPHY_K, get_alpha_bound ());
              /* restore g */
              mpz_poly_mul_si (poly->pols[RAT_SIDE], poly->pols[RAT_SIDE], mod);
              /* this can only occur for one thread, thus no need to put
                 #pragma omp critical */
              gmp_printf ("u=%ld v=%ld w=%ld lognorm=%.2f est_alpha_aff=%.2f E=%.2e [original]\n",
                          u, v, mod * w + modw, lognorm, A[j], E);
              fflush (stdout);
              /* restore the original polynomial (w=0) and skewness */
              poly->skew = skew;
              rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], w, 0, 0);
            }
          if (A[j] < best_alpha + guard_alpha)
            {
              w = wcur + j;
              /* compute E */
              rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0, w, 0);
              double skew = poly->skew; /* save skewness */
              poly->skew = L2_skewness (poly->pols[ALG_SIDE]);
              double lognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew);
              /* to compute E, we need to divide g by mod */
              mpz_poly_divexact_ui (poly->pols[RAT_SIDE], poly->pols[RAT_SIDE], mod);
              double E = MurphyE (poly, Bf, Bg, area, MURPHY_K, get_alpha_bound ());
              /* restore g */
              mpz_poly_mul_si (poly->pols[RAT_SIDE], poly->pols[RAT_SIDE], mod);
              /* restore the original polynomial (w=0) and skewness */
              poly->skew = skew;
              rotate_aux (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], w, 0, 0);

              if (optimizeE == 0 || (optimizeE == 1 && E > best_E))
#pragma omp critical
                {
                  bestu = u;
                  bestv = v;
                  mpz_set_si (bestw, mod * w + modw);
                  best_alpha = (double) A[j];
                  best_E = E;
                  gmp_printf ("u=%ld v=%ld w=%Zd lognorm=%.2f est_alpha_aff=%.2f E=%.2e\n",
                              u, v, bestw, lognorm, best_alpha, E);
                  fflush (stdout);
                }
            }
        }

#if defined(TRACE_V) && defined(TRACE_W)
      if (v == TRACE_V && (mod * wcur + modw <= TRACE_W &&
                           TRACE_W < mod * (wcur + LEN) + modw))
        printf ("A[%d] = %f\n", TRACE_W, A[(TRACE_W - modw) / mod - wcur]);
#endif

      wcur += LEN;
    }

  free (sieve_d);
  free (L);
  free (roots);
  mpz_clear (ump);

  free (A);

 end:
  cado_poly_clear (poly);
  mpz_clear (wminz);
  mpz_clear (wmaxz);

  /* accumulate the number of polynomials and the alpha values */
#pragma omp critical
  {
    tot_pols += tot_pols_local;
    tot_alpha += tot_alpha_local;
  }
}

static void
rotate (cado_poly cpoly, long B, double maxlognorm, double Bf, double Bg,
        double area, long u, gmp_randstate_ptr rstate)
{
  /* determine range [vmin,vmax] */
  rotation_space r;
  expected_growth (&r, cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 1,
                   maxlognorm, cpoly->skew);
  long vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  long vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
  if (verbose)
    {
      printf ("u=%ld: vmin=%ld vmax=%ld\n", u, vmin, vmax);
      fflush (stdout);
    }

  class *c = NULL;
  int n;
  c = best_classes (cpoly, mod, keep, vmin, vmax, &n, u, rstate);
  ASSERT_ALWAYS (n <= keep);
  printf ("u=%ld: kept %d class(es)", u, n);
  if (n == 0)
    printf ("\n");
  else
    printf (" with alpha from %.2f to %.2f\n", c[0].alpha, c[n-1].alpha);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++)
    {
      if (verbose > 1)
#pragma omp critical
        printf ("u=%ld, class i=%d: -mod %ld -modv %ld -modw %ld %.2f\n",
                u0 + u, i, mod, get_mod (v0 + c[i].vmod, mod),
                get_mod (w0 + c[i].wmod, mod), c[i].alpha);

      long vmin1 = vmin;

      /* first ensure that vmin1 = modv % mod */
      long t = get_mod (c[i].vmod - vmin1, mod);
      ASSERT_ALWAYS(0 <= t && t < mod);
      vmin1 += t;
      ASSERT_ALWAYS(get_mod (vmin1, mod) == c[i].vmod);
      for (long v = vmin1; v <= vmax; v += mod)
        rotate_v (cpoly, v, B, maxlognorm, Bf, Bg, area, u, c[i].wmod);
    }

  free (c);
}

/* don't modify poly, which is the size-optimized polynomial
   (poly0 is the initial polynomial) */
static void
print_transformation (cado_poly_ptr poly0, cado_poly_srcptr cpoly)
{
  mpz_t k;
  int d = poly0->pols[ALG_SIDE]->deg;

  mpz_init (k);
  /* first compute the translation k: g(x+k) = g1*x + g1*k + g0 */
  mpz_sub (k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 0), mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 0));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1));
  gmp_printf ("translation %Zd, ", k);

  mpz_poly_translation(poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], k);
  mpz_poly_translation(poly0->pols[RAT_SIDE], poly0->pols[RAT_SIDE], k);

  /* size_optimization might multiply f0 by some integer t */
  ASSERT_ALWAYS(mpz_divisible_p (mpz_poly_coeff_const(poly0->pols[ALG_SIDE], d),
				 mpz_poly_coeff_const(poly0->pols[ALG_SIDE], d)));
  mpz_divexact (k, mpz_poly_coeff_const(poly0->pols[ALG_SIDE], d),
		mpz_poly_coeff_const(poly0->pols[ALG_SIDE], d));
  if (mpz_cmp_ui (k, 1) != 0)
    {
      gmp_printf ("multiplier %Zd, ", k);
      mpz_poly_mul_mpz (poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], k);
    }
  /* now compute rotation by x^2 */
  mpz_sub (k, mpz_poly_coeff_const(cpoly->pols[ALG_SIDE], 3), mpz_poly_coeff_const(poly0->pols[ALG_SIDE], 3));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1));
  gmp_printf ("rotation [%Zd,", k);
  ASSERT (mpz_fits_slong_p (k));
  u0 = mpz_get_si (k);
  mpz_poly_rotation(poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], poly0->pols[RAT_SIDE], k, 2);
  mpz_sub (k, mpz_poly_coeff_const(cpoly->pols[ALG_SIDE], 2), mpz_poly_coeff_const(poly0->pols[ALG_SIDE], 2));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1));
  gmp_printf ("%Zd,", k);
  ASSERT (mpz_fits_slong_p (k));
  v0 = mpz_get_si (k);
  mpz_poly_rotation(poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], poly0->pols[RAT_SIDE], k, 1);
  mpz_sub (k, mpz_poly_coeff_const(cpoly->pols[ALG_SIDE], 1), mpz_poly_coeff_const(poly0->pols[ALG_SIDE], 1));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0->pols[RAT_SIDE], 1));
  gmp_printf ("%Zd]\n", k);
  ASSERT (mpz_fits_slong_p (k));
  w0 = mpz_get_si (k);
  mpz_poly_rotation(poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], poly0->pols[RAT_SIDE], k, 0);
  ASSERT_ALWAYS(mpz_cmp (mpz_poly_coeff_const(poly0->pols[ALG_SIDE], 0),
                         mpz_poly_coeff_const(cpoly->pols[ALG_SIDE], 0)) == 0);
  mpz_clear (k);
}

double
rotate_area_v (cado_poly_srcptr poly0, double maxlognorm, long v)
{
  cado_poly cpoly;
  double area;

  cado_poly_init (cpoly);
  cado_poly_set (cpoly, poly0);
  rotate_aux (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 0, v, 1);
  rotation_space r;
  expected_growth (&r, cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 0,
                   maxlognorm, cpoly->skew);
  area = r.kmax - r.kmin;
  cado_poly_clear (cpoly);
  return area;
}

/* estimate the rootsieve area for a given u */
double
rotate_area_u (cado_poly_srcptr poly0, double maxlognorm, long u)
{
  double area, sum = 0.0;
  long h, vmin, vmax;
  cado_poly cpoly;

  cado_poly_init (cpoly);
  cado_poly_set (cpoly, poly0);
  rotate_aux (cpoly->pols[ALG_SIDE],
	      cpoly->pols[RAT_SIDE], 0, u, 2);
  rotation_space r;
  expected_growth (&r, cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 1,
                   maxlognorm, cpoly->skew);
  vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
#define SAMPLE 100
  if (vmax / SAMPLE - vmin / SAMPLE > 1)
    h = vmax / SAMPLE - vmin / SAMPLE;
  else
    h = 1;

  for (long v = vmin; v <= vmax; v += h)
    {
      area = rotate_area_v (cpoly, maxlognorm, v);
      sum += area;
    }
  cado_poly_clear (cpoly);
  return sum * h;
}

/* estimate the rootsieve area for umin <= u <= umax */
double
rotate_area (cado_poly_srcptr cpoly, double maxlognorm, long umin, long umax)
{
  double area, sum = 0.0;

  for (long u = umin; u <= umax; u++)
    {
      area = rotate_area_u (cpoly, maxlognorm, u);
      sum += area;
    }
  return sum;
}

/* Given a sieving area, a maximal effort, and a value of keep,
   compute the best 'mod' value. */
long
best_mod (double area, double maxeffort, double keep)
{
  long l[] = {1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360, 720720, 12252240,
              232792560, 5354228880, 26771144400, 80313433200, 2329089562800};
  int i = 0;
  double e;
  do {
    mod = l[i];
    /* The number of polynomials sieved is approximately area/mod^2*keep.
       Note: there is a bias when vmax-vmin is smaller than mod, since we
       only keep classes that contain at least an element in [vmin, vmax],
       thus the probability is larger than (vmax-vmin)/mod. */
    e = area / (double) mod / (double) mod * (double) keep;
    if (e <= maxeffort)
      break;
    i += 1;
  }
  while (i < 17);
  printf ("using mod = %ld, effort = %.2e\n", mod, e);
  return mod;
}

int main(int argc, char const * argv[])
{
    int argc0 = argc;
    char const **argv0 = argv;
    cado_poly cpoly;
    int I = 0;
    double margin = NORM_MARGIN;
    long umin = LONG_MIN, umax = LONG_MAX;
    int sopt = 0;
    double time = seconds ();
    long B = ALPHA_BOUND;
    gmp_randstate_t rstate;

    gmp_randinit_default(rstate);

    while (argc >= 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-area") == 0)
          {
            area = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-I") == 0)
          {
            I = atoi (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bf") == 0)
          {
            bound_f = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bg") == 0)
          {
            bound_g = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-margin") == 0)
          {
            margin = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-effort") == 0)
          {
            effort = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-v") == 0)
          {
            verbose ++;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-E") == 0)
          {
            optimizeE = 1;
            guard_alpha = GUARD_ALPHA;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-sopt") == 0)
          {
            sopt = 1;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-B") == 0)
          {
            unsigned long B = strtol (argv [2], NULL, 10);
            set_alpha_bound (B);
            argv += 2;
            argc -= 2;
          }
        /* use http://oeis.org/A051451 for -mod:
           1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360, 720720, 12252240,
           232792560, 5354228880, 26771144400, 80313433200, 2329089562800,
           72201776446800, 144403552893600 */
        else if (strcmp (argv[1], "-mod") == 0)
          {
            mod = strtol (argv [2], NULL, 10);
            ASSERT_ALWAYS(mod >= 1);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-umin") == 0)
          {
            umin = strtol (argv [2], NULL, 10);
            ASSERT_ALWAYS(umin > LONG_MIN); /* LONG_MIN is reserved */
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-umax") == 0)
          {
            umax = strtol (argv [2], NULL, 10);
            ASSERT_ALWAYS(umax < LONG_MAX); /* LONG_MAX is reserved */
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-keep") == 0)
          {
            keep = atoi (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else
          break;
      }
    if (argc != 2)
        usage_and_die (argv[0]);

#pragma omp parallel
#pragma omp master
    printf ("# Using %d thread(s)\n", omp_get_num_threads ());

    ASSERT_ALWAYS(B <= 65536);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    cado_poly_init (cpoly);
    if (!cado_poly_read (cpoly, argv[1]))
      {
        fprintf (stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die (argv[0]);
      }

    if (cpoly->skew == 0.0)
      cpoly->skew = L2_skewness (cpoly->pols[ALG_SIDE]);

    /* if -sopt, size-optimize */
    if (sopt)
      {
        cado_poly c;
        cado_poly_init (c);
        cado_poly_set (c, cpoly);
        size_optimization (c->pols[ALG_SIDE], c->pols[RAT_SIDE],
                           cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE],
                           SOPT_DEFAULT_EFFORT, verbose);
        printf ("# initial polynomial:\n");
        cado_poly_fprintf (stdout, "", cpoly);
        print_transformation (cpoly, c);
        cado_poly_set (cpoly, c);
        printf ("# size-optimized polynomial:\n");
        cado_poly_fprintf (stdout, "", cpoly);
        cado_poly_clear (c);
      }

    nprimes = initPrimes (B);

    /* compute the skewness */
    cpoly->skew = L2_skewness (cpoly->pols[ALG_SIDE]);
    double lognorm = L2_lognorm (cpoly->pols[ALG_SIDE], cpoly->skew);
    double maxlognorm = lognorm + margin;
    printf ("initial lognorm %.2f, maxlognorm %.2f\n", lognorm, maxlognorm);

    /* determine range [umin,umax] */
    rotation_space r;
    expected_growth (&r, cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 2,
                     maxlognorm, cpoly->skew);
    if (umin == LONG_MIN) /* umin was not given by the user */
      umin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
    if (umax == LONG_MAX) /* umax was not given by the user */
      umax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
    if (verbose)
      printf ("umin=%ld umax=%ld\n", umin, umax);

    if (mod == 0) /* compute best 'mod' for given effort */
      {
        double sieving_area = rotate_area (cpoly, maxlognorm, umin, umax);
        /* print total sieving area */
        printf ("sieving area %.2e\n", sieving_area);
        mod = best_mod (sieving_area, effort, keep);
      }

    mpz_init (bestw);
    long u0 = 0; /* current translation in u */
    for (long u = umin; u <= umax; u++)
      {
        rotate_aux (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], u0, u, 2);
        u0 = u;

        rotate (cpoly, B, maxlognorm, bound_f, bound_g, area, u, rstate);
      }

    /* restore original polynomial */
    rotate_aux (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], u0, 0, 2);

    /* perform the best rotation */
    gmp_printf ("best rotation: u=%ld v=%ld w=%Zd alpha=%1.2f\n",
            bestu, bestv, bestw, best_alpha);

    /* perform the best rotation */
    rotate_aux (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 0, bestu, 2);
    rotate_aux (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 0, bestv, 1);
    mpz_poly_rotation(cpoly->pols[ALG_SIDE], cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], bestw, 0);

    /* recompute the skewness of the best polynomial */
    cpoly->skew = L2_combined_skewness2 (cpoly->pols[0], cpoly->pols[1]);

    print_cadopoly_extra (stdout, cpoly, argc0, argv0, 0);

    time = seconds () - time;
    printf ("# Sieved %.2e polynomials in %.2f seconds (%.2es/p)\n",
            tot_pols, time, time / tot_pols);
    printf ("# Average alpha %.2f\n", tot_alpha / tot_pols);

    free (Primes);
    free (Q);
    cado_poly_clear (cpoly);
    mpz_clear (bestw);

    gmp_randclear(rstate);

    return 0;
}
