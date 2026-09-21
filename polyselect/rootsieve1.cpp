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

#include <cfloat>
#include <climits>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <gmp.h>

#include "area.hpp"
#include "auxiliary.hpp" /* for common routines with polyselect_old.c */
#include "cado_poly.hpp"
#include "gcd.h"
#include "gmp_aux.h"
#include "macros.h"
#include "mpz_poly.h"
#include "murphyE.hpp"
#include "omp_proxy.h"
#include "polyselect_alpha.h"
#include "polyselect_norms.hpp"
#include "size_optimization.hpp"
#include "timing.h"
#include "macros.h"
#include "polyselect_alpha.h"

#include <memory>

#include "utils_cxx.hpp"
#include "verbose.hpp"
#include "cado_main.hpp"


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
std::unique_ptr<long[]> Primes; /* primes less than B */
long nprimes;
std::unique_ptr<long[]> Q;      /* largest p^k < B */
long bestu = 0, bestv = 0;      /* current best rotation */
mpz_t bestw;                    /* current best rotation in w */
double best_alpha = DBL_MAX;    /* alpha of best rotation */
double best_E = 0;              /* E of best rotation (with -E) */
double tot_pols = 0;            /* number of sieved polynomial */
long u0 = 0, v0 = 0, w0 = 0;    /* initial translation */
int optimizeE = 0;              /* if not zero, optimize E instead of alpha */
double guard_alpha = 0.0;       /* guard when -E */
long mod = 0;                   /* consider congruence class of u,v,w % mod, 0 = undef */
double tot_alpha = 0;           /* sum of alpha's */
int keep = 10;                  /* number of best congruences kept */
double effort = DBL_MAX;        /* total effort */

typedef struct sieve_data {
  uint16_t q;
  uint16_t s;
  float nu;
} sieve_data;

static void
declare_usage (cxx_param_list & pl)
{
  pl.declare_usage("poly", "filename of polynomial (may be given positionally)");
  pl.declare_usage("area", fmt::format("area parameter for Murphy-E computation (default {:.2e})", AREA));
  pl.declare_usage("I", "I-value for Murphy-E computation (overrides area)");
  pl.declare_usage("Bf", fmt::format("Bf bound for Murphy-E computation (default {:.2e})", BOUND_F));
  pl.declare_usage("Bg", fmt::format("Bg bound for Murphy-E computation (default {:.2e})", BOUND_G));
  pl.declare_usage("margin", fmt::format("allows a lognorm increase of x (default {:.2f})", NORM_MARGIN));
  pl.declare_usage("effort", "total effort");
  pl.declare_usage("B", fmt::format("parameter for alpha computation (default {})", ALPHA_BOUND));
  /* use http://oeis.org/A051451 for -mod:
     1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360, 720720, 12252240,
     232792560, 5354228880, 26771144400, 80313433200, 2329089562800,
     72201776446800, 144403552893600 */
  pl.declare_usage("mod", "consider congruence classes of u,v,w mod this");
  pl.declare_usage("umin", "lower end of the rotation range");
  pl.declare_usage("umax", "upper end of the rotation range");
  pl.declare_usage("keep", "number of best congruences kept");
  pl.declare_usage("v", "(switch) verbose mode, repeat for more");
  pl.declare_usage("sopt", "(switch) first size-optimize the given polynomial");
  pl.declare_usage("E", "(switch) optimize E instead of alpha");
  verbose_decl_usage(pl);
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

  /* Count first, so that the array comes out at exactly the right size
   * and needs no shrinking afterwards. B is at most 65536 here, so the
   * extra primality sweep costs nothing worth measuring. */
  for (p = 2; p < B; p += 1 + (p > 2))
    if (ulong_isprime (p))
      nprimes++;

  Primes = std::make_unique<long[]>(nprimes);
  {
    unsigned long k = 0;
    for (p = 2; p < B; p += 1 + (p > 2))
      if (ulong_isprime (p))
        Primes[k++] = p;
    ASSERT_ALWAYS(k == nprimes);
  }

  /* compute prime powers */
  Q = std::make_unique<long[]>(nprimes);
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

/* Return the average value of alpha in the congruence (v,w) = (modv,modw) % mod.
   Assume q is a prime power. */
static double
average_alpha (cxx_cado_poly const & poly0, long modv, long modw, long q, gmp_randstate_ptr rstate)
{
  cxx_cado_poly poly;
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
  poly = poly0;

  /* first rotate by modv*x+modw */
  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, modv, 1);
  v0 = modv;
  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, modw, 0);
  w0 = modw;

  for (long j = 0; j < TRIES; j++)
    {
      rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], v0, q * j + modv, 1);
      v0 = q * j + modv;
      for (long k = 0; k < TRIES; k++)
        {
          rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], w0, q*k + modw, 0);
          w0 = q * k + modw;
          alpha = get_alpha_affine_p (poly[ALG_SIDE], p, rstate);
          s += alpha;
        }
    }

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
} congruence;

/* Insert alpha into c, where c has already n entries (maximum is keep).
   Return the new value of n. */
static int
insert_congruence (congruence *c, int n, int keep, double alpha, long v, long w,
	      long vmin, long vmax, long mod)
{
  int i;

  /* check if this congruence has at least one representative in [vmin,vmax] */
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

/* Return the (at most keep) best congruences (v,w) mod 'mod'.
   Put in *nc the number of returned congruences. */
static congruence*
best_congruences (cxx_cado_poly const & poly0, long mod, int keep, long vmin, long vmax,
              int *nc, long u, gmp_randstate_ptr rstate)
{
  int nfactors = 0;
  unsigned long *factors = NULL, p;
  congruence *c, *d, *e;
  int nd, ne;
  int i;
  cxx_cado_poly poly;
  long q, Q = 1;

  if (mod == 1)
    {
      c = new congruence[1];
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
          checked_realloc(factors, nfactors);
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
  poly = poly0;

  c = new congruence[keep + 1];
  d = new congruence[keep + 1];
  e = new congruence[keep + 1];

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
              nd = insert_congruence (d, nd, keep, alpha, v, w, vmin, vmax, q);
            }
        }
      if (i == 0) /* copy d into c */
        {
          memcpy (c, d, nd * sizeof (congruence));
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
                ne = insert_congruence (e, ne, keep, alpha, v, w, vmin, vmax, Q * q);
              }
          /* copy back e to c */
          memcpy (c, e, ne * sizeof (congruence));
          *nc = ne;
        }
      Q *= q;
    }

  /* if u = -u0, check the congruence of the initial polynomial (-v0,-w0) */
  if (u == -u0)
    {
      int included = -1;
      for (i = 0; i < *nc; i++)
        if (get_mod (-v0, mod) == c[i].vmod && get_mod (-w0, mod) == c[i].wmod)
          included = i;
      if (included >= 0)
        printf ("congruence of initial polynomial has rank %d (%.2f)\n",
                included, c[included].alpha);
      else
        {
          double alpha = 0;
          for (int j = 0; j < nfactors; j++)
            alpha += average_alpha (poly, get_mod (-v0, factors[j]),
                                    get_mod (-w0, factors[j]), factors[j], rstate);
          printf ("congruence of initial polynomial is not included");
          if (*nc > 0)
            printf (" (last %.2f wrt %.2f)\n", c[*nc - 1].alpha, alpha);
          else
            printf (" (%.2f)\n", alpha);
        }
    }

  free (factors);
  delete[] d;
  delete[] e;
  return c;
}

/* rotation for a fixed value of v */
static void
rotate_v (cxx_cado_poly const & poly0, long v, long B,
          double maxlognorm, double Bf, double Bg, double area, long u,
          long modw)
{
  long w, wmin, wmax;
  cxx_cado_poly poly;
  long l;
  mpz_t wminz, wmaxz;
  double tot_pols_local = 0;
  double tot_alpha_local = 0;

  mpz_init (wminz);
  mpz_init (wmaxz);

  /* first make a local copy of the original polynomial */
  poly = poly0;

  /* compute f + (v*x)*g */
  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, v, 1);

  rotation_space r;
  expected_growth (&r, poly[ALG_SIDE], poly[RAT_SIDE], 0,
                   maxlognorm, poly.skew);
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

  if (mpz_cmp (wminz, wmaxz) < 0) {

      /* if mod != 1, we have f + (k*mod+modw)*g = (f+modw*g) + k*(mod*g) */
      if (mod > 1)
      {
          rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, modw, 0); /* f <- f+modw*g */
          mpz_poly_mul_si (poly[RAT_SIDE], poly[RAT_SIDE], mod);
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
      auto * roots = new unsigned long[B];
      double nu;
      auto * L = new float[B];

      /* sieve data: sieve_q contains the largest p^k < B for each prime p,
         it thus fits in an uint16_t if B does;
         sieve_s contains the first index i multiple of q where the contribution
         sieve_nu should be added, it is thus smaller than q and thus fits too. */
      auto *sieve_d = new sieve_data[sum_of_prime_powers];
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
                  mpz_poly_eval_ui (ump, poly[ALG_SIDE], x);
                  fx = mpz_fdiv_ui (ump, q);
                  mpz_poly_eval_ui (ump, poly[RAT_SIDE], x);
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

      auto * A = new float[LEN];

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
                  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, w, 0);
                  double skew = poly.skew; /* save skewness */
                  poly.skew = L2_skewness (poly[ALG_SIDE]);
                  double lognorm = L2_lognorm (poly[ALG_SIDE], poly.skew);
                  /* to compute E, we need to divide g by mod */
                  mpz_poly_divexact_ui (poly[RAT_SIDE], poly[RAT_SIDE], mod);
                  double E = MurphyE (poly, Bf, Bg, area, MURPHY_K, get_alpha_bound ());
                  /* restore g */
                  mpz_poly_mul_si (poly[RAT_SIDE], poly[RAT_SIDE], mod);
                  /* this can only occur for one thread, thus no need to put
#pragma omp critical */
                  gmp_printf ("u=%ld v=%ld w=%ld lognorm=%.2f est_alpha_aff=%.2f E=%.2e [original]\n",
                          u, v, mod * w + modw, lognorm, A[j], E);
                  fflush (stdout);
                  /* restore the original polynomial (w=0) and skewness */
                  poly.skew = skew;
                  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], w, 0, 0);
              }
              if (A[j] < best_alpha + guard_alpha)
              {
                  w = wcur + j;
                  /* compute E */
                  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], 0, w, 0);
                  double skew = poly.skew; /* save skewness */
                  poly.skew = L2_skewness (poly[ALG_SIDE]);
                  double lognorm = L2_lognorm (poly[ALG_SIDE], poly.skew);
                  /* to compute E, we need to divide g by mod */
                  mpz_poly_divexact_ui (poly[RAT_SIDE], poly[RAT_SIDE], mod);
                  double E = MurphyE (poly, Bf, Bg, area, MURPHY_K, get_alpha_bound ());
                  /* restore g */
                  mpz_poly_mul_si (poly[RAT_SIDE], poly[RAT_SIDE], mod);
                  /* restore the original polynomial (w=0) and skewness */
                  poly.skew = skew;
                  rotate_aux (poly[ALG_SIDE], poly[RAT_SIDE], w, 0, 0);

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

      delete[] A;

      delete[] sieve_d;
      delete[] L;
      delete[] roots;
      mpz_clear (ump);
  }

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
rotate (cxx_cado_poly & cpoly, long B, double maxlognorm, double Bf, double Bg,
        double area, long u, gmp_randstate_ptr rstate)
{
  /* determine range [vmin,vmax] */
  rotation_space r;
  expected_growth (&r, cpoly[ALG_SIDE], cpoly[RAT_SIDE], 1,
                   maxlognorm, cpoly.skew);
  long vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  long vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
  if (verbose)
    {
      printf ("u=%ld: vmin=%ld vmax=%ld\n", u, vmin, vmax);
      fflush (stdout);
    }

  int n;
  congruence * c = best_congruences (cpoly, mod, keep, vmin, vmax, &n, u, rstate);
  ASSERT_ALWAYS (n <= keep);
  printf ("u=%ld: kept %d congruence(s)", u, n);
  if (n == 0)
    printf ("\n");
  else
    printf (" with alpha from %.2f to %.2f\n", c[0].alpha, c[n-1].alpha);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++)
    {
      if (verbose > 1)
#pragma omp critical
        printf ("u=%ld, congruence i=%d: -mod %ld -modv %ld -modw %ld %.2f\n",
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

  delete[] c;
}

/* don't modify poly, which is the size-optimized polynomial
   (poly0 is the initial polynomial) */
static void
print_transformation (cxx_cado_poly & poly0, cxx_cado_poly const & cpoly)
{
  mpz_t k;
  int d = poly0[ALG_SIDE]->deg;

  mpz_init (k);
  /* first compute the translation k: g(x+k) = g1*x + g1*k + g0 */
  mpz_sub (k, mpz_poly_coeff_const(cpoly[RAT_SIDE], 0), mpz_poly_coeff_const(poly0[RAT_SIDE], 0));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1));
  gmp_printf ("translation %Zd, ", k);

  mpz_poly_translation(poly0[ALG_SIDE], poly0[ALG_SIDE], k);
  mpz_poly_translation(poly0[RAT_SIDE], poly0[RAT_SIDE], k);

  /* size_optimization might multiply f0 by some integer t */
  ASSERT_ALWAYS(mpz_divisible_p (mpz_poly_coeff_const(cpoly[ALG_SIDE], d),
				 mpz_poly_coeff_const(poly0[ALG_SIDE], d)));
  mpz_divexact (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], d),
		mpz_poly_coeff_const(poly0[ALG_SIDE], d));
  if (mpz_cmp_ui (k, 1) != 0)
    {
      gmp_printf ("multiplier %Zd, ", k);
      mpz_poly_mul_mpz (poly0[ALG_SIDE], poly0[ALG_SIDE], k);
    }
  /* now compute rotation by x^2 */
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 3), mpz_poly_coeff_const(poly0[ALG_SIDE], 3));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1));
  gmp_printf ("rotation [%Zd,", k);
  ASSERT (mpz_fits_slong_p (k));
  u0 = mpz_get_si (k);
  mpz_poly_rotation(poly0[ALG_SIDE], poly0[ALG_SIDE], poly0[RAT_SIDE], k, 2);
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 2), mpz_poly_coeff_const(poly0[ALG_SIDE], 2));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1));
  gmp_printf ("%Zd,", k);
  ASSERT (mpz_fits_slong_p (k));
  v0 = mpz_get_si (k);
  mpz_poly_rotation(poly0[ALG_SIDE], poly0[ALG_SIDE], poly0[RAT_SIDE], k, 1);
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 1), mpz_poly_coeff_const(poly0[ALG_SIDE], 1));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly0[RAT_SIDE], 1));
  gmp_printf ("%Zd]\n", k);
  ASSERT (mpz_fits_slong_p (k));
  w0 = mpz_get_si (k);
  mpz_poly_rotation(poly0[ALG_SIDE], poly0[ALG_SIDE], poly0[RAT_SIDE], k, 0);
  ASSERT_ALWAYS(mpz_cmp (mpz_poly_coeff_const(poly0[ALG_SIDE], 0),
                         mpz_poly_coeff_const(cpoly[ALG_SIDE], 0)) == 0);
  mpz_clear (k);
}

double
rotate_area_v (cxx_cado_poly const & poly0, double maxlognorm, long v)
{
  double area;

  cxx_cado_poly cpoly = poly0;
  rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, v, 1);
  rotation_space r;
  expected_growth (&r, cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0,
                   maxlognorm, cpoly.skew);
  area = r.kmax - r.kmin;
  return area;
}

/* estimate the rootsieve area for a given u */
double
rotate_area_u (cxx_cado_poly const & poly0, double maxlognorm, long u)
{
  double area, sum = 0.0;
  long h, vmin, vmax;

  cxx_cado_poly cpoly = poly0;

  rotate_aux (cpoly[ALG_SIDE],
	      cpoly[RAT_SIDE], 0, u, 2);
  rotation_space r;
  expected_growth (&r, cpoly[ALG_SIDE], cpoly[RAT_SIDE], 1,
                   maxlognorm, cpoly.skew);
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
  return sum * h;
}

/* estimate the rootsieve area for umin <= u <= umax */
double
rotate_area (cxx_cado_poly const & cpoly, double maxlognorm, long umin, long umax)
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
       only keep congruences that contain at least an element in [vmin, vmax],
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

static int main_(int argc, char const * argv[]);

int main(int argc, char const * argv[])
{
    return cado::main_wrapper(main_, argc, argv);
}

static int main_(int argc, char const * argv[])
{
    int argc0 = argc;
    char const **argv0 = argv;
    cxx_cado_poly cpoly;
    int I = 0;
    double margin = NORM_MARGIN;
    long umin = LONG_MIN, umax = LONG_MAX;
    int sopt = 0;
    double time = seconds ();
    long B = ALPHA_BOUND;
    gmp_randstate_t rstate;

    gmp_randinit_default(rstate);

    cxx_param_list pl;
    const char * polyfilename = NULL;

    declare_usage(pl);
    pl.configure_switch_old("-v", &verbose);
    pl.configure_switch_old("-sopt", &sopt);
    pl.configure_switch_old("-E", &optimizeE);

    if (argc == 1)
        pl.fail("Error, a polynomial file is mandatory");

    argv++, argc--;
    for (int wild = 0 ; argc ; ) {
        if (pl.update_cmdline(argc, argv)) continue;
        if (wild == 0 && argv[0][0] != '-') {
            polyfilename = argv[0];
            argc--, argv++, wild++;
            continue;
        }
        pl.fail("Unhandled parameter {}", argv[0]);
    }

    pl.parse("area", area);
    pl.parse("I", I);
    pl.parse("Bf", bound_f);
    pl.parse("Bg", bound_g);
    pl.parse("margin", margin);
    pl.parse("effort", effort);
    pl.parse("keep", keep);
    if (pl.parse("mod", mod) && mod < 1)
        pl.fail("Error, -mod must be at least 1");
    /* LONG_MIN and LONG_MAX are reserved to mean "the user said nothing",
     * so they are only rejected when the option is actually given. */
    if (pl.parse("umin", umin) && umin == LONG_MIN)
        pl.fail("Error, -umin is out of range");
    if (pl.parse("umax", umax) && umax == LONG_MAX)
        pl.fail("Error, -umax is out of range");
    if (pl.parse("B", B))
        set_alpha_bound(B);
    if (!polyfilename)
        polyfilename = pl.lookup_old("poly");

    if (optimizeE)
        guard_alpha = GUARD_ALPHA;

    if (pl.warn_unused())
        pl.fail("unexpected parameter(s) on the command line");
    verbose_interpret_parameters(pl);
    pl.print_command_line(stdout);

    if (!polyfilename)
        pl.fail("Error, a polynomial file is mandatory");

#pragma omp parallel
#if defined(_OPENMP) && _OPENMP >= 202011
    #pragma omp masked
#else
    #pragma omp master
#endif
    printf ("# Using %d thread(s)\n", omp_get_num_threads ());

    ASSERT_ALWAYS(B <= 65536);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    if (!cpoly.read(polyfilename))
        pl.fail("Problem when reading file {}", polyfilename);

    if (cpoly.skew == 0.0)
      cpoly.skew = L2_skewness (cpoly[ALG_SIDE]);

    /* if -sopt, size-optimize */
    if (sopt)
      {
        cxx_cado_poly c = cpoly;
        size_optimization (c[ALG_SIDE], c[RAT_SIDE],
                           cpoly[ALG_SIDE], cpoly[RAT_SIDE],
                           SOPT_DEFAULT_EFFORT, verbose);
        printf ("# initial polynomial:\n");
        cpoly.fprintf(stdout);
        print_transformation (cpoly, c);
        cpoly = c;
        printf ("# size-optimized polynomial:\n");
        cpoly.fprintf(stdout);
      }

    nprimes = initPrimes (B);

    /* compute the skewness */
    cpoly.skew = L2_skewness (cpoly[ALG_SIDE]);
    double lognorm = L2_lognorm (cpoly[ALG_SIDE], cpoly.skew);
    double maxlognorm = lognorm + margin;
    printf ("initial lognorm %.2f, maxlognorm %.2f\n", lognorm, maxlognorm);

    /* determine range [umin,umax] */
    rotation_space r;
    expected_growth (&r, cpoly[ALG_SIDE], cpoly[RAT_SIDE], 2,
                     maxlognorm, cpoly.skew);
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
        rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], u0, u, 2);
        u0 = u;

        rotate (cpoly, B, maxlognorm, bound_f, bound_g, area, u, rstate);
      }

    /* restore original polynomial */
    rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], u0, 0, 2);

    /* perform the best rotation */
    gmp_printf ("best rotation: u=%ld v=%ld w=%Zd alpha=%1.2f\n",
            bestu, bestv, bestw, best_alpha);

    /* perform the best rotation */
    rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, bestu, 2);
    rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, bestv, 1);
    mpz_poly_rotation(cpoly[ALG_SIDE], cpoly[ALG_SIDE], cpoly[RAT_SIDE], bestw, 0);

    /* recompute the skewness of the best polynomial */
    cpoly.skew = L2_combined_skewness2 (cpoly[0], cpoly[1]);

    print_cadopoly_extra (stdout, cpoly, argc0, argv0, 0);

    time = seconds () - time;
    printf ("# Sieved %.2e polynomials in %.2f seconds (%.2es/p)\n",
            tot_pols, time, time / tot_pols);
    printf ("# Average alpha %.2f\n", tot_alpha / tot_pols);

    mpz_clear (bestw);

    gmp_randclear(rstate);

    return 0;
}
