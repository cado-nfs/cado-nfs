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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <vector>

#include <gmp.h>

#include "fmt/base.h"

#include "area.hpp"
#include "auxiliary.hpp" /* for common routines with polyselect_old.c */
#include "cado_main.hpp"
#include "cado_poly.hpp"
#include "cxx_mpz.hpp"
#include "gcd.h"
#include "gmp_aux.h"
#include "macros.h"
#include "mpz_poly.h"
#include "murphyE.hpp"
#include "omp_proxy.h"
#include "params.hpp"
#include "polyselect_alpha.h"
#include "polyselect_norms.hpp"
#include "size_optimization.hpp"
#include "timing.h"
#include "verbose.hpp"


/* define ORIGINAL if you want the original algorithm from the paper */
// #define ORIGINAL

// #define TRACE_V 7
// #define TRACE_W 3

/* The algorithm is very sensitive to GUARD_ALPHA: with GUARD_ALPHA=1.0,
   almost 87% of the time is spent checking potential records.
   With GUARD_ALPHA=0.5, only about 20% of the time is spent for that. */
constexpr double GUARD_ALPHA = 0.5;

/* number of (v,w) pairs tried per congruence in average_alpha() */
constexpr long TRIES = 10;

/* number of samples used to estimate the rootsieve area */
constexpr long SAMPLE = 100;

/* length of the sieve array */
constexpr long LEN = 1 << 14;

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

struct congruence {
  long vmod, wmod;
  double alpha;
};

/* sieve_q contains the largest p^k < B for each prime p, it thus fits in
   an uint16_t if B does; sieve_s contains the first index i multiple of q
   where the contribution sieve_nu should be added, it is thus smaller
   than q and thus fits too. */
struct sieve_data {
  uint16_t q;
  uint16_t s;
  float nu;
};

/* The root sieve: the parameters it was given, the tables it precomputes
   from them, and the record it keeps while it runs. Only the record is
   written to during the sieve, and only under omp critical. */
struct rootsieve {
  cxx_cado_poly poly;          /* rotated in place as u varies */
  int verbose = 0;             /* verbosity level */
  int optimizeE = 0;           /* if not zero, optimize E instead of alpha */
  int keep = 10;               /* number of best congruences kept */
  long B = ALPHA_BOUND;        /* alpha is computed with the primes < B */
  long mod = 0;                /* congruence class of u,v,w % mod, 0 = undef */
  double effort = DBL_MAX;     /* total effort */
  double guard_alpha = 0.0;    /* guard when -E */
  double maxlognorm = 0;       /* largest lognorm we accept */
  double Bf = 0, Bg = 0;       /* smoothness bounds, for Murphy-E */
  double area = 0;             /* sieving area, for Murphy-E */
  long u0 = 0, v0 = 0, w0 = 0; /* translation performed by -sopt */

  std::vector<long> primes;        /* the primes less than B */
  std::vector<long> prime_powers;  /* for each, its largest power < B */
  cxx_gmp_randstate rstate;

  /* best rotation seen so far, and the running totals */
  long bestu = 0, bestv = 0;
  cxx_mpz bestw;
  double best_alpha = DBL_MAX;
  double best_E = 0;
  double tot_pols = 0;         /* number of sieved polynomials */
  double tot_alpha = 0;        /* sum of alpha's */

  void init_primes ();
  void insert_congruence (std::vector<congruence> & c, double alpha, long v,
                          long w, long vmin, long vmax, long modulus) const;
  std::vector<congruence> best_congruences (long vmin, long vmax, long u);
  std::vector<sieve_data> prepare_sieve (cxx_cado_poly const & cpoly,
                                         long wmin) const;
  double murphyE_at (cxx_cado_poly & cpoly, long w, double & lognorm) const;
  void rotate_v (long v, long u, long modw);
  void rotate (long u);
  void print_transformation (cxx_cado_poly const & cpoly);
  double rotate_area (long umin, long umax) const;
  long best_mod (double sieving_area);
};

/* fill primes[] with the primes less than B, and prime_powers[] with the
   largest power of each that is still less than B */
void
rootsieve::init_primes ()
{
  for (unsigned long p = 2; p < (unsigned long) B; p += 1 + (p > 2))
    if (ulong_isprime (p))
      primes.push_back (p);

  prime_powers.reserve (primes.size ());
  for (long p : primes)
    {
      long q;
      for (q = p; q * p < B; q *= p);
      prime_powers.push_back (q);
    }
}

/* Insert alpha into c, which holds at most keep entries sorted by
   increasing alpha. */
void
rootsieve::insert_congruence (std::vector<congruence> & c, double alpha,
                              long v, long w, long vmin, long vmax,
                              long modulus) const
{
  /* check if this congruence has at least one representative in [vmin,vmax] */
  long t = get_mod (v - vmin, modulus);
  if (vmin + t > vmax)
    return; /* no representative in [vmin,vmax] */

  /* if alpha exceeds the best alpha value + guard_alpha, then it cannot
     yield A[j] < best_alpha + guard_alpha in rotate_v(). Note that this
     remains true after crt: if mod=mod1*mod2, and alpha1 > best_alpha1
     + guard_alpha, then since alpha = alpha1 + alpha2, then
     alpha > best_alpha1 + best_alpha2 + guard_alpha */
  if (!c.empty () && alpha > c.front ().alpha + guard_alpha)
    return;

  auto it = std::upper_bound (c.begin (), c.end (), alpha,
                              [](double a, congruence const & x) {
                                return a < x.alpha;
                              });
  if (it - c.begin () >= keep)
    return;
  c.insert (it, { v, w, alpha });
  if ((int) c.size () > keep)
    c.pop_back ();
}

/* Return the (at most keep) best congruences (v,w) mod 'mod'. */
std::vector<congruence>
rootsieve::best_congruences (long vmin, long vmax, long u)
{
  long q, Q = 1;

  if (mod == 1)
    return { { 0, 0, 0 /* value does not matter */ } };

  /* first determine the prime factors of mod */
  std::vector<long> factors;
  for (long t = mod, p = 2; t != 1; p += 1 + (p & 1))
    {
      if ((t % p) == 0)
        {
          q = 1;
          while ((t % p) == 0)
            {
              t /= p;
              q *= p;
            }
          factors.push_back (q);
        }
    }

  std::vector<congruence> c, d, e;

  for (size_t i = 0; i < factors.size (); i++)
    {
      d.clear ();
      q = factors[i];
      for (long v = 0; v < q; v++)
        {
          for (long w = 0; w < q; w++)
            {
              double alpha = average_alpha (poly, v, w, q, rstate);
              insert_congruence (d, alpha, v, w, vmin, vmax, q);
            }
        }
      if (i == 0)
        c = d;
      else /* merge c and d into e */
        {
          /* Q and q are coprime, being powers of distinct primes */
          long const inv = (long) invert_ul ((unsigned long) (Q % q),
                                             (unsigned long) q);
          e.clear ();
          for (auto const & ci : c)
            for (auto const & di : d)
              {
                double alpha = ci.alpha + di.alpha;
                /* if alpha is larger (i.e., worse) than the last element,
                   since d[] is sorted by increasing values of alpha, we
                   assume all further values will be worse */
                if ((int) e.size () == keep && e.back ().alpha < alpha)
                  break;
                long v = crt (ci.vmod, di.vmod, Q, q, inv);
                long w = crt (ci.wmod, di.wmod, Q, q, inv);
                insert_congruence (e, alpha, v, w, vmin, vmax, Q * q);
              }
          c = e;
        }
      Q *= q;
    }

  /* if u = -u0, check the congruence of the initial polynomial (-v0,-w0) */
  if (u == -u0)
    {
      int included = -1;
      for (size_t i = 0; i < c.size (); i++)
        if (get_mod (-v0, mod) == c[i].vmod && get_mod (-w0, mod) == c[i].wmod)
          included = i;
      if (included >= 0)
        printf ("congruence of initial polynomial has rank %d (%.2f)\n",
                included, c[included].alpha);
      else
        {
          double alpha = 0;
          for (long f : factors)
            alpha += average_alpha (poly, get_mod (-v0, f), get_mod (-w0, f),
                                    f, rstate);
          printf ("congruence of initial polynomial is not included");
          if (!c.empty ())
            printf (" (last %.2f wrt %.2f)\n", c.back ().alpha, alpha);
          else
            printf (" (%.2f)\n", alpha);
        }
    }

  return c;
}

/* Compute the contribution of each prime power to alpha, for the
   polynomial cpoly. An entry (q, s, nu) means that nu is to be added to
   every cell of index s, s+q, s+2q, ... of a sieve array that starts at
   w = wmin. */
std::vector<sieve_data>
rootsieve::prepare_sieve (cxx_cado_poly const & cpoly, long wmin) const
{
    cxx_mpz ump;
    std::vector<unsigned long> roots (B);
    std::vector<float> L (B);

    /* the sieve holds at most one entry per residue class modulo the
       largest power of each prime */
    unsigned long sum_of_prime_powers = 0;
    for (long q : prime_powers)
        sum_of_prime_powers += q;

    std::vector<sieve_data> sieve_d;
    sieve_d.reserve (sum_of_prime_powers);

    for (size_t l = 0; l < primes.size (); l++)
    {
        long const p = primes[l];
        long const qmax = prime_powers[l];
        double const logp = log ((double) p);

        std::fill (L.begin (), L.end (), 0);
        for (long q = p; q <= qmax; q *= p)
        {
            /* the contribution is log(p)/p^(k-1)/(p+1) when the exponent k
               is not the largest one, and log(p)/p^(k-1)/(p+1)*p/(p-1) for
               the largest exponent k */
            double nu = logp / (double) q * (double) p / (double) (p + 1);
#ifndef ORIGINAL
            if (q == qmax)
                nu *= (double) p / (double) (p - 1);
#endif
            for (long x = 0; x < q; x++)
            {
                /* compute f(x) and g(x) mod p^k, where q = p^k */
                mpz_poly_eval_ui (ump, cpoly[ALG_SIDE], x);
                unsigned long const fx = mpz_fdiv_ui (ump, q);
                mpz_poly_eval_ui (ump, cpoly[RAT_SIDE], x);
                unsigned long const gx = mpz_fdiv_ui (ump, q);
                /* search roots w of fx + w*gx = 0 mod q */
                unsigned long const nroots = get_roots (roots.data (), fx, gx, q);
                for (unsigned long i = 0; i < nroots; i++)
                {
                    long const w = roots[i];
                    /* update for w+t*q */
                    for (long t = 0; t < qmax / q; t++)
                        L[w + t * q] += nu;
                }
            }
        }

        /* prepare data for the sieve */
        for (long w = 0; w < qmax; w++)
        {
            double const nu = L[w];
            if (nu == 0.0)
                continue;
            /* compute s = w+t*qmax-wmin such that s - qmax < 0 <= s, i.e.,
               (t-1)*qmax < wmin-w <= t*qmax: t = ceil((wmin-w)/qmax) */
            long t;
            if (wmin - w < 0)
                t = (wmin - w) / qmax;
            else
                t = (wmin - w + qmax - 1) / qmax;
            long const s = w + t * qmax - wmin;
            sieve_d.push_back ({ (uint16_t) qmax, (uint16_t) s, (float) nu });
        }
    }

    ASSERT_ALWAYS(sieve_d.size () <= sum_of_prime_powers);

    return sieve_d;
}

/* Rotate cpoly by w -- the local value, the global rotation being
   mod*w+modw -- and return the Murphy-E value of the result, with its
   lognorm in *lognorm. cpoly comes back as it was. */
double
rootsieve::murphyE_at (cxx_cado_poly & cpoly, long w, double & lognorm) const
{
  rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, w, 0);
  double const skew = cpoly.skew; /* save skewness */
  cpoly.skew = L2_skewness (cpoly[ALG_SIDE]);
  lognorm = L2_lognorm (cpoly[ALG_SIDE], cpoly.skew);
  /* to compute E, we need to divide g by mod */
  mpz_poly_divexact_ui (cpoly[RAT_SIDE], cpoly[RAT_SIDE], mod);
  double const E = MurphyE (cpoly, Bf, Bg, area, MURPHY_K, get_alpha_bound ());
  /* restore g, the skewness, and the polynomial we were given */
  mpz_poly_mul_si (cpoly[RAT_SIDE], cpoly[RAT_SIDE], mod);
  cpoly.skew = skew;
  rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], w, 0, 0);
  return E;
}

/* rotation for a fixed value of v */
void
rootsieve::rotate_v (long v, long u, long modw)
{
  long wmin, wmax;
  cxx_mpz wminz, wmaxz;
  double tot_pols_local = 0;
  double tot_alpha_local = 0;

  /* first make a local copy of the original polynomial */
  cxx_cado_poly cpoly = poly;

  /* compute f + (v*x)*g */
  rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, v, 1);

  rotation_space r;
  expected_growth (&r, cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0,
                   maxlognorm, cpoly.skew);
  mpz_set_d (wminz, r.kmin);
  mpz_set_d (wmaxz, r.kmax);

  if (verbose > 1)
#pragma omp critical
    {
      fmt::print ("v={}: wmin={} wmax={}\n", v, wminz, wmaxz);
      fflush (stdout);
    }

  /* Ensure wminz % mod = modw. Since mod <= MAX_LONG, we have
     s := wminz % mod < MAX_LONG thus modw - s fits in a long and
     there is no underflow below. */
  long t = get_mod (modw - mpz_fdiv_ui (wminz, mod), mod);
  ASSERT_ALWAYS(0 <= t && t < mod);
  mpz_add_ui (wminz, wminz, t);

  /* nothing to sieve for this v */
  if (mpz_cmp (wminz, wmaxz) >= 0)
    return;

  /* if mod != 1, we have f + (k*mod+modw)*g = (f+modw*g) + k*(mod*g) */
  if (mod > 1)
  {
      rotate_aux (cpoly[ALG_SIDE], cpoly[RAT_SIDE], 0, modw, 0); /* f <- f+modw*g */
      mpz_poly_mul_si (cpoly[RAT_SIDE], cpoly[RAT_SIDE], mod);
      /* wmin -> (wmin - modw) / mod */
      mpz_sub_ui (wminz, wminz, modw);
      ASSERT_ALWAYS(mpz_divisible_ui_p (wminz, mod));
      mpz_divexact_ui (wminz, wminz, mod);
      /* wmax -> (wmax - modw) / mod */
      mpz_sub_ui (wmaxz, wmaxz, modw);
      mpz_cdiv_q_ui (wmaxz, wmaxz, mod);
  }

  /* compute the expected value sum(log(p)/(p-1), p < B) */
  double expected = 0.0;
  for (long p : primes)
      expected += log ((double) p) / (double) (p - 1);

  ASSERT_ALWAYS (mpz_fits_slong_p (wmaxz));
  ASSERT_ALWAYS (mpz_fits_slong_p (wminz));

  wmax = mpz_get_si (wmaxz);
  wmin = mpz_get_si (wminz);

  /* the sieve carries its offsets from one chunk to the next, so
     this is modified as we go */
  std::vector<sieve_data> sieve_d = prepare_sieve (cpoly, wmin);

  std::vector<float> A (LEN);

  /* we sieve by chunks of LEN cells at a time */
  long wcur = wmin;
  while (wcur < wmax)
  {
      /* A[j] corresponds to w = wcur + j */
      std::fill (A.begin (), A.end (), (float) expected);

#if defined(TRACE_V) && defined(TRACE_W)
      if (v == TRACE_V && (mod * wcur + modw <= TRACE_W &&
                  TRACE_W < mod * (wcur + LEN) + modw))
          printf ("initialized A[%d] to %f\n", TRACE_W, A[(TRACE_W - modw) / mod - wcur]);
#endif

      /* now perform the sieve */
      for (auto & sd : sieve_d)
      {
          long q = sd.q;
          long s = sd.s;
          float nu = sd.nu;
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
          sd.s = (uint16_t) (s - LEN);
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
              /* local value of w, the global one is mod * w + modw */
              long const w = wcur + j;
              double lognorm;
              double const E = murphyE_at (cpoly, w, lognorm);
              /* this can only occur for one thread, thus no need to put
#pragma omp critical */
              printf ("u=%ld v=%ld w=%ld lognorm=%.2f est_alpha_aff=%.2f E=%.2e [original]\n",
                      u, v, mod * w + modw, lognorm, (double) A[j], E);
              fflush (stdout);
          }
          if (A[j] < best_alpha + guard_alpha)
          {
              long const w = wcur + j;
              double lognorm;
              double const E = murphyE_at (cpoly, w, lognorm);

              if (optimizeE == 0 || (optimizeE == 1 && E > best_E))
#pragma omp critical
              {
                  bestu = u;
                  bestv = v;
                  mpz_set_si (bestw, mod * w + modw);
                  best_alpha = (double) A[j];
                  best_E = E;
                  fmt::print ("u={} v={} w={} lognorm={:.2f} est_alpha_aff={:.2f} E={:.2e}\n",
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

  /* accumulate the number of polynomials and the alpha values */
#pragma omp critical
  {
    tot_pols += tot_pols_local;
    tot_alpha += tot_alpha_local;
  }
}

void
rootsieve::rotate (long u)
{
  /* determine range [vmin,vmax] */
  rotation_space r;
  expected_growth (&r, poly[ALG_SIDE], poly[RAT_SIDE], 1,
                   maxlognorm, poly.skew);
  long vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  long vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
  if (verbose)
    {
      printf ("u=%ld: vmin=%ld vmax=%ld\n", u, vmin, vmax);
      fflush (stdout);
    }

  auto const c = best_congruences (vmin, vmax, u);
  int const n = c.size ();
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
        rotate_v (v, u, c[i].wmod);
    }
}

/* Recover the transformation that took our own polynomial to cpoly, the
   size-optimized one, print it, and apply it to ours. cpoly is left
   alone. */
void
rootsieve::print_transformation (cxx_cado_poly const & cpoly)
{
  cxx_mpz k;
  int d = poly[ALG_SIDE]->deg;

  /* first compute the translation k: g(x+k) = g1*x + g1*k + g0 */
  mpz_sub (k, mpz_poly_coeff_const(cpoly[RAT_SIDE], 0), mpz_poly_coeff_const(poly[RAT_SIDE], 0));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly[RAT_SIDE], 1));
  fmt::print ("translation {}, ", k);

  mpz_poly_translation(poly[ALG_SIDE], poly[ALG_SIDE], k);
  mpz_poly_translation(poly[RAT_SIDE], poly[RAT_SIDE], k);

  /* size_optimization might multiply f0 by some integer t */
  ASSERT_ALWAYS(mpz_divisible_p (mpz_poly_coeff_const(cpoly[ALG_SIDE], d),
				 mpz_poly_coeff_const(poly[ALG_SIDE], d)));
  mpz_divexact (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], d),
		mpz_poly_coeff_const(poly[ALG_SIDE], d));
  if (mpz_cmp_ui (k, 1) != 0)
    {
      fmt::print ("multiplier {}, ", k);
      mpz_poly_mul_mpz (poly[ALG_SIDE], poly[ALG_SIDE], k);
    }
  /* now compute rotation by x^2 */
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 3), mpz_poly_coeff_const(poly[ALG_SIDE], 3));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly[RAT_SIDE], 1));
  fmt::print ("rotation [{},", k);
  ASSERT (mpz_fits_slong_p (k));
  u0 = mpz_get_si (k);
  mpz_poly_rotation(poly[ALG_SIDE], poly[ALG_SIDE], poly[RAT_SIDE], k, 2);
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 2), mpz_poly_coeff_const(poly[ALG_SIDE], 2));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly[RAT_SIDE], 1));
  fmt::print ("{},", k);
  ASSERT (mpz_fits_slong_p (k));
  v0 = mpz_get_si (k);
  mpz_poly_rotation(poly[ALG_SIDE], poly[ALG_SIDE], poly[RAT_SIDE], k, 1);
  mpz_sub (k, mpz_poly_coeff_const(cpoly[ALG_SIDE], 1), mpz_poly_coeff_const(poly[ALG_SIDE], 1));
  ASSERT_ALWAYS(mpz_divisible_p (k, mpz_poly_coeff_const(poly[RAT_SIDE], 1)));
  mpz_divexact (k, k, mpz_poly_coeff_const(poly[RAT_SIDE], 1));
  fmt::print ("{}]\n", k);
  ASSERT (mpz_fits_slong_p (k));
  w0 = mpz_get_si (k);
  mpz_poly_rotation(poly[ALG_SIDE], poly[ALG_SIDE], poly[RAT_SIDE], k, 0);
  ASSERT_ALWAYS(mpz_cmp (mpz_poly_coeff_const(poly[ALG_SIDE], 0),
                         mpz_poly_coeff_const(cpoly[ALG_SIDE], 0)) == 0);
}

static double
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
static double
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
rootsieve::rotate_area (long umin, long umax) const
{
  double sum = 0.0;

  for (long u = umin; u <= umax; u++)
    sum += rotate_area_u (poly, maxlognorm, u);
  return sum;
}

/* Given a sieving area, our maximal effort, and our value of keep,
   compute the best 'mod' value. */
long
rootsieve::best_mod (double sieving_area)
{
  /* http://oeis.org/A051451 */
  static constexpr long l[] = {1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360,
    720720, 12252240, 232792560, 5354228880, 26771144400, 80313433200,
    2329089562800};
  double e = 0;

  for (long m : l)
    {
      mod = m;
      /* The number of polynomials sieved is approximately area/mod^2*keep.
         Note: there is a bias when vmax-vmin is smaller than mod, since we
         only keep congruences that contain at least an element in
         [vmin, vmax], thus the probability is larger than (vmax-vmin)/mod. */
      e = sieving_area / (double) mod / (double) mod * (double) keep;
      if (e <= effort)
        break;
    }
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
    rootsieve rs;
    int I = 0;
    double margin = NORM_MARGIN;
    long umin = LONG_MIN, umax = LONG_MAX;
    int sopt = 0;
    double time = seconds ();

    cxx_param_list pl;
    const char * polyfilename = NULL;

    declare_usage(pl);
    pl.configure_switch_old("-v", &rs.verbose);
    pl.configure_switch_old("-sopt", &sopt);
    pl.configure_switch_old("-E", &rs.optimizeE);

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
    pl.parse("effort", rs.effort);
    pl.parse("keep", rs.keep);
    if (pl.parse("mod", rs.mod) && rs.mod < 1)
        pl.fail("Error, -mod must be at least 1");
    /* LONG_MIN and LONG_MAX are reserved to mean "the user said nothing",
     * so they are only rejected when the option is actually given. */
    if (pl.parse("umin", umin) && umin == LONG_MIN)
        pl.fail("Error, -umin is out of range");
    if (pl.parse("umax", umax) && umax == LONG_MAX)
        pl.fail("Error, -umax is out of range");
    if (pl.parse("B", rs.B))
        set_alpha_bound(rs.B);
    if (!polyfilename)
        polyfilename = pl.lookup_old("poly");

    if (rs.optimizeE)
        rs.guard_alpha = GUARD_ALPHA;

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

    ASSERT_ALWAYS(rs.B <= 65536);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    rs.Bf = bound_f;
    rs.Bg = bound_g;
    rs.area = area;

    if (!rs.poly.read(polyfilename))
        pl.fail("Problem when reading file {}", polyfilename);

    if (rs.poly.skew == 0.0)
      rs.poly.skew = L2_skewness (rs.poly[ALG_SIDE]);

    /* if -sopt, size-optimize */
    if (sopt)
      {
        cxx_cado_poly c = rs.poly;
        size_optimization (c[ALG_SIDE], c[RAT_SIDE],
                           rs.poly[ALG_SIDE], rs.poly[RAT_SIDE],
                           SOPT_DEFAULT_EFFORT, rs.verbose);
        printf ("# initial polynomial:\n");
        rs.poly.fprintf(stdout);
        rs.print_transformation (c);
        rs.poly = c;
        printf ("# size-optimized polynomial:\n");
        rs.poly.fprintf(stdout);
      }

    rs.init_primes ();

    /* compute the skewness */
    rs.poly.skew = L2_skewness (rs.poly[ALG_SIDE]);
    double lognorm = L2_lognorm (rs.poly[ALG_SIDE], rs.poly.skew);
    rs.maxlognorm = lognorm + margin;
    printf ("initial lognorm %.2f, maxlognorm %.2f\n", lognorm, rs.maxlognorm);

    /* determine range [umin,umax] */
    rotation_space r;
    expected_growth (&r, rs.poly[ALG_SIDE], rs.poly[RAT_SIDE], 2,
                     rs.maxlognorm, rs.poly.skew);
    if (umin == LONG_MIN) /* umin was not given by the user */
      umin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
    if (umax == LONG_MAX) /* umax was not given by the user */
      umax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
    if (rs.verbose)
      printf ("umin=%ld umax=%ld\n", umin, umax);

    if (rs.mod == 0) /* compute best 'mod' for given effort */
      {
        double sieving_area = rs.rotate_area (umin, umax);
        /* print total sieving area */
        printf ("sieving area %.2e\n", sieving_area);
        rs.best_mod (sieving_area);
      }

    long ucur = 0; /* current translation in u */
    for (long u = umin; u <= umax; u++)
      {
        rotate_aux (rs.poly[ALG_SIDE], rs.poly[RAT_SIDE], ucur, u, 2);
        ucur = u;

        rs.rotate (u);
      }

    /* restore original polynomial */
    rotate_aux (rs.poly[ALG_SIDE], rs.poly[RAT_SIDE], ucur, 0, 2);

    fmt::print ("best rotation: u={} v={} w={} alpha={:.2f}\n",
            rs.bestu, rs.bestv, rs.bestw, rs.best_alpha);

    /* perform the best rotation */
    rotate_aux (rs.poly[ALG_SIDE], rs.poly[RAT_SIDE], 0, rs.bestu, 2);
    rotate_aux (rs.poly[ALG_SIDE], rs.poly[RAT_SIDE], 0, rs.bestv, 1);
    mpz_poly_rotation(rs.poly[ALG_SIDE], rs.poly[ALG_SIDE], rs.poly[RAT_SIDE],
                      rs.bestw, 0);

    /* recompute the skewness of the best polynomial */
    rs.poly.skew = L2_combined_skewness2 (rs.poly[0], rs.poly[1]);

    print_cadopoly_extra (stdout, rs.poly, argc0, argv0, 0);

    time = seconds () - time;
    printf ("# Sieved %.2e polynomials in %.2f seconds (%.2es/p)\n",
            rs.tot_pols, time, time / rs.tot_pols);
    printf ("# Average alpha %.2f\n", rs.tot_alpha / rs.tot_pols);

    return 0;
}
