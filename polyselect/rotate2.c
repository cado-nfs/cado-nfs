#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "area.h"
#include "murphyE.h"
#include "size_optimization.h"
#include "omp.h"

/* global variables */
int verbose = 0;                /* verbosity level */
unsigned long *Primes, nprimes; /* primes less than B */
unsigned long *Q;               /* largest p^k < B */
long bestu = 0, bestv = 0;      /* current best rotation */
double best_alphaprime;         /* alpha' value of best rotation */
double best_alpha, best_var;    /* (alpha,var) of best rotation */
double best_E;                  /* best MurphyE_chi2 value */
double best_Eprime;             /* lognorm + alphaprime of best rotation */
mpz_t bestw;                    /* current best rotation in w */

static void
usage_and_die (char *argv0)
{
  fprintf (stderr, "usage: %s [-area a] [-I nnn] [-Bf b] [-Bg c] [-skew s] poly\n", argv0);
  fprintf (stderr, "  apply rotation f += (v*x+w)*g to poly.\n");
  fprintf (stderr, "  poly: filename of polynomial\n");
  exit (1);
}

/* returns x mod m, with 0 <= x < m */
static long
umod (long x, long m)
{
  x = x % m;
  return (x >= 0) ? x : x + m;
}

/* get 0.8981 quantile of N(alpha,var) */
static double
quantile (double alpha, double var)
{
  return alpha - 1.2707994625519516 * sqrt (var);
}

static unsigned long
initPrimes (unsigned long B)
{
  unsigned long nprimes = 0, p, q, l;

  Primes = malloc (B * sizeof (unsigned long));
  ASSERT_ALWAYS(Primes != NULL);
  for (p = 2; p < 2000; p += 1 + (p > 2))
    if (ulong_isprime (p))
      Primes[nprimes++] = p;
  Primes = realloc (Primes, nprimes * sizeof (unsigned long));
  ASSERT_ALWAYS(Primes != NULL);

  /* compute prime powers */
  Q = malloc (nprimes * sizeof (unsigned long));
  ASSERT_ALWAYS(Q != NULL);
  for (l = 0; l < nprimes; l++)
    {
      p = Primes[l];
      for (q = p; q * p < B; q *= p);
      Q[l] = q;
    }

  return nprimes;
}

/* determine the best value of 'mod' <= maxmod:
   (1) prime factors of 'mod' should divide both f[d] and f[d-1]
   (2) for two prime factors p and q, if p^e divides exactly mod and
       q^(f+1) does not divide mod, we should have p^e < q^(f+1).
*/
long
get_mod (mpz_poly_srcptr f, long maxmod)
{
  int d = f->deg, nprimes = 0, i, i0;
  mpz_t g;
  long P[8], Q[8] = {1,1,1,1,1,1,1,1}, e[8] = {0,0,0,0,0,0,0,0}, p, mod;
  long primes[8] = {2, 3, 5, 7, 11, 13, 17, 19};

  mpz_init (g);
  mpz_gcd (g, f->coeff[d], f->coeff[d - 1]);
  /* first determine primes in {2, 3, 5, 7, 11, 13, 17, 19} dividing g */
  for (mod = 1, i = 0; i < 8; i++)
    {
      p = primes[i];
      if (mpz_divisible_ui_p (g, p))
	{
	  P[nprimes++] = p;
	  mod *= p;
	}
    }
  /* Ensure nprimes large enough to avoid q too large in rotate_v.
     As a rule of thumb we want each q <= 1024 thus
     nprimes >= log(maxmod)/log(1024) */
  int min_nprimes = ceil (log ((double) maxmod) / log (1024.0));
  for (i = 0; i < 8 && nprimes < min_nprimes; i++)
    if (mod % primes[i] != 0)
      P[nprimes++] = primes[i];
  mod = 1;
  while (1)
    {
      /* determine the smallest Q[i]*P[i] such that mod*P[i] <= maxmod */
      for (i0 = -1, i = 0; i < nprimes; i++)
	if (mod * P[i] <= maxmod && (i0 == -1 || Q[i] * P[i] < Q[i0] * P[i0]))
	  i0 = i;
      if (i0 == -1)
	break;
      Q[i0] *= P[i0];
      e[i0] ++;
      mod *= P[i0];
      ASSERT_ALWAYS(mod <= maxmod);
    }
  printf ("mod=%ld=", mod);
  for (i = 0; i < nprimes; i++)
    {
      if (i > 0)
	printf ("*");
      if (e[i] > 0)
	{
	  printf ("%ld", P[i]);
	  if (e[i] > 1)
	    printf ("^%ld", e[i]);
	}
    }
  printf ("\n");
  mpz_clear (g);
  return mod;
}

static void
print_disc (mpz_poly_srcptr f)
{
  mpz_t d;
  int first = 1;
  mpz_init (d);
  mpz_poly_discriminant (d, f);
  for (unsigned long p = 2; p < 20; p += 1 + (p > 2))
    {
      unsigned long e = 0;
      while (mpz_divisible_ui_p (d, p))
	{
	  e ++;
	  mpz_divexact_ui (d, d, p);
	}
      if (e != 0)
	{
	  if (!first)
	    printf ("*");
	  first = 0;
	  printf ("%lu", p);
	  if (e > 1)
	    printf ("^%lu", e);
	}
    }
  printf ("\n");
  mpz_clear (d);
}

/* rotation for a fixed value of v.
   lognorm0 is the lognorm of the initial polynomial (for v=w=0)
   skew is the skewness of the initial polynomial (for v=w=0) */
static void
rotate_v (cado_poly_srcptr poly0, long v, long B,
          double maxlognorm, double Bf, double Bg, double area,
          long mod, long wmod, long u)
{
  double s1, s2;
  double lognorm, alphaprime, Eprime, E, Alphaprime = 0.0, V_alphaprime = 0.0;
  long w0, w, ww, wmin, wmax;
  cado_poly poly;
  float *alpha, *var, **A, **V;
  unsigned long l;
  mpz_t g1, g0, wminz, wmaxz;

  /* first make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  mpz_init_set (g1, poly->pols[RAT_SIDE]->coeff[1]);
  mpz_init_set (g0, poly->pols[RAT_SIDE]->coeff[0]);
  mpz_init (wminz);
  mpz_init (wmaxz);

  /* compute f + (v*x)*g */
  rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, 0, v, 1);

  /* compute lognorm of current polynomial at (v,0) */
  double skew = L2_skewness (poly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);
  lognorm = L2_lognorm (poly->pols[ALG_SIDE], skew);

  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0,
                   maxlognorm, poly->skew);
  mpz_set_d (wminz, r.kmin);
  mpz_set_d (wmaxz, r.kmax);
  mpz_add_ui (wmaxz, wmaxz, 1); /* wmax is not included in the search range */

  /* ensure wminz = wmod mod 'mod' */
  ASSERT_ALWAYS(0 <= wmod && wmod < mod);
  long tmp = wmod - mpz_fdiv_ui (wminz, mod);
  if (tmp < 0)
    tmp += mod;
  if (tmp >= mod)
    tmp -= mod;
  mpz_add_ui (wminz, wminz, tmp);
  ASSERT_ALWAYS(mpz_fdiv_ui (wminz, mod) == (unsigned long) wmod);

  /* ensure wmaxz = wmod mod 'mod' */
  tmp = wmod - mpz_fdiv_ui (wmaxz, mod);
  if (tmp < 0)
    tmp += mod;
  if (tmp >= mod)
    tmp -= mod;
  mpz_add_ui (wmaxz, wmaxz, tmp);
  ASSERT_ALWAYS(mpz_fdiv_ui (wmaxz, mod) == (unsigned long) wmod);

  if (mpz_cmp (wminz, wmaxz) >= 0)
    goto end;

  if (verbose)
#pragma omp critical
    {
      gmp_printf ("v=%ld: wmin=%Zd wmax=%Zd ", v, wminz, wmaxz);
      print_disc (poly->pols[ALG_SIDE]);
      fflush (stdout);
    }

  /* rotating by w * g where w = t*mod + wmod is equivalent to first add
     wmod*g, then rotating by t * (mod*g) */
  rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, 0, wmod, 0);
  w0 = 0; /* the initial polynomial is f + wmod * g */
  /* multiply g1 and g0 by mod */
  mpz_mul_si (g1, g1, mod);
  mpz_mul_si (g0, g0, mod);
  /* divide wmin and wmax by mod */
  mpz_sub_ui (wminz, wminz, wmod);
  ASSERT_ALWAYS(mpz_divisible_ui_p (wminz, mod));
  mpz_divexact_ui (wminz, wminz, mod);
  mpz_sub_ui (wmaxz, wmaxz, wmod);
  ASSERT_ALWAYS(mpz_divisible_ui_p (wmaxz, mod));
  mpz_divexact_ui (wmaxz, wmaxz, mod);
  ASSERT_ALWAYS(mpz_fits_slong_p (wminz));
  wmin = mpz_get_si (wminz);
  ASSERT_ALWAYS(mpz_fits_slong_p (wmaxz));
  wmax = mpz_get_si (wmaxz);

  A = malloc (nprimes * sizeof (float*));
  ASSERT_ALWAYS(A != NULL);
  V = malloc (nprimes * sizeof (float*));
  ASSERT_ALWAYS(V != NULL);
  for (l = 0; l < nprimes; l++)
    {
      A[l] = malloc (B * sizeof (float));
      ASSERT_ALWAYS(A[l] != NULL);
      V[l] = malloc (B * sizeof (float));
      ASSERT_ALWAYS(V[l] != NULL);
    }

  /* first compute A[l][] and V[l][] */
  s1 = seconds_thread ();
  for (l = 0; l < nprimes; l++)
    {
      long p = Primes[l], q = Q[l];
      for (w = 0; w < q; w++)
        {
          double t;
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w0, w, 0);
          w0 = w;
          A[l][w] = dist_alpha_p (poly->pols[ALG_SIDE], p, &t);
          V[l][w] = t;
        }
      /* reset initial polynomial */
      rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w0, 0, 0);
      w0 = 0;
    }
  s1 = seconds_thread () - s1;

#define BUCKET (1<<16)
  /* fill the alpha[] and var[] arrays: alpha[w - start] is the alpha value for
     f+w*g, wmin <= w < wmax, and var[w - wmin] the corresponding variance */
  s2 = seconds_thread ();
  alpha = malloc (BUCKET * sizeof (float));
  ASSERT_ALWAYS(alpha != NULL);
  var = malloc (BUCKET * sizeof (float));
  ASSERT_ALWAYS(var != NULL);
  long start, end, len, q, wq;
  long *Wq = malloc (nprimes * sizeof (long));
  ASSERT_ALWAYS(Wq != NULL);
  for (l = 0; l < nprimes; l++)
    {
      q = Q[l];
      wq = wmin % q;
      if (wq < 0)
        wq += q;
      Wq[l] = wq;
    }
  start = wmin;
  while (start < wmax)
    {
      end = start + BUCKET;
      if (end > wmax)
        end = wmax;
      len = end - start;
      memset (alpha, 0, len * sizeof (float));
      memset (var, 0, len * sizeof (float));
      for (l = 0; l < nprimes; l++)
        {
          q = Q[l];
          wq = Wq[l];
          /* process the first q-wq elements */
          for (ww = 0; ww < len && wq < q; ww++, wq++)
            {
              alpha[ww] += A[l][wq];
              var[ww] += V[l][wq];
            }
          if (wq == q)
            wq = 0;
          /* now either ww=len, and the last loops do nothing;
             or ww < len, which implies wq = 0 */
          /* now process batches of q elements */
          for (; ww + q <= len; ww += q)
            for (long i = 0; i < q; i++)
              {
                alpha[ww + i] += A[l][i];
                var[ww + i] += V[l][i];
            }
          /* process the last <q elements */
          for (; ww < len; ww++)
            {
              alpha[ww] += A[l][wq];
              var[ww] += V[l][wq];
              if (++wq == q)
                wq = 0;
            }
          Wq[l] = wq;
        }

      /* now scan for best alpha_prime */
      for (w = start; w < end; w++)
        {
          ww = w - start;
          alphaprime = quantile (alpha[ww], var[ww]);
          Alphaprime += alphaprime;
          V_alphaprime += alphaprime * alphaprime;
          Eprime = lognorm + alphaprime;
#define MARGIN 1.0 /* 2.0 increases the time by about 32% */
          if (Eprime < best_Eprime + MARGIN)
            {
              rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w0, w, 0);
              w0 = w;
              alphaprime = get_alpha_prime (poly->pols[ALG_SIDE], B);
              Eprime = lognorm + alphaprime;
              if (Eprime < best_Eprime + MARGIN)
                {
                  poly->skew = L2_combined_skewness2 (poly->pols[0],
                                                      poly->pols[1],
                                                      SKEWNESS_DEFAULT_PREC);
                  E = MurphyE_chi2 (poly, Bf, Bg, area, MURPHY_K);
                  if (E > best_E)
#pragma omp critical
                    {
                      best_alphaprime = alphaprime;
                      best_Eprime = Eprime;
                      best_alpha = alpha[ww];
                      best_var = var[ww];
                      bestv = v;
                      bestu = u;
		      mpz_set_si (bestw, mod);
		      mpz_mul_si (bestw, bestw, w);
		      mpz_add_si (bestw, bestw, wmod);
                      best_E = E;
                      gmp_printf ("u=%ld v=%ld w=%Zd alpha'=%1.2f (%1.2f,%1.2f) E=%.2e\n",
                              u, v, bestw, alphaprime,
                              best_alpha, best_var, best_E);
                      fflush (stdout);
                    }
                }
            }
        }
      start = end;
    }
  free (Wq);
  free (alpha);
  free (var);
  for (l = 0; l < nprimes; l++)
    {
      free (A[l]);
      free (V[l]);
    }
  free (A);
  free (V);

  Alphaprime /= (wmax - wmin);
  V_alphaprime /= (wmax - wmin);
  V_alphaprime -= Alphaprime * Alphaprime;
  s2 = seconds_thread () - s2;

  if (verbose)
#pragma omp critical
    {
      printf ("v=%ld: prepare %.2fs, sieve+scan %.2fs, avg %.2f, var %.2f\n",
              v, s1, s2, Alphaprime, V_alphaprime);
      fflush (stdout);
    }

 end:
  cado_poly_clear (poly);
  mpz_clear (g1);
  mpz_clear (g0);
  mpz_clear (wminz);
  mpz_clear (wmaxz);
}

static void
rotate (cado_poly poly, long B, double maxlognorm, double Bf, double Bg,
        double area, long vmin, long vmax, long mod, long vmod, long wmod,
	long u)
{
  /* ensure vmin = vmod mod mod */
  long tmp = (vmod - vmin) % mod;
  vmin += (tmp < 0) ? tmp + mod : tmp;
  ASSERT_ALWAYS((vmin % mod) == vmod || (vmin % mod) == vmod - mod);

  if (mod == 1)
#pragma omp parallel for schedule(dynamic)
    for (long v = vmin; v <= vmax; v += mod)
      rotate_v (poly, v, B, maxlognorm, Bf, Bg, area, mod, wmod, u);
  else
    for (long v = vmin; v <= vmax; v += mod)
      rotate_v (poly, v, B, maxlognorm, Bf, Bg, area, mod, wmod, u);
}

typedef struct
{
  long vmod, wmod;
  double alpha, var, alphaprime;
} class;

/* Insert (alpha,var,alphaprime) into c, where c has already n entries
   (maximum is keep). Return the new value of n. */
static int
insert_class (class *c, int n, int keep, double alpha, double var,
              double alphaprime, long v, long w,
	      long vmin, long vmax, long mod)
{
  int i;

  /* check if this class has at least one representative in [vmin,vmax] */
  long t = (v - vmin) % mod; /* vmin + t = v mod 'mod' */
  t = (t < 0) ? t + mod : t; /* vmin + t = v mod 'mod' and 0 <= t < mod */
  if (vmin + t > vmax)
    return n; /* no representative in [vmin,vmax] */

  for (i = n; i > 0 && alphaprime < c[i-1].alphaprime; i--)
    c[i] = c[i-1];
  /* now i = 0 or alphaprime >= c[i-1].alphaprime */
  if (i < keep)
    {
      c[i].vmod = v;
      c[i].wmod = w;
      c[i].alpha = alpha;
      c[i].var = var;
      c[i].alphaprime = alphaprime;
    }
  n += (n < keep);
  return n;
}

#if 0
static double
alphaprime_class (cado_poly poly0, long w, unsigned long *primes,
                  int nprimes, double *alpha, double *var)
{
  cado_poly poly;
  double tmp, alphaprime;

  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[0]->coeff[1],
              poly->pols[0]->coeff[0], 0, w, 0);

  *alpha = *var = 0.0;
  for (int i = 0; i < nprimes; i++)
    {
      *alpha += dist_alpha_p (poly->pols[ALG_SIDE], primes[i], &tmp);
      *var += tmp;
    }
  alphaprime = quantile (*alpha, *var);
  cado_poly_clear (poly);
  return alphaprime;
}

/* return the 'keep' best classes (v,w) mod 'mod' */
static class*
best_classes (cado_poly poly0, long mod, int keep, long vmin, long vmax)
{
  class *c, *T;
  int n = 0; /* number of classes so far */
  int nprimes = 0;
  long v, w, t, v0 = 0;
  unsigned long *primes = NULL, p;
  cado_poly poly;

  /* first determine the prime factors of mod */
  for (t = mod, p = 2; t != 1; p += 1 + (p & 1))
    {
      if ((t % p) == 0)
        {
          nprimes ++;
          primes = realloc (primes, nprimes * sizeof (unsigned long));
          ASSERT_ALWAYS(primes != NULL);
          primes[nprimes - 1] = p;
          while ((t % p) == 0)
            t /= p;
        }
    }

  /* make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  c = malloc ((keep + 1) * sizeof (class));
  ASSERT_ALWAYS(c != NULL);
  T = malloc (mod * sizeof (class));
  ASSERT_ALWAYS(T != NULL);
  for (v = 0; v < mod; v++)
    {
      /* check if this class has at least one representative in [vmin,vmax] */
      long t = (v - vmin) % mod; /* vmin + t = v mod 'mod' */
      t = (t < 0) ? t + mod : t; /* vmin + t = v mod 'mod' and 0 <= t < mod */
      if (vmin + t > vmax)
        continue;

      rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[0]->coeff[1],
                  poly->pols[0]->coeff[0], v0, v, 1);
      v0 = v;
#pragma omp parallel for schedule(dynamic)
      for (w = 0; w < mod; w++)
        T[w].alphaprime = alphaprime_class (poly, w, primes, nprimes,
                                            &(T[w].alpha), &(T[w].var));
      for (w = 0; w < mod; w++)
	n = insert_class (c, n, keep, T[w].alpha, T[w].var, T[w].alphaprime,
			  v, w, vmin, vmax, mod);
    }

  cado_poly_clear (poly);
  free (T);
  free (primes);
  return c;
}
#else
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

/* return the 'keep' best classes (v,w) mod 'mod' */
static class*
best_classes (cado_poly poly0, long mod, int keep, long vmin, long vmax)
{
  int nfactors = 0;
  unsigned long *factors = NULL, p;
  class *c, *d, *e;
  int nc = 0; /* number of classes so far in c */
  int nd, ne;
  int i;
  cado_poly poly;
  long q, Q = 1, v0 = 0, w0 = 0;

  /* first determine the prime factors of mod */
  for (long t = mod, p = 2; t != 1; p += 1 + (p & 1))
    {
      if ((t % p) == 0)
        {
          nfactors ++;
          factors = realloc (factors, nfactors * sizeof (unsigned long));
          ASSERT_ALWAYS(factors != NULL);
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
          rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[0]->coeff[1],
                      poly->pols[0]->coeff[0], v0, v, 1);
          v0 = v;
          for (long w = 0; w < q; w++)
            {
              double alphaprime, alpha, var;
              rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[0]->coeff[1],
                          poly->pols[0]->coeff[0], w0, w, 0);
              w0 = w;
              alpha = dist_alpha_p (poly->pols[ALG_SIDE], p, &var);
              alphaprime = quantile (alpha, var);
              nd = insert_class (d, nd, keep, alpha, var, alphaprime, v, w,
				 vmin, vmax, q);
            }
        }
      if (i == 0) /* copy d into c */
        {
          memcpy (c, d, nd * sizeof (class));
          nc = nd;
        }
      else /* merge c and d into e */
        {
          long inv = invert (Q, q);
          ne = 0;
          for (int ic = 0; ic < nc; ic++)
            for (int id = 0; id < nd; id++)
              {
                double alpha = c[ic].alpha + d[id].alpha;
                double var = c[ic].var + d[id].var;
                double alphaprime = quantile (alpha, var);
		/* if alphaprime is larger (i.e., worse) than the last element,
		   since d[] is sorted by increasing values of alphaprime, we
		   assume all further values will be worse */
		if (ne == keep && e[ne-1].alphaprime < alphaprime)
		  break;
                long v = crt (c[ic].vmod, d[id].vmod, Q, q, inv);
                long w = crt (c[ic].wmod, d[id].wmod, Q, q, inv);
                ne = insert_class (e, ne, keep, alpha, var, alphaprime, v, w,
				   vmin, vmax, Q * q);
              }
          /* copy back e to c */
          memcpy (c, e, ne * sizeof (class));
          nc = ne;
        }
      Q *= q;
    }

  cado_poly_clear (poly);
  free (factors);
  free (d);
  free (e);
  return c;
}
#endif

/* don't modify poly, which is the size-optimized polynomial
   (poly0 is the initial polynomial) */
static void
print_transformation (cado_poly_ptr poly0, cado_poly_srcptr poly)
{
  mpz_t k;
  int d = poly0->pols[ALG_SIDE]->deg;

  mpz_init (k);
  /* first compute the translation k: g(x+k) = g1*x + g1*k + g0 */
  mpz_sub (k, poly->pols[RAT_SIDE]->coeff[0], poly0->pols[RAT_SIDE]->coeff[0]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("translation %Zd, ", k);
  do_translate_z (poly0->pols[ALG_SIDE], poly0->pols[RAT_SIDE]->coeff, k);
  /* size_optimization might multiply f0 by some integer t */
  ASSERT_ALWAYS(mpz_divisible_p (poly->pols[ALG_SIDE]->coeff[d],
				 poly0->pols[ALG_SIDE]->coeff[d]));
  mpz_divexact (k, poly->pols[ALG_SIDE]->coeff[d],
		poly0->pols[ALG_SIDE]->coeff[d]);
  if (mpz_cmp_ui (k, 1) != 0)
    {
      gmp_printf ("multiplier %Zd, ", k);
      mpz_poly_mul_mpz (poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], k);
    }
  /* now compute rotation by x^2 */
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[3], poly0->pols[ALG_SIDE]->coeff[3]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("rotation [%Zd,", k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 2);
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[2], poly0->pols[ALG_SIDE]->coeff[2]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("%Zd,", k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 1);
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[1], poly0->pols[ALG_SIDE]->coeff[1]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("%Zd]\n", k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 0);
  ASSERT_ALWAYS(mpz_cmp (poly0->pols[ALG_SIDE]->coeff[0],
                         poly->pols[ALG_SIDE]->coeff[0]) == 0);
  mpz_clear (k);
}

double
rotate_area_v (cado_poly_srcptr poly0, double maxlognorm, long v)
{
  cado_poly poly;
  double area;

  cado_poly_init (poly);
  cado_poly_set (poly, poly0);
  rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
	      poly->pols[RAT_SIDE]->coeff[0], 0, v, 1);
  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0,
                   maxlognorm, poly->skew);
  area = r.kmax - r.kmin;
  cado_poly_clear (poly);
  return area;
}
  
/* estimate the rootsieve area for a given u */
double
rotate_area_u (cado_poly_srcptr poly0, double maxlognorm, long u)
{
  double area, sum = 0.0;
  long h, vmin, vmax;
  cado_poly poly;

  cado_poly_init (poly);
  cado_poly_set (poly, poly0);
  rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
	      poly->pols[RAT_SIDE]->coeff[0], 0, u, 2);
  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 1,
                   maxlognorm, poly->skew);
  vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
#define SAMPLE 100
  if (vmax / SAMPLE - vmin / SAMPLE > 1)
    h = vmax / SAMPLE - vmin / SAMPLE;
  else
    h = 1;

  for (long v = vmin; v <= vmax; v += h)
    {
      area = rotate_area_v (poly, maxlognorm, v);
      sum += area;
    }
  cado_poly_clear (poly);
  return sum * h;
}

/* estimate the rootsieve area for umin <= u <= umax */
double
rotate_area (cado_poly_srcptr poly, double maxlognorm, long umin, long umax)
{
  double area, sum = 0.0;

  for (long u = umin; u <= umax; u++)
    {
      area = rotate_area_u (poly, maxlognorm, u);
      sum += area;
    }
  return sum;
}

int
main (int argc, char **argv)
{
    int argc0 = argc;
    char **argv0 = argv;
    cado_poly poly;
    int I = 0;
    double skew = 0.0;
    double margin = NORM_MARGIN;
    long mod = 1;
    long maxmod = LONG_MAX;            /* compute 'mod' automatically for each
					  polynomial, with mod <= maxmod */
    class *c = NULL;
    int keep = 1, nbest = 0;
    long umin, umax;
    int sopt = 0;
    long tracev = 0, tracew = 0;
    double effort = 0;

    while (argc >= 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-area") == 0)
          {
            area = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bf") == 0)
          {
            bound_f = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-I") == 0)
          {
            I = atoi (argv [2]);
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
        /* For mod we should give numbers from http://oeis.org/A051451:
           1, 2, 6, 12, 60, 420, 840, 2520, 27720, 360360, 720720, 12252240,
           232792560, 5354228880, 26771144400, 80313433200, 2329089562800.
	   If we restrict to prime factors 2, 3, 5, 7:
	   1, 2, 6, 12, 60, 420, 840, 2520, 5040, 25200, 75600, 151200,
	   1058400, 2116800, 6350400, 31752000, 63504000, 190512000,
	   381024000, 2667168000, 5334336000, 26671680000, 80015040000. */
        else if (strcmp (argv[1], "-mod") == 0)
          {
            mod = strtol (argv [2], NULL, 10);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-maxmod") == 0)
          {
            maxmod = strtol (argv [2], NULL, 10);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-tracev") == 0)
          {
            tracev = strtol (argv [2], NULL, 10);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-tracew") == 0)
          {
            tracew = strtol (argv [2], NULL, 10);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-keep") == 0)
          {
            keep = atoi (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-v") == 0)
          {
            verbose ++;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-sopt") == 0)
          {
            sopt = 1;
            argv ++;
            argc --;
          }
        else
          break;
      }
    if (argc != 2)
        usage_and_die (argv[0]);

#pragma omp parallel
#pragma omp master
  printf ("Using %d thread(s)\n", omp_get_num_threads ());

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    cado_poly_init (poly);
    if (!cado_poly_read (poly, argv[1]))
      {
        fprintf (stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die (argv[0]);
      }

    if (poly->skew == 0.0)
      poly->skew = L2_combined_skewness2 (poly->pols[0], poly->pols[1],
					  SKEWNESS_DEFAULT_PREC);

    /* if -sopt, size-optimize */
    if (sopt)
      {
        cado_poly c;
        cado_poly_init (c);
        cado_poly_set (c, poly);
        size_optimization (c->pols[ALG_SIDE], c->pols[RAT_SIDE],
                           poly->pols[ALG_SIDE], poly->pols[RAT_SIDE],
                           SOPT_DEFAULT_EFFORT, verbose);
        printf ("initial polynomial:\n");
        cado_poly_fprintf (stdout, poly, "");
        print_transformation (poly, c);
        cado_poly_set (poly, c);
        printf ("size-optimized polynomial:\n");
        cado_poly_fprintf (stdout, poly, "");
        cado_poly_clear (c);
      }

    /* compute the skewness */
    skew = L2_skewness (poly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);

    nprimes = initPrimes (ALPHA_BOUND);

    /* first compute alpha' of the original polynomial */
    best_alpha = dist_alpha (poly->pols[ALG_SIDE], ALPHA_BOUND, &best_var);
    best_alphaprime = get_alpha_prime (poly->pols[ALG_SIDE], ALPHA_BOUND);
    poly->skew = L2_combined_skewness2 (poly->pols[0], poly->pols[1],
                                      SKEWNESS_DEFAULT_PREC);
    best_E = MurphyE_chi2 (poly, bound_f, bound_g, area, MURPHY_K);
    printf ("u=0 v=0 w=0 alpha'=%1.2f (%1.2f,%1.2f) E=%.2e\n",
            best_alphaprime, best_alpha, best_var, best_E);

    best_E = 0.0;
    best_Eprime = DBL_MAX;

    /* determine range [umin,umax] */
    rotation_space r;
    double maxlognorm = L2_lognorm (poly->pols[ALG_SIDE], skew) + margin;
    expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 2,
                     maxlognorm, poly->skew);
    umin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
    umax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
    printf ("umin=%ld umax=%ld\n", umin, umax);

    double tot_area = rotate_area (poly, maxlognorm, umin, umax);
    printf ("total rootsieve area: %.2e\n", tot_area);

    /* if -effort, compute -maxmod from wanted effort */
    if (effort != 0.0)
      {
	long local_maxmod;
	/* we want keep/mod^2*tot_area ~ effort,
	   thus maxmod = sqrt(tot_area*keep/effort) */
	local_maxmod = (long) sqrt (tot_area * (double) keep / effort);
	if (local_maxmod < maxmod)
	  maxmod = local_maxmod;
	printf ("using maxmod=%ld\n", maxmod);
      }

    /* if -maxmod, compute -mod automatically from given polynomial */
    if (maxmod != 0)
      {
	if (mod != 1)
	  {
	    fprintf (stderr, "Error, -mod and -maxmod are incompatible\n");
	    exit (1);
	  }
	/* now mod = 1 */
	mod = get_mod (poly->pols[ALG_SIDE], maxmod);
      }

    mpz_init (bestw);

    long u0 = 0;
    for (long u = umin; u <= umax; u++)
      {
        long vmin, vmax;
	double area_u = rotate_area_u (poly, maxlognorm, u);

        rotate_aux (poly->pols[ALG_SIDE]->coeff,
                    poly->pols[RAT_SIDE]->coeff[1],
                    poly->pols[RAT_SIDE]->coeff[0], u0, u, 2);
        u0 = u;

        /* recompute the skewness */
        skew = L2_skewness (poly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);
        expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 1,
                         maxlognorm, poly->skew);
        vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
        vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
	if (verbose)
	  printf ("u=%ld: vmin=%ld vmax=%ld (area %.2e)\n", u, vmin, vmax,
		  area_u);

        if (mod != 1)
          {
            /* determine best classes (v,w) mod 'mod' */
            c = best_classes (poly, mod, keep, vmin, vmax);
            nbest = (mod < keep && mod * mod < keep) ? mod * mod : keep;

            /* print best classes */
	    if (verbose)
	      {
		printf ("Best %d/%ld^2 classes mod %ld:\n", nbest, mod, mod);
		for (int i = 0; i < nbest; i++)
		  if (i < 3 || nbest - 3 <= i)
		    printf ("v=%ld w=%ld alpha=%.2f var=%.2f alpha'=%.2f\n",
			    c[i].vmod, c[i].wmod, c[i].alpha, c[i].var, c[i].alphaprime);
		  else if (i == 3 && nbest > 6)
		    printf ("...\n");
	      }

	    if (tracev != 0 || tracew != 0)
	      {
		tracev = umod (tracev, mod);
		tracew = umod (tracew, mod);
		for (int i = 0; i < nbest; i++)
		  if (c[i].vmod == tracev && c[i].wmod == tracew)
		    printf ("TRACE: v=%ld w=%ld alpha=%.2f var=%.2f alpha'=%.2f (rank %d)\n",
			    c[i].vmod, c[i].wmod, c[i].alpha, c[i].var, c[i].alphaprime, i);
	      }
          }

        if (mod == 1)
          rotate (poly, ALPHA_BOUND, maxlognorm, bound_f, bound_g, area,
                  vmin, vmax, 1, 0, 0, u);
        else /* try only the 'keep' best classes modulo 'mod' */
#pragma omp parallel for schedule(dynamic)
          for (int i = 0; i < nbest; i++)
	    rotate (poly, ALPHA_BOUND, maxlognorm, bound_f, bound_g, area,
		    vmin, vmax, mod, c[i].vmod, c[i].wmod, u);
        if (mod != 1)
          free (c);
      }

    /* restore original polynomial */
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], u0, 0, 2);

    /* perform the best rotation */
    gmp_printf ("best rotation: u=%ld v=%ld w=%Zd alpha'=%1.2f (%1.2f,%1.2f) E=%.2e\n",
            bestu, bestv, bestw, best_alphaprime,
            best_alpha, best_var, best_E);

    /* perform the best rotation */
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], 0, bestu, 2);
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], 0, bestv, 1);
    rotate_auxg_z (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], bestw, 0);

    /* recompute the skewness of the best polynomial */
    poly->skew = L2_combined_skewness2 (poly->pols[0], poly->pols[1],
                                        SKEWNESS_DEFAULT_PREC);

    print_cadopoly_extra (stdout, poly, argc0, argv0, 0);

    free (Primes);
    free (Q);
    cado_poly_clear (poly);
    mpz_clear (bestw);

    return 0;
}
