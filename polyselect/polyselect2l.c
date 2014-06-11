/*
  Polynomial selection using Kleinjung's algorithm (cf slides presented
  at the CADO Workshop in October 2008, Nancy, France).

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-nq xxx" denotes the number of special-q's trials for each ad;

  Please report bugs to the Bug Tracking System on:
  https://gforge.inria.fr/tracker/?atid=7442&group_id=2065&func=browse
*/

#define EMIT_ADDRESSABLE_shash_add

#include "cado.h"
#include "polyselect2l.h"
#include "portability.h"
#include "mpz_poly.h"
#include "area.h"

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */
#define NEW_ROOTSIEVE
#define INIT_FACTOR 8UL
#define PREFIX_HASH
//#define DEBUG_POLYSELECT2L

#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#endif

#ifdef PREFIX_HASH
char *phash = "# ";
#else
char *phash = "";
#endif

#define BATCH_SIZE 20 /* number of special (q, r) per batch */
#define KEEP 10       /* number of best raw polynomials kept */

/* Read-Only */
uint32_t *Primes = NULL;
unsigned long lenPrimes = 1; // length of Primes[]
int nq = INT_MAX;
int keep = KEEP;
const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
static int verbose = 0;
int rseffort = DEFAULT_RSEFFORT; /* sieving effort, among 1-5 */
static unsigned long incr = DEFAULT_INCR;
const char *out = NULL; /* output file for msieve input (msieve.dat.m) */
cado_poly best_poly, curr_poly;
double best_E = 0.0; /* Murphy's E (the larger the better) */

/* read-write global variables */
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */
int tot_found = 0; /* total number of polynomials */
int opt_found = 0; /* number of size-optimized polynomials */
int ros_found = 0; /* number of rootsieved polynomials */
double potential_collisions = 0.0, aver_opt_lognorm = 0.0,
  aver_raw_lognorm = 0.0, aver_lognorm_ratio = 0.0,
  var_opt_lognorm = 0.0, var_raw_lognorm = 0.0;
double min_raw_lognorm = 999.99, max_raw_lognorm = 0.0;
double min_opt_lognorm = 999.99, max_opt_lognorm = 0.0;
unsigned long collisions = 0;
unsigned long collisions_good = 0;
double total_adminus2 = 0.0;
double best_raw_logmu[KEEP], best_opt_logmu[KEEP], best_logmu[KEEP + 1];
double rootsieve_time = 0.0;
int raw = 0;

static void
mutex_lock(pthread_mutex_t *lock)
{
#ifdef MAX_THREADS
  pthread_mutex_lock (lock);
#endif
}

static void
mutex_unlock(pthread_mutex_t *lock)
{
#ifdef MAX_THREADS
  pthread_mutex_unlock (lock);
#endif
}

/* inline function */
extern void shash_add (shash_t, uint64_t);

/* -- functions starts here -- */

/* crt, set r and qqz */
void
crt_sq ( mpz_t qqz,
         mpz_t r,
         unsigned long *q,
         unsigned long *rq,
         unsigned long lq )
{
  mpz_t prod, pprod, mod, inv, sum;
  unsigned long i;
  unsigned long qq[lq];

  mpz_init_set_ui (prod, 1);
  mpz_init (pprod);
  mpz_init (mod);
  mpz_init (inv);
  mpz_init_set_ui (sum, 0);

  for (i = 0; i < lq; i ++) {
    qq[i] = q[i] * q[i]; // q small
    mpz_mul_ui (prod, prod, qq[i]);
  }

  for (i = 0; i < lq; i ++) {
    mpz_divexact_ui (pprod, prod, qq[i]);
    mpz_set_ui (mod, qq[i]);
    mpz_invert (inv, pprod, mod);
    mpz_mul_ui (inv, inv, rq[i]);
    mpz_mul (inv, inv, pprod);
    mpz_add (sum, sum, inv);
  }

  mpz_mod (sum, sum, prod);
  mpz_set (r, sum);
  mpz_set (qqz, prod);

  mpz_clear (prod);
  mpz_clear (pprod);
  mpz_clear (mod);
  mpz_clear (inv);
  mpz_clear (sum);
}

/* check that l/2 <= d*m0/P^2, where l = p1 * p2 * q with P <= p1, p2 <= 2P
   q is the product of special-q primes. It suffices to check that
   q <= d*m0/(2P^4). */
static int
check_parameters (mpz_t m0, unsigned long d, unsigned long lq)
{
  double maxq = 1.0, maxP;
  int k = lq;

  while (k > 0)
    maxq *= (double) SPECIAL_Q[LEN_SPECIAL_Q - 1 - (k--)];

  maxP = (double) Primes[lenPrimes - 1];
  if (2.0 * pow (maxP, 4.0) * maxq >= (double) d * mpz_get_d (m0))
    return 0;

  if (maxq > pow (maxP, 2.0))
    return 0;

  return 1;
}


/* print poly info */
void
print_poly_info ( mpz_t *f,
                  const unsigned int d,
                  mpz_t g[2],
                  const mpz_t n,
                  const int raw,
                  const char *prefix )
{
  unsigned int i, nroots;
  double skew, skew2, logmu, alpha, alpha_proj, exp_E;
  mpz_t k0, k1, k2;
  mpz_poly_t F;
  mpz_init (k0);
  mpz_init (k1);
  mpz_init (k2);
  F->coeff = f;
  F->deg = d;

  gmp_printf ("%sn: %Zd\n", prefix, n);
  gmp_printf ("%sY1: %Zd\n%sY0: %Zd\n", prefix, g[1], prefix, g[0]);
  for (i = d + 1; i -- != 0; )
    gmp_printf ("%sc%u: %Zd\n", prefix, i, f[i]);

  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);

  if (d == 6) {
    mpz_set_d (k2, skew);
    rotate_auxg_z (f, g[1], g[0], k2, 2);
    mpz_mul_ui (k1, k2, (unsigned long) skew);
    rotate_auxg_z (f, g[1], g[0], k1, 1);
    mpz_mul_ui (k0, k1, (unsigned long) skew);
    rotate_auxg_z (f, g[1], g[0], k0, 0);
    skew2 = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu = L2_lognorm (F, skew2);
    exp_E = logmu - 0.824 * sqrt (2.0 * exp_rot[d] * log (skew)),
    mpz_neg (k2, k2);
    mpz_neg (k1, k1);
    mpz_neg (k0, k0);
    rotate_auxg_z (f, g[1], g[0], k2, 2);
    rotate_auxg_z (f, g[1], g[0], k1, 1);
    rotate_auxg_z (f, g[1], g[0], k0, 0);
  }
  else {
    mpz_set_d (k1, skew);
    rotate_auxg_z (f, g[1], g[0], k1, 1);
    mpz_mul_ui (k0, k1, (unsigned long) skew);
    rotate_auxg_z (f, g[1], g[0], k0, 0);
    skew2 = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu = L2_lognorm (F, skew2);
    exp_E = logmu - 0.824 * sqrt (2.0 * exp_rot[d] * log (skew)),
    mpz_neg (k1, k1);
    mpz_neg (k0, k0);
    rotate_auxg_z (f, g[1], g[0], k1, 1);
    rotate_auxg_z (f, g[1], g[0], k0, 0);
  }
  
  nroots = numberOfRealRoots (f, d, 0, 0, NULL);
  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm (F, skew);
  alpha = get_alpha (F, ALPHA_BOUND);
  alpha_proj = get_biased_alpha_projective (F, ALPHA_BOUND);

  if (raw == 1)
    printf ("# raw lognorm ");
  else
    printf ("# lognorm ");
  printf ("%1.2f, skew %1.2f, alpha %1.2f (proj: %1.2f), E %1.2f"
          ", exp_E %1.2f, %u rroots\n",
          logmu, skew, alpha, alpha_proj, logmu + alpha,
          exp_E, nroots);

  mpz_clear (k0);
  mpz_clear (k1);
  mpz_clear (k2);
}


/* the number of expected collisions is 8*lenPrimes^2/2/(2P)^2 */
static double
expected_collisions (uint32_t twoP)
{
  double m = (lenPrimes << 1) / (double) twoP;
  return m * m;
}

static void
check_divexact_ui(mpz_t r, const mpz_t d, const char *d_name MAYBE_UNUSED,
                  const unsigned long q, const char *q_name MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_ui_p (d, q) == 0)
  {
    gmp_fprintf (stderr, "Error: %s=%Zd not divisible by %s=%lu\n",
                 d_name, d, q_name, q);
    exit (1);
  }
#endif
  mpz_divexact_ui (r, d, q);
}

static void
check_divexact(mpz_t r, const mpz_t d, const char *d_name MAYBE_UNUSED, const mpz_t q,
               const char *q_name MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (d, q) == 0)
  {
    gmp_fprintf (stderr, "Error: %s=%Zd not divisible by %s=%Zd\n",
                 d_name, d, q_name, q);
    exit (1);
  }
#endif
  mpz_divexact (r, d, q);
}


void rootsieve_poly(mpz_t *g, const unsigned long d,
    const mpz_t N, mpz_poly_t F)
{
  ros_found ++;
  /* root sieve */
#ifndef NEW_ROOTSIEVE
  unsigned long alim = 2000;
  long jmin, kmin;
#endif

  mutex_lock (&lock);
  rootsieve_time -= seconds_thread ();
  mutex_unlock (&lock);

#ifdef NEW_ROOTSIEVE
  if (d > 3) {
    /* verbose = 2 to see details */
    ropt_polyselect (F->coeff, d, g[0], g[1], N, rseffort, 0);
  }
  else {
    unsigned long alim = 2000;
    long jmin, kmin;
    mpz_t m;
    mpz_init(m);
    mpz_neg (m, g[0]);
    rotate (F, alim, m, g[1], &jmin, &kmin, 0, verbose);
    mpz_neg (g[0], m);
    mpz_clear(m);
    /* optimize again, but only translation */
    optimize_aux (F, g, 0, 0);
  }
#else
  mpz_t m;
  mpz_init(m);
  mpz_neg (m, g[0]);
  rotate (F, alim, m, g[1], &jmin, &kmin, 0, verbose);
  mpz_neg (g[0], m);
  mpz_clear(m);
  /* optimize again, but only translation */
  optimize_aux (F, g, 0, 0);
#endif

  mutex_lock (&lock);
  rootsieve_time += seconds_thread ();
  mutex_unlock (&lock);
}

void
output_polynomials(mpz_t *fold, const unsigned long d, mpz_t *gold,
    const mpz_t N, const double logmu0c3, const double logmu0c4, mpz_t *f, mpz_t *g,
    const double E)
{
  mutex_lock (&lock);
  if (fold != NULL && gold != NULL) {
    printf ("# Raw polynomial:\n");
    print_poly_info (fold, d, gold, N, 1, phash);
  }
  if (d == 6 && verbose >= 1)
    gmp_printf ("# noc4/noc3: %.2f/%.2f (%.2f)\n",
                logmu0c4, logmu0c3, logmu0c4/logmu0c3);
  if (raw)
    gmp_printf ("# Size-optimized polynomial:\n");
  else
    gmp_printf ("# Optimized polynomial:\n");
  print_poly_info (f, d, g, N, 0, raw ? "" : phash);
  printf ("# Murphy's E(Bf=%.2e,Bg=%.2e,area=%.2e)=%.2e (best so far %.2e)\n",
          bound_f, bound_g, area, E, best_E);
  printf ("\n");
  fflush (stdout);
  mutex_unlock (&lock);
}

void
output_skipped_poly(const unsigned long d, const double logmu, const uint64_t ad,
    const mpz_t l, const mpz_t g0, const double logmu0c3, const double logmu0c4)
{
  mpz_t m;
  mpz_init(m);
  mpz_neg(m, g0); 
  mutex_lock (&lock);
  if (d == 6)
    gmp_printf ("# Skip polynomial: %.2f, ad: %llu, l: %Zd, m: %Zd, noc4/noc3: %.2f/%.2f (%.2f)\n",
                logmu, (unsigned long long) ad, l, m, logmu0c4, logmu0c3, logmu0c4/logmu0c3);
  else
    gmp_printf ("# Skip polynomial: %.2f, ad: %llu, l: %Zd, m: %Zd\n",
                logmu, (unsigned long long) ad, l, m);
  mutex_unlock (&lock);
  mpz_clear(m);
}


void
output_msieve(const char *out, const unsigned long d, mpz_t *f, mpz_t *g)
{
  FILE *fp;
  unsigned long j;
  mpz_t m;
  mutex_lock (&lock);
  fp = fopen (out, (ros_found == 0) ? "w" : "a");
  if (fp == NULL)
  {
    fprintf (stderr, "Error, cannot open file %s\n", out);
    exit (1);
  }
  fprintf (fp, "0");
  for (j = d + 1; j -- != 0; )
    gmp_fprintf (fp, " %Zd", f[j]);
  mpz_init(m);
  mpz_neg (m, g[0]);
  gmp_fprintf (fp, " %Zd %Zd\n", g[1], m);
  mpz_clear(m);
  fclose (fp);
  mutex_unlock (&lock);
}

/* Insert a value into a sorted array of length len.
   Returns 1 if element was inserted, 0 if it was too big */
int
sorted_insert_double(double *array, const size_t len, const double value)
{
  size_t k;
  int result = 0;
  if (len == 0)
    return 0;
  mutex_lock (&lock);
  if (value < array[len - 1]) {
    for (k = len - 1; k > 0 && value < array[k-1]; k--)
      array[k] = array[k-1];
    array[k] = value;
    result = 1;
  }
  mutex_unlock (&lock);
  return result;
}

/* return 1 if the polynomial is ok and among the best ones,
   otherwise return 0 */
static int
optimize_raw_poly (double *logmu, mpz_poly_t F, mpz_t *g,
                   const unsigned long d, mpz_t N, double *E)
{
  unsigned long j;
  double skew;
  mpz_t t;

  /* check that the algebraic polynomial has content 1, otherwise skip it */
  mpz_init (t);
  mpz_poly_content (t, F);
  if (mpz_cmp_ui (t, 1) != 0)
    {
      mpz_clear (t);
      return 0;
    }
  mpz_clear (t);

  /* optimize size */
  mutex_lock (&lock);
  opt_found ++;
  mutex_unlock (&lock);
  optimize (F, g, verbose, 1);

  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
  *logmu = L2_lognorm (F, skew);

  if (!sorted_insert_double (best_opt_logmu, keep, *logmu))
    return 0; /* not among the best 'keep' ones */

  if (!raw)
    {
      rootsieve_poly (g, d, N, F);
      skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
      *logmu = L2_lognorm (F, skew);
    }

  sorted_insert_double (best_logmu, keep, *logmu);
  mutex_lock (&lock);
  
  collisions_good ++;
  aver_opt_lognorm += *logmu;
  var_opt_lognorm += *logmu * *logmu;
  if (*logmu < min_opt_lognorm)
    min_opt_lognorm = *logmu;
  if (*logmu > max_opt_lognorm)
    max_opt_lognorm = *logmu;

  /* MurphyE */
  mpz_set (curr_poly->n, N);
  curr_poly->rat->deg = 1;
  mpz_set (curr_poly->rat->coeff[0], g[0]);
  mpz_set (curr_poly->rat->coeff[1], g[1]);
  curr_poly->alg->deg = F->deg;
  for (j = d + 1; j -- != 0; )
    mpz_set (curr_poly->alg->coeff[j], F->coeff[j]);
  curr_poly->skew = skew;
  *E = MurphyE (curr_poly, bound_f, bound_g, area, MURPHY_K);

  if (*E > best_E)
  {
    best_E = *E;
    cado_poly_set (best_poly, curr_poly);
  }
  mutex_unlock (&lock);

  return 1;
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
match (unsigned long p1, unsigned long p2, const int64_t i, mpz_t m0,
       uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
       mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], adz;
  int cmp;
  double skew, logmu;
  mpz_poly_t F;

  /* the expected rotation space is S^5 for degree 6 */
#ifdef DEBUG_POLYSELECT2L
  gmp_printf ("Found match: (%lu,%lld) (%lu,%lld) for "
	      "ad=%llu, q=%llu, rq=%Zd\n",
              p1, (long long) i, p2, (long long) i, (unsigned long long) ad,
              (unsigned long long) q, rq);
  gmp_printf ("m0=%Zd\n", m0);
#endif

  mpz_init (l);
  mpz_init (m);
  mpz_init (t);
  mpz_init (k);
  mpz_init (adm1);
  mpz_init (adz);
  mpz_init (mtilde);
  mpz_init (g[0]);
  mpz_init (g[1]);
  mpz_init (gold[0]);
  mpz_init (gold[1]);
  mpz_poly_init (F, d);
  F->deg = d;
  f = F->coeff;
  fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (fold == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in match\n");
    exit (1);
  }
  for (unsigned long j = 0; j <= d; j++)
    mpz_init (fold[j]);
  /* we have l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_mul_ui (l, l, q);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set_si (mtilde, i);
  mpz_mul_ui (mtilde, mtilde, q);
  mpz_mul_ui (mtilde, mtilde, q);
  mpz_add (mtilde, mtilde, rq);
  mpz_add (mtilde, mtilde, m0);
  /* we want mtilde = d*ad*m + a_{d-1}*l with 0 <= a_{d-1} < d*ad.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_set_uint64 (m, ad);
  mpz_mul_ui (m, m, d);
  if (mpz_invert (adm1, l, m) == 0)
  {
    fprintf (stderr, "Error in 1/l mod (d*ad)\n");
    exit (1);
  }
  mpz_mul (adm1, adm1, mtilde);
  mpz_mod (adm1, adm1, m); /* m is d*ad here */

  /* we make -d*ad/2 <= adm1 < d*ad/2 */
  mpz_mul_2exp (t, adm1, 1);
  if (mpz_cmp (t, m) >= 0)
    mpz_sub (adm1, adm1, m);

  mpz_mul (m, adm1, l);
  mpz_sub (m, mtilde, m);
  check_divexact_ui (m, m, "m-a_{d-1}*l", d, "d");

  mpz_set_uint64 (adz, ad);

  check_divexact (m, m, "(m-a_{d-1}*l)/d", adz, "ad");
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set (f[d], adz);
  mpz_pow_ui (t, m, d);
  mpz_mul (t, t, adz);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
  check_divexact (t, t, "t", l, "l");
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (unsigned long j = d - 2; j > 0; j--)
  {
    check_divexact (t, t, "t", l, "l");
    /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
    mpz_pow_ui (mtilde, m, j);
    /* fdiv rounds toward -infinity: adm1 = floor(t/mtilde) */
    mpz_fdiv_q (adm1, t, mtilde); /* t -> adm1 * mtilde + t */
    mpz_invert (k, mtilde, l); /* search adm1 + k such that
                                  t = (adm1 + k) * m^j mod l */
    mpz_mul (k, k, t);
    mpz_sub (k, k, adm1);
    mpz_mod (k, k, l);

    mpz_mul_2exp (k, k, 1);
    cmp = mpz_cmp (k, l);
    mpz_div_2exp (k, k, 1);
    if (cmp >= 0)
      mpz_sub (k, k, l);
    mpz_add (adm1, adm1, k);
    mpz_set (f[j], adm1);
    /* subtract adm1*m^j */
    mpz_submul (t, mtilde, adm1);
  }
  check_divexact (t, t, "t", l, "l");
  mpz_set (f[0], t);

  /* save unoptimized polynomial to fold */
  for (unsigned long j = d + 1; j -- != 0; )
    mpz_set (fold[j], f[j]);
  mpz_set (gold[1], g[1]);
  mpz_set (gold[0], g[0]);

  /* old lognorm */
  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm (F, skew);

  /* for degree 6 polynomials, find bottleneck coefficient */
  double skewtmp = 0.0, logmu0c4 = 0.0, logmu0c3 = 0.0;
  if (d == 6) {
    mpz_set (adz, f[3]);
    mpz_set_ui (f[3], 0);
    skewtmp = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu0c3 = L2_lognorm (F, skewtmp);
    mpz_set (f[3], adz);
    mpz_set (adz, f[4]);
    mpz_set_ui (f[4], 0);
    skewtmp = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu0c4 = L2_lognorm (F, skewtmp);
    mpz_set (f[4], adz);
  }

  double g0 = mpz_get_d (g[0]);
  g0 /= mpz_get_d (f[d-2]);
  g0 = (g0 > 0)? g0 : -g0;

  mutex_lock (&lock);
  /* information on all polynomials */
  total_adminus2 += g0;
  collisions ++;
  tot_found ++;
  aver_raw_lognorm += logmu;
  var_raw_lognorm += logmu * logmu;
  if (d == 6) {
    aver_lognorm_ratio += logmu0c4/logmu0c3;
  }
  if (logmu < min_raw_lognorm)
      min_raw_lognorm = logmu;
  if (logmu > max_raw_lognorm)
    max_raw_lognorm = logmu;
  mutex_unlock (&lock);

  /* if the polynomial has small norm, we optimize it */
  int did_optimize = 0;
  double E;
  if (sorted_insert_double(best_raw_logmu, keep, logmu)) {
    did_optimize = optimize_raw_poly(&logmu, F, g, d, N, &E);
  }

  if (did_optimize && out != NULL)
    output_msieve(out, d, F->coeff, g);
  
  /* print optimized (maybe size- or size-root- optimized) polynomial */
  if (did_optimize && verbose >= 0)
    output_polynomials(fold, d, gold, N, logmu0c3, logmu0c4, F->coeff, g, E);
  
  if (!did_optimize && verbose >= 1)
    output_skipped_poly(d, logmu, ad, l, g[0], logmu0c3, logmu0c4);

  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (t);
  mpz_clear (k);
  mpz_clear (adm1);
  mpz_clear (adz);
  mpz_clear (mtilde);
  mpz_clear (g[0]);
  mpz_clear (g[1]);
  mpz_clear (gold[0]);
  mpz_clear (gold[1]);
  mpz_poly_clear (F);
  for (unsigned long j = 0; j <= d; j++)
    mpz_clear (fold[j]);
  free (fold);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
gmp_match (uint32_t p1, uint32_t p2, int64_t i, mpz_t m0,
	   uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
	   mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], qq, adz, tmp;
  int cmp;
  double skew, logmu;
  mpz_poly_t F;

#ifdef DEBUG_POLYSELECT2L
  gmp_printf ("Found match: (%" PRIu32 ",%lld) (%" PRIu32 ",%lld) for "
	      "ad=%llu, q=%llu, rq=%Zd\n",
              p1, (long long) i, p2, (long long) i, (unsigned long long) ad,
              (unsigned long long) q, rq);
  gmp_printf ("m0=%Zd\n", m0);
#endif
  mpz_init (tmp);
  mpz_init (l);
  mpz_init (m);
  mpz_init (t);
  mpz_init (k);
  mpz_init (qq);
  mpz_init (adm1);
  mpz_init (adz);
  mpz_init (mtilde);
  mpz_init (g[0]);
  mpz_init (g[1]);
  mpz_init (gold[0]);
  mpz_init (gold[1]);
  mpz_poly_init (F, d);
  F->deg = d;
  f = F->coeff;
  fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (fold == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in match\n");
    exit (1);
  }
  for (unsigned long j = 0; j <= d; j++)
    mpz_init (fold[j]);
  /* we have l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_set_uint64 (tmp, q);
  mpz_mul (l, l, tmp);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set (qq, tmp); // qq = q
  mpz_mul (qq, qq, tmp); // qq = q^2
  if (i >= 0)
    mpz_set_uint64 (tmp, (uint64_t) i);
  else {
    mpz_set_uint64 (tmp, (uint64_t) (-i));
    mpz_neg (tmp, tmp);
  }
  mpz_set (mtilde, tmp);
  mpz_mul (mtilde, mtilde, qq);
  mpz_add (mtilde, mtilde, rq);
  mpz_add (mtilde, mtilde, m0);
  /* we want mtilde = d*ad*m + a_{d-1}*l with 0 <= a_{d-1} < d*ad.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_set_uint64 (m, ad);
  mpz_mul_ui (m, m, d);
  if (mpz_invert (adm1, l, m) == 0)
  {
    fprintf (stderr, "Error in 1/l mod (d*ad)\n");
    exit (1);
  }
  mpz_mul (adm1, adm1, mtilde);
  mpz_mod (adm1, adm1, m); /* m is d*ad here */
  /* we make -d*ad/2 <= adm1 < d*ad/2 */
  mpz_mul_2exp (t, adm1, 1);
  if (mpz_cmp (t, m) >= 0)
    mpz_sub (adm1, adm1, m);
  mpz_mul (m, adm1, l);
  mpz_sub (m, mtilde, m);
  check_divexact_ui (m, m, "m-a_{d-1}*l", d, "d");
  mpz_set_uint64 (adz, ad);
  check_divexact (m, m, "(m-a_{d-1}*l)/d", adz, "ad");
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set (f[d], adz);
  mpz_pow_ui (t, m, d);
  mpz_mul (t, t, adz);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
  check_divexact (t, t, "t", l, "l");
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (unsigned long j = d - 2; j > 0; j--)
  {
    check_divexact (t, t, "t", l, "l");
    /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
    mpz_pow_ui (mtilde, m, j);
    mpz_fdiv_q (adm1, t, mtilde); /* t -> adm1 * mtilde + t */
    mpz_invert (k, mtilde, l); /* search adm1 + k such that
                                  t = (adm1 + k) * m^j mod l */
    mpz_mul (k, k, t);
    mpz_sub (k, k, adm1);
    mpz_mod (k, k, l);
    mpz_mul_2exp (k, k, 1);
    cmp = mpz_cmp (k, l);
    mpz_div_2exp (k, k, 1);
    if (cmp >= 0)
      mpz_sub (k, k, l);
    mpz_add (adm1, adm1, k);
    mpz_set (f[j], adm1);
    /* subtract adm1*m^j */
    mpz_submul (t, mtilde, adm1);
  }

  check_divexact (t, t, "t", l, "l");
  mpz_set (f[0], t);

  /* save unoptimized polynomial to fold */
  for (i = d + 1; i -- != 0; )
    mpz_set (fold[i], f[i]);
  mpz_set (gold[1], g[1]);
  mpz_set (gold[0], g[0]);

  /* old lognorm */
  skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm (F, skew);

  /* for degree 6 polynomials, find bottleneck coefficient */
  double skewtmp = 0.0, logmu0c4 = 0.0, logmu0c3 = 0.0;
  if (d == 6) {
    mpz_set (tmp, f[3]);
    mpz_set_ui (f[3], 0);
    skewtmp = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu0c3 = L2_lognorm (F, skewtmp);
    mpz_set (f[3], tmp);
    mpz_set (tmp, f[4]);
    mpz_set_ui (f[4], 0);
    skewtmp = L2_skewness (F, SKEWNESS_DEFAULT_PREC);
    logmu0c4 = L2_lognorm (F, skewtmp);
    mpz_set (f[4], tmp);
  }

  double g0 = mpz_get_d (g[0]);
  g0 /= mpz_get_d (f[d-2]);
  g0 = (g0 > 0)? g0 : -g0;

  mutex_lock (&lock);
  /* information on all polynomials */
  total_adminus2 += g0;
  collisions ++;
  tot_found ++;
  aver_raw_lognorm += logmu;
  var_raw_lognorm += logmu * logmu;
  if (d == 6) {
    aver_lognorm_ratio += logmu0c4/logmu0c3;
  }
  if (logmu < min_raw_lognorm)
      min_raw_lognorm = logmu;
  if (logmu > max_raw_lognorm)
    max_raw_lognorm = logmu;
  mutex_unlock (&lock);

  /* if the polynomial has small norm, we optimize it */
  int did_optimize = 0;
  double E;
  if (sorted_insert_double(best_raw_logmu, keep, logmu)) {
    did_optimize = optimize_raw_poly(&logmu, F, g, d, N, &E);
  }

  if (did_optimize && out != NULL)
    output_msieve(out, d, F->coeff, g);
  
  /* print optimized (maybe size- or size-root- optimized) polynomial */
  if (did_optimize && verbose >= 0)
    output_polynomials(fold, d, gold, N, logmu0c3, logmu0c4, F->coeff, g, E);
  
  if (!did_optimize && verbose >= 1)
    output_skipped_poly(d, logmu, ad, l, g[0], logmu0c3, logmu0c4);

  mpz_clear (tmp);
  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (t);
  mpz_clear (k);
  mpz_clear (qq);
  mpz_clear (adm1);
  mpz_clear (adz);
  mpz_clear (mtilde);
  mpz_clear (g[0]);
  mpz_clear (g[1]);
  mpz_clear (gold[0]);
  mpz_clear (gold[1]);
  mpz_poly_clear (F);
  for (unsigned long j = 0; j <= d; j++)
    mpz_clear (fold[j]);
  free (fold);
}


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long j, nprimes, p, nrp, c = 0, tot_roots = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  double pc1;
  mpz_t zero;
  int found = 0;
  shash_t H;
  int st = 0;

  /* init zero */
  mpz_init_set_ui (zero, 0);

  rp = (uint64_t*) malloc (header->d * sizeof (uint64_t));
  if (rp == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }

  shash_init (H, 4 * lenPrimes);
  umax = (int64_t) Primes[lenPrimes - 1] * (int64_t) Primes[lenPrimes - 1];
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++)
    {
      p = Primes[nprimes];
      ppl = (int64_t) p * (int64_t) p;

      /* add fake roots to keep indices */
      if ((header->d * header->ad) % p == 0)
        {
          R->nr[nprimes] = 0; // nr = 0.
          R->roots[nprimes] = NULL;
          continue;
        }

      st -= milliseconds ();
      nrp = roots_mod_uint64 (rp, mpz_fdiv_ui (header->Ntilde, p), header->d,
                              p);
      st += milliseconds ();
      tot_roots += nrp;
      roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
      proots_add (R, nrp, rp, nprimes);
      for (j = 0; j < nrp; j++, c++)
            {
              for (u = (int64_t) rp[j]; u < umax; u += ppl)
                shash_add (H, u);
              for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
                shash_add (H, -u);
            }
        }
  found = shash_find_collision (H);
  shash_clear (H);
  free (rp);

  if (verbose > 2)
    fprintf (stderr, "# computing %lu p-roots took %dms\n", tot_roots, st);

  if (found) /* do the real work */
    {
      hash_t H;

      hash_init (H, INIT_FACTOR * lenPrimes);
      for (nprimes = 0; nprimes < lenPrimes; nprimes ++)
        {
          nrp = R->nr[nprimes];
          if (nrp == 0)
            continue;
          p = Primes[nprimes];
          ppl = (int64_t) p * (int64_t) p;
          rp = R->roots[nprimes];

          for (j = 0; j < nrp; j++)
            {
              for (u = (int64_t) rp[j]; u < umax; u += ppl)
                hash_add (H, p, u, header->m0, header->ad, header->d,
                          header->N, 1, zero);
              for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
                hash_add (H, p, -u, header->m0, header->ad,
                          header->d, header->N, 1, zero);
            }
        }
#ifdef DEBUG_POLYSELECT2L
      fprintf (stderr, "# collision_on_p took %lums\n", milliseconds () - st);
      fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

#ifdef DEBUG_HASH_TABLE
      fprintf (stderr, "# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
      fprintf (stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif
      hash_clear (H);
    }

  mpz_clear (zero);

  pc1 = expected_collisions (Primes[lenPrimes - 1]);
  pthread_mutex_lock (&lock);
  potential_collisions += pc1;
  pthread_mutex_unlock (&lock);
  return c;
}


/* collision on each special-q, call collision_on_batch_p() */
static inline void
collision_on_each_sq ( header_t header,
                       proots_t R,
                       unsigned long q,
                       mpz_t rqqz,
                       unsigned long *inv_qq )
{
  shash_t H;
  uint64_t **cur1, **cur2, *ccur1, *ccur2;
  long *pc, *epc;
  double pc2;
  uint64_t pp;
  int64_t ppl, neg_umax, umax, v1, v2, nv;
  unsigned long p, nprimes, c;
  uint8_t vpnr, *pnr, nr, j;
  uint32_t *pprimes, i;
  int found;

#ifdef DEBUG_POLYSELECT2L
  int st = milliseconds();
#endif
#if SHASH_NBUCKETS == 256
#define CURRENT(V) (H->current + (uint8_t) (V))
#else
#define CURRENT(V) (H->current + ((V) & (SHASH_NBUCKETS - 1)))
#endif

  shash_init (H, 4 * lenPrimes);

  pc = (long *) inv_qq;
  nv = *pc;
  pprimes = Primes - 1;
  pnr = R->nr;
  R->nr[R->size] = 0xff; /* I use guard to end */
  umax = Primes[lenPrimes - 1];
  umax *= umax;
  neg_umax = -umax;

  /* This define inserts 2 values v1 and v2 with a interlace.
     The goal is to have a little time to read ccurX from L0
     cache before to use it. The best seems a
     three read interlacing in fact, two seems too short. */
#define INSERT_2I(I1,I2)                                                \
  do {                                                                  \
    cur1 = CURRENT(I1); ccur1 = *cur1;					\
    cur2 = CURRENT(I2); ccur2 = *cur2;					\
    *ccur1++ = I1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;	\
    *ccur2++ = I2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;	\
  } while (0)
  /* This version is slow because ccur1 is used immediatly after
     it has been read from L0 cache -> 3 ticks of latency on P4 Nehalem. */
#define INSERT_I(I)						\
  do {								\
    cur1 = CURRENT(I); ccur1 = *cur1; *ccur1++ = I;		\
    __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;		\
  } while (0)

  int64_t b;
  b = (int64_t) ((double) umax * 0.3333333333333333);
  do {
    do {
      vpnr = *pnr++;
      pprimes++;
    } while (!vpnr);
    if (UNLIKELY(vpnr == 0xff)) goto bend;
    ppl = *pprimes;
    __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
    __builtin_prefetch(((void *) pprimes) + 0x80, 0, 3);
    __builtin_prefetch(((void *) pc) + 0x100, 0, 3);
    ppl *= ppl;
    epc = pc + vpnr;
    if (UNLIKELY(ppl > b)) { b = umax >> 1; goto iter2; }
    do {
      v1 = nv;                    cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 = v1 - ppl;              cur2 = CURRENT(v2); ccur2 = *cur2;
      nv = *++pc;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl;                  cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 -= ppl;                  cur2 = CURRENT(v2); ccur2 = *cur2;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl;                  cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 -= ppl;                  cur2 = CURRENT(v2); ccur2 = *cur2;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl; v2 -= ppl;
      if (LIKELY (v1 > umax)) {
	if (UNLIKELY (v2 >= neg_umax)) INSERT_I(v2);
      } else if (UNLIKELY (v2 >= neg_umax)) INSERT_2I(v1, v2);
      else INSERT_I(v1);
    } while (pc != epc);
  } while (1);

  do {
    do {
      vpnr = *pnr++;
      pprimes++;
    } while (!vpnr);
    if (UNLIKELY(vpnr == 0xff)) goto bend;
    ppl = *pprimes;
    __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
    __builtin_prefetch(((void *) pprimes) + 0x100, 0, 3);
    __builtin_prefetch(((void *) pc) + 0x280, 0, 3);
    ppl *= ppl;
    epc = pc + vpnr;
  iter2:
    if (UNLIKELY(ppl > b)) goto iter1;
    do {
      v1 = nv;                    cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 = v1 - ppl;              cur2 = CURRENT(v2); ccur2 = *cur2;
      nv = *++pc;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl;                  cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 -= ppl;                  cur2 = CURRENT(v2); ccur2 = *cur2;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl; v2 -= ppl;
      if (LIKELY (v1 > umax)) {
	if (UNLIKELY (v2 >= neg_umax)) INSERT_I(v2);
      } else if (UNLIKELY (v2 >= neg_umax)) INSERT_2I(v1, v2);
      else INSERT_I(v1);
    } while (pc != epc);
  } while (1);

  do {
    do {
      vpnr = *pnr++;
      pprimes++;
    } while (!vpnr);
    if (UNLIKELY(vpnr == 0xff)) goto bend;
    ppl = *pprimes;
    __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
    __builtin_prefetch(((void *) pprimes) + 0x100, 0, 3);
    __builtin_prefetch(((void *) pc) + 0x280, 0, 3);
    ppl *= ppl;
    epc = pc + vpnr;
  iter1:
    do {
      v1 = nv;                    cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 = v1 - ppl;              cur2 = CURRENT(v2); ccur2 = *cur2;
      nv = *++pc;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl; v2 -= ppl;
      if (LIKELY (v1 > umax)) {
	if (UNLIKELY (v2 >= neg_umax)) INSERT_I(v2);
      } else if (UNLIKELY (v2 >= neg_umax)) INSERT_2I(v1, v2);
      else INSERT_I(v1);
    } while (pc != epc);
  } while (1);

 bend:
#undef INSERT_2I
#undef INSERT_I

  for (i = 0; i < SHASH_NBUCKETS; i++) assert (H->current[i] <= H->base[i+1]);

  found = shash_find_collision (H);
  shash_clear (H);

  if (found) /* do the real work */
    {
      hash_t H;

      hash_init (H, INIT_FACTOR * lenPrimes);

      umax = (int64_t) Primes[lenPrimes - 1] * (int64_t) Primes[lenPrimes - 1];
      for (nprimes = c = 0; nprimes < lenPrimes; nprimes ++)
        {
          p = Primes[nprimes];
          if ((header->d * header->ad) % p == 0)
            continue;
          pp = p * p;
          ppl = (long) pp;
          nr = R->nr[nprimes];
          for (j = 0; j < nr; j++, c++)
            {
              v1 = (long) inv_qq[c];
              for (v2 = v1; v2 < umax; v2 += ppl)
                hash_add (H, p, v2, header->m0, header->ad, header->d,
                          header->N, q, rqqz);
              for (v2 = ppl - v1; v2 < umax; v2 += ppl)
                hash_add (H, p, -v2, header->m0, header->ad, header->d,
                          header->N, q, rqqz);
            }
        }
      hash_clear (H);
    }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %lums\n", milliseconds () - st);
  fprintf (stderr, "# - q hash_alloc (q=%lu): %u\n", q, H->alloc);
#endif

#ifdef DEBUG_HASH_TABLE
  fprintf (stderr, "# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
  fprintf (stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif

  pc2 = expected_collisions (Primes[lenPrimes - 1]);
  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* Given p, rp, q, invqq[], for each rq of q, compute (rp - rq) / q^2 */
static inline void
collision_on_each_sq_r ( header_t header,
                         proots_t R,
                         unsigned long q,
                         mpz_t *rqqz,
                         unsigned long *inv_qq,
                         unsigned long number_pr,
                         int count )
{
  if (count == 0)
    return;

  uint8_t i, nr, *pnr;
  unsigned long nprimes, p, c = 0, rp, rqi;
  int k;
  uint64_t pp;
  unsigned long **tinv_qq = malloc (count * sizeof (unsigned long*));

  if (!tinv_qq)
  {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }
  for (k = 0; k < count; k++) {
    /* number_pr + 1 for guard for pre-load in collision_on_each_sq (nv) */
    tinv_qq[k] = malloc ((number_pr + 1) * sizeof (unsigned long));
    tinv_qq[k][number_pr] = 0;
  }

  int st = milliseconds();
  pnr = R->nr;

  /* for each rp, compute (rp-rq)*1/q^2 (mod p^2) */
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++)
  {
    if (!pnr[nprimes]) continue;
    nr = pnr[nprimes];
    p = Primes[nprimes];
    pp = p*p;

    modulusredcul_t modpp;
    residueredcul_t res_rqi, res_rp, res_tmp;
    modredcul_initmod_ul_raw (modpp, pp);
    modredcul_init (res_rqi, modpp);
    modredcul_init (res_rp, modpp);
    modredcul_init (res_tmp, modpp);

    for (k = 0; k < count; k ++)
    {
      rqi = mpz_fdiv_ui (rqqz[k], pp);
      modredcul_intset_ul (res_rqi, rqi);
      modredcul_intset_ul (res_tmp, inv_qq[nprimes]);
      for (i = 0; i < nr; i ++, c++)
      {
        rp = R->roots[nprimes][i];
        modredcul_intset_ul (res_rp, rp);
        /* rp - rq */
        modredcul_sub (res_rp, res_rp, res_rqi, modpp);
        /* res_rp = (rp - rq) / q[i]^2 */
        modredcul_mul (res_rp, res_rp, res_tmp, modpp);
        tinv_qq[k][c] = modredcul_intget_ul (res_rp);
      }
      c -= nr;
    }
    c += nr;

    modredcul_clear (res_rp, modpp);
    modredcul_clear (res_rqi, modpp);
    modredcul_clear (res_tmp, modpp);
    modredcul_clearmod (modpp);
  }

  if (verbose > 2) {
    fprintf (stderr, "#  substage: batch %d many (rp-rq)*1/q^2 took %lums\n",
             count, milliseconds () - st);
    st = milliseconds();
  }

  /* core function to find collisions */
  for (k = 0; k < count; k ++) {
    collision_on_each_sq (header, R, q, rqqz[k], tinv_qq[k]);
  }

  if (verbose > 2)
    fprintf (stderr, "#  substage: collision-detection %d many rq took %lums\n",
             count, milliseconds () - st);

  for (k = 0; k < count; k++)
    free (tinv_qq[k]);
  free (tinv_qq);
}


/* Next combination */
static inline unsigned int
aux_nextcomb ( unsigned int *ind,
               unsigned int len_q,
               unsigned int *len_nr )
{
  unsigned int i;

  /* bottom change first */
  for (i = len_q - 1; ; i--) {
    if (ind[i] < (len_nr[i] - 1)) {
      ind[i]++;
      return 1;
    }
    else {
      if (i == 0)
        break;
      ind[i] = 0;
    }
  }
  return 0;
}


/* Compute crted rq */
static inline void
aux_return_rq ( qroots_t SQ_R,
                unsigned long *idx_q,
                unsigned int *idx_nr,
                unsigned long k,
                mpz_t qqz,
                mpz_t rqqz,
                unsigned long lq )
{
  unsigned long i, q[k], rq[k];

  /* q and roots */
  for (i = 0; i < k; i ++) {
    q[i] = SQ_R->q[idx_q[i]];
    rq[i] = SQ_R->roots[idx_q[i]][idx_nr[i]];
  }

  /* crt roots */
  crt_sq (qqz, rqqz, q, rq, lq);

  return;
}


/* Consider each rq */
static inline void
collision_on_batch_sq_r ( header_t header,
                          proots_t R,
                          qroots_t SQ_R,
                          unsigned long q,
                          unsigned long *idx_q,
                          unsigned long *inv_qq,
                          unsigned long number_pr,
                          int *curr_nq,
                          unsigned long lq )
{
  int count;
  unsigned int ind_qr[lq]; /* indices of roots for each small q */
  unsigned int len_qnr[lq]; /* for each small q, number of roots */
  unsigned long i;
  mpz_t qqz, rqqz[BATCH_SIZE];

  mpz_init (qqz);
  for (i = 0; i < BATCH_SIZE; i ++)
    mpz_init (rqqz[i]);

  /* initialization indices */
  for (i = 0; i < lq; i ++) {
    ind_qr[i] = 0;
    len_qnr[i] = SQ_R->nr[idx_q[i]];
  }

#if 0
  fprintf (stderr, "q: %lu, ", q);
  for (i = 0; i < lq; i ++)
    fprintf (stderr, "%u ", SQ_R->q[idx_q[i]]);
  fprintf (stderr, ", ");
  for (i = 0; i < lq; i ++)
    fprintf (stderr, "%u ", SQ_R->nr[idx_q[i]]);
  fprintf (stderr, "\n");
#endif

  /* we proceed with BATCH_SIZE many rq for each time */
  i = count = 0;
  int re = 1, num_rq;
  while (re) {
    /* compute BATCH_SIZE such many rqqz[] */
    num_rq = 0;
    for (count = 0; count < BATCH_SIZE; count ++)
    {
        aux_return_rq (SQ_R, idx_q, ind_qr, lq, qqz, rqqz[count], lq);
        re = aux_nextcomb (ind_qr, lq, len_qnr);
        (*curr_nq)++;
        num_rq ++;
        if ((*curr_nq) >= nq)
          re = 0;
        if (!re)
          break;
    }

    /* core function for a fixed qq and several rqqz[] */
    collision_on_each_sq_r (header, R, q, rqqz, inv_qq, number_pr, num_rq);
  }

  mpz_clear (qqz);
  for (i = 0; i < BATCH_SIZE; i ++)
    mpz_clear (rqqz[i]);
}


/* SQ inversion, write 1/q^2 (mod p_i^2) to invqq[i] */
static inline void
collision_on_batch_sq ( header_t header,
                        proots_t R,
                        qroots_t SQ_R,
                        unsigned long q,
                        unsigned long *idx_q,
                        unsigned long number_pr,
                        unsigned long lq )
{
  unsigned nr;
  int curr_nq = 0;
  uint64_t pp;
  unsigned long nprimes, p;
  unsigned long *invqq = malloc (lenPrimes * sizeof (unsigned long));
  if (!invqq) {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

  int st = milliseconds();

  /* Step 1: inversion */
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;
    nr = R->nr[nprimes];
    if (nr == 0)
      continue;
    pp = p * p;

    modulusredcul_t modpp;
    residueredcul_t qq, tmp;
    modredcul_initmod_ul (modpp, pp);
    modredcul_init (qq, modpp);
    modredcul_init (tmp, modpp);

    /* q^2/B (mod pp) */
    modredcul_intset_ul (tmp, q);
    modredcul_sqr (qq, tmp, modpp);
    /* B/q^2 (mod pp) */
    modredcul_intinv (tmp, qq, modpp);
    invqq[nprimes] = modredcul_intget_ul (tmp);

    modredcul_clear (tmp, modpp);
    modredcul_clear (qq, modpp);
    modredcul_clearmod (modpp);
  }

  if (verbose > 2)
    fprintf (stderr, "# stage (1/q^2 inversion) for %lu primes took %lums\n",
             lenPrimes, milliseconds () - st);

  /* Step 2: find collisions on q. */
  int st2 = milliseconds();

  collision_on_batch_sq_r ( header, R, SQ_R, q, idx_q, invqq, number_pr,
                            &curr_nq, lq );
  if (verbose > 2)
    fprintf (stderr, "#  stage (special-q) for %d special-q's took %lums\n",
             curr_nq, milliseconds() - st2);

  free (invqq);
}


/* collision on special-q, call collision_on_batch_sq */
static inline void
collision_on_sq ( header_t header,
                  proots_t R,
                  unsigned long c )
{
  int prod = 1;
  unsigned int i;
  unsigned long j, lq = 0UL;
  qroots_t SQ_R;

  /* init special-q roots */
  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  //qroots_print (SQ_R);

  /* find a suitable lq */
  for (i = 0; i < SQ_R->size; i++) {
    if (prod < nq) {
      if (!check_parameters (header->m0, header->d, lq))
        break;
      prod *= SQ_R->nr[i];
      lq ++;
    }
  }

  /* lq < 8 for the moment */
  if (lq > 7)
    lq = 7;
  if (lq < 1)
    lq = 1;

  unsigned long q, idx_q[lq];
  mpz_t qqz;
  mpz_init (qqz);

  for (j = 0; j < lq; j ++)
    idx_q[j] = j;
  q = return_q_norq (SQ_R, idx_q, lq, qqz);

  /* collision batch */
  collision_on_batch_sq (header, R, SQ_R, q, idx_q, c, lq);

  /* clean */
  mpz_clear (qqz);
  qroots_clear (SQ_R);
  return;
}


// separator between modredc_ul and gmp


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
gmp_collision_on_p ( header_t header,
		     proots_t R )
{
  unsigned long j, nprimes, p, nrp, c = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  double pc1;
  mpz_t zero;

  /* init zero */
  mpz_init_set_ui (zero, 0);

  rp = (uint64_t*) malloc (header->d * sizeof (uint64_t));
  if (rp == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }

  hash_t H;

  hash_init (H, INIT_FACTOR * lenPrimes);

#ifdef DEBUG_POLYSELECT2L
  int st = milliseconds();
#endif

  umax = (int64_t) Primes[lenPrimes - 1] * (int64_t) Primes[lenPrimes - 1];
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    ppl = (int64_t) p * (int64_t) p;

    /* add fake roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    nrp = roots_mod_uint64 (rp, mpz_fdiv_ui (header->Ntilde, p), header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);
    for (j = 0; j < nrp; j++, c++) {
      for (u = (int64_t) rp[j]; u < umax; u += ppl)
        gmp_hash_add (H, p, u, header->m0, header->ad,
                      header->d, header->N, 1, zero);
      for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
        gmp_hash_add (H, p, -u, header->m0, header->ad,
                      header->d, header->N, 1, zero);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %lums\n", milliseconds () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

  hash_clear (H);

  free (rp);
  mpz_clear (zero);

  pc1 = expected_collisions (Primes[lenPrimes - 1]);
  pthread_mutex_lock (&lock);
  potential_collisions += pc1;
  pthread_mutex_unlock (&lock);
  return c;
}


/* collision on each special-q, call collision_on_batch_p() */
static inline void
gmp_collision_on_each_sq ( header_t header,
			   proots_t R,
			   uint64_t q,
			   mpz_t rqqz,
			   uint64_t *inv_qq )
{
  unsigned int nr, j;
  unsigned long nprimes, p, c = 0;
  uint64_t pp;
  int64_t ppl, u, v, umax;
  double pc2;

#ifdef DEBUG_POLYSELECT2L
  int st = milliseconds();
#endif

  hash_t H;

  hash_init (H, INIT_FACTOR * lenPrimes);

  umax = (int64_t) Primes[lenPrimes - 1] * (int64_t) Primes[lenPrimes - 1];
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;

    /* set p, p^2, ppl */
    pp = (uint64_t) p;
    pp *= (uint64_t) p;
    ppl = (int64_t) pp;
    nr = R->nr[nprimes];

    for (j = 0; j < nr; j++, c++)
    {
      u = (int64_t) inv_qq[c];

      for (v = u; v < umax; v += ppl)
        gmp_hash_add (H, p, v, header->m0, header->ad, header->d,
                      header->N, q, rqqz);
      for (v = ppl - u; v < umax; v += ppl)
        gmp_hash_add (H, p, -v, header->m0, header->ad, header->d,
                      header->N, q, rqqz);

    }  // next rp
  } // next p

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %lums\n",
	   milliseconds () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  hash_clear (H);

  pc2 = expected_collisions (Primes[lenPrimes - 1]);
  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* Batch SQ mode */
static inline void
gmp_collision_on_batch_sq ( header_t header,
			    proots_t R,
			    uint64_t *q,
			    mpz_t *qqz,
			    mpz_t *rqqz,
			    unsigned long size,
          unsigned long number_pr )
{
  if (size == 0)
    return;

  unsigned int i, j, nr;
  unsigned long nprimes, p, c = 0;
  uint64_t pp, **invqq, rp;
  mpz_t qprod[size], modpp, tmp, tmp1, tmp2, rpmp;

  mpz_init (tmp);
  mpz_init (tmp1);
  mpz_init (tmp2);
  mpz_init (rpmp);
  mpz_init (modpp);
  for (i = 0; i < size; i++)
    mpz_init (qprod[i]);

  invqq = (uint64_t **) malloc (size * sizeof (uint64_t *));
  if (invqq) {
    for (i = 0; i < size; i++)
      invqq[i] = malloc (number_pr * sizeof (uint64_t));
  }
  else {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

  mpz_set (qprod[0], qqz[0]);
  for (i = 1; i < size; i ++)
    mpz_mul(qprod[i], qqz[i], qprod[i-1]);

  /* Step 1: batch inversion */
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    pp = (uint64_t) p;
    pp *= (uint64_t) p;
    if ((header->d * header->ad) % p == 0)
      continue;
    if (R->nr[nprimes] == 0)
      continue;
    nr = R->nr[nprimes];

    /* inversion */
    mpz_set_uint64 (modpp, pp);
    mpz_invert (tmp1 ,qprod[size-1], modpp);

    /* for each (q, r) \in a batch */
    for (i = size-1; i > 0; i --) {
      mpz_mul(tmp, qprod[i-1], tmp1);
      /* for each rp, compute (rp-rq)*1/q^2 (mod p^2) */
      for (j = 0; j < nr; j ++, c++) {
        rp = R->roots[nprimes][j];
	mpz_set_uint64 (rpmp, rp);
	mpz_sub (rpmp, rpmp, rqqz[i]);
	mpz_mul (rpmp, rpmp, tmp);
	mpz_mod (rpmp, rpmp, modpp);
        invqq[i][c] = mpz_get_uint64 (rpmp);
      }
      mpz_mul (tmp, qqz[i], tmp1);
      mpz_mod (tmp1, tmp, modpp);
      c -= nr;
    }
    /* last q in the batch is in tmp_modul */
    for (j = 0; j < nr; j ++, c ++) {
      rp = R->roots[nprimes][j];
      mpz_set_uint64 (rpmp, rp);
      mpz_sub (rpmp, rpmp, rqqz[0]);
      mpz_fdiv_r (rpmp, rpmp, modpp);
      mpz_mul (tmp2, rpmp, tmp1);
      mpz_mod (tmp2, tmp2, modpp);
      invqq[0][c] = mpz_get_uint64 (tmp2);
    }
  } // next prime p

  /* Step 2: find collisions on q. */
  for (i = 0; i < size; i ++) {
    //int st2 = milliseconds();
    gmp_collision_on_each_sq (header, R, q[i], rqqz[i], invqq[i]);
    //printf ("# outer collision_on_each_sq took %dms\n", milliseconds () - st2);
  }

  for (i = 0; i < size; i++)
    free (invqq[i]);
  free (invqq);
  mpz_clear (tmp);
  mpz_clear (tmp1);
  mpz_clear (tmp2);
  mpz_clear (rpmp);
  mpz_clear (modpp);
  for (i = 0; i < size; i++)
    mpz_clear (qprod[i]);
}


/* collision on special-q, call collisio_on_batch_sq */
static inline void
gmp_collision_on_sq ( header_t header,
		      proots_t R,
		      unsigned long c )
{
  // init special-q roots
  int lq = 2; // fixed for the moment
  qroots_t SQ_R;
  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  // qroots_print (SQ_R);

  unsigned long K = lq, N = SQ_R->size, tot, i, l, idx_q[K];
  uint64_t q[BATCH_SIZE];
  mpz_t *qqz, *rqqz;

  qqz = (mpz_t*) malloc (BATCH_SIZE * sizeof (mpz_t));
  rqqz = (mpz_t*) malloc (BATCH_SIZE * sizeof (mpz_t));
  if (!qqz || !rqqz) {
    fprintf (stderr, "Error, cannot allocate memory "
	     "in gmp_collision_on_sq \n");
    exit (1);
  }
  for (l = 0; l < BATCH_SIZE; l++) {
    mpz_init (qqz[l]);
    mpz_init (rqqz[l]);
  }

  // less than lq special primes having roots for this ad
  if (N == 0 || N < K) {
    fprintf (stderr, "# Info: binomial(%lu, %lu) error in "
             "collision_on_sq(). ad=%" PRIu64 ".\n", N, K, header->ad);
    return;
  }

  tot =  binom (N, K);

  if (tot > (unsigned long) nq)
    tot = (unsigned long) nq;

  if (tot < BATCH_SIZE)
    tot = BATCH_SIZE;

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# Info: n=%lu, k=%lu, (n,k)=%lu"
	   ", maxnq=%d, nq=%lu\n", N, K, binom(N, K), nq, tot);
#endif

  i = 0;
  while ( i <= (tot-BATCH_SIZE) ) {

    l = i; // why do I use an extra l here?
    if (l == 0) {

      // enumerate first combination
      first_comb (K, idx_q);
      //print_comb (K, idx_q);
      q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l], lq);

      for (l = 1; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l], lq);
      }
    }
    else {
      for (l = 0; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l], lq);
      }
    }

#ifdef DEBUG_POLYSELECT2L
    unsigned long j;
    for (j = 0; j < BATCH_SIZE; j++)
      gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n",
		   q[j], qqz[j], rqqz[j]);
#endif

    // collision batch
    gmp_collision_on_batch_sq (header, R, q, qqz, rqqz, BATCH_SIZE, c);
    i += BATCH_SIZE;
  }

  // tail batch
  for (l = 0; l < (tot % BATCH_SIZE); l++) {
    next_comb (N, K, idx_q);
    q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l], lq);

#ifdef DEBUG_POLYSELECT2L
    gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n",
		 q[l], qqz[l], rqqz[l]);
#endif

  }

  gmp_collision_on_batch_sq (header, R, q, qqz, rqqz, tot % BATCH_SIZE, c);

  for (l = 0; l < BATCH_SIZE; l++) {
    mpz_clear (qqz[l]);
    mpz_clear (rqqz[l]);
  }
  free (qqz);
  free (rqqz);
  qroots_clear (SQ_R);
}

static void
newAlgo (mpz_t N, unsigned long d, uint64_t ad)
{
  unsigned long c = 0;
  header_t header;
  proots_t R;

  header_init (header, N, d, ad);
  proots_init (R, lenPrimes);

  if (sizeof (unsigned long int) == 8) {
    c = collision_on_p (header, R);
    if (nq > 0)
      collision_on_sq (header, R, c);
  }
  else {
    c = gmp_collision_on_p (header, R);
    if (nq > 0)
      gmp_collision_on_sq (header, R, c);
  }

  proots_clear (R, lenPrimes);
  header_clear (header);
}

int
read_mpz(mpz_t result, const char *line, const char prefix, const unsigned int index)
{
  char pattern[64];
  snprintf(pattern, sizeof(pattern)/sizeof(char), "%c%u: %%Zd\n", prefix, index);
  return gmp_sscanf(line, pattern, result);
}

/* Read polynomials from a file, separated by empty lines,
   optimize them (size and roots) and print them */
void
read_raw_poly_file(const char *filename)
{
  char line[MAX_LINE_LENGTH];
  mpz_t g[2];
  mpz_poly_t F;
  mpz_t N;
  FILE *file;
  const int max_degree = 10;
  
  file = fopen(filename, "ra");
  if (file == NULL) {
    perror("Could not open file");
    exit(EXIT_FAILURE);
  }
  mpz_init(g[0]);
  mpz_init(g[1]);
  mpz_init(N);
  mpz_poly_init (F, max_degree);
  F->deg = 0;
  memset(line, 0, MAX_LINE_LENGTH); /* For nicer gdb output */

  while (!feof(file)) {
    line[0] = 0;
    if (fgets(line, MAX_LINE_LENGTH, file) == NULL) {
      /* Don't break here, in fact, as we want to output the polynomial
         when we hit the EOF on the input */
    }
    /*  Empty line separates polynomials. At EOF, we also output */
    if (line[0] == 0 || strcmp(line, "\n") == 0) {
      /* Did we get all required fields? */
      if (mpz_sgn(N) && mpz_sgn(g[0]) && mpz_sgn(g[1]) && F->deg > 0) {
        double logmu, E;
        /* Optimize and, if good enough, print */
        if (optimize_raw_poly(&logmu, F, g, F->deg, N, &E)) {
          output_polynomials(NULL, F->deg, NULL, N, 0., 0., F->coeff, g, E);
        }
        mpz_set_ui(N, 0);
        for (int i = 0; i < 2; i++)
          mpz_set_ui(g[i], 0);
        F->deg = 0;
      }
      continue;
    }
    gmp_sscanf(line, "n: %Zd\n", N);
    for (int i = 0; i < 2; i++)
      read_mpz(g[i], line, 'Y', i);
    for (int i = 0; i < max_degree; i++) {
      if (read_mpz(F->coeff[i], line, 'c', i))
        F->deg = MAX(F->deg, i);
    }
  }
  
  mpz_clear(g[0]);
  mpz_clear(g[1]);
  mpz_clear(N);
  mpz_poly_clear (F);
  fclose(file); 
}


void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
  return NULL;
}

static void
declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "degree", "(required, alias d) polynomial degree (max=6)");
  param_list_decl_usage(pl, "n", "(required, alias N) input number");
  param_list_decl_usage(pl, "P", "(required) deg-1 coeff of g(x) has two prime factors in [P,2P]\n");

  param_list_decl_usage(pl, "admax", "max value for ad");
  param_list_decl_usage(pl, "admin", "min value for ad (default 0)");
  param_list_decl_usage(pl, "incr", "(alias i) forced factor of ad (default 60)");
  param_list_decl_usage(pl, "maxtime", "stop the search after maxtime seconds");

  char str[200];
  snprintf (str, 200, "maximum number of special-q's considered\n"
            "               for each ad (default %d)", INT_MAX);
  param_list_decl_usage(pl, "nq", str);
  param_list_decl_usage(pl, "keep", "number of polynomials kept (default 10)");
  param_list_decl_usage(pl, "out", "filename for msieve-format output");
  param_list_decl_usage(pl, "r", "(switch) size-optimize polynomial only (skip root-optimization)");
  param_list_decl_usage(pl, "resume", "resume state from given file");
  param_list_decl_usage(pl, "rootsieve", "root-sieve the size-optimized polynomials in given file");
  snprintf (str, 200, "root-sieve effort ranging from 1 to 5 (default %d)",
            DEFAULT_RSEFFORT);
  param_list_decl_usage(pl, "rseffort", str);
  snprintf(str, 200, "time interval (seconds) for printing statistics (default %d)", TARGET_TIME / 1000);
  param_list_decl_usage(pl, "s", str);
  param_list_decl_usage(pl, "save", "save state in given file");
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
  param_list_decl_usage(pl, "q", "(switch) quiet mode");
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
}

static void
usage (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
        missing);
  }
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  int argc0 = argc;
  char **argv0 = argv;
  const char *save = NULL, *resume = NULL, *rootsieve_filename = NULL;
  double st0 = seconds (), maxtime = DBL_MAX;
  mpz_t N;
  unsigned int d = 0;
  unsigned long P = 0, admin, admax;
  double admin_d, admax_d;
  int quiet = 0, tries = 0, i, nthreads = 1, st,
    target_time = TARGET_TIME, incr_target_time = TARGET_TIME;
  tab_t *T;
  FILE *fp;
#ifdef MAX_THREADS
  pthread_t tid[MAX_THREADS];
#endif

  mpz_init (N);
  cado_poly_init (best_poly);
  cado_poly_init (curr_poly);

  /* read params */
  param_list pl;
  param_list_init (pl);

  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-r", &raw);
  param_list_configure_switch (pl, "-q", &quiet);
  param_list_configure_alias(pl, "degree", "-d");
  param_list_configure_alias(pl, "incr", "-i");
  param_list_configure_alias(pl, "n", "-N");

  if (argc == 1)
    usage (argv0[0], NULL, pl);

  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0[0], NULL, pl);
  }

  /* initialize best norms */
  for (i = 0; i < keep; i++)
    {
      best_raw_logmu[i] = 999.99; /* best logmu before size optimization */
      best_opt_logmu[i] = 999.99;   /* best logmu after size optimization */
      best_logmu[i] = 999.99;       /* best logmu after rootsieve */
    }

  /* These parameters need to be parsed even for -rootsieve */
  if (param_list_parse_double (pl, "area", &area) == 0) /* no -area */
    area = AREA;
  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;

  /* sieving effort that passed to ropt */
  param_list_parse_int (pl, "rseffort", &rseffort);
  if (rseffort < 1 || rseffort > 5)
  {
    fprintf (stderr, "Error, -rseffort should be in [1,5]\n");
    exit (1);
  }

  /* filename for doing rootsieve only */
  rootsieve_filename = param_list_lookup_string (pl, "rootsieve");
  if (rootsieve_filename != NULL) {
    read_raw_poly_file(rootsieve_filename);
    goto print_statistics;
  }

  /* parse and check N in the first place */
  int have_n = param_list_parse_mpz(pl, "n", N);

  if (!have_n) {
    fprintf(stderr, "# Reading n from stdin\n");
    param_list_read_stream(pl, stdin);
    have_n = param_list_parse_mpz(pl, "n", N);
  }

  if (!have_n) {
      fprintf(stderr, "No n defined ; sorry.\n");
      exit(1);
  }

  if (mpz_cmp_ui (N, 0) <= 0) usage(argv0[0], "n", pl);

  param_list_parse_ulong(pl, "P", &P);
  if (P == 0) usage(argv0[0], "P", pl);

  param_list_parse_int (pl, "t", &nthreads);
  param_list_parse_int (pl, "nq", &nq);
  param_list_parse_int (pl, "keep", &keep);
  if (keep <= 0 || keep > KEEP)
    {
      fprintf (stderr, "Error, keep should be in [1,%d]\n", KEEP);
      exit (1);
    }
  param_list_parse_int (pl, "s", &target_time);
  incr_target_time = target_time;
  param_list_parse_uint (pl, "degree", &d);
  if (param_list_parse_double (pl, "admin", &admin_d) == 0) /* no -admin */
    admin = 0;
  else
    admin = (unsigned long) admin_d;
  if (param_list_parse_double (pl, "admax", &admax_d) == 0) /* no -admax */
    admax = ULONG_MAX;
  else
    admax = (unsigned long) admax_d;
  param_list_parse_ulong (pl, "incr", &incr);
  param_list_parse_double (pl, "maxtime", &maxtime);
  save = param_list_lookup_string (pl, "save");
  resume = param_list_lookup_string (pl, "resume");
  out = param_list_lookup_string (pl, "out");

  if (param_list_warn_unused(pl))
    usage (argv0[0], NULL, pl);

  /* print command line */
  param_list_print_command_line (stdout, pl);

  /* check degree */
  if (d <= 0) usage(argv0[0], "degree", pl);

  /* check lq and nq */
  if (nq < 0) {
    fprintf (stderr, "Error, number of special-q's should >= 0\n");
    exit (1);
  }

  /* check nthreads */
#ifdef MAX_THREADS
  if (nthreads > MAX_THREADS) {
    fprintf (stderr, "Error, nthreads should be <= %d\n", MAX_THREADS);
    exit (1);
  }
#endif

  /* quiet mode */
  if (quiet == 1)
    verbose = -1;

  /* set cpoly */
  mpz_set (best_poly->n, N);
  mpz_set (curr_poly->n, N);
  best_poly->alg->deg = d;
  best_poly->rat->deg = 1;
  curr_poly->alg->deg = d;
  curr_poly->rat->deg = 1;

  /* if resume, read admin */
  if (resume != NULL)
  {
    fp = fopen (resume, "r");
    if (fp == NULL)
    {
      fprintf (stderr, "Cannot open resume file %s\n", resume);
      exit (1);
    }
    if (fscanf (fp, "%lu", &admin) != 1)
    {
      fprintf (stderr, "Cannot read ad value from resume file %s\n",
               resume);
      exit (1);
    }
    fclose (fp);
  }

  /* init primes */
  double Pd;
  Pd = (double) P;
  if (Pd > (double) UINT_MAX) {
    fprintf (stderr, "Error, too large value of P\n");
    exit (1);
  }
  if (P <= (unsigned long) SPECIAL_Q[LEN_SPECIAL_Q - 2]) {
    fprintf (stderr, "Error, too small value of P\n");
    exit (1);
  }
  /* since for each prime p in [P,2P], we convert p^2 to int64_t, we need
     (2P)^2 < 2^63, thus P < 2^30.5 */
  if (P >= 1518500250UL)
    {
      fprintf (stderr, "Error, too large value of P\n");
      exit (1);
    }

  /* detect L1 cache size */
  ropt_L1_cachesize ();

  st = milliseconds ();
  lenPrimes = initPrimes (P, &Primes);

  printf ( "# Info: initializing %lu P primes took %lums,"
           " rawonly=%d, nq=%d, target_time=%d\n",
           lenPrimes,
           milliseconds () - st,
           raw,
           nq,
           target_time / 1000 );


  printf ( "# Info: estimated peak memory=%.2fMB (%d thread(s),"
           " batch %d inversions on SQ)\n",
           (double) (nthreads * (BATCH_SIZE * 2 + INIT_FACTOR) * lenPrimes
           * (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           BATCH_SIZE );

  //printPrimes (Primes, lenPrimes);

  /* init tabs_t for threads */
  T = malloc (nthreads * sizeof (tab_t));
  if (T == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in main\n");
    exit (1);
  }
  for (i = 0; i < nthreads ; i++)
  {
    mpz_init_set (T[i]->N, N);
    T[i]->d = d;
    T[i]->thread = i;
  }

  if (incr <= 0)
  {
    fprintf (stderr, "Error, incr should be positive\n");
    exit (1);
  }

  /* force admin to be divisible by incr */
  if (admin == 0)
    admin = incr;
  admin = ((admin + incr - 1) / incr) * incr; /* incr * ceil (admin/incr) */

  while (admin <= admax && seconds () - st0 <= maxtime)
  {
    for (i = 0; i < nthreads ; i++)
    {
      tries ++;
      if (verbose >= 1)
      {
        gmp_printf ("# %d ad=%lu\n", tries, admin);
      }
      T[i]->ad = admin;
#ifndef MAX_THREADS
      newAlgo (N, d, admin);
#else
      pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
#endif
      admin += incr;
    }
#ifdef MAX_THREADS
    for (i = 0 ; i < nthreads ; i++)
      pthread_join(tid[i], NULL);
#endif

    if (save != NULL)
    {
      fp = fopen (save, "w");
      if (fp == NULL)
      {
        fprintf (stderr, "Cannot open save file %s\n", save);
        exit (1);
      }
      if (fprintf (fp, "%lu\n", admin) < 0)
      {
        fprintf (stderr, "Cannot print ad value to save file %s\n",
                 save);
        exit (1);
      }
      fclose (fp);
    }

    if (milliseconds () > (unsigned long) target_time || verbose > 0)
    {
      double mean = aver_opt_lognorm / collisions_good;
      double rawmean = aver_raw_lognorm / collisions;

      printf ("# Stat: ad=%lu, exp. coll.=%1.2f (%0.2e/s), got %lu with %lu good ones, av. lognorm=%1.2f (min=%1.2f,std=%1.2f), av. raw. lognorm=%1.2f (min=%1.2f,std=%1.2f), time=%lums\n",
              admin,
              potential_collisions,
              1000.0 * (double) potential_collisions / milliseconds (),
              collisions,
              collisions_good,
              mean, min_opt_lognorm,
              sqrt (var_opt_lognorm / collisions_good - mean * mean),
              rawmean, min_raw_lognorm,
              sqrt (var_raw_lognorm / collisions - rawmean * rawmean),
              milliseconds () );
      fflush (stdout);
      target_time += incr_target_time;
    }
  }

  /* finishing up statistics */
  if (verbose >= 0)
    {
      printf ("# Stat: potential collisions=%1.2f (%1.2e/s)\n",
              potential_collisions, 1000.0 * potential_collisions
              / (double) milliseconds ());
      if (collisions > 0)
        {
          double mean = aver_opt_lognorm / collisions_good;
          double rawmean = aver_raw_lognorm / collisions;

          printf ("# Stat: raw lognorm (nr/min/av/max/std): %lu/%1.2f/%1.2f/%1.2f/%1.2f\n",
                  collisions, min_raw_lognorm, rawmean, max_raw_lognorm,
                  sqrt (var_raw_lognorm / collisions - rawmean * rawmean));
          if (collisions_good > 0)
            printf ("# Stat: optimized lognorm (nr/min/av/max/std): %lu/%1.2f/%1.2f/%1.2f/%1.2f\n",
                    collisions_good, min_opt_lognorm, mean, max_opt_lognorm,
                    sqrt (var_opt_lognorm / collisions_good - mean * mean));
          printf ("# Stat: av. g0/adm2 ratio: %.3e\n",
                  total_adminus2 / (double) collisions);
          if (d == 6)
            printf ("# Stat: av. logmu noc4/noc3 ratio: %.3f\n",
                    aver_lognorm_ratio / (double) collisions);
        }
    }

  printf ("# Stat: tried %d ad-value(s), found %d polynomial(s), %d size-optimized, %d rootsieved\n",
          tries, tot_found, opt_found, ros_found);

  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  clearPrimes (&Primes);

print_statistics:
  /* print best keep values of logmu */
  if (collisions_good > 0)
    {
      printf ("# Stat: best raw logmu:");
      for (i = 0; i < keep; i++)
        printf (" %1.2f", best_raw_logmu[i]);
      printf ("\n");
      printf ("# Stat: best opt logmu:");
      for (i = 0; i < keep; i++)
        printf (" %1.2f", best_opt_logmu[i]);
      printf ("\n");
      printf ("# Stat: best logmu:");
      for (i = 0; i < keep; i++)
        printf (" %1.2f", best_logmu[i]);
      printf ("\n");
    }

  /* print total time (this gets parsed by the scripts) */
  printf ("# Stat: total phase took %.2fs\n", seconds () - st0);
#ifndef HAVE_RUSAGE_THREAD /* rootsieve_time is correct only if RUSAGE_THREAD
                              works or in mono-thread mode */
  if (nthreads == 1)
#endif
    printf ("# Stat: rootsieve took %.2fs\n", rootsieve_time);

  if (raw) {
    /* We did not optimize roots, so Murphy_E is almost meaningless.
       Don't print the "best" polynomial because it probably isn't */
  } else if (best_E == 0.0) {
    /* This line is required by the script: */
    printf ("# No polynomial found, please increase the ad range or decrease P\n");
  } else {
    /* This line is required by the script: */
    printf ("# Best polynomial found:\n");
    print_cadopoly_extra (stdout, best_poly, argc0, argv0, st0);
  }

  mpz_clear (N);
  cado_poly_clear (best_poly);
  cado_poly_clear (curr_poly);
  param_list_clear (pl);

  return 0;
}
