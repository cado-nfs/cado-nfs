/*
  polyselect2l.c is a variant of Paul Zimmermann's polyselect2.c

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-nq xxx" denotes the number of special-q's trials for each ad;

  "-maxnorm xxx" only optimize raw polynomials with size <= xxx.
  If the raw polynomial is not good enough, we will still stream
  it to STDERR for further reference.

  Please report bugs to shi.bai AT anu.edu.au.
*/

#define EMIT_ADDRESSABLE_shash_add

#include "polyselect2l.h"

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */
#define NEW_ROOTSIEVE
#define MAX_THREADS 16
#define INIT_FACTOR 7
// #define DEBUG_POLYSELECT2L

#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#endif

#define BATCH_SIZE 20 /* number of special (q, r) per batch */

/* Read-Only */
uint32_t *Primes = NULL;
unsigned long lenPrimes = 1; // length of Primes[]
int nq = INT_MAX;
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
static int verbose = 0;
static unsigned long incr = DEFAULT_INCR;
const char *out = NULL; /* output file for msieve input (msieve.dat.m) */
cado_poly best_poly, curr_poly;
double best_E = 0.0; /* Murphy's E (the larger the better) */

int argc0;
char ** argv0;
double st0;

/* read-write global variables */
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */
int tot_found = 0; /* total number of polynomials */
int found = 0; /* number of polynomials below maxnorm */
double potential_collisions = 0.0, aver_opt_lognorm = 0.0,
  aver_raw_lognorm = 0.0, aver_lognorm_ratio = 0.0,
  var_opt_lognorm = 0.0, var_raw_lognorm = 0.0;
double min_raw_lognorm = DBL_MAX, max_raw_lognorm = 0.0;
double min_opt_lognorm = DBL_MAX, max_opt_lognorm = 0.0;
unsigned long collisions = 0;
unsigned long collisions_good = 0;
double total_adminus2 = 0.0;
double best_logmu[11];
double rootsieve_time = 0.0;
int raw = 0;

static inline uint64_t cputicks()
{
        uint64_t r;
        __asm__ __volatile__(
                "rdtsc\n\t"
                "shlq $32, %%rdx\n\t"
                "orq %%rdx, %%rax\n\t"
                : "=a"(r)
                :
                : "rdx");
        return r;
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

  return 1;
}

/* print poly info */
void
print_poly_info ( mpz_t *f,
                  unsigned int d,
                  mpz_t g[2],
                  int raw )
{
  unsigned int i, nroots;
  double skew, logmu, alpha;

  gmp_printf ("Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
  for (i = d + 1; i -- != 0; )
    gmp_printf ("c%u: %Zd\n", i, f[i]);

  nroots = numberOfRealRoots (f, d, 0, 0, NULL);
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  if (raw == 1)
    printf ("# raw lognorm ");
  else
    printf ("# lognorm ");
  printf ("%1.2f, skew %1.2f, alpha %1.2f, E %1.2f,  exp_E %1.2f, %u rroots\n",
          logmu, skew, alpha, logmu + alpha,
          logmu - 0.824 * sqrt (2.0 * exp_rot[d] * log (skew)),
          nroots);
}

/* the number of expected collisions is 8*lenPrimes^2/2/(2P)^2 */
static double
expected_collisions (uint32_t twoP)
{
  double m = (lenPrimes << 1) / (double) twoP;
  return m * m;
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
       mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], adz;
  unsigned long j;
  int cmp;
  double skew, logmu, E;
  /* the expected rotation space is S^5 for degree 6 */
#ifdef DEBUG_POLYSELECT2L
  gmp_printf ("Found match: (%lu,%"PRId64") (%lu,%"PRId64") for "
	      "ad=%"PRIu64", q=%"PRIu64", rq=%Zd\n",
              p1, i, p2, i, ad, q, rq);
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
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (f == NULL || fold == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in match\n");
    exit (1);
  }
  for (j = 0; j <= d; j++) {
    mpz_init (f[j]);
    mpz_init (fold[j]);
  }
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
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_ui_p (m, d) == 0)
  {
    fprintf (stderr, "Error: m-a_{d-1}*l not divisible by d\n");
    exit (1);
  }
#endif
  mpz_divexact_ui (m, m, d);

  mpz_set_uint64 (adz, ad);

#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (m, adz) == 0)
  {
    fprintf (stderr, "Error: (m-a_{d-1}*l)/d not divisible by ad\n");
    exit (1);
  }
#endif
  mpz_divexact (m, m, adz);
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set (f[d], adz);
  mpz_pow_ui (t, m, d);
  mpz_mul (t, t, adz);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif

  mpz_divexact (t, t, l);
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (j = d - 2; j > 0; j--)
  {
#ifdef DEBUG_POLYSELECT2L
    if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
#endif
    mpz_divexact (t, t, l);
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
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif
  mpz_divexact (t, t, l);
  mpz_set (f[0], t);

  /* save unoptimized polynomial to fold */
  for (i = d + 1; i -- != 0; )
    mpz_set (fold[i], f[i]);
  mpz_set (gold[1], g[1]);
  mpz_set (gold[0], g[0]);

  /* old lognorm */
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

  /* raw poly */
  if (logmu <= max_norm)
    {
      printf ("# Raw polynomial:\n");
      gmp_printf ("n: %Zd\n", N);
      print_poly_info (fold, d, gold, 1);
    }

  /* for degree 6 polynomials, find bottleneck coefficient */
  double skewtmp = 0.0, logmu0c4 = 0.0, logmu0c3 = 0.0;
  if (d == 6) {
    mpz_set (adz, f[3]);
    mpz_set_ui (f[3], 0);
    skewtmp = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu0c3 = L2_lognorm (f, d, skewtmp, DEFAULT_L2_METHOD);
    mpz_set (f[3], adz);
    mpz_set (adz, f[4]);
    mpz_set_ui (f[4], 0);
    skewtmp = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu0c4 = L2_lognorm (f, d, skewtmp, DEFAULT_L2_METHOD);
    mpz_set (f[4], adz);
  }

  double g0 = mpz_get_d (g[0]);
  g0 /= mpz_get_d (f[d-2]);
  g0 = (g0 > 0)? g0 : -g0;

#ifdef MAX_THREADS
  pthread_mutex_lock (&lock);
#endif
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
#ifdef MAX_THREADS
  pthread_mutex_unlock (&lock);
#endif

  fprintf(stderr, "# Found match with norm = %.1f\n", logmu);

  double preopt_logmu = logmu;

  /* if the polynomial has norm < "-maxnorm", we optimize it */
  if (logmu <= max_norm)
  {

    /* optimize size */
    optimize (f, d, g, 0, 1);

    if (!raw) {

/* root sieve */
#ifndef NEW_ROOTSIEVE
      unsigned long alim = 2000;
      long jmin, kmin;
#endif
      mpz_neg (m, g[0]);

      rootsieve_time -= seconds_thread ();

#ifdef NEW_ROOTSIEVE
      if (d > 3) {
        ropt_polyselect (f, d, m, g[1], N, 0); // verbose = 2 to see details.
        mpz_neg (g[0], m);
      }
      else {
        unsigned long alim = 2000;
        long jmin, kmin;
        rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
        mpz_neg (g[0], m);
        /* optimize again, but only translation */
        optimize_aux (f, d, g, 0, 0, CIRCULAR);
      }
#else
      rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
      mpz_neg (g[0], m);
      /* optimize again, but only translation */
      optimize_aux (f, d, g, 0, 0, CIRCULAR);
#endif

      rootsieve_time += seconds_thread ();

    } // raw and sopt only ?

    skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

    for (i = 10; i > 0 && logmu < best_logmu[i-1]; i--)
      best_logmu[i] = best_logmu[i-1];
    best_logmu[i] = logmu;
#ifdef MAX_THREADS
    pthread_mutex_lock (&lock);
#endif
    collisions_good ++;
    aver_opt_lognorm += logmu;
    var_opt_lognorm += logmu * logmu;
    if (logmu < min_opt_lognorm)
      min_opt_lognorm = logmu;
    if (logmu > max_opt_lognorm)
      max_opt_lognorm = logmu;
#ifdef MAX_THREADS
    pthread_mutex_unlock (&lock);
#endif

    /* MurphyE */
    mpz_set (curr_poly->rat->f[0], g[0]);
    mpz_set (curr_poly->rat->f[1], g[1]);
    for (j = d + 1; j -- != 0; )
      mpz_set (curr_poly->alg->f[j], f[j]);
    curr_poly->skew = skew;
    E =  MurphyE (curr_poly, BOUND_F, BOUND_G, AREA, MURPHY_K);

    mpz_neg (m, g[0]);

#ifdef MAX_THREADS
		  pthread_mutex_lock (&lock);
#endif
    if (E > best_E)
    {
      best_E = E;
      cado_poly_set (best_poly, curr_poly);
      double alpha = get_alpha (curr_poly->alg->f, curr_poly->alg->degree, ALPHA_BOUND);
      fprintf(stderr, "# Best polynomial found so far, logmu=%.2f (was %.2f), alpha=%.2f, logmu+alpha=%.2f\n", logmu, preopt_logmu, alpha, logmu+alpha);
      print_cadopoly_extra (stdout, best_poly, argc0, argv0, st0, 1 /* raw */);
    }
    if (out != NULL) /* msieve output */
    {
      FILE *fp;
      fp = fopen (out, (found == 0) ? "w" : "a");
      if (fp == NULL)
      {
        fprintf (stderr, "Error, cannot open file %s\n", out);
        exit (1);
      }
      fprintf (fp, "0");
      for (j = d + 1; j -- != 0; )
        gmp_fprintf (fp, " %Zd", f[j]);
      mpz_neg (m, g[0]);
      gmp_fprintf (fp, " %Zd %Zd\n", g[1], m);
      fclose (fp);
    }
    found ++;
#ifdef MAX_THREADS
		  pthread_mutex_unlock (&lock);
#endif

    /* print optimized (maybe size- or size-root- optimized) polynomial */
    if (verbose >= 0)
      {
        if (d == 6 && verbose >= 1)
          gmp_printf ("# noc4/noc3: %.2f/%.2f (%.2f)\n",
                      logmu0c4, logmu0c3, logmu0c4/logmu0c3);
        gmp_printf ("# Optimized polynomial:\n");
        gmp_printf ("n: %Zd\n", N);
        print_poly_info (f, d, g, 0);
        printf ("# Murphy's E(Bf=%.0f,Bg=%.0f,area=%.2e)=%1.2e (best so far %1.2e)\n",
                BOUND_F, BOUND_G, AREA, E, best_E);
        printf ("\n");
        fflush (stdout);
      }
  }
  else {
    if (verbose >= 1) {
      if (d == 6)
        gmp_printf ("# Skip polynomial: %.2f, ad: %"PRIu64", l: %Zd, m: %Zd, noc4/noc3: %.2f/%.2f (%.2f)\n",
                    logmu, ad, l, m, logmu0c4, logmu0c3, logmu0c4/logmu0c3);
      else
        gmp_printf ("# Skip polynomial: %.2f, ad: %"PRIu64", l: %Zd, m: %Zd\n",
                    logmu, ad, l, m);
    }
  }

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
  for (j = 0; j <= d; j++) {
    mpz_clear (f[j]);
    mpz_clear (fold[j]);
  }
  free (f);
  free (fold);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
gmp_match (uint32_t p1, uint32_t p2, int64_t i, mpz_t m0,
	   uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
	   mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], qq, adz, tmp;
  unsigned int j;
  int cmp;
  double skew, logmu, E;


#ifdef DEBUG_POLYSELECT2L
  gmp_printf ("Found match: (%"PRIu32",%"PRId64") (%"PRIu32",%"PRId64") for "
	      "ad=%"PRIu64", q=%"PRIu64", rq=%Zd\n",
              p1, i, p2, i, ad, q, rq);
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
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (f == NULL || fold == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in match\n");
    exit (1);
  }
  for (j = 0; j <= d; j++) {
    mpz_init (f[j]);
    mpz_init (fold[j]);
  }
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
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_ui_p (m, d) == 0)
  {
    fprintf (stderr, "Error: m-a_{d-1}*l not divisible by d\n");
    exit (1);
  }
#endif
  mpz_divexact_ui (m, m, d);
  mpz_set_uint64 (adz, ad);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (m, adz) == 0)
  {
    fprintf (stderr, "Error: (m-a_{d-1}*l)/d not divisible by ad\n");
    exit (1);
  }
#endif
  mpz_divexact (m, m, adz);
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set (f[d], adz);
  mpz_pow_ui (t, m, d);
  mpz_mul (t, t, adz);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif
  mpz_divexact (t, t, l);
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (j = d - 2; j > 0; j--)
  {
#ifdef DEBUG_POLYSELECT2L
    if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
#endif
    mpz_divexact (t, t, l);
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

#ifdef DEBUG_POLYSELECT2L
  if (mpz_divisible_p (t, l) == 0)
  {
    fprintf (stderr, "Error: t not divisible by l\n");
    exit (1);
  }
#endif
  mpz_divexact (t, t, l);
  mpz_set (f[0], t);

  /* save unoptimized polynomial to fold */
  for (i = d + 1; i -- != 0; )
    mpz_set (fold[i], f[i]);
  mpz_set (gold[1], g[1]);
  mpz_set (gold[0], g[0]);

  /* old lognorm */
  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

  /* for degree 6 polynomials, find bottleneck coefficient */
  double skewtmp = 0.0, logmu0c4 = 0.0, logmu0c3 = 0.0;
  if (d == 6) {
    mpz_set (tmp, f[3]);
    mpz_set_ui (f[3], 0);
    skewtmp = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu0c3 = L2_lognorm (f, d, skewtmp, DEFAULT_L2_METHOD);
    mpz_set (f[3], tmp);
    mpz_set (tmp, f[4]);
    mpz_set_ui (f[4], 0);
    skewtmp = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu0c4 = L2_lognorm (f, d, skewtmp, DEFAULT_L2_METHOD);
    mpz_set (f[4], tmp);
  }

  double g0 = mpz_get_d (g[0]);
  g0 /= mpz_get_d (f[d-2]);
  g0 = (g0 > 0)? g0 : -g0;

#ifdef MAX_THREADS
  pthread_mutex_lock (&lock);
#endif
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
#ifdef MAX_THREADS
  pthread_mutex_unlock (&lock);
#endif

  /* if the polynomial has norm < "-maxnorm", we optimize it */
  if (logmu <= max_norm)
  {
    /* optimize size */
    optimize (f, d, g, 0, 1);

    if (!raw) {

/* root sieve */
#ifndef NEW_ROOTSIEVE
      unsigned long alim = 2000;
      long jmin, kmin;
#endif
      mpz_neg (m, g[0]);

      rootsieve_time -= seconds_thread ();

#ifdef NEW_ROOTSIEVE
      if (d > 3) {
        ropt_polyselect (f, d, m, g[1], N, 0); // verbose = 2 to see details.
        mpz_neg (g[0], m);
      }
      else {
        unsigned long alim = 2000;
        long jmin, kmin;
        rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
        mpz_neg (g[0], m);
        /* optimize again, but only translation */
        optimize_aux (f, d, g, 0, 0, CIRCULAR);
      }
#else
      rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
      mpz_neg (g[0], m);
      /* optimize again, but only translation */
      optimize_aux (f, d, g, 0, 0, CIRCULAR);
#endif

      rootsieve_time += seconds_thread ();

    } // raw and sopt only ?

    skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

    for (i = 10; i > 0 && logmu < best_logmu[i-1]; i--)
      best_logmu[i] = best_logmu[i-1];
    best_logmu[i] = logmu;

#ifdef MAX_THREADS
    pthread_mutex_lock (&lock);
#endif
    collisions_good ++;
    aver_opt_lognorm += logmu;
    var_opt_lognorm += logmu * logmu;
    if (logmu < min_opt_lognorm)
      min_opt_lognorm = logmu;
    if (logmu > max_opt_lognorm)
      max_opt_lognorm = logmu;
#ifdef MAX_THREADS
    pthread_mutex_unlock (&lock);
#endif

    /* MurphyE */
    mpz_set (curr_poly->rat->f[0], g[0]);
    mpz_set (curr_poly->rat->f[1], g[1]);
    for (j = d + 1; j -- != 0; )
      mpz_set (curr_poly->alg->f[j], f[j]);
    curr_poly->skew = skew;
    E =  MurphyE (curr_poly, BOUND_F, BOUND_G, AREA, MURPHY_K);

    mpz_neg (m, g[0]);

#ifdef MAX_THREADS
		  pthread_mutex_lock (&lock);
#endif
    if (E > best_E)
    {
      best_E = E;
      cado_poly_set (best_poly, curr_poly);
    }
    if (out != NULL) /* msieve output */
    {
      FILE *fp;
      fp = fopen (out, (found == 0) ? "w" : "a");
      if (fp == NULL)
      {
        fprintf (stderr, "Error, cannot open file %s\n", out);
        exit (1);
      }
      fprintf (fp, "0");
      for (j = d + 1; j -- != 0; )
        gmp_fprintf (fp, " %Zd", f[j]);
      mpz_neg (m, g[0]);
      gmp_fprintf (fp, " %Zd %Zd\n", g[1], m);
      fclose (fp);
    }
    found ++;
#ifdef MAX_THREADS
		  pthread_mutex_unlock (&lock);
#endif

    /* raw poly */
    if (raw) {
      printf ("# Raw polynomial:\n");
      gmp_printf ("n: %Zd\n", N);
      print_poly_info (fold, d, gold, 1);
      if (d == 6)
        gmp_printf ("# noc4/noc3: %.2f/%.2f (%.2f)\n",
                    logmu0c4, logmu0c3, logmu0c4/logmu0c3);
      gmp_printf ("# Optimized polynomial:\n");
    }

    /* print optimized (maybe size- or size-root- optimized) polynomial */
    if (verbose >= 0)
      {
        gmp_printf ("n: %Zd\n", N);
        print_poly_info (f, d, g, 0);
        printf ("# Murphy's E(Bf=%.0f,Bg=%.0f,area=%.2e)=%1.2e (best so far %1.2e)\n",
                BOUND_F, BOUND_G, AREA, E, best_E);
        printf ("\n");
        fflush (stdout);
      }
  }
  else {
    if (verbose >= 0) {
      if (d == 6)
        gmp_printf ("# Skip polynomial: %.2f, ad: %"PRIu64", l: %Zd, m: %Zd, noc3: %.2f, noc4: %.2f\n",
                    logmu, ad, l, m, logmu0c3, logmu0c4);
      else
        gmp_printf ("# Skip polynomial: %.2f, ad: %"PRIu64", l: %Zd, m: %Zd\n",
                    logmu, ad, l, m);
    }
  }

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
  for (j = 0; j <= d; j++) {
    mpz_clear (f[j]);
    mpz_clear (fold[j]);
  }
  free (f);
  free (fold);
}


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long i, j, nprimes, p, nrp, c = 0, tot_roots = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  double pc1;
  mpz_t *f, tmp;
  int found = 0;
  shash_t H;
  int st = 0;

  /* init f for roots computation */
  mpz_init_set_ui (tmp, 0);
  f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
  if (f == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }
  for (i = 0; i <= header->d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[header->d], 1);

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
                          header->N, 1, tmp);
              for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
                hash_add (H, p, -u, header->m0, header->ad,
                          header->d, header->N, 1, tmp);
            }
        }
#ifdef DEBUG_POLYSELECT2L
      fprintf (stderr, "# collision_on_p took %dms\n", milliseconds () - st);
      fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

#ifdef DEBUG_HASH_TABLE
      fprintf (stderr, "# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
      fprintf (stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif
      hash_clear (H);
    }

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  mpz_clear (tmp);

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
  MAYBE_UNUSED uint64_t cpt = 0, cpt1 = 0, cpt2 = 0, cpt3 = 0;

#ifdef DEBUG_POLYSELECT2L
  int st = milliseconds();
#endif
#if SHASH_NBUCKETS == 256
#define CURRENT(V) (H->current + (uint8_t) (V))
#else
#define CURRENT(V) (H->current + ((V) & (SHASH_NBUCKETS - 1)))
#endif
  /*
  uint64_t t1, t2;
  static uint64_t sum1 = 0, sum2 = 0;
  */

  shash_init (H, 4 * lenPrimes);

  /*
  t1 = cputicks();
  */

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

#if 0 /* Old algo */

  for (;;) {
    do {
      vpnr = *pnr++;
      pprimes++;
    } while (!vpnr);
    if (UNLIKELY(vpnr == 0xff)) break;
    ppl = *pprimes;
    __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
    __builtin_prefetch(((void *) pprimes) + 0x100, 0, 3);
    __builtin_prefetch(((void *) pc) + 0x280, 0, 3);
    ppl *= ppl;
    epc = pc + vpnr;
    for (;;) {
      /* v2 = nv [- ppl] insertions == nv + [ppl] insertions -/+ 1;
	 max inserts == 4+4, min == 1+1. */
      v1 = nv;       cur1 = CURRENT(v1); ccur1 = *cur1;
      v2 = v1 - ppl; cur2 = CURRENT(v2); ccur2 = *cur2;
      nv = *++pc;
      *ccur1++ = v1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;
      *ccur2++ = v2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;
      v1 += ppl; v2 -= ppl;
      if (v1 > umax)     goto l1v2;
      if (v2 < neg_umax) goto l1v1;
      INSERT_2I(v1,v2); v1 += ppl; v2 -= ppl;
      if (v1 > umax)     goto l1v2;
      if (v2 < neg_umax) goto l1v1;
      INSERT_2I(v1,v2); v1 += ppl; v2 -= ppl;
      if (LIKELY(v1 > umax))
      l1v2: if (LIKELY(v2 < neg_umax))
	  if (LIKELY(pc != epc)) continue; else break;
	else {
	  INSERT_I(v2); if (LIKELY(pc != epc)) continue; else break;
	}
      else if (LIKELY(v2 < neg_umax)) {
      l1v1: INSERT_I(v1); if (LIKELY(pc != epc)) continue; else break;
      }
      else {
	INSERT_2I(v1,v2); if (LIKELY(pc != epc)) continue; else break;
      }
    }
  }
    
#else /* New Algo */

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
#endif
#undef INSERT_2I
#undef INSERT_I

  /* fprintf (stderr, "%lu %lu %lu %lu\n", cpt, cpt1, cpt2, cpt3); */
  for (i = 0; i < SHASH_NBUCKETS; i++) assert (H->current[i] <= H->base[i+1]);

  /*
  t2 = cputicks();
  sum1 += t2 - t1;
  */

  found = shash_find_collision (H);
  shash_clear (H);

  /*
  sum2 += cputicks() - t2;
  */

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
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n", milliseconds () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

#ifdef DEBUG_HASH_TABLE
  fprintf (stderr, "# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
  fprintf (stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif

  pc2 = expected_collisions (Primes[lenPrimes - 1]);
  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);

  /*
  fprintf (stderr, "%lu %lu\n", sum1 / 3000000, sum2 / 3000000);
  */
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
        tinv_qq[k][c] = modredcul_intget_ul (res_rp, modpp);
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
    invqq[nprimes] = modredcul_intget_ul (tmp, modpp);

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
  unsigned long i, j, nprimes, p, nrp, c = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  double pc1;
  mpz_t *f, tmp;

  /* init f for roots computation */
  mpz_init_set_ui (tmp, 0);
  f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
  if (f == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
    exit (1);
  }
  for (i = 0; i <= header->d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[header->d], 1);

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
                      header->d, header->N, 1, tmp);
      for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
        gmp_hash_add (H, p, -u, header->m0, header->ad,
                      header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", milliseconds () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

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
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n",
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
             "collision_on_sq(). ad=%"PRIu64".\n", N, K, header->ad);
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

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
  return NULL;
}


static void
usage (const char *argv, const char * missing)
{
  fprintf (stderr, "Usage: %s [options] P\n", argv);
  fprintf (stderr, "Required parameters and options:\n");
  fprintf (stderr, "P            --- degree-1 coefficient of g(x) has\n");
  fprintf (stderr, "                 two prime factors in [P,2P]\n");
  fprintf (stderr, "-v           --- verbose mode\n");
  fprintf (stderr, "-q           --- quiet mode\n");
  fprintf (stderr, "-r           --- size-optimize polynomial only (skip ropt)\n");
  fprintf (stderr, "-t nnn       --- use n threads (default 1)\n");
  fprintf (stderr, "-admin nnn   --- start from ad=nnn (default 0)\n");
  fprintf (stderr, "-admax nnn   --- stop at ad=nnn\n");
  fprintf (stderr, "-incr nnn    --- forced factor of ad (default 60)\n");
  fprintf (stderr, "-N nnn       --- input number\n");
  fprintf (stderr, "-degree nnn  --- wanted polynomial degree\n");
  fprintf (stderr, "-nq nnn      --- maximum number of special-q's considered\n");
  fprintf (stderr, "                 for each ad (default %d)\n", INT_MAX);
  fprintf (stderr, "-save xxx    --- save state in file xxx\n");
  fprintf (stderr, "-resume xxx  --- resume state from file xxx\n");
  fprintf (stderr, "-maxnorm xxx --- only optimize polynomials with norm <= xxx\n");
  fprintf (stderr, "-maxtime xxx --- stop the search after xxx seconds\n");
  fprintf (stderr, "-out xxx     --- for msieve-format output\n");
  fprintf (stderr, "-s xxx       --- time intervals (seconds) for printing\n");
  fprintf (stderr, "                 out statistics (default %d)\n", TARGET_TIME / 1000);
  if (missing) {
      fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
              missing);
  }
  exit (1);
}


int
main (int argc, char *argv[])
{
  argc0 = argc;
  argv0 = argv;
  const char *save = NULL, *resume = NULL;
  st0 = seconds ();
  double maxtime = DBL_MAX;
  mpz_t N;
  unsigned int d = 0;
  unsigned long P, admin = 0, admax = ULONG_MAX;
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

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-r", &raw);
  param_list_configure_switch (pl, "-q", &quiet);
  param_list_configure_alias(pl, "degree", "-d");
  param_list_configure_alias(pl, "incr", "-i");
  param_list_configure_alias(pl, "N", "-n");

  argv++, argc--;
  for ( ; argc-1; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0[0],NULL);
  }

  /* parse and check N in the first place */
  int have_n = param_list_parse_mpz(pl, "N", N);

  if (!have_n) {
    fprintf(stderr, "Reading n from stdin\n");
    param_list_read_stream(pl, stdin);
    have_n = param_list_parse_mpz(pl, "n", N);
  }

  if (!have_n) {
      fprintf(stderr, "No n defined ; sorry.\n");
      exit(1);
  }

  if (mpz_cmp_ui (N, 0) <= 0) usage(argv0[0], "N");

  param_list_parse_int (pl, "t", &nthreads);
  param_list_parse_int (pl, "nq", &nq);
  param_list_parse_int (pl, "s", &target_time);
  incr_target_time = target_time;
  param_list_parse_uint (pl, "degree", &d);
  param_list_parse_ulong (pl, "admin", &admin);
  param_list_parse_ulong (pl, "admax", &admax);
  param_list_parse_ulong (pl, "incr", &incr);
  param_list_parse_double (pl, "maxnorm", &max_norm);
  param_list_parse_double (pl, "maxtime", &maxtime);
  save = param_list_lookup_string (pl, "save");
  resume = param_list_lookup_string (pl, "resume");
  out = param_list_lookup_string (pl, "out");

  /* print command line */
  param_list_print_command_line (stdout, pl);

  /* check degree */
  if (d <= 0) usage(argv0[0], "degree");

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
  best_poly->alg->degree = d;
  best_poly->rat->degree = 1;
  curr_poly->alg->degree = d;
  curr_poly->rat->degree = 1;

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

  for (i = 0; i <= 10; i++)
    best_logmu[i] = 999.9;

  /* init primes */
  double Pd;
  Pd = strtod (argv[argc-1], NULL);
  if (Pd > (double) UINT_MAX) {
    fprintf (stderr, "Error, too large value of P\n");
    exit (1);
  }
  P = (unsigned long) Pd;
  if (P <= (unsigned long) SPECIAL_Q[LEN_SPECIAL_Q - 2]) {
    fprintf (stderr, "Error, too small value of P\n");
    exit (1);
  }

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
        gmp_printf ("%d ad=%lu\n", tries, admin);
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

          printf ("# Stat: raw lognorm (min/av/max/std): %1.2f/%1.2f/%1.2f/%1.2f\n",
                  min_raw_lognorm, rawmean, max_raw_lognorm,
                  sqrt (var_raw_lognorm / collisions - rawmean * rawmean));
          if (collisions_good > 0)
            printf ("# Stat: optimized lognorm (min/av/max/std): %1.2f/%1.2f/%1.2f/%1.2f\n",
                    min_opt_lognorm, mean, max_opt_lognorm,
                    sqrt (var_opt_lognorm / collisions_good - mean * mean));
          printf ("# Stat: av. g0/adm2 ratio: %.3e\n",
                  total_adminus2 / (double) collisions);
          if (d == 6)
            printf ("# Stat: av. logmu noc4/noc3 ratio: %.3f\n",
                    aver_lognorm_ratio / (double) collisions);
        }
    }

  printf ("# Stat: tried %d ad-value(s), found %d polynomial(s), %d below maxnorm\n",
          tries, tot_found, found);

  /* print best 10 values of logmu */
  if (collisions_good > 0)
    {
      printf ("# Stat: best logmu:");
      for (i = 0; i < 10; i++)
        printf (" %1.2f", best_logmu[i]);
      printf ("\n");
    }

  /* print total time (format for cpu_time.sh) */
  printf ("# Stat: total phase took %.2fs\n", seconds () - st0);
#ifndef HAVE_RUSAGE_THREAD /* rootsieve_time is correct only if RUSAGE_THREAD
                              works or in mono-thread mode */
  if (nthreads == 1)
#endif
    printf ("# Stat: rootsieve took %.2fs\n", rootsieve_time);

  if (best_E == 0.0)
    printf ("No polynomial found, please increase the ad range or decrease P\n");
  else
    print_cadopoly_extra (stdout, best_poly, argc0, argv0, st0, 1 /* raw */);

  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  mpz_clear (N);
  clearPrimes (&Primes);
  cado_poly_clear (best_poly);
  cado_poly_clear (curr_poly);
  param_list_clear (pl);

  return 0;
}
