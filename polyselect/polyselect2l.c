/*
  polyselect2l.c is a variant of Paul Zimmermann's polyselect2.c

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-nq xxx" denotes the number of special-q's trials for each ad;

  "-lq xxx" denotes the number of small factors (<= 251) in the special-q
  (see SPECIAL_Q[] in polyselect2l_str.c);

  "-maxnorm xxx" only optimize raw polynomials with size <= xxx.
  If the raw polynomial is not good enough, we will still stream
  it to STDERR for further reference.

  Please report bugs to shi.bai AT anu.edu.au.
*/

#include "polyselect2l.h"

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */
#define NEW_ROOTSIEVE
#define MAX_THREADS 16
#define INIT_FACTOR 4
// #define DEBUG_POLYSELECT2L

#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#endif

/* Two modes: batch P or batch SQ. The latter
   seems faster, but need more memory. Use the latter by default. */
//#define BATCH_P
#ifdef BATCH_P
#define BATCH_SIZE 8 /* batch 8 p */
#else
#define BATCH_SIZE 10 /* batch 10 sq */
#endif

/* consider only two roots or more */
#define CONSIDER_ONLY_TWO_ROOTS
#ifndef CONSIDER_ONLY_TWO_ROOTS
#define NUMBER_CONSIDERED_ROOTS 16
#endif

#define LQ_DEFAULT 1 /* default number of factors in special-q part */

/* Read-Only */
uint32_t *Primes = NULL;
unsigned long lenPrimes = 1; // length of Primes[]
int nq = INT_MAX;
int lq = LQ_DEFAULT;
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
static int verbose = 0;
static unsigned long incr = DEFAULT_INCR;
const char *out = NULL; /* output file for msieve input (msieve.dat.m) */
cado_poly best_poly, curr_poly;
double best_E = 0.0; /* Murphy's E (the larger the better) */
int seed = 0; /* seed */

int argc0;
char ** argv0;
double st0;

/* read-write global variables */
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */
int tot_found = 0; /* total number of polynomials */
int found = 0; /* number of polynomials below maxnorm */
double potential_collisions = 0.0, aver_opt_lognorm = 0.0,
  aver_raw_lognorm = 0.0, aver_lognorm_ratio = 0.0;
double min_raw_lognorm = DBL_MAX, max_raw_lognorm = 0.0;
double min_opt_lognorm = DBL_MAX, max_opt_lognorm = 0.0;
unsigned long collisions = 0;
unsigned long collisions_good = 0;
double total_adminus2 = 0.0;
double best_logmu[11];
double rootsieve_time = 0.0;
int raw = 0;

/* -- functions starts here -- */

/* crt, set r and qqz */
void
crt_sq ( mpz_t qqz,
         mpz_t r,
         unsigned long *q,
         unsigned long *rq )
{
  mpz_t prod, pprod, mod, inv, sum;
  int i;
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
static void
check_parameters (mpz_t m0, unsigned long d)
{
  double maxq = 1.0, maxP;
  int k = lq;
  
  while (k > 0)
    maxq *= (double) SPECIAL_Q[LEN_SPECIAL_Q - 1 - (k--)];

  maxP = (double) Primes[lenPrimes - 1];
  if (2.0 * pow (maxP, 4.0) * maxq >= (double) d * mpz_get_d (m0))
    {
      fprintf (stderr, "Error, too large value of -lq parameter\n");
      exit (1);
    }
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

  nroots = numberOfRealRoots (f, d, 0, 0);
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


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       uint64_t ad, unsigned long d, mpz_t N, uint64_t q,
       mpz_t rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], qq, adz;
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
  mpz_mul_ui (l, l, q);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set_si (mtilde, i);
  mpz_set_ui (qq, q);
  mpz_mul_ui (qq, qq, q);
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


// Do we batch p or batch q?
#ifdef BATCH_P


/* find collisions between "P" primes */
static inline void
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long i, j, nprimes, p, nrp;
  uint64_t *rp;
  int64_t ppl = 0;
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
#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H, (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H, (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				 * (double) lenPrimes));
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    ppl = (int64_t) p;
    ppl *= (int64_t) p;
    printf ("p: %lu, ppl %"PRId64"\n", p, ppl);

    /* add faked roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    mpz_mod_ui (f[0], header->Ntilde, p);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrp = poly_roots_uint64 (rp, f, header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);
    for (j = 0; j < nrp; j++) {
      /* only consider r[j] and r[j] - pp */
      hash_add (H, p, (int64_t) rp[j], header->m0, header->ad,
		header->d, header->N, 1, tmp);
      hash_add (H, p, (int64_t) (rp[j] - ppl), header->m0,
		header->ad, header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", cputime () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %"PRIu64"\n", H->size,
	   header->ad);
#endif

  /* if the hash table contains n entries, each one smaller than (2P)^2,
     the number of potential collisions is about 0.5*n^2/(2P)^2 */
  pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);
  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

  pthread_mutex_lock (&lock);
  potential_collisions += pc1;
  pthread_mutex_unlock (&lock);
}


/* collision on each special-q, called by collision_on_batch_p() */
static inline void
collision_on_each_sq ( header_t header,
                       proots_t R,
                       unsigned long q,
                       mpz_t rqqz,
                       unsigned long *inv_qq ) // inv_qq is of length #Primes
{
  unsigned int nr, j;
  unsigned long nprimes, p, rp, pp;
  long ppl, u;
  mpz_t rppz;
  double pc2;

#ifndef CONSIDER_ONLY_TWO_ROOTS
  int i;
  long h;
#endif  

  mpz_init (rppz);

#ifdef DEBUG_POLYSELECT2L
  mpz_t tmp_debug, tmp_debug2, qqz;
  mpz_init_set (tmp_debug, header->Ntilde);
  mpz_init_set_ui (qqz, q);
  mpz_mul_ui (qqz, qqz, q);
  mpz_mod (tmp_debug, tmp_debug, qqz);
  mpz_init_set (tmp_debug2, header->m0);
  mpz_add (tmp_debug2, tmp_debug2, rqqz);
  mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
  mpz_mod (tmp_debug2, tmp_debug2, qqz);
  if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
    fprintf (stderr, "Error: crt root is wrong in collision_on_each_sq\n");
    exit (1);
  }
  mpz_clear (tmp_debug);
  mpz_clear (tmp_debug2);
  mpz_clear (qqz);
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  hash_t H;

#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS * (double) lenPrimes));
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;
    /* set p, p^2, ppl */
    pp = p * p;
    ppl = (long) pp;
    nr = R->nr[nprimes];

    modulusul_t modul_pp;
    residueul_t res_rp, res_tmp;
    modul_initmod_ul (modul_pp, pp);
    modul_init (res_rp, modul_pp);
    modul_init (res_tmp, modul_pp);
    modul_set_ul (res_tmp, inv_qq[nprimes], modul_pp);

    for (j = 0; j < nr; j++) {
      rp = R->roots[nprimes][j];

      modul_set_ul (res_rp, rp, modul_pp);
      modul_sub_ul (res_rp, res_rp, mpz_fdiv_ui (rqqz, pp), modul_pp);
      modul_mul (res_rp, res_rp, res_tmp, modul_pp);
      u = (long) modul_get_ul (res_rp, modul_pp);

#ifdef DEBUG_POLYSELECT2L // this is perhaps wrong.
      mpz_t Ntmp, m0tmp, tmp, qqz;
      mpz_init_set_ui (qqz, q);
      mpz_mul_ui (qqz, qqz, q);
      mpz_init_set (Ntmp, header->Ntilde);
      mpz_mod_ui (Ntmp, Ntmp, ppl); // tildeN (mod p^2)
      mpz_init_set (m0tmp, header->m0);
      mpz_add (m0tmp, m0tmp, rqqz); // m0 + rq
      mpz_init_set_ui (tmp, u);
      mpz_mul (tmp, tmp, qqz); // i*q^2
      mpz_add (m0tmp, m0tmp, tmp); // m0 + rq + i*q^2
      mpz_pow_ui (m0tmp, m0tmp, header->d);
      mpz_mod_ui (m0tmp, m0tmp, ppl);
      if (mpz_cmp (m0tmp, Ntmp) != 0) {
        fprintf (stderr, "Error: i computation is wrong in "
                 "collision_on_each_sq\n");
        gmp_printf ("Details: (p=%lu, i(mod p^2)=%ld) for "
                    "ad=%"PRIu64", q=%lu, rq=%Zd, rp=%lu, invqq=%lu\n",
                    p, u, header->ad, q, rqqz, rp, inv_qq[nprimes]);
        gmp_printf ("m0=%Zd\n", header->m0);
        gmp_printf ("Ntilde=%Zd\n", header->Ntilde);
        exit (1);
      }
      mpz_clear (Ntmp);
      mpz_clear (m0tmp);
      mpz_clear (tmp);
      mpz_clear (qqz);
#endif

#ifdef CONSIDER_ONLY_TWO_ROOTS
      hash_add (H, p, u, header->m0, header->ad, header->d, header->N, q, rqqz);
      hash_add (H, p, u - ppl, header->m0, header->ad, header->d, header->N, q, rqqz);
#else
      h = u - ppl;
      for (i = 0; i < NUMBER_CONSIDERED_ROOTS; i ++) {
        hash_add (H, p, u, header->m0, header->ad, header->d, header->N, q, rqqz);
        hash_add (H, p, h, header->m0, header->ad, header->d, header->N, q, rqqz);
        u += ppl;
        h -= ppl;
      }

#endif

    }  // next rp

  modul_clear (res_rp, modul_pp);
  modul_clear (res_tmp, modul_pp);
  modul_clearmod (modul_pp);

  } // next p

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n", cputime () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  pc2 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);

  mpz_clear (rppz);
  hash_clear (H);

  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* Batch P mode */
static inline void
collision_on_batch_p ( header_t header,
                       proots_t R,
                       uint64_t q,
                       mpz_t qqz,
                       mpz_t rqqz )
{
  unsigned int size = BATCH_SIZE;
  if (size == 0)
    return;

  unsigned int i, j;
  unsigned long nprimes, pul, ppul[size], index[size];
  unsigned long *invqq = malloc (lenPrimes * sizeof (unsigned long));
  if (!invqq) {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

  mpz_t tmp, p_prod;

  mpz_init (tmp);
  mpz_init (p_prod);

#ifdef DEBUG_POLYSELECT2L  // check roots
  mpz_t tmp_debug, tmp_debug2;
  mpz_init_set (tmp_debug, header->Ntilde);
  mpz_mod (tmp_debug, tmp_debug, qqz);
  mpz_init_set (tmp_debug2, header->m0);
  mpz_add (tmp_debug2, tmp_debug2, rqqz);
  mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
  mpz_mod (tmp_debug2, tmp_debug2, qqz);
  if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
    fprintf (stderr, "Error: crt root is wrong in %s\n", __FUNCTION__);
    exit (1);
  }
  mpz_clear (tmp_debug);
  mpz_clear (tmp_debug2);
#endif

  //int st = cputime();

  /* Step 1: batch inversion */
  i = 0;
  mpz_set_ui (p_prod, 1);

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    pul = Primes[nprimes];
    if ((header->d * header->ad) % pul == 0)
      continue;

    // batch or reset
    if ( (i == size - 1) || (nprimes == (lenPrimes-1)) ) {

      ppul[i] = pul*pul;
      mpz_mul_ui (p_prod, p_prod, ppul[i]);
      index[i] = nprimes;

      // the inversion
      mpz_invert (tmp, qqz, p_prod);

      // individual inversions
      for (j = 0; j <= i; j ++) {
        pul = Primes[index[j]];
        invqq[index[j]] = mpz_fdiv_ui (tmp, ppul[j]);
      }

#ifdef DEBUG_POLYSELECT2L // check inversions for p in this batch
      for (j = 0; j < i; j ++) {
        mpz_t tmp_debug, tmp_debug2;
        mpz_init_set_ui (tmp_debug, invqq[index[j]]);
        mpz_init (tmp_debug2);
        mpz_mul (tmp_debug2, tmp_debug, qqz);
        unsigned long re = mpz_fdiv_ui (tmp_debug2, ppul[j]);
        if (re != 1) {
          fprintf (stderr, "Error, batch inversion is wrong in %s, iter: %d\n", __FUNCTION__, i);
          exit (1);
        }
        mpz_clear (tmp_debug);
        mpz_clear (tmp_debug2);
      }
#endif
      i = 0;
      mpz_set_ui (p_prod, 1);
    }
    else {
      ppul[i] = pul*pul;
      mpz_mul_ui (p_prod, p_prod, ppul[i]);
      index[i++] = nprimes;
    }

  } // next prime

  //fprintf (stderr, "# one batch P inversion took %dms\n", cputime () - st);

#ifndef DEBUG_POLYSELECT2L // check inversions for all p
  unsigned long pp;
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    mpz_t tmp_debug, tmp_debug2;
    pul = Primes[nprimes];
    if ((header->d * header->ad) % pul == 0)
      continue;
    pp = pul * pul;
    mpz_init_set_ui (tmp_debug, invqq[nprimes]);
    mpz_init (tmp_debug2);
    mpz_mul (tmp_debug2, tmp_debug, qqz);
    unsigned long re = mpz_fdiv_ui (tmp_debug2, pp);
    if (re != 1) {
      fprintf (stderr, "Error, batch inversion is wrong in %s, iter: %lu\n", __FUNCTION__, nprimes);
      exit (1);
    }
    mpz_clear (tmp_debug);
    mpz_clear (tmp_debug2);
  }
#endif

  /* Step 2: find collisions on sq. */
  //st = cputime();
  collision_on_each_sq ( header,
                         R,
                         q,
                         rqqz,
                         invqq );
  //printf ("# collision_on_each_sq took %dms\n", cputime () - st);
  mpz_clear (p_prod);
  mpz_clear (tmp);
  free (invqq);
}


/* collision on special-q, call collision_on_batch_p */
static void
collision_on_sq ( header_t header,
                  proots_t R )
{
  mpz_t qqz, rqqz;
  qroots_t SQ_R;

  mpz_init (qqz);
  mpz_init (rqqz);

  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  unsigned long k = lq, n = SQ_R->size;
  uint64_t q;
  // roots_print (SQ_R);

  /* less than lq special primes having roots for this ad */
  if (n == 0 || n < k)
    return;

  unsigned long idx_q[n], tot, i;
  /* if tot (n, k) < wanted, use tot as wanted */
  tot =  binom (n, k);

  if (tot > (unsigned long) nq)
    tot = (unsigned long) nq;
  /*
   fprintf (stderr, "n=%"PRIu64", k=%"PRIu64", (n,k)=%"PRIu64",
   nq:%d\n", n, k, tot, nq);
  */

  /* enumerate each combination */
  first_comb (n, idx_q);
  q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
  collision_on_batch_p (header, R, q, qqz, rqqz);

  //fprintf (stderr, "# - q ");
  for (i = 1; i < tot; i ++) {
    next_comb (n, k, idx_q); // print_comb (k, idx_q);
    q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
    collision_on_batch_p (header, R, q, qqz, rqqz);
  }

  mpz_clear (qqz);
  mpz_clear (rqqz);
  qroots_clear (SQ_R);
}


// separator between modredc_ul and gmp


/* find collisions between "P" primes */
static inline void
gmp_collision_on_p ( header_t header,
		     proots_t R )
{
  unsigned long i, j, nprimes, p, nrp;
  uint64_t *rp;
  int64_t ppl = 0;
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
#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H, (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H, (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				 * (double) lenPrimes));
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    ppl = (int64_t) p;
    ppl *= (int64_t) p;

    /* add faked roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    mpz_mod_ui (f[0], header->Ntilde, p);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrp = poly_roots_uint64 (rp, f, header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);
    for (j = 0; j < nrp; j++) {
      /* only consider r[j] and r[j] - pp */
      gmp_hash_add (H, p, (int64_t) rp[j], header->m0, header->ad,
		    header->d, header->N, 1, tmp);
      gmp_hash_add (H, p, (int64_t)((int64_t) (rp[j])-ppl),  header->m0,
		    header->ad, header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", cputime () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %"PRIu64"\n", H->size,
	   header->ad);
#endif

  /* if the hash table contains n entries, each one smaller than (2P)^2,
     the number of potential collisions is about 0.5*n^2/(2P)^2 */
  pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);
  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

  pthread_mutex_lock (&lock);
  potential_collisions += pc1;
  pthread_mutex_unlock (&lock);
}


/* collision on each special-q, called by collision_on_batch_p() */
static inline void
gmp_collision_on_each_sq ( header_t header,
			   proots_t R,
			   uint64_t q,
			   mpz_t rqqz,
			   uint64_t *inv_qq )
{
  unsigned int nr, j;
  unsigned long nprimes, p;
  mpz_t rppz, ppmp, tmp;
  double pc2;
  uint64_t pp, rp;
  int64_t ppl, u;

#ifndef CONSIDER_ONLY_TWO_ROOTS
  int i;
  int64_t h;
#endif  

  mpz_init (rppz);
  mpz_init (ppmp);
  mpz_init (tmp);
  /*
  gmp_printf ("ad: %"PRIu64", rq: %Zd, q: %"PRIu64"\n",
	      header->ad, rqqz, q);
  */
#ifdef DEBUG_POLYSELECT2L
  mpz_t tmp_debug, tmp_debug2, qqz;
  mpz_init_set (tmp_debug, header->Ntilde);
  mpz_init_set_ui (qqz, q);
  mpz_mul_ui (qqz, qqz, q);
  mpz_mod (tmp_debug, tmp_debug, qqz);
  mpz_init_set (tmp_debug2, header->m0);
  mpz_add (tmp_debug2, tmp_debug2, rqqz);
  mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
  mpz_mod (tmp_debug2, tmp_debug2, qqz);
  if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
    fprintf (stderr, "Error: crt root is wrong in collision_on_each_sq\n");
    exit (1);
  }
  mpz_clear (tmp_debug);
  mpz_clear (tmp_debug2);
  mpz_clear (qqz);
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  hash_t H;

#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				  * (double) lenPrimes));
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;
    /* set p, p^2, ppl */
    pp = (uint64_t) p;
    pp *= (uint64_t) p;
    ppl = (int64_t) pp;
    nr = R->nr[nprimes];

    for (j = 0; j < nr; j++) {
      rp = R->roots[nprimes][j];
      mpz_set_uint64 (rppz, rp);
      mpz_sub (rppz, rppz, rqqz);
      mpz_mul_ui (rppz, rppz, inv_qq[nprimes]);
      mpz_set_uint64 (ppmp, pp);
      mpz_fdiv_r (tmp, rppz, ppmp);
      u = (int64_t) mpz_get_uint64 (tmp);

#ifdef DEBUG_POLYSELECT2L
      mpz_t Ntmp, m0tmp, qqz, tmp2, tmp3;
      mpz_init (Ntmp);
      mpz_init (m0tmp);
      mpz_init (qqz);
      mpz_init (tmp2);
      mpz_init (tmp3);
      /* qqz = q^2 */
      mpz_set_ui (qqz, q);
      mpz_mul_ui (qqz, qqz, q);
      /* tmp3 = p^2 */
      mpz_set_uint64 (tmp3, pp);
      /* tildeN (mod p^2) */
      mpz_set (Ntmp, header->Ntilde);
      mpz_fdiv_r (Ntmp, Ntmp, tmp3);
      /* m0 + rq */
      mpz_set (m0tmp, header->m0);
      mpz_add (m0tmp, m0tmp, rqqz);
      /* i*q^2 */
      if (u >= 0)
	mpz_set_uint64 (tmp2, (uint64_t) u);
      else
	mpz_set_uint64 (tmp2, (uint64_t) (-u));
      mpz_mul (tmp2, tmp2, qqz);
      mpz_add (m0tmp, m0tmp, tmp2); // m0 + rq + i*q^2
      mpz_pow_ui (m0tmp, m0tmp, header->d);
      mpz_mod (m0tmp, m0tmp, tmp3);
      if (mpz_cmp (m0tmp, Ntmp) != 0) {
        fprintf (stderr, "Error: i computation is wrong in "
		 "collision_on_each_sq\n");
        gmp_printf ("Details: (p=%lu, i(mod p^2)=%"PRId64") "
		    "for ad=%"PRIu64", q=%lu, rq=%Zd, rp=%"PRIu64", "
		    "invqq=%"PRIu64"\n",
                    p, u, header->ad, q, rqqz, rp, inv_qq[nprimes]);
        gmp_printf ("m0=%Zd\n", header->m0);
        gmp_printf ("Ntilde=%Zd\n", header->Ntilde);
        exit (1);
      }
      mpz_clear (Ntmp);
      mpz_clear (m0tmp);
      mpz_clear (tmp2);
      mpz_clear (tmp3);
      mpz_clear (qqz);
#endif

#ifdef CONSIDER_ONLY_TWO_ROOTS
      gmp_hash_add (H, p, u, header->m0, header->ad, header->d, header->N,
		q, rqqz);
      gmp_hash_add (H, p, u - ppl, header->m0, header->ad, header->d,
		header->N, q, rqqz);
#else
      h = u - ppl;
      for (i = 0; i < NUMBER_CONSIDERED_ROOTS; i ++) {
        gmp_hash_add (H, p, u, header->m0, header->ad, header->d, header->N,
		  q, rqqz);
        gmp_hash_add (H, p, h, header->m0, header->ad, header->d, header->N,
		  q, rqqz);
        u += ppl;
        h -= ppl;
      }

#endif

    }  // next rp
  } // next p


#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n", cputime () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  pc2 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);

  /* clear */
  mpz_clear (rppz);
  mpz_clear (ppmp);
  mpz_clear (tmp);
  hash_clear (H);

  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* Batch P mode */
static inline void
gmp_collision_on_batch_p ( header_t header,
			   proots_t R,
			   uint64_t q,
			   mpz_t qqz,
			   mpz_t rqqz )
{
  unsigned int size = BATCH_SIZE;
  if (size == 0)
    return;
  unsigned int i, j;
  unsigned long nprimes, pul, index[size];
  uint64_t ppul[size], *invqq;
  mpz_t tmp, tmp1, p_prod, *ppmp;

  mpz_init (tmp);
  mpz_init (tmp1);
  mpz_init (p_prod);
  invqq = (uint64_t*) malloc (lenPrimes * sizeof (uint64_t));
  ppmp = (mpz_t*) malloc (size * sizeof(mpz_t));
  if (!invqq || !ppmp) {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

  for (j = 0; j < size; j ++)
    mpz_init (ppmp[j]);

  //int st = cputime();

  /* Step 1: batch inversion */
  i = 0;
  mpz_set_ui (p_prod, 1);

  //int st = cputime();

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    pul = Primes[nprimes];
    if ((header->d * header->ad) % pul == 0)
      continue;

    // batch or reset
    if ( (i == size - 1) || (nprimes == (lenPrimes-1)) ) {

      ppul[i] = (uint64_t) pul;
      ppul[i] *= (uint64_t) pul;
      mpz_set_uint64 (ppmp[i], ppul[i]);
      mpz_mul (p_prod, p_prod, ppmp[i]);
      index[i] = nprimes;

      // the inversion
      mpz_invert (tmp, qqz, p_prod);

      // individual inversions
      for (j = 0; j <= i; j ++) {
	mpz_fdiv_r (tmp1, tmp, ppmp[j]);
	invqq[index[j]] = mpz_get_uint64 (tmp1);
      }

#ifdef DEBUG_POLYSELECT2L // check inversions for p in this batch
      for (j = 0; j < i; j ++) {
        mpz_t tmp_debug, tmp_debug2;
        mpz_init (tmp_debug);
	mpz_set_uint64 (tmp_debug, invqq[index[j]]);
        mpz_init (tmp_debug2);
        mpz_mul (tmp_debug2, tmp_debug, qqz);
	mpz_set_uint64 (tmp_debug, ppul[j]);
	mpz_fdiv_r (tmp_debug, tmp_debug2, tmp_debug);
        if (mpz_cmp_ui (tmp_debug, 1) != 0) {
          fprintf (stderr, "Error, batch inversion is wrong "
		   "in %s, iter: %d\n", __FUNCTION__, i);
          exit (1);
        }
        mpz_clear (tmp_debug);
        mpz_clear (tmp_debug2);
      }
#endif
      i = 0;
      mpz_set_ui (p_prod, 1);
    }
    else {
      ppul[i] = (uint64_t) pul;
      ppul[i] *= (uint64_t) pul;
      mpz_set_uint64 (ppmp[i], ppul[i]);
      mpz_mul (p_prod, p_prod, ppmp[i]);
      index[i++] = nprimes;
    }

  } // next prime

  //fprintf (stderr, "# one batch P inversion took %dms\n", cputime () - st);

  /* Step 2: find collisions on sq. */

  //st = cputime();
  gmp_collision_on_each_sq (header, R, q, rqqz, invqq);
  //printf ("# collision_on_each_sq took %dms\n", cputime () - st);

  mpz_clear (p_prod);
  mpz_clear (tmp);
  mpz_clear (tmp1);
  for (j = 0; j < size; j ++)
    mpz_clear (ppmp[j]);
  free (ppmp);
  free (invqq);
}


/* collision on special-q, call gmp_collision_on_batch_p */
static void
gmp_collision_on_sq ( header_t header,
		      proots_t R )
{
  mpz_t qqz, rqqz;
  qroots_t SQ_R;

  mpz_init (qqz);
  mpz_init (rqqz);

  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  unsigned long k = lq, n = SQ_R->size;
  uint64_t q;
  // roots_print (SQ_R);

  /* less than lq special primes having roots for this ad */
  if (n == 0 || n < k)
    return;

  unsigned long idx_q[n], tot, i;
  /* if tot (n, k) < wanted, use tot as wanted */
  tot =  binom (n, k);

  if (tot > (unsigned long) nq)
    tot = (unsigned long) nq;
  /*
   fprintf (stderr, "n=%"PRIu64", k=%"PRIu64", (n,k)=%"PRIu64",
   nq:%d\n", n, k, tot, nq);
  */

  /* enumerate each combination */
  first_comb (n, idx_q);
  q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
  gmp_collision_on_batch_p (header, R, q, qqz, rqqz);

  //fprintf (stderr, "# - q ");
  for (i = 1; i < tot; i ++) {
    next_comb (n, k, idx_q); // print_comb (k, idx_q);
    q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
    gmp_collision_on_batch_p (header, R, q, qqz, rqqz);
  }

  mpz_clear (qqz);
  mpz_clear (rqqz);
  qroots_clear (SQ_R);
}


static void
newAlgo (mpz_t N, unsigned long d, uint64_t ad)
{
  header_t header;
  header_init (header, N, d, ad);
  check_parameters (header->m0, d);

  proots_t R;
  proots_init (R, lenPrimes);

  if (sizeof (unsigned long int) == 8) {
    collision_on_p (header, R);
    collision_on_sq (header, R);
  }
  else {
    gmp_collision_on_p (header, R);
    gmp_collision_on_sq (header, R);
  }

  proots_clear (R, lenPrimes);
  header_clear (header);
}


#else // TAG


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long i, j, nprimes, p, nrp, c = 0;
  uint64_t *rp;
  int64_t ppl = 0;
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
  
#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS * (double) lenPrimes));
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    ppl = (int64_t) p;
    ppl *= (int64_t) p;

    /* add faked roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    mpz_mod_ui (f[0], header->Ntilde, p);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrp = poly_roots_uint64 (rp, f, header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);
    for (j = 0; j < nrp; j++, c++) {
      /* only consider r[j] and r[j] - pp */
      hash_add (H, p, (int64_t) rp[j], header->m0, header->ad, header->d,
		header->N, 1, tmp);
      hash_add (H, p, (int64_t) (rp[j] - ppl), header->m0, header->ad,
		header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", cputime () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

  /* if the hash table contains n entries, each one smaller than (2P)^2,
     the number of potential collisions is about 0.5*n^2/(2P)^2 */
  pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);
  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

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
  unsigned int nr, j;
  unsigned long nprimes, p, c = 0;
  uint64_t pp;
  int64_t ppl, u;
  double pc2;

#ifndef CONSIDER_ONLY_TWO_ROOTS
  int i;
  long h;
#endif  

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  hash_t H;

#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				  * (double) lenPrimes));
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    if ((header->d * header->ad) % p == 0)
      continue;

    /* set p, p^2, ppl */
    pp = p * p;
    ppl = (long) pp;
    nr = R->nr[nprimes];

    for (j = 0; j < nr; j++, c++)
    {
      u = (long) inv_qq[c];

#ifdef CONSIDER_ONLY_TWO_ROOTS
      hash_add (H, p, u, header->m0, header->ad, header->d, 
		header->N, q, rqqz);
      hash_add (H, p, u - ppl, header->m0, header->ad, header->d,
		header->N, q, rqqz);
#else
      h = u - ppl;
      for (i = 0; i < NUMBER_CONSIDERED_ROOTS; i ++) {
        hash_add (H, p, u, header->m0, header->ad, header->d,
		  header->N, q, rqqz);
        hash_add (H, p, h, header->m0, header->ad, header->d, 
		  header->N, q, rqqz);
        u += ppl;
        h -= ppl;
      }

#endif

    }  // next rp
  } // next p

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n", cputime () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  pc2 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);

  hash_clear (H);

  pthread_mutex_lock (&lock);
  potential_collisions += pc2;
  pthread_mutex_unlock (&lock);
}


/* Batch SQ mode */
static inline void
collision_on_batch_sq ( header_t header,
                        proots_t R,
                        unsigned long *q,
                        mpz_t *rqqz,
                        unsigned long size,
                        unsigned long number_pr )
{
  if (size == 0)
    return;

  unsigned int i, j, nr;
  unsigned long nprimes, p, c = 0, rp;
  uint64_t pp;
  unsigned long **invqq = malloc (size * sizeof (unsigned long *));

  if (invqq) {
    for (i = 0; i < size; i++)
      invqq[i] = malloc (number_pr * sizeof (unsigned long));
  }
  else {
    fprintf (stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
    exit (1);
  }

#ifdef DEBUG_POLYSELECT2L  /* check if crt roots are correct */
  for (i = 0; i < size; i ++) {
    mpz_t tmp_debug, tmp_debug2, qqz;
    mpz_init_set_ui (qqz, q[i]);
    mpz_mul_ui (qqz, qqz, q[i]);
    mpz_init_set (tmp_debug, header->Ntilde);
    mpz_mod (tmp_debug, tmp_debug, qqz);
    mpz_init_set (tmp_debug2, header->m0);
    mpz_add (tmp_debug2, tmp_debug2, rqqz[i]);
    mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
    mpz_mod (tmp_debug2, tmp_debug2, qqz);
    if (mpz_cmp (tmp_debug, tmp_debug2) != 0) {
      fprintf (stderr, "Error: crt root is wrong in %s, iter: %d, "
               "ad: %"PRIu64"\n",
               __FUNCTION__, i, header->ad);
      fprintf (stderr, "We should have Ntilde = (m0+r)^d mod q^2\n");
      gmp_fprintf (stderr, "Ntilde=%Zd m0=%Zd r=%Zd d=%d q=%lu\n",
                   header->Ntilde, header->m0, rqqz[i], header->d, q[i]);
      gmp_fprintf (stderr, "Ntilde mod q^2=%Zd\n", tmp_debug);
      gmp_fprintf (stderr, "(m0+r)^d mod q^2=%Zd\n", tmp_debug2);
      exit (1);
    }
    mpz_clear (tmp_debug);
    mpz_clear (tmp_debug2);
    mpz_clear (qqz);
    fprintf (stderr, "i: %u OK\n", i);
  }
#endif

  //int st = cputime();
   
  /* Step 1: batch inversion */
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {

    p = Primes[nprimes];
    pp = p*p;
    if ((header->d * header->ad) % p == 0)
      continue;
    if (R->nr[nprimes] == 0)
      continue;
    nr = R->nr[nprimes];
    
    modulusredcul_t modpp;
    residueredcul_t qprod[size], tmp_modul, tmp2_modul, tmp3_modul;
    residueredcul_t res_rp, res_tmp;

    modredcul_initmod_ul (modpp, pp);
    modredcul_init (tmp_modul, modpp);
    modredcul_init (tmp2_modul, modpp);
    modredcul_init (tmp3_modul, modpp);
    modredcul_init (res_rp, modpp);
    modredcul_init (res_tmp, modpp);
    for (i = 0; i < size; i++) {
      modredcul_init (qprod[i], modpp);
    }

    // (size -1) multiplications
    modredcul_set_ul (qprod[0], q[0], modpp);
    modredcul_sqr (qprod[0], qprod[0], modpp);
    for (i = 1; i < size; i ++) {
      modredcul_set_ul (tmp_modul, q[i], modpp);
      modredcul_sqr (tmp_modul, tmp_modul, modpp);
      modredcul_mul(qprod[i], tmp_modul, qprod[i-1], modpp);
    }
    // the inversion, also reduce qprod[size-1] (mod pp).
    modredcul_inv (tmp_modul, qprod[size-1], modpp);

    // for each q in a batch
    for (i = size-1; i > 0; i --) {

      modredcul_mul (tmp2_modul, qprod[i-1], tmp_modul, modpp);
      
#ifdef DEBUG_POLYSELECT2L /* check if batch inversions correct */
      unsigned long tmp_1 = modredcul_get_ul (tmp2_modul, modpp);
      mpz_t tmp_debug, tmp_debug2, qqz;
      mpz_init_set_ui (qqz, q[i]);
      mpz_mul_ui (qqz, qqz, q[i]);
      mpz_init_set_ui (tmp_debug, tmp_1);
      mpz_init (tmp_debug2);
      mpz_mul (tmp_debug2, tmp_debug, qqz);
      tmp_1 = mpz_fdiv_ui (tmp_debug2, pp);
      if (tmp_1 != 1) {
        fprintf (stderr, "Error, batch inversion is wrong in %s, iter: %d\n",
                 __FUNCTION__, i);
        exit (1);
      }
      mpz_clear (qqz);
      mpz_clear (tmp_debug);
      mpz_clear (tmp_debug2);
#endif

      // for each rp, compute (rp-rq)*1/q^2 (mod p^2)
      for (j = 0; j < nr; j ++, c++) {

        rp = R->roots[nprimes][j];
        modredcul_set_ul (res_rp, rp, modpp);
        modredcul_sub_ul (res_rp, res_rp, mpz_fdiv_ui (rqqz[i], pp), modpp);
        modredcul_mul (res_rp, res_rp, tmp2_modul, modpp);
        invqq[i][c] = modredcul_get_ul (res_rp, modpp);

#ifdef DEBUG_POLYSELECT2L  /* check if u is correct */
        mpz_t Ntmp, m0tmp, tmp, qqz;
        mpz_init_set_ui (qqz, q[i]);
        mpz_mul_ui (qqz, qqz, q[i]);
        mpz_init_set (Ntmp, header->Ntilde);
        mpz_mod_ui (Ntmp, Ntmp, pp); // tildeN (mod p^2)
        mpz_init_set (m0tmp, header->m0);
        mpz_add (m0tmp, m0tmp, rqqz[i]); // m0 + rq
        mpz_init_set_ui (tmp, invqq[i][c]);
        mpz_mul (tmp, tmp, qqz); // i*q^2
        mpz_add (m0tmp, m0tmp, tmp); // m0 + rq + i*q^2
        mpz_pow_ui (m0tmp, m0tmp, header->d);
        mpz_mod_ui (m0tmp, m0tmp, pp);
        if (mpz_cmp (m0tmp, Ntmp) != 0) {
          fprintf (stderr, "Error: u computation is wrong in"
                   "collision_on_each_sq, i: %d\n", i);
          gmp_printf ("Details: (p=%lu, i(mod p^2)=%ld) for "
                      "ad=%"PRIu64", q=%lu, rq=%Zd, rp=%lu, invqq=%lu\n",
                      p, invqq[i][c], header->ad, q[i], rqqz[i], rp, tmp2_modul[0]);
          gmp_printf ("m0=%Zd\n", header->m0);
          gmp_printf ("Ntilde=%Zd\n", header->Ntilde);
          exit (1);
        }
        mpz_clear (Ntmp);
        mpz_clear (m0tmp);
        mpz_clear (tmp);
        mpz_clear (qqz);
#endif
      }
      
      modredcul_set_ul (tmp3_modul, q[i], modpp);
      modredcul_sqr (tmp3_modul, tmp3_modul, modpp);
      modredcul_mul (tmp_modul, tmp3_modul, tmp_modul, modpp);
      c -= nr;
    }

    // last q in the batch is in tmp_modul
    // invqq[0][nprimes] = modredcul_get_ul (tmp_modul, modpp);
    for (j = 0; j < nr; j ++, c ++) {
      rp = R->roots[nprimes][j];
      modredcul_set_ul (res_rp, rp, modpp);
      modredcul_sub_ul (res_rp, res_rp, mpz_fdiv_ui (rqqz[0], pp), modpp);
      modredcul_mul (tmp2_modul, res_rp, tmp_modul, modpp); // tmp_modul should be retained!
      invqq[0][c] = modredcul_get_ul (tmp2_modul, modpp);

#ifdef DEBUG_POLYSELECT2L /* check if batch inversions correct */
      unsigned long tmp_1 = modredcul_get_ul (tmp_modul, modpp);
      mpz_t tmp_debug, tmp_debug2, qqz1;
      mpz_init_set_ui (qqz1, q[0]);
      mpz_mul_ui (qqz1, qqz1, q[0]);
      mpz_init_set_ui (tmp_debug, tmp_1);
      mpz_init (tmp_debug2);
      mpz_mul (tmp_debug2, tmp_debug, qqz1);
      tmp_1 = mpz_fdiv_ui (tmp_debug2, pp);
      if (tmp_1 != 1) {
        fprintf (stderr, "Error, batch inversion is wrong in %s, iter: %d\n",
                 __FUNCTION__, i);
        exit (1);
      }
      mpz_clear (qqz1);
      mpz_clear (tmp_debug);
      mpz_clear (tmp_debug2);
#endif

#ifdef DEBUG_POLYSELECT2L  /* check if u is correct */
      mpz_t Ntmp, m0tmp, tmp, qqz;
      mpz_init_set_ui (qqz, q[0]);
      mpz_mul_ui (qqz, qqz, q[0]);
      mpz_init_set (Ntmp, header->Ntilde);
      mpz_mod_ui (Ntmp, Ntmp, pp); // tildeN (mod p^2)
      mpz_init_set (m0tmp, header->m0);
      mpz_add (m0tmp, m0tmp, rqqz[0]); // m0 + rq
      mpz_init_set_ui (tmp, invqq[0][c]);
      mpz_mul (tmp, tmp, qqz); // i*q^2
      mpz_add (m0tmp, m0tmp, tmp); // m0 + rq + i*q^2
      mpz_pow_ui (m0tmp, m0tmp, header->d);
      mpz_mod_ui (m0tmp, m0tmp, pp);
      if (mpz_cmp (m0tmp, Ntmp) != 0) {
        fprintf (stderr, "Error: u computation is wrong in"
                 "collision_on_each_sq, i: %d\n", i);
        gmp_printf ("Details: (p=%lu, i(mod p^2)=%ld) for "
                    "ad=%"PRIu64", q=%lu, rq=%Zd, rp=%lu, invqq=%lu\n",
                    p, invqq[0][c], header->ad, q[0], rqqz[0], rp, tmp_modul[0]);
        gmp_printf ("m0=%Zd\n", header->m0);
        gmp_printf ("Ntilde=%Zd\n", header->Ntilde);
        exit (1);
      }
      mpz_clear (Ntmp);
      mpz_clear (m0tmp);
      mpz_clear (tmp);
      mpz_clear (qqz);
#endif
    }

    modredcul_clear (res_rp, modpp);
    modredcul_clear (res_tmp, modpp);
    modredcul_clear (tmp_modul, modpp);
    modredcul_clear (tmp2_modul, modpp);
    modredcul_clear (tmp3_modul, modpp);
    for (i = 0; i < size; i++)
      modredcul_clear (qprod[i], modpp);
    modredcul_clearmod (modpp);

  } // next prime p

  //fprintf (stderr, "# one batch SQ inversion took %dms\n", cputime () - st);

  /* Step 2: find collisions on q. */
  for (i = 0; i < size; i ++) {
    //int st2 = cputime();
    collision_on_each_sq ( header,
                           R,
                           q[i],
                           rqqz[i],
                           invqq[i] );

    //printf ("# outer collision_on_each_sq took %dms\n", cputime () - st2);
  }

  for (i = 0; i < size; i++)
    free (invqq[i]);
  free (invqq);
}


/* collision on special-q, call collision_on_batch_sq */
static inline void
collision_on_sq ( header_t header,
                  proots_t R,
                  unsigned long c )
{
  // init special-q roots
  qroots_t SQ_R;
  qroots_init (SQ_R);
  comp_sq_roots (header, SQ_R);
  // qroots_print (SQ_R);

  unsigned long K = lq, N = SQ_R->size, tot, i, l;
  unsigned long idx_q[K], q[BATCH_SIZE];
  mpz_t qqz[BATCH_SIZE], rqqz[BATCH_SIZE];

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
      q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

      for (l = 1; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
      }
    }
    else {
      for (l = 0; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
      }
    }

#ifdef DEBUG_POLYSELECT2L
    unsigned long j;
    for (j = 0; j < BATCH_SIZE; j++)
      gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n", q[j], qqz[j], rqqz[j]);
#endif

    // collision batch
    collision_on_batch_sq ( header,
                            R,
                            q,
                            rqqz,
                            BATCH_SIZE,
                            c );
    i += BATCH_SIZE;
  }

  // tail batch
  for (l = 0; l < (tot % BATCH_SIZE); l++) {
    next_comb (N, K, idx_q);
    q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

#ifdef DEBUG_POLYSELECT2L
    gmp_fprintf (stderr, "q: %lu, qq: %Zd, rqq: %Zd\n",
		 q[l], qqz[l], rqqz[l]);
#endif

  }

  collision_on_batch_sq ( header,
                          R,
                          q,
                          rqqz,
                          tot % BATCH_SIZE,
                          c );

  for (l = 0; l < BATCH_SIZE; l++) {
    mpz_clear (qqz[l]);
    mpz_clear (rqqz[l]);
  }
  qroots_clear (SQ_R);
}


// separator between modredc_ul and gmp


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
gmp_collision_on_p ( header_t header,
		     proots_t R )
{
  unsigned long i, j, nprimes, p, nrp, c = 0;
  uint64_t *rp;
  int64_t ppl = 0;
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
  
#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				  * (double) lenPrimes));
#endif

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  for (nprimes = 0; nprimes < lenPrimes; nprimes ++) {
    p = Primes[nprimes];
    ppl = (int64_t) p;
    ppl *= (int64_t) p;

    /* add faked roots to keep indices */
    if ((header->d * header->ad) % p == 0) {
      R->nr[nprimes] = 0; // nr = 0.
      R->roots[nprimes] = NULL;
      continue;
    }

    /* we want p^2 | N - (m0 + i)^d, thus
       (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
    mpz_mod_ui (f[0], header->Ntilde, p);
    mpz_neg (f[0], f[0]); /* f = x^d - N */
    nrp = poly_roots_uint64 (rp, f, header->d, p);
    roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
    proots_add (R, nrp, rp, nprimes);
    for (j = 0; j < nrp; j++, c++) {
      /* only consider r[j] and r[j] - pp */
      gmp_hash_add (H, p, (int64_t) rp[j], header->m0, header->ad,
		    header->d, header->N, 1, tmp);
      gmp_hash_add (H, p, (int64_t) (rp[j] - ppl), header->m0, header->ad,
		header->d, header->N, 1, tmp);
    }
  }

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# collision_on_p took %dms\n", cputime () - st);
  fprintf (stderr, "# p hash_size: %u for ad = %lu\n", H->size, header->ad);
#endif

  /* if the hash table contains n entries, each one smaller than (2P)^2,
     the number of potential collisions is about 0.5*n^2/(2P)^2 */
  pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);
  hash_clear (H);

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  free (rp);
  mpz_clear (tmp);

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
  int64_t ppl, u;
  double pc2;

#ifndef CONSIDER_ONLY_TWO_ROOTS
  int i;
  int64_t h;
#endif  

#ifdef DEBUG_POLYSELECT2L
  int st = cputime();
#endif

  hash_t H;

#ifdef CONSIDER_ONLY_TWO_ROOTS
  hash_init (H,  (unsigned long) (INIT_FACTOR * (double) lenPrimes));
#else
  hash_init (H,  (unsigned long) (INIT_FACTOR * NUMBER_CONSIDERED_ROOTS
				  * (double) lenPrimes));
#endif

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

#ifdef CONSIDER_ONLY_TWO_ROOTS
      gmp_hash_add (H, p, u, header->m0, header->ad, header->d, 
		    header->N, q, rqqz);
      gmp_hash_add (H, p, u - ppl, header->m0, header->ad, header->d,
		    header->N, q, rqqz);
#else
      h = u - ppl;
      for (i = 0; i < NUMBER_CONSIDERED_ROOTS; i ++) {
        gmp_hash_add (H, p, u, header->m0, header->ad, header->d,
		      header->N, q, rqqz);
        gmp_hash_add (H, p, h, header->m0, header->ad, header->d,
		      header->N, q, rqqz);
        u += ppl;
        h -= ppl;
      }

#endif

    }  // next rp
  } // next p

#ifdef DEBUG_POLYSELECT2L
  fprintf (stderr, "# inner collision_on_each_sq took %dms\n",
	   cputime () - st);
  fprintf (stderr, "# - q hash_size (q=%lu): %u\n", q, H->size);
#endif

  pc2 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 1], 2.0);

  hash_clear (H);

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
    //int st2 = cputime();
    gmp_collision_on_each_sq (header, R, q[i], rqqz[i], invqq[i]);
    //printf ("# outer collision_on_each_sq took %dms\n", cputime () - st2);
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
      q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

      for (l = 1; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
      }
    }
    else {
      for (l = 0; l < BATCH_SIZE; l++) {
        next_comb (N, K, idx_q);
        q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);
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
    q[l] = return_q_rq (SQ_R, idx_q, K, qqz[l], rqqz[l]);

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
  check_parameters (header->m0, d);
  proots_init (R, lenPrimes);

  if (sizeof (unsigned long int) == 8) {
    c = collision_on_p (header, R);
    collision_on_sq (header, R, c);
  }
  else {
    c = gmp_collision_on_p (header, R);
    gmp_collision_on_sq (header,R, c);
  }

  proots_clear (R, lenPrimes);
  header_clear (header);
}


#endif // end collision_on_sq(), whatever batch mode is used.


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
  fprintf (stderr, "-lq nnn      --- number of factors in the special-q"
           " (default %d)\n", LQ_DEFAULT);
  fprintf (stderr, "-seed nnn    --- seed for srand (default by time(NULL))\n");
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
  param_list_parse_int (pl, "lq", &lq);
  param_list_parse_int (pl, "seed", &seed);
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
  if (lq < 1 || nq < 1) {
    fprintf (stderr, "Error, number of factors in special-q "
             "should >= 1 and/or number of special-q's should >=1\n");
    exit (1);
  }

  /* check seed */
  if (seed < 1) {
    seed = time(NULL);
    srand(seed);
  }
  else
    srand(seed);

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

  st = cputime ();
  lenPrimes = initPrimes (P, &Primes);

  printf ( "# Info: initializing %lu P primes took %dms, seed=%d,"
           " rawonly=%d, nq=%d, target_time=%d\n",
           lenPrimes,
           cputime () - st,
           seed,
           raw,
           nq,
           target_time / 1000 );

#ifdef CONSIDER_ONLY_TWO_ROOTS

#ifdef BATCH_P
  printf ( "# Info: estimated peak memory=%.2fMB (%d thread(s),"
           " batch %d inversions on P)\n",
           (double) (nthreads * INIT_FACTOR * lenPrimes *
           (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           BATCH_SIZE );
#else
  printf ( "# Info: estimated peak memory=%.2fMB (%d thread(s),"
           " batch %d inversions on SQ)\n",
           (double) (nthreads * (BATCH_SIZE + INIT_FACTOR) * lenPrimes
           * (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           BATCH_SIZE );
#endif

#else

#ifdef BATCH_P
  printf ( "# Info: estimated peak memory=%.2fMB (%d threads,"
           " batch %d inversions on P)\n",
           (double) (nthreads * INIT_FACTOR * NUMBER_CONSIDERED_ROOTS *
           lenPrimes * (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           BATCH_SIZE );
#else
  printf ( "# Info: estimated peak memory=%.2fMB (%d threads,"
           " batch %d inversions on SQ)\n",
           (double) (nthreads * (BATCH_SIZE + INIT_FACTOR * 
           NUMBER_CONSIDERED_ROOTS) * lenPrimes * (sizeof(uint32_t)
            + sizeof(uint64_t)) / 1024 / 1024),
           nthreads,
           BATCH_SIZE );
#endif


#endif

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

    if (cputime () > target_time || verbose > 0)
    {
      printf ("# Stat: ad=%lu, exp. coll.=%1.1f (%0.3f/s), got %lu with %lu good ones, av. lognorm=%1.2f, av. raw. lognorm=%1.2f, time=%dms\n",
              admin,
              potential_collisions,
              1000.0 * (double) potential_collisions / cputime (),
              collisions,
              collisions_good,
              aver_opt_lognorm / collisions_good,
              aver_raw_lognorm / collisions,
              cputime () );
      fflush (stdout);
      target_time += incr_target_time;
    }
  }

  /* finishing up statistics */
  if (verbose >= 0)
    {
      printf ("# Stat: tried %d ad-value(s), found %d polynomial(s), %d below maxnorm\n",
              tries, tot_found, found);
      printf ("# Stat: raw lognorm (min/av/max): %1.2f/%1.2f/%1.2f\n",
              min_raw_lognorm, aver_raw_lognorm / collisions, max_raw_lognorm);
      printf ("# Stat: optimized lognorm (min/av/max): %1.2f/%1.2f/%1.2f\n",
              min_opt_lognorm, aver_opt_lognorm / collisions_good,
              max_opt_lognorm);
      printf ("# Stat: potential collisions=%1.2e (%1.2e/s)\n",
              potential_collisions, 1000.0 * potential_collisions
              / (double) cputime ());
      printf ("# Stat: av. g0/adm2 ratio: %.3e\n",
              total_adminus2 / (double) collisions);
      if (d == 6)
        printf ("# Stat: av. logmu noc4/noc3 ratio: %.3f\n",
                aver_lognorm_ratio / (double) collisions);
    }

  printf ("# Stat: tried %d ad-value(s), found %d polynomial(s), %d below maxnorm\n",
          tries, tot_found, found);

  /* print best 10 values of logmu */
  printf ("# Stat: best logmu:");
  for (i = 0; i < 10; i++)
    printf (" %1.2f", best_logmu[i]);
  printf ("\n");

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
