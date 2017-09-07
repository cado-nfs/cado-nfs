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
#include "twocubics.h"
#include "portability.h"
#include "mpz_poly.h"
#include "area.h"

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */
#define NEW_ROOTSIEVE
#define INIT_FACTOR 8UL
#define PREFIX_HASH
//#define DEBUG_POLYSELECT

#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#endif

#ifdef PREFIX_HASH
char *phash = "# ";
#else
char *phash = "";
#endif

#define BATCH_SIZE 20 /* number of special (q, r) per batch */

/* Read-Only */
uint32_t *Primes = NULL;
unsigned long lenPrimes = 1; // length of Primes[]
int nq = INT_MAX;
static int verbose = 0;
static unsigned long incr = DEFAULT_INCR;
cado_poly best_poly, curr_poly;
double best_E = DBL_MAX; /* combined score E (the smaller the better) */
double aver_E = 0.0;
unsigned long found = 0; /* number of polynomials found so far */

/* read-write global variables */
static pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                   lock for those variables */

mpz_t maxS; /* maximun skewness. O for default max */

/* inline function */
extern void shash_add (shash_t, uint64_t);

/* -- functions starts here -- */

static void
mutex_lock(pthread_mutex_t *lock)
{
  pthread_mutex_lock (lock);
}

static void
mutex_unlock(pthread_mutex_t *lock)
{
  pthread_mutex_unlock (lock);
}

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


/* Compute maximun skewness, which in floor(N^(1/d^2)) */
void
compute_default_max_skew (mpz_t skew, mpz_t N, int d)
{
  mpz_root (skew, N, (unsigned long) d*d);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       mpz_t ad, unsigned long d, mpz_t N, unsigned long q, mpz_t rq)
{
#if 0
  unsigned long j;
  int cmp;
  double skew, logmu, E;
  mpz_poly F;
#endif

  mpz_t l, r, k, mprime, Nprime, C, l2, tmp, r1, r0, t, adm1, m, skew, root;
  mpz_vector_t a, b, reduced_a, reduced_b;
  mpz_poly f, g;

  mpz_init (root);
  mpz_init (m);
  mpz_init (adm1);
  mpz_init (l);
  mpz_init (l2);
  mpz_init (r);
  mpz_init (k);
  mpz_init (mprime);
  mpz_init (Nprime);
  mpz_init (C);
  mpz_init (tmp);
  mpz_init (r1);
  mpz_init (t);
  mpz_init (r0);
  mpz_init (skew);

  mpz_vector_init (a, d+1);
  mpz_vector_init (b, d+1);
  mpz_vector_init (reduced_a, d+1);
  mpz_vector_init (reduced_b, d+1);

  mpz_poly_init (f, d);
  mpz_poly_init (g, d);

#ifdef DEBUG_POLYSELECT
  gmp_printf ("#### MATCH ######\nN = %Zd\nd = %d\nad = %Zd\n"
              "p1 = %lu\np2 = %lu\nq = %lu\ni = %" PRId64 "\n"
              "rq = %Zd\n", N, d, ad, p1, p2, q, i, rq);
#endif

  /* l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_mul_ui (l, l, q);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("l = p1 * p2 * q\nl == %Zd # has %lu bits\n", l,
              mpz_sizeinbase(l, 2));
#endif
  mpz_mul (l2, l, l); /* l2 = l^2 */

  /* r = rq + i*q^2 */
  mpz_set_si (r, i);
  mpz_mul_ui (r, r, q);
  mpz_mul_ui (r, r, q);
  mpz_add (r, r, rq);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("r = rq + i * q^2 \nr == %Zd # has %lu bits\n", r,
              mpz_sizeinbase(r, 2));
#endif

  /* k = d^d * ad^(d-1) */
  mpz_mul_ui (k, ad, d);
  mpz_pow_ui (k, k, d-1);
  mpz_mul_ui (k, k, d);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("k = d^d * ad^(d-1)\nk == %Zd\n", k);
#endif

  /* Nprime = k * N */
  mpz_mul (Nprime, k, N);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("Nprime = k * N\nNprime == %Zd\n", Nprime);
#endif

  /* mprime = m0 + r */
  mpz_add (mprime, m0, r);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("m0 = %Zd\nmprime = m0 + r\nmprime == %Zd\n", m0, mprime);
#endif

  /* C = mprime^d - Nprime */
  mpz_pow_ui (C, mprime, d);
  mpz_sub (C, C, Nprime);
  ASSERT_ALWAYS (mpz_divisible_p (C, l2));
#ifdef DEBUG_POLYSELECT
  gmp_printf ("(mprime^d - Nprime) %% l^2 == 0\n");
#endif

  /* adm1 is such that mprime = d*ad*m + adm1*l and -d*ad/2 <= adm1 < d*ad/2
     We have adm1 = mprime/l mod (d*ad). */
  mpz_mul_ui (tmp, ad, d); /* tmp = d*ad */
  if (mpz_invert (adm1, l, tmp) == 0)
  {
    fprintf (stderr, "Error in 1/l mod (d*ad)\n");
    abort();
  }
  mpz_mul (adm1, adm1, mprime);
  mpz_mod (adm1, adm1, tmp);
  mpz_mul_2exp (t, adm1, 1);
  if (mpz_cmp (t, tmp) >= 0)
    mpz_sub (adm1, adm1, tmp);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("adm1 = %Zd\n", adm1);
#endif

  /* m = (mprime - adm1 * l)/ (d * ad) */
  mpz_mul (m, adm1, l);
  mpz_sub (m, mprime, m);
  ASSERT_ALWAYS (mpz_divisible_ui_p (m, d));
  mpz_divexact_ui (m, m, d);
  ASSERT_ALWAYS (mpz_divisible_p (m, ad));
  mpz_divexact (m, m, ad);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("m = (mprime - adm1*l) / (d*ad)\nm == %Zd\n", m);
#endif

  /* Set vector a = (-m, l, 0, ..., 0) */
  mpz_neg (tmp, m);
  mpz_vector_setcoordinate (a, 0, tmp); /* a[0] = -m */
  mpz_vector_setcoordinate (a, 1, l); /* a[1] = -l */
  for (unsigned int j = 2; j <= d; j++)
    mpz_vector_setcoordinate_ui (a, j, 0); /* a[j] = 0 */

  /* Set vector b = (a0, a1, ..., ai, ..., adm1, ad) */
  mpz_vector_setcoordinate (b, d, ad);  /* b[d] = ad */
  mpz_vector_setcoordinate (b, d-1, adm1); /* b[d-1] = adm1 */




  mpz_pow_ui (t, m, d);
  mpz_mul (t, t, ad);
  mpz_sub (t, N, t);
  ASSERT_ALWAYS (mpz_divisible_p (t, l));

  mpz_divexact (t, t, l);
  mpz_pow_ui (tmp, m, d-1);
  mpz_mul (tmp, tmp, adm1);
  mpz_sub (t, t, tmp);
  for (int j = d - 2; j > 0; j--)
  {
    ASSERT_ALWAYS (mpz_divisible_p (t, l));
    mpz_divexact (t, t, l);
    /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
    mpz_pow_ui (tmp, m, j);
    /* fdiv rounds toward -infinity: r1 = floor(t/tmp) */
    mpz_fdiv_q (r1, t, tmp); /* t -> r1 * tmp + t */
    mpz_invert (k, tmp, l); /* search r1 + k such that */

    mpz_mul (k, k, t);
    mpz_sub (k, k, r1);
    mpz_mod (k, k, l);

    mpz_mul_2exp (k, k, 1);
    int cmp = mpz_cmp (k, l);
    mpz_div_2exp (k, k, 1);
    if (cmp >= 0)
      mpz_sub (k, k, l);
    mpz_add (r1, r1, k);
    mpz_vector_setcoordinate (b, j, r1);
    /* subtract r1*m^j */
    mpz_submul (t, tmp, r1);
  }
  ASSERT_ALWAYS (mpz_divisible_p (t, l));
  mpz_divexact (t, t, l);
  mpz_vector_setcoordinate (b, 0, t);


  mpz_vector_get_mpz_poly(f, a);
  mpz_vector_get_mpz_poly(g, b);
#ifdef DEBUG_POLYSELECT
  printf ("a = ");
  mpz_poly_fprintf (stdout, f);
  printf ("b = ");
  mpz_poly_fprintf (stdout, g);
#endif

  mpz_vector_reduce_with_max_skew (reduced_a, reduced_b, skew, a, b, maxS, d);

  mpz_vector_get_mpz_poly(f, reduced_a);
  mpz_vector_get_mpz_poly(g, reduced_b);

#ifdef DEBUG_POLYSELECT
  gmp_printf ("skew = %Zd\nf = ", skew);
  mpz_poly_fprintf (stdout, f);
  printf ("g = ");
  mpz_poly_fprintf (stdout, g);
#endif

  mpz_invert (root, l, N);
  mpz_mul (root, root, m);
  mpz_mod (root, root, N);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("root = (m / l) %% N\nroot == %Zd\n", root);

  gmp_printf ("## Begin poly file for ad = %Zd and l = %Zd\n", ad, l);
#endif

  double skewness, logmu[2], alpha[2], E;

  skewness = L2_combined_skewness2 (f, g, SKEWNESS_DEFAULT_PREC);
  logmu[0] = L2_lognorm (g, skewness);
  logmu[1] = L2_lognorm (f, skewness);
  alpha[0] = get_alpha (g, ALPHA_BOUND_SMALL);
  alpha[1] = get_alpha (f, ALPHA_BOUND_SMALL);
  E = logmu[1] + alpha[1] + logmu[0] + alpha[0];

  found ++;
  aver_E += E;

  if (E < best_E)
    {
      mutex_lock (&lock);
      best_E = E;
      gmp_printf("n: %Zd\n", N);
      for (int i = 0; i <= f->deg; i++)
        gmp_printf ("c%d: %Zd\n", i, f->coeff[i]);
      for (int i = 0; i <= g->deg; i++)
        gmp_printf ("Y%d: %Zd\n", i, g->coeff[i]);
      printf ("skew: %1.2f\n", skewness);
      printf ("# f lognorm %1.2f, alpha %1.2f, score %1.2f\n",
              logmu[1], alpha[1], logmu[1] + alpha[1]);
      printf ("# g lognorm %1.2f, alpha %1.2f, score %1.2f\n",
              logmu[0], alpha[0], logmu[0] + alpha[0]);
      printf ("# f+g score %1.2f\n", E);
      printf ("# found %lu polynomial(s) so far, aver. E = %1.2f\n",
              found, aver_E / (double) found);
#ifdef DEBUG_POLYSELECT
      printf ("## End poly file\n");
#else
      printf ("\n");
#endif
      mutex_unlock (&lock);
    }
  
  mpz_clear (root);
  mpz_clear (skew);
  mpz_clear (adm1);
  mpz_clear (l);
  mpz_clear (l2);
  mpz_clear (C);
  mpz_clear (r);
  mpz_clear (k);
  mpz_clear (mprime);
  mpz_clear (m);
  mpz_clear (Nprime);
  mpz_clear (tmp);
  mpz_clear (r1);
  mpz_clear (r0);
  mpz_clear (t);

  mpz_poly_clear (f);
  mpz_poly_clear (g);

  mpz_vector_clear (a);
  mpz_vector_clear (b);
  mpz_vector_clear (reduced_a);
  mpz_vector_clear (reduced_b);
}

void
gmp_match (uint32_t p1, uint32_t p2, int64_t i, mpz_t m0,
	   mpz_t ad, unsigned long d, mpz_t N, uint64_t q, mpz_t rq)
{
  match (p1, p2, i, m0, ad, d, N, q, rq);
}


/* find collisions between "P" primes, return number of loops */
static inline unsigned long
collision_on_p ( header_t header,
                 proots_t R )
{
  unsigned long i, j, nprimes, p, nrp, c = 0, tot_roots = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  mpz_t *f, tmp;
  int found = 0;
  shash_t H;
  int st = milliseconds ();

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
  shash_reset (H);
  umax = (int64_t) Primes[lenPrimes - 1] * (int64_t) Primes[lenPrimes - 1];
  for (nprimes = 0; nprimes < lenPrimes; nprimes ++)
    {
      p = Primes[nprimes];
      ppl = (int64_t) p * (int64_t) p;

      /* add fake roots to keep indices */
      if (header_skip (header, p))
        {
          R->nr[nprimes] = 0; // nr = 0.
          R->roots[nprimes] = NULL;
          continue;
        }

      nrp = roots_mod_uint64 (rp, mpz_fdiv_ui (header->Ntilde, p), header->d,
                              p);
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
    printf ("# computing %lu p-roots took %dms\n", tot_roots, st);

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
#ifdef DEBUG_POLYSELECT
      printf ("# collision_on_p took %lums\n", milliseconds () - st);
      gmp_printf ("# p hash_size: %u for ad = %Zd\n", H->size, header->ad);
#endif

#ifdef DEBUG_HASH_TABLE
      printf ("# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
      printf ("# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif
      hash_clear (H);
    }

  for (i = 0; i <= header->d; i++)
    mpz_clear (f[i]);
  free (f);
  mpz_clear (tmp);

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
  uint64_t pp;
  int64_t ppl, neg_umax, umax, v1, v2, nv;
  unsigned long p, nprimes, c;
  uint8_t vpnr, *pnr, nr, j;
  uint32_t *pprimes, i;
  int found;

#ifdef DEBUG_POLYSELECT
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
  shash_reset (H);

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
          if (header_skip (header, p))
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

#ifdef DEBUG_POLYSELECT
  printf ("# inner collision_on_each_sq took %lums\n", milliseconds () - st);
  printf ("# - q hash_alloc (q=%lu): %u\n", q, H->alloc);
#endif

#ifdef DEBUG_HASH_TABLE
  printf ("# p hash_size: %u, hash_alloc: %u\n", H->size, H->alloc);
  printf ("# hash table coll: %lu, all_coll: %lu\n", H->coll, H->coll_all);
#endif
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
    printf ("#  substage: batch %d many (rp-rq)*1/q^2 took %lums\n",
            count, milliseconds () - st);
    st = milliseconds();
  }

  /* core function to find collisions */
  for (k = 0; k < count; k ++) {
    collision_on_each_sq (header, R, q, rqqz[k], tinv_qq[k]);
  }

  if (verbose > 2)
    printf ("#  substage: collision-detection %d many rq took %lums\n",
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
  printf ("q: %lu, ", q);
  for (i = 0; i < lq; i ++)
    printf ("%u ", SQ_R->q[idx_q[i]]);
  printf (", ");
  for (i = 0; i < lq; i ++)
    printf ("%u ", SQ_R->nr[idx_q[i]]);
  printf ("\n");
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
    if (header_skip (header, p))
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
    printf ("# stage (1/q^2 inversion) for %lu primes took %lums\n",
            lenPrimes, milliseconds () - st);

  /* Step 2: find collisions on q. */
  int st2 = milliseconds();

  collision_on_batch_sq_r ( header, R, SQ_R, q, idx_q, invqq, number_pr,
                            &curr_nq, lq );
  if (verbose > 2)
    printf ("#  stage (special-q) for %d special-q's took %lums\n",
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

  for (j = 0; j < lq; j ++)
    idx_q[j] = j;
  q = return_q_norq (SQ_R, idx_q, lq);

  /* collision batch */
  collision_on_batch_sq (header, R, SQ_R, q, idx_q, c, lq);

  /* clean */
  qroots_clear (SQ_R);
  return;
}


static void
newAlgo (mpz_t N, unsigned long d, mpz_t ad)
{
  unsigned long c = 0;
  header_t header;
  proots_t R;

  header_init (header, N, d, ad);
  proots_init (R, lenPrimes);

  c = collision_on_p (header, R);
  if (nq > 0)
    collision_on_sq (header, R, c);

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
declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "degree", "(alias d) polynomial degree (2 or 3, "
                                      "default is 3)");
  param_list_decl_usage(pl, "n", "(required, alias N) input number");
  param_list_decl_usage(pl, "P", "(required) deg-1 coeff of g(x) has two prime factors in [P,2P]\n");

  param_list_decl_usage(pl, "admax", "max value for ad (+1)");
  param_list_decl_usage(pl, "admin", "min value for ad (default 0)");
  param_list_decl_usage(pl, "incr", "(alias i) increment of ad (default 60)");
  param_list_decl_usage(pl, "skewness", "maximun skewness possible "
                                        "(default N^(1/9))");
  param_list_decl_usage(pl, "maxtime", "stop the search after maxtime seconds");

  char str[200];
  snprintf (str, 200, "maximum number of special-q's considered\n"
            "               for each ad (default %d)", INT_MAX);
  param_list_decl_usage(pl, "nq", str);
  param_list_decl_usage(pl, "keep", "number of polynomials kept (default 10)");
  snprintf(str, 200, "time interval (seconds) for printing statistics (default %d)", TARGET_TIME / 1000);
  param_list_decl_usage(pl, "s", str);
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
  param_list_decl_usage(pl, "q", "(switch) quiet mode");
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
  verbose_decl_usage(pl);
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
  char **argv0 = argv;
  double st0 = seconds (), maxtime = DBL_MAX;
  mpz_t N, admin, admax;
  unsigned int d = 3;
  unsigned long P = 0;
  int quiet = 0, tries = 0, i, nthreads = 1, st,
    target_time = TARGET_TIME;
  tab_t *T;
  pthread_t *tid;

  mpz_init (N);
  mpz_init (maxS);
  cado_poly_init (best_poly);
  cado_poly_init (curr_poly);

  /* read params */
  param_list pl;
  param_list_init (pl);

  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-q", &quiet);
  param_list_configure_alias(pl, "incr", "-i");
  param_list_configure_alias(pl, "n", "-N");
  param_list_configure_alias(pl, "degree", "-d");

  if (argc == 1)
    usage (argv0[0], NULL, pl);

  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0[0], NULL, pl);
  }

  /* parse and check N in the first place */
  int have_n = param_list_parse_mpz(pl, "n", N);

  if (!have_n) {
    printf ("# Reading n from stdin\n");
    param_list_read_stream(pl, stdin, 0);
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
  param_list_parse_int (pl, "s", &target_time);
  param_list_parse_uint (pl, "degree", &d);
  ASSERT_ALWAYS (2 <= d && d <= 3);
  if (param_list_parse_double (pl, "area", &area) == 0) /* no -area */
    area = AREA;
  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;

  mpz_init (admin);
  param_list_parse_mpz (pl, "admin", admin);

  mpz_init (admax);
  if (param_list_parse_mpz (pl, "admax", admax) == 0) /* no -admax */
    mpz_set (admax, N);

  param_list_parse_ulong (pl, "incr", &incr);
  param_list_parse_double (pl, "maxtime", &maxtime);

  if (!param_list_parse_mpz(pl, "skewness", maxS))
    mpz_set_ui (maxS, 0);
  else if (mpz_cmp_ui (maxS, 1) < 0)
  {
    gmp_fprintf(stderr, "Error, skewness (%Zd) should be greater or equal "
                        "to 1\n", maxS);
    abort();
  }

  if (param_list_warn_unused(pl))
    usage (argv0[0], NULL, pl);

  /* print command line */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* check lq and nq */
  if (nq < 0) {
    fprintf (stderr, "Error, number of special-q's should >= 0\n");
    exit (1);
  }

  /* allocate threads */
  tid = malloc (nthreads * sizeof (pthread_t));
  ASSERT_ALWAYS(tid != NULL);

  /* quiet mode */
  if (quiet == 1)
    verbose = -1;

  /* set cpoly */
  mpz_set (best_poly->n, N);
  mpz_set (curr_poly->n, N);
  best_poly->pols[ALG_SIDE]->deg = d;
  best_poly->pols[RAT_SIDE]->deg = d;
  curr_poly->pols[ALG_SIDE]->deg = d;
  curr_poly->pols[RAT_SIDE]->deg = d;

  /* Compute maxS: use maxS if maxS argument is greater than 0 and lesser than
     default value */
  mpz_t tmp;
  mpz_init (tmp);
  compute_default_max_skew (tmp, N, 2);
  if (mpz_cmp_ui(maxS, 0) == 0 || mpz_cmp(maxS, tmp) > 0)
    mpz_set (maxS, tmp);

  /* Initialize primes. Since we use modredc to perform computations mod p^2
     for p < 2P, we need that 4P^2 fits in an unsigned long. */
  double Pd;
  Pd = (double) P;
  if (4.0 * Pd * Pd >= (double) ULONG_MAX)
    {
      mpz_set_ui (tmp, ULONG_MAX >> 2);
      mpz_sqrt (tmp, tmp);
      gmp_fprintf (stderr, "Error, too large value of P, maximum is %Zd\n",
                   tmp);
      exit (1);
    }
  mpz_clear (tmp);
  if (P <= (unsigned long) SPECIAL_Q[LEN_SPECIAL_Q - 2]) {
    fprintf (stderr, "Error, too small value of P, need P > %u\n",
             SPECIAL_Q[LEN_SPECIAL_Q - 2]);
    exit (1);
  }

  /* detect L1 cache size */
  ropt_L1_cachesize ();

  st = milliseconds ();
  lenPrimes = initPrimes (P, &Primes);

  printf ( "# Info: initializing %lu P primes took %lums,"
           " nq=%d, target_time=%d\n",
           lenPrimes, milliseconds () - st, nq, target_time / 1000 );


  printf ( "# Info: estimated peak memory=%.2fMB (%d thread(s),"
           " batch %d inversions on SQ)\n",
           (double) (nthreads * (BATCH_SIZE * 2 + INIT_FACTOR) * lenPrimes
           * (sizeof(uint32_t) + sizeof(uint64_t)) / 1024 / 1024),
           nthreads, BATCH_SIZE );

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
    mpz_init (T[i]->ad);
    T[i]->d = d;
    T[i]->thread = i;
  }

  if (incr <= 0)
  {
    fprintf (stderr, "Error, incr should be positive\n");
    exit (1);
  }

  /* admin should be nonnegative */
  if (mpz_cmp_ui (admin, 0) < 0)
    {
      fprintf (stderr, "Error, admin should be nonnegative\n");
      exit (1);
    }

  /* ensure admin > 0 and admin is a multiple of incr */
  if (mpz_cmp_ui (admin, 0) == 0)
    mpz_set_ui (admin, incr);
  else if (mpz_fdiv_ui (admin, incr) != 0)
    mpz_add_ui (admin, admin, incr - mpz_fdiv_ui (admin, incr));

  while (mpz_cmp (admin, admax) < 0 && seconds () - st0 <= maxtime)
  {
    for (i = 0; i < nthreads ; i++)
    {
      tries ++;
      if (verbose >= 1)
        gmp_printf ("# %d ad=%Zd\n", tries, admin);
      mpz_set (T[i]->ad, admin);
      pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
      mpz_add_ui (admin, admin, incr);
    }
    for (i = 0 ; i < nthreads ; i++)
      pthread_join(tid[i], NULL);
  }

  printf ("# Stat: tried %d ad-value(s), found %lu polynomial(s)\n", tries,
          found);

  /* print total time (this gets parsed by the scripts) */
  printf ("# Stat: total phase took %.2fs\n", seconds () - st0);

  if (best_E == DBL_MAX)
    /* This line is required by the script: */
    printf ("# No polynomial found, please increase the ad range or decrease P\n");

  for (i = 0; i < nthreads ; i++)
    {
      mpz_clear (T[i]->N);
      mpz_clear (T[i]->ad);
    }
  free (T);
  mpz_clear (N);
  mpz_clear (admin);
  mpz_clear (admax);
  mpz_clear (maxS);
  clearPrimes (&Primes);
  cado_poly_clear (best_poly);
  cado_poly_clear (curr_poly);
  param_list_clear (pl);
  free (tid);

  return 0;
}
