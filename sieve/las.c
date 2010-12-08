#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h> /* for PRIx64 macro and strtoumax */
#include <math.h>   // for ceiling, floor in cfrac
#include "fb.h"
#include "utils.h"           /* lots of stuff */
#include "basicnt.h"         /* ctzl bin_gcd */
#include "ecm/facul.h"
#include "bucket.h"
#include "trialdiv.h"
#include <pthread.h>

#ifdef SSE_NORM_INIT
#include <emmintrin.h>
#endif

#ifndef HAVE_LOG2
static double
log2 (double x)
{
  return log (x) / log (2.0);
}
#endif

/* As its name says, this is a ugly hack that initializes all lognorms to the
   maximal value (255) on the rational side. But it seems to work well, and to
   miss only about 7% to 8% relations wrt a more accurate estimation. */
//#define UGLY_HACK

/* number of bits used to estimate the norms */
#define NORM_BITS 10

/* define PROFILE to keep certain function from being inlined, in order to
   make them show up on profiler output */
//#define PROFILE

/* Trick to discard lognorms that will probably lead to non L-smooth
   cofactors. Disabled for now since it requires accurate sieving
   (in particular by prime powers and bad primes). */
// #define COFACTOR_TRICK

/* default sieve region side is 2^DEFAULT_I */
#define DEFAULT_I 12

/* default bucket region: 2^16 = 64K == close to L1 size, but this is the
   (current) largest possible value, otherwise bucket.h must be changed,
   since it stores positions on 16 bits */
#ifndef LOG_BUCKET_REGION
#define LOG_BUCKET_REGION 16
#endif

#if LOG_BUCKET_REGION > 16
#error "Too large LOG_BUCKET_REGION, please adapt bucket.h first"
#endif

/* Define SKIP_GCD3 to skip updates where 3|gcd(i,j) in the
   bucket sieving phase. Slightly slower than not skipping them
   in single-thread mode, but might be useful for multi-threading,
   or when memory is tight */
// #define SKIP_GCD3

/* These parameters control the size of the buckets. 
 * The number of updates that a bucket can accumulate is estimated as
 *   (loglog(factor base bound) - loglog(bucket sieving threshold)) 
 *     * BUCKET_LIMIT_FACTOR * I * J + BUCKET_LIMIT_ADD 
 * We don't store updates where 2|gcd(i,j) which reduces the number 
 * of updates by about 1/4, so 0.8 should be safe.
 * If we don't store updates where 3>gcd(i,j) either, their numher
 * is reduced by another factor 0f 0.11, then 0.7 should be ok
 */

#ifndef BUCKET_LIMIT_FACTOR
#ifdef SKIP_GCD3
#define BUCKET_LIMIT_FACTOR 0.7
#else
#define BUCKET_LIMIT_FACTOR 0.8
#endif
#endif

#ifndef BUCKET_LIMIT_ADD
#define BUCKET_LIMIT_ADD 1024
#endif

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */

/* Guard for the logarithms of norms, so that the value does not wrap around
   zero due to roundoff errors. */
#define GUARD 4

/* GUARD+LOG_MAX should be as near as possible from 256, to get more accuracy
   in the norm computations, but not too much, otherwise a norm might be
   rounded to zero. */
#define LOG_MAX (255.9 - (double) GUARD)

/* define TRACE_K to -I/2+i+I*j to trace the sieve array entry corresponding
   to (i,j), i.e., a = a0*i+a1*j, b = b0*i+b1*j */
// #define TRACE_K 5418470

/* Define CHECK_UNDERFLOW to check for underflow when subtracting
   the rounded log(p) from sieve array locations */
//#define CHECK_UNDERFLOW

/* Some functions should not be inlined when we profile or it's hard or
   impossible to tell them apart from the rest in the profiler output */
#ifdef PROFILE
#define NOPROFILE_INLINE
#define NOPROFILE_STATIC
#else
#define NOPROFILE_INLINE static inline
#define NOPROFILE_STATIC static
#endif

/* uintmax_t is guaranteed to be larger or equal to uint64_t */
#define strtouint64(nptr,endptr,base) (uint64_t) strtoumax(nptr,endptr,base)



/* This global mutex should be locked in multithreaded parts when a
 * thread does a read / write, especially on stdout, stderr...
 */
pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER; 


// General information about the siever
typedef struct {
    int bench;
    // multithreading info
    int nb_threads;
    // sieving area
    uint32_t I;
    uint32_t J;
    int logI; // such that I = 1<<logI
    // description of the q-lattice
    int ratq;   // 0 means special q on alg side, otherwise on rat side
    uint64_t q;
    uint64_t rho;
    int32_t a0, b0, a1, b1;
    // parameters for bucket sieving
    int bucket_thresh;    // bucket sieve primes >= bucket_thresh
    int bucket_region;    // should be around L1 cache size, a power of 2
                          // and a multiple of I.
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
    unsigned int degree;   /* polynomial degree */
    double scale_alg;      /* norm scale used on the algebraic side */
    double scale_rat;      /* norm scale used on the rational side */
    double logmax_alg;     /* norms on the alg. side are < 2^logmax_alg */
    double logmax_rat;     /* norms on the rat. side are < 2^logmax_rat */
    mpz_t *fij, gij[2];    /* coefficients of F,G(a0*i+a1*j, b0*i+b1*j)  */
    double *fijd;     /* coefficients of F_q/q */
    double B;         /* bound for the norm computation */
    unsigned char alg_Bound[256]; /* zero for good lognorms, 127 otherwise */
    unsigned char rat_Bound[256]; /* zero for good lognorms, 127 otherwise */
    fbprime_t *trialdiv_primes_alg;
    fbprime_t *trialdiv_primes_rat;
    trialdiv_divisor_t *trialdiv_data_alg;
    trialdiv_divisor_t *trialdiv_data_rat;
    unsigned char S_rat[1 << NORM_BITS];
    unsigned char S_alg[1 << NORM_BITS];
    unsigned int *lpf; /* lpf[i] is largest prime factor of i, for i < I */
    facul_strategy_t *strategy;
} sieve_info_t;

/* Forward declarations of functions */
MAYBE_UNUSED 
static void NxToIJ(int *, unsigned int *, const unsigned int, 
                   const unsigned int, const sieve_info_t *);
MAYBE_UNUSED 
static void IJToAB(int64_t *, uint64_t *, const int, const unsigned int, 
                   const sieve_info_t *);
static void xToAB(int64_t *a, uint64_t *b, const unsigned int x,
                  const sieve_info_t *si);
/* Test if entry x in bucket region n is divisible by p */
void test_divisible_x (const fbprime_t p, const unsigned long x, const int n,
		       const sieve_info_t *si, const char side);
int factor_leftover_norm (mpz_t n, unsigned int b, mpz_array_t* const factors,
			  uint32_array_t* const multis,
			  facul_strategy_t *strategy);
NOPROFILE_STATIC void
eval_fij (mpz_t v, const mpz_t *f, const unsigned int d, const long i,
	  const unsigned long j);


/************************** sieve info stuff *********************************/

static double get_maxnorm (cado_poly, sieve_info_t *, uint64_t);

/* initialize array C[0..255]: C[i] is non-zero whenever the log-norm i
   is considered as a potential report.
   The large prime bound is L = 2^l.
*/
static void
sieve_info_init_lognorm (unsigned char *C, unsigned char threshold,
                         unsigned long B MAYBE_UNUSED,
                         unsigned long l MAYBE_UNUSED,
                         double scale MAYBE_UNUSED)
{
  unsigned long k;

  for (k = 0; k < 256; k++)
    C[k] = (k <= threshold) ? 0 : 127;

#ifdef COFACTOR_TRICK
  {
    unsigned char k0, k1;
    /* for L < R <= B^2, R cannot be L-smooth (see lattice.tex) */
    k0 = (unsigned char) ((double) l * scale + 0.5) + GUARD;
    k1 = (unsigned char) (2.0 * log2 ((double) B) * scale + 0.5)
      + GUARD;
    for (k = k0 + GUARD; k + GUARD <= k1 && k < 256; k++)
      C[k] = 0;
  }
#endif
}

static void
sieve_info_init (sieve_info_t *si, cado_poly cpoly, int logI, uint64_t q0,
		 const int bucket_thresh, FILE *output, int nb_threads)
{
  unsigned int d = cpoly->degree;
  unsigned int k;
  double r, scale;
  unsigned char alg_bound, rat_bound;

  si->nb_threads = nb_threads;
  si->degree = d;
  si->fij = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  si->fijd = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    mpz_init (si->fij[k]);
  mpz_init (si->gij[0]);
  mpz_init (si->gij[1]);
  si->logI = logI;
  si->I = 1 << si->logI;
  si->J = 1 << (si->logI - 1);

  /* initialize bounds for the norm computation, see lattice.tex */
  si->B = sqrt (2.0 * (double) q0 / (cpoly->skew * sqrt (3.0)));
  si->logmax_alg = get_maxnorm (cpoly, si, q0); /* log2(max norm) */

  /* We want some margin, (see below), so that we can set 255 to discard
   * non-survivors.*/
  scale = si->logmax_alg + cpoly->alambda * (double) cpoly->lpba;

  fprintf (output, "# Alg. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           scale, exp2 (scale / LOG_MAX));
  // second guard, due to the 255 trick!
  scale = (LOG_MAX - GUARD) / scale;
  alg_bound = (unsigned char) (cpoly->alambda * (double) cpoly->lpba *  scale)
            + GUARD;
  fprintf (output, " bound=%u\n", alg_bound);
  sieve_info_init_lognorm (si->alg_Bound, alg_bound, cpoly->alim, cpoly->lpba,
                           scale);
  si->scale_alg = scale;

  /* similar bound on the rational size: |a| <= s*I*B and |b| <= I*B */
  scale = fabs (mpz_get_d (cpoly->g[1])) * cpoly->skew
        + fabs (mpz_get_d (cpoly->g[0]));
  scale *= si->B * (double) si->I;
  if (si->ratq)
      scale /= (double) q0;
  si->logmax_rat = scale = log2 (scale);
  /* on the rational side, we want that the non-reports on the algebraic
     side, which are set to 255, remain over the report bound R, even if
     the rational norm is totally smooth. For this, we simply add R to the
     maximal lognorm to compute the log base */
  r = cpoly->rlambda * (double) cpoly->lpbr; /* base-2 logarithm of the
                                                report bound */
  fprintf (output, "# Rat. side: log2(maxnorm)=%1.2f ", scale);
  fprintf (output, "logbase=%1.6f", exp2 (scale / LOG_MAX ));
  /* we subtract again GUARD to avoid that non-reports overlap the report
     region due to roundoff errors */
  si->scale_rat = LOG_MAX / scale;
  rat_bound = (unsigned char) (r * si->scale_rat) + GUARD;
  fprintf (output, " bound=%u\n", rat_bound);
  sieve_info_init_lognorm (si->rat_Bound, rat_bound, cpoly->rlim, cpoly->lpbr,
                           si->scale_rat);

  // bucket info
  // TODO: be more clever, here.
  si->bucket_thresh = si->I; /* Default value */
  if (si->bucket_thresh < bucket_thresh)
    si->bucket_thresh = bucket_thresh;
  ASSERT_ALWAYS(LOG_BUCKET_REGION > si->logI); /* so that we sieve an even
                                                      number of rows */
  si->bucket_region = 1<<LOG_BUCKET_REGION;
  si->nb_buckets = 1 + (si->I * si->J - 1) / si->bucket_region;

  double limit_factor = log(log(MAX(cpoly->rlim, cpoly->alim))) - 
                        log(log(si->bucket_thresh));
  si->bucket_limit = limit_factor * si->bucket_region * BUCKET_LIMIT_FACTOR 
                     + BUCKET_LIMIT_ADD; 

  /* Store largest prime factor of k in si->lpf[k], 0 for k=0, 1 for k=1 */
  si->lpf = (unsigned int *) malloc (si->I * sizeof (unsigned int));
  ASSERT_ALWAYS (si->lpf != NULL);
  si->lpf[0] = 0U;
  si->lpf[1] = 1U;
  for (k = 2U; k < si->I; k++)
    {
      unsigned int p, c = k;
      for (p = 2U; p * p <= c; p += 1U + p % 2U)
        {
          while (c % p == 0U)
            c /= p;
          if (c == 1U)
            break;
        }
      si->lpf[k] = (c == 1U) ? p : c;
    }

  fprintf(output, "# bucket_region = %u\n", si->bucket_region);
  fprintf(output, "# nb_buckets = %u\n", si->nb_buckets);
  fprintf(output, "# bucket_limit = %u\n", si->bucket_limit);
}


/* Finds prime factors p < lim of n and returns a pointer to a zero-terminated
   list of those factors. Repeated factors are stored only once. */
static fbprime_t *
factor_small (mpz_t n, fbprime_t lim)
{
  unsigned long p;
  unsigned long l; /* number of prime factors */
  fbprime_t *f;

  l = 0;
  f = (fbprime_t*) malloc (sizeof (fbprime_t));
  for (p = 2; p <= lim; p = getprime (p))
    {
      if (mpz_divisible_ui_p (n, p))
        {
          l ++;
          f = (fbprime_t*) realloc (f, (l + 1) * sizeof (fbprime_t));
          f[l - 1] = p;
        }
    }
  f[l] = 0; /* end of list marker */
  getprime (0);
  return f;
}

static void
sieve_info_update (sieve_info_t *si, const int verbose, FILE *output)
{
  unsigned int k;
  double invq;
  
  if (verbose)
    fprintf (output, "# I=%u; J=%u\n", si->I, si->J);

  /* update number of buckets */
  si->nb_buckets = 1 + (si->I * si->J - 1) / si->bucket_region;
  
  /* Update floating point version of algebraic poly */
  if (!si->ratq)
      invq = 1.0 / (double) si->q;
  else 
      invq = 1.0;
  for (k = 0; k <= si->degree; k++)
    si->fijd[k] = mpz_get_d (si->fij[k]) * invq;
}

static void
sieve_info_clear (sieve_info_t *si)
{
  unsigned int d = si->degree;
  unsigned int k;

  for (k = 0; k <= d; k++)
    mpz_clear (si->fij[k]);
  free (si->fij);
  mpz_clear (si->gij[0]);
  mpz_clear (si->gij[1]);
  free (si->fijd);
  free (si->lpf);
}

/*****************************************************************************/

// Redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32. FIXME: p < 2^31? (see below)
//   * invp is -1/p mod 2^32.
//   * x is some integer in [0, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = x*invp;                            /* t = x * invp mod 2^32 */
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32; /* x + t*p is bounded
            by 2^32*p-1+(2^32-1)*p < 2*2^32*p: we might want p < 2^31 so that
            there is no overflow */
    if (y >= p)
        y -= p;
    ASSERT(y<p);
    return y;
}

// Signed redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some signed integer in ]-2^32*p, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = ((uint64_t)x)*invp;
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32;
    // might be too large by p, or too small by p.
    if ((int32_t)y < 0)
        y += p;
    else if (y >= p)
        y -= p;
#ifndef NDEBUG
    if (y >= p) {
        fprintf(stderr, "BUG in redc_32. x = %" PRId64
                " p = %u, invp = %u, y = %d\n", x, p, invp, y);
        ASSERT(0);
    }
#endif

    return y;
}

NOPROFILE_INLINE int
invmod (uint64_t *pa, uint64_t b)
{
  modulusul_t m;
  residueul_t r;
  int rc;
  modul_initmod_ul (m, b);
  modul_init (r, m);
  modul_set_ul (r, *pa, m); /* With mod reduction */
  if ((rc = modul_inv(r, r, m)))
    *pa = modul_get_ul (r, m);
  modul_clear (r, m);
  modul_clearmod (m);
  return rc;
}


/* Define to 0 or 1. Table lookup seems to be slightly faster than
   ctzl() on Opteron and slightly slower on Core2, but doesn't really
   make much difference in either case. */
#define LOOKUP_TRAILING_ZEROS 1

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// a is modified in place
// return 1 on succes, 0 on failure
NOPROFILE_INLINE int
invmod_REDC(uint64_t *pa, uint64_t b) {
#if LOOKUP_TRAILING_ZEROS
  static const unsigned char trailing_zeros[256] =
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
#endif
  uint64_t a, u, v, fix, p = b;
  int t, lsh;
  a = *pa;

  if (UNLIKELY(*pa == 0))
    return 0; /* or we get infinite loop */

  if (UNLIKELY(b % 2UL == 0))
    return invmod(pa, b);

  fix = (b+1)>>1;

  ASSERT (a < b);

  u = 1; v = 0; t = 0;

  // make a odd
  lsh = ctzl(a);
  a >>= lsh;
  t += lsh;
  /* v <<= lsh; ??? v is 0 here */

  // Here a and b are odd, and a < b
  do {
    do {
      b -= a; v += u;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b];
	b >>= lsh;
	t += lsh;
	u <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(b);
      b >>= lsh;
      t += lsh;
      u <<= lsh;
#endif
    } while (a<b);
    if (UNLIKELY(a == b))
      break;
    do {
      a -= b; u += v;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a];
	a >>= lsh;
	t += lsh;
	v <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
#endif
    } while (b < a);
  } while (a != b);
  if (a != 1)
    return 0;

  // Here, the inverse of a is u/2^t mod b.
  while (t>32) {
    unsigned long sig = u & 1UL;
    u >>= 1;
    if (sig)
      u += fix;
    --t;
  }
  while (t < 32)
    {
      u <<= 1;
      if (u >= p)
        u -= p;
      t ++;
    }
  *pa = u;
  return 1;
}

MAYBE_UNUSED
static inline fbprime_t
invmod_po2 (fbprime_t n)
{
  fbprime_t r;

  ASSERT (n % 2 != 0);

  r = (3 * n) ^ 2;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  if (sizeof (fbprime_t) > 4)
    r = 2 * r - r * r * n;

  return r;
}

/* Return if a is divisible by 3 */
static inline int 
is_divisible_3_u32 (uint32_t a)
{
  return (a * (- (UINT32_MAX / 3)) <= UINT32_MAX / 3);
}

// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-R*b0) mod p
// Assumes p < 2^32

/* General version of the lattice transform function. Allows projective
   roots in input and output, and handles prime powers.
   In input, if the root is projective, say s/t (mod p) and t is
   non-invertible (mod p), then we expect R = p + (t/s mod p).
   On output, if the root is projective, say u/v (mod p) and v is
   non-invertible (mod p), then return value r = p + (v/u mod p).
   So projective roots are stored as their reciprocal, and have p added
   to signal the fact that it's a reciprocal value.
*/

NOPROFILE_INLINE fbprime_t
fb_root_in_qlattice (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, const sieve_info_t * si)
{
  int64_t aux1, aux2;
  uint64_t u, v;
  fbprime_t add;

    /* Handle powers of 2 separately, REDC doesn't like them */
    if (UNLIKELY(p % 2 == 0))
      {
	fbprime_t u, v;
	ASSERT(p == (p & -p)); /* Test that p is power of 2 */
	if (R < p) /* Root in a,b-plane is non-projective */
	  {
	    u = R * si->b1 - si->a1;
	    v = si->a0 - R * si->b0;
	  }
	else /* Root in a,b-plane is projective */
	  {
	    u = si->b1 - (R - p) * si->a1;
	    v = (R - p) * si->a0 - si->b0;
	  }

	if (v % 2 != 0)
	  {
	    /* root in i,j-plane is non-projective */
	    v = invmod_po2 (v);
	    return (u * v) & (p - 1);
	  }
	else
	  {
	    /* root in i,j-plane is projective */
	    u = invmod_po2 (u);
	    return p + ((u * v) & (p - 1));
	  }
      }

    // Use Signed Redc for the computation:
    // Numerator and denominator will get divided by 2^32, but this does
    // not matter, since we take their quotient.

    if (LIKELY(R < p)) /* Root in a,b-plane is affine */
      {
	aux1 = ((int64_t)R)*((int64_t)si->b1) - ((int64_t)si->a1);
	aux2 = ((int64_t)si->a0) - ((int64_t)R)*((int64_t)si->b0);
      }
    else /* Root in a,b-plane is projective */
      {
	aux1 = ((int64_t)si->b1) - ((int64_t)(R - p))*((int64_t)si->a1);
	aux2 = ((int64_t)(R - p))*((int64_t)si->a0) - ((int64_t)si->b0);
      }
    u = redc_32(aux1, p, invp); /* 0 <= u < p */
    v = redc_32(aux2, p, invp); /* 0 <= den < p */

    add = 0;
    if (UNLIKELY(!invmod_REDC(&v, p)))
      {
	/* root in i,j-plane is projective */
	if (UNLIKELY(!invmod_REDC(&u, p)))
          {
            fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
            exit (1); /* Should never happen! */
          }
	add = p;
      }
    u *= v;
    return (fbprime_t) redc_u32 (u, p, invp) + add;
}



/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 */

typedef struct {
    int32_t alpha, beta, gamma, delta;  // coordinates of the basis
    uint32_t a, c;       // sieving offsets in x-coordinate
    uint32_t b0, b1;     // thresholds for the branch inside siever
} plattice_info_t;

// Proposition 1 of [FrKl05]:
// Compute a basis <(alpha, beta), (gamma, delta)> of the p-lattice
// inside the q-lattice, such that
//    beta, delta > 0
//    -I < alpha <= 0 <= gamma < I
//    gamma-alpha >= I
//
// Sizes:
//    p is less than 32 bits and I fits easily in 32 bits.
//    So, alpha and beta fit easily in 32 bits, since they are less than I
//    Now, gamma and delta are also bounded by p, so 32 bits is enough
//    However: a and c can be as large as p*I (not both ?).
//    We still store them in 32 bits, since if they are larger, it means
//    that as soon as they are added to the offset for S, the index will
//    be out of range for S and the loop stops. Hence, this is safe to
//    replace a and c by a large value within 32 bits, when they are
//    larger than 32 bits.
//
// Return value:
// * non-zero if everything worked ok
// * zero when the algorithm failed. This can happen when p is a prime power,
//   and g, gcd(p,r) >= I, since then the subtractive Euclidean algorithm will
//   yield (a0=g, b0=0) at some point --- or the converse --- and the loop
//   while (|a0| >= I) a0 += b0 will loop forever.
//
NOPROFILE_INLINE int
reduce_plattice (plattice_info_t *pli, const fbprime_t p, const fbprime_t r,
                 const sieve_info_t * si)
{
    int32_t a0, a1, b0, b1, I, IJ, k;
    int logI;

    I = si->I;
    a0 = -((int32_t)p); a1 = 0;
    b0 = r;  b1 = 1;

    /* subtractive variant of Euclid's algorithm */
    while (b0 >= I)
    {
      /* a0 < 0 < b0 < -a0 */
        do {
            a0 += b0;
            a1 += b1;
        } while (a0 + b0 <= 0);
        /* -b0 < a0 <= 0 < b0 */
        if (-a0 < I)
          {
            if (UNLIKELY(a0 == 0))
              return 0;
            /* Now that |a0| < I, we switch to classical division, since
               if say |a0|=1 and b0 is large, the subtractive variant
               will be very expensive.
               We want b0 + k*a0 < I, i.e., b0 - I + 1 <= k*(-a0),
               i.e., k = ceil((b0-I+1)/a0). */
            k = 1 + (b0 - I) / (-a0);
            b0 += k * a0;
            b1 += k * a1;
            goto case_even;
          }
        /* -b0 < a0 < 0 < b0 */
        do {
            b0 += a0;
            b1 += a1;
        } while (b0 + a0 >= 0);
        /* a0 < 0 <= b0 < -a0 */
    }
    if (UNLIKELY(b0 == 0))
      return 0;
    /* we switch to the classical algorithm here too */
    k = 1 + (-a0 - I) / b0;
    a0 += k * b0;
    a1 += k * b1;

 case_even:
    pli->alpha = a0;
    pli->beta =  a1;
    pli->gamma = b0;
    pli->delta = b1;

    ASSERT (pli->beta > 0);
    ASSERT (pli->delta > 0);
    ASSERT ((pli->alpha <= 0) && (pli->alpha > -I));
    ASSERT ((pli->gamma >= 0) && (pli->gamma < I));
    ASSERT (pli->gamma - pli->alpha >= I);

    // WARNING: Here, we assume a lot on a bound on I,J
    // TODO: clean these bound problems
    logI = si->logI;
    IJ = si->J << logI;
    int64_t aa = ((int64_t)pli->beta << logI) + (int64_t)(pli->alpha);
    if (aa > IJ)
        pli->a = (uint32_t)(INT32_MAX/2);
    else
        pli->a = (uint32_t)aa;
    int64_t cc = ((int64_t)pli->delta << logI) + (int64_t)(pli->gamma);
    if (cc > IJ)
        pli->c = (uint32_t)(INT32_MAX/2);
    else
        pli->c = (uint32_t)cc;
    pli->b0 = -pli->alpha;
    pli->b1 = I - pli->gamma;

    return 1; /* algorithm was ok */
}

/***************************************************************************/

/********        Main bucket sieving functions                    **********/

void
fill_in_buckets(bucket_array_t *BA_param, factorbase_degn_t *fb,
        fbprime_t pmax, const sieve_info_t * si) 
{
    bucket_array_t BA = *BA_param; /* local copy */
    // Loop over all primes in the factor base >= bucket_thresh
    // and up to FB_END or pmax (pmax included)

    /* Skip forward to first FB prime p >= si->bucket_thresh */
    while (fb->p != FB_END && fb->p < (fbprime_t) si->bucket_thresh)
        fb = fb_next (fb); // cannot do fb++, due to variable size !

    /* Init first logp value */
    if (fb->p != FB_END && fb->p <= pmax)
      bucket_new_logp (&BA, fb->plog);

    while (fb->p != FB_END && fb->p <= pmax) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;
        ASSERT_ALWAYS (p % 2 == 1);

        /* Write new set of pointers if the logp value changed */
        bucket_new_logp (&BA, logp);

        /* If we sieve for special-q's smaller than the algebraic factor
           base bound, the prime p might equal the special-q prime q.
           Note that usually, this doesn't happen on the rational side, since
           the prime q cannot divide both sides, unless q divides Res(f,g). */
        if (UNLIKELY(p == si->q))
          goto next_fb;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            const uint32_t I = si->I;
            const int logI = si->logI;
            const uint32_t even_mask = (1U << logI) | 1U;
            const uint32_t maskI = I-1;
            const uint32_t maskbucket = (1U << LOG_BUCKET_REGION) - 1;
            const int shiftbucket = LOG_BUCKET_REGION;
            const uint32_t IJ = si->I * si->J;
            fbprime_t r, R;

            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, fb->invp, si);
            // TODO: should be line sieved in the non-bucket phase?
            // Or should we have a bucket line siever?
            if (UNLIKELY(r == 0))
              {
                /* If r == 0 (mod p), this prime hits for i == 0 (mod p),
                   but since p > I, this implies i = 0 or i > I. We don't
                   sieve i > I. Since gcd(i,j) | gcd(a,b), for i = 0 we
                   only need to sieve j = 1 */
                /* x = j*I + (i + I/2) = I + I/2 */
                bucket_update_t update;
                uint32_t x = I + I / 2;
                update.x = (uint16_t) (x & maskbucket);
                update.p = bucket_encode_prime (p);
                push_bucket_update(BA, x >> shiftbucket, update);
                continue;
              }
            if (UNLIKELY(r == p))
              {
                /* r == p means root at infinity, which hits for
                   j == 0 (mod p). Since q > I > J, this implies j = 0
                   or j > J. This means we sieve only (i,j) = (1,0) here.
                   Since I < bucket_region, this always goes in bucket 0.
                   FIXME: what about (-1,0)? It's the same (a,b) as (1,0)
                   but which of these two (if any) do we sieve? */
                bucket_update_t update;
                update.x = (uint16_t) I / 2 + 1;
                update.p = bucket_encode_prime (p);
                push_bucket_update(BA, 0, update);
                continue;
              }
            if (UNLIKELY(r > p))
              {
                continue;
              }

            plattice_info_t pli;
            if (reduce_plattice(&pli, p, r, si) == 0)
              {
                  pthread_mutex_lock(&io_mutex);
                  fprintf (stderr, "# fill_in_buckets: reduce_plattice() "
                          "returned 0 for p = " FBPRIME_FORMAT ", r = "
                          FBPRIME_FORMAT "\n", p, r);
                  pthread_mutex_unlock(&io_mutex);
                  continue; /* Simply don't consider that (p,r) for now.
                             FIXME: can we find the locations to sieve? */
              }

            // Start sieving from (0,0) which is I/2 in x-coordinate
            uint32_t x;
            x = (I>>1);
            // Skip (0,0), since this can not be a valid report.
            {
                uint32_t i = x;
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            // TODO: check the generated assembly, in particular, the
            // push function should be reduced to a very simple step.
            bucket_update_t update;
            update.p = bucket_encode_prime (p);
            __asm__("## Inner bucket sieving loop starts here!!!\n");
             while (x < IJ) {
                uint32_t i;
                i = x & maskI;   // x mod I
                /* if both i = x % I and j = x / I are even, then
                   both a, b are even, thus we can't have a valid relation */
                /* i-coordinate = (x % I) - I/2
                   (I/2) % 3 == (-I) % 3, hence
                   3|i-coordinate iff (x%I+I) % 3 == 0 */
                if ((x & even_mask) 
#ifdef SKIP_GCD3
                    && (!is_divisible_3_u32 (i + I) ||
                        !is_divisible_3_u32 (x >> logI))
#endif
                   )
                  {
#if defined(__x86_64__) && defined(__GNUC__)
                    /* The x value in update can be set by a write to 
                       the low word of the register, but gcc does not 
                       do so - it writes the word to memory, then reads 
                       the dword back again. */
                    __asm__ (
                      "movw %1, %w0\n\t"
                      : "+r" (update)
                      : "r" ((uint16_t) (x & maskbucket))
                    );
#else
                    update.x = (uint16_t) (x & maskbucket);
#endif
#ifdef PROFILE
		    /* To make it visible in profiler */
		    *(BA.bucket_write[x >> shiftbucket])++ = update;
#else
		    push_bucket_update(BA, x >> shiftbucket, update);
#endif
                  }
#ifdef TRACE_K
                if (x == TRACE_K)
                  fprintf (stderr, "Pushed (%u, %u) (%u) to BA[%u]\n",
                           x & maskbucket, logp, p, x >> shiftbucket);
#endif
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            __asm__("## Inner bucket sieving loop stops here!!!\n");
        }
        int i;
next_fb:
        for (i = 0; i < si->nb_threads; ++i) {
            fb = fb_next (fb);
            if (fb->p == FB_END || fb->p > pmax) {
                break;
            }
        }
    }
  *BA_param = BA; /* Write back BA so the nr_logp etc get copied to caller */
}


typedef struct {
    bucket_array_t *alg_BA;
    factorbase_degn_t *fb_alg;
    fbprime_t alg_pmax;
    bucket_array_t *rat_BA;
    factorbase_degn_t *fb_rat;
    fbprime_t rat_pmax;
    sieve_info_t *si;
} fill_in_buckets_mt_arg_t;

void *
call_fill_in_buckets(fill_in_buckets_mt_arg_t *arg) {
    fill_in_buckets(arg->alg_BA, arg->fb_alg, arg->alg_pmax, arg->si);
    fill_in_buckets(arg->rat_BA, arg->fb_rat, arg->rat_pmax, arg->si);
    return NULL;
}

void
fill_in_buckets_mt(bucket_array_t *alg_BA, factorbase_degn_t **fb_alg,
        fbprime_t *alg_pmax,
        bucket_array_t *rat_BA, factorbase_degn_t **fb_rat,
        fbprime_t *rat_pmax,
        sieve_info_t * si){
    int i;
    pthread_t *th;
    fill_in_buckets_mt_arg_t *my_arg;
    th = malloc(si->nb_threads*sizeof(pthread_t)); 
    ASSERT_ALWAYS(th);
    my_arg = malloc(si->nb_threads*sizeof(fill_in_buckets_mt_arg_t)); 
    ASSERT_ALWAYS(my_arg);
    for (i = 0; i < si->nb_threads; ++i) {
        int ret;
        my_arg[i].alg_BA = &(alg_BA[i]);
        my_arg[i].fb_alg = fb_alg[i];
        my_arg[i].alg_pmax = alg_pmax[i];
        my_arg[i].rat_BA = &(rat_BA[i]);
        my_arg[i].fb_rat = fb_rat[i];
        my_arg[i].rat_pmax = rat_pmax[i];
        my_arg[i].si = si;
        ret = pthread_create(&(th[i]), NULL, 
		(void * (*)(void *)) &call_fill_in_buckets,
                (void *)(&(my_arg[i])));
        ASSERT_ALWAYS(!ret);
    }
    for (i = 0; i < si->nb_threads; ++i) {
        int ret;
        ret = pthread_join(th[i], NULL);
        ASSERT_ALWAYS(!ret);
    }
    free(my_arg);
    free(th);
}


/* Decrease the sieve array entry *S by logp, with underflow checking 
   and tracing if desired. Variables x, bucket_nr, p, si, and caller 
   are used only for trace test and output */

static inline void
sieve_decrease (unsigned char *S, const unsigned char logp, 
                MAYBE_UNUSED const uint16_t x, 
                MAYBE_UNUSED const unsigned int bucket_nr,
                MAYBE_UNUSED const fbprime_t p,
                MAYBE_UNUSED const sieve_info_t *si, 
                MAYBE_UNUSED const int caller)
{
#ifdef TRACE_K
  if ((x + bucket_nr * si->bucket_region) == TRACE_K)
    fprintf(stderr, "# (%d) Subtract log(" FBPRIME_FORMAT ") = %u from "
            "S[%u] = %hhu, from BA[%u], ", 
            caller, p, logp, TRACE_K, *S, bucket_nr);
#endif

#ifdef CHECK_UNDERFLOW
  if (*S < logp)
    {
      int i;
      unsigned int j;
      int64_t a;
      uint64_t b;
      NxToIJ (&i, &j, bucket_nr, x, si);
      IJToAB (&a, &b, i, j, si);
      if (i % 2 == 0 && j % 2 == 0)
        fprintf (stderr, "# Error, underflow (%d) at (N,x)=(%u, %u), "
                 "(i,j)=(%d, %u), both even!\n", caller, bucket_nr, x, i, j);
        
      else
        fprintf (stderr, "# Error, underflow (%d) at (N,x)=(%u, %u), "
                 "(i,j)=(%d, %u), (a,b)=(%ld, %lu), S[x] = %hhu, log(" 
                 FBPRIME_FORMAT ") = %hhu\n", 
                 caller, bucket_nr, x, i, j, a, b, S, p, logp);
      *S = 0;
    }
  else
    *S -= logp;
#else
  *S -= logp;
#endif

#ifdef TRACE_K
  if ((x + bucket_nr * si->bucket_region) == TRACE_K)
    fprintf(stderr, "new value is %hhu\n", *S);
#endif
}


NOPROFILE_STATIC void
apply_one_bucket (unsigned char *S, bucket_array_t BA, const int i,
                     sieve_info_t * si MAYBE_UNUSED)
{
    int j = nb_of_updates(BA, i);
    int next_logp_j = 0;
    unsigned char logp = 0;
    bucket_update_t *next_logp_change, *read_ptr;

    /* Having the read_ptr here defeats the whole idea of having 
       nice inline functions to handle all the BA stuff, but I yet 
       need to figure out how to keep gcc from writing 
       BA.bucket_read[i] back to memory and reading it from memory again
       on every get_next_bucket_update()  */
    read_ptr = BA.bucket_read[i];

    /* Init so that first access fetches logp */
    next_logp_j = 0;
    next_logp_change = read_ptr;

    for (; j > 0; --j) {
       uint16_t x;

       /* Do we need a new logp ? */
       if (read_ptr >= next_logp_change)
         {
           ASSERT_ALWAYS (next_logp_j < BA.nr_logp);
           ASSERT_ALWAYS (BA.logp_idx[next_logp_j * BA.n_bucket + i] 
                           == next_logp_change);
           logp = BA.logp_val[next_logp_j++];
           /* Get pointer telling when to fetch new logp next time */
           if (next_logp_j < BA.nr_logp)
             next_logp_change = BA.logp_idx[next_logp_j * BA.n_bucket + i];
           else
             next_logp_change = BA.bucket_write[i]; /* effectively: never */
         }
       
       x = (read_ptr++)->x;
       sieve_decrease (S + x, logp, x, i, 0, si, 1);
    }
}

/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
NOPROFILE_STATIC void
init_rat_norms (sieve_info_t *si)
{
  double e, m, norm, h;
  int i, j, k, l, K, L;
  unsigned char *S = si->S_rat;

  k = (int) ceil (log2 (si->logmax_rat));
  K = 1 << k;
  l = NORM_BITS - k;
  L = 1 << l;

  /* extract k bits from the exponent, and l bits from the mantissa */
  h = 1.0 / (double) L;
  for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
    {
      /* e = 2^i for 0 <= i < 2^k */
      for (m = 1.0, j = 0; j < L; j++, m += h)
        {
          /* m = 1 + j/2^l for 0 <= j < 2^l */
          norm = m * e;
          /* Warning: since si->logmax_rat does not usually correspond to
             a power a two, and we consider full binades here, we have to
             take care that values > si->logmax_rat do not wrap around to 0 */
          norm = log2 (norm);
          if (norm >= si->logmax_rat)
            S[(i << l) + j] = 255;
          else
            {
              norm = norm * si->scale_rat;
              S[(i << l) + j] = GUARD + (unsigned char) norm;
            }
        }
    }
}

/* idem as init_rat_norms, but for the algebraic side */
NOPROFILE_STATIC void
init_alg_norms (sieve_info_t *si)
{
  double e, m, norm, h;
  int i, j, k, l, K, L;
  unsigned char *S = si->S_alg;

  k = (int) ceil (log2 (si->logmax_alg));
  K = 1 << k;
  ASSERT_ALWAYS(NORM_BITS >= k);
  l = NORM_BITS - k;
  L = 1 << l;

  /* extract k bits from the exponent, and l bits from the mantissa */
  h = 1.0 / (double) L;
  for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
    {
      for (m = 1.0, j = 0; j < L; j++, m += h)
        {
          norm = m * e;
          norm = log2 (norm);
          if (norm >= si->logmax_alg)
            S[(i << l) + j] = 255;
          else
            {
              norm = norm * si->scale_alg;
              S[(i << l) + j] = GUARD + (unsigned char) norm;
            }
        }
    }
}

/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
NOPROFILE_STATIC void
init_rat_norms_bucket_region (unsigned char *S, int N, cado_poly cpoly,
                              sieve_info_t *si)
{
    double g1, g0, gi, gj, norm MAYBE_UNUSED, gjj;
    int i MAYBE_UNUSED, halfI MAYBE_UNUSED, l;
    unsigned int j, lastj;
    uint64_t mask = (1 << NORM_BITS) - 1;
    union { double z; uint64_t x; } zx[1];

    l = NORM_BITS - (int) ceil (log2 (si->logmax_rat));

    /* G_q(i,j) = g1 * (a0*i+a1*j) + g0 * (b0*i+b1*j)
       = (g1*a0+g0*b0) * i + (g1*a1+g0*b1) * j
       = gi * i + gj * j
    */

    g1 = mpz_get_d (cpoly->g[1]);
    g0 = mpz_get_d (cpoly->g[0]);
    gi = g1 * (double) si->a0 + g0 * (double) si->b0;
    gj = g1 * (double) si->a1 + g0 * (double) si->b1;

    if (si->ratq) {
        double invq = 1.0 / (double) si->q;
        gi *= invq;
        gj *= invq;
    }

    /* bucket_region is a multiple of I. Let's find the starting j
     * corresponding to N and the last j.
     */
    halfI = si->I / 2;
    j = (N*si->bucket_region) >> si->logI;
    lastj = j + (si->bucket_region >> si->logI);
#ifdef UGLY_HACK
    memset (S, 255, (lastj - j) * si->I);
#else
    /* if j = 0, it will be the first value */
    if (j == 0)
      {
        // compute only the norm for i = 1. Everybody else is 255.
        norm = log2 (fabs (gi));
        norm = norm * si->scale_rat;
        for (i = -halfI; i < halfI; i++)
          *S++ = UCHAR_MAX;
        S[-halfI + 1] = GUARD + (unsigned char) (norm);
        ++j;
      }
    for (; j < lastj; ++j)
      {
        unsigned char n;
        gjj = gj * (double) j;
        zx->z = gjj - gi * (double) halfI;
        __asm__("### Begin rational norm loop\n");
#ifndef SSE_NORM_INIT
        uint64_t y;
        for (i = 0; i < (int) si->I; i++) {
          /* the double precision number 1.0 has high bit 0 (sign),
             then 11-bit biased exponent 1023, and 52-bit mantissa 0 */
          /* the magic constant here is simply 1023*2^52, where
             1023 is the exponent bias in binary64 */
          y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
          n = si->S_rat[y & mask];
          ASSERT (n > 0);
          *S++ = n;
          zx->z += gi;
        }
#else
        union { __v2df dble;
            __v2di intg;
            struct {uint64_t y0; uint64_t y1; } intpair;
        } y_vec;

        __v2df gi_vec = { 2*gi, 2*gi };
        __v2di mask_vec = { mask, mask };
        __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
            (uint64_t) 0x3FF0000000000000 };

        // in spite of the appearance, only the low word gives the shift
        // count. The high word is ignored.
        __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
        __v2df z_vec = { zx->z, zx->z+gi };

        for (i = 0; i < halfI; ++i) {
            y_vec.dble = z_vec;
            y_vec.intg -= cst_vec;
            y_vec.intg = _mm_srl_epi64(y_vec.intg, shift_value);
            y_vec.intg &= mask_vec;
//            *S++ = si->S_rat[((uint32_t *)(&y_vec.intg))[0]];
//             *S++ = si->S_rat[((uint32_t *)(&y_vec.intg))[2]];
            *S++ = si->S_rat[y_vec.intpair.y0];
            *S++ = si->S_rat[y_vec.intpair.y1];
            z_vec += gi_vec;
        }
#endif
        __asm__("### End rational norm loop\n");
      }
#endif
}

#ifdef SSE_NORM_INIT
static inline void
init_fpoly_v2df(__v2df *F, const double *f, const int deg)
{
    int i;
    for (i = 0; i <= deg; ++i) {
        __v2df tmp = { f[i], f[i] };
        F[i] = tmp;
    }
}


static inline __v2df
fpoly_eval_v2df_deg5(const __v2df *f, const __v2df x)
{
    __v2df r;
    r = f[5];
    r = r * x + f[4];
    r = r * x + f[3];
    r = r * x + f[2];
    r = r * x + f[1];
    r = r * x + f[0];
    return r;
}
#endif

/* Initialize lognorms on the algebraic side for the bucket_region
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 */
NOPROFILE_STATIC int
init_alg_norms_bucket_region (unsigned char *alg_S, 
			      const unsigned char *rat_S, const int N, 
			      cado_poly cpoly, sieve_info_t *si)
{
  unsigned int j, lastj, k, d = cpoly->degree;
  double *t, *u, powj, i, halfI;
  int report = 0, l;
  uint64_t mask = (1 << NORM_BITS) - 1;
  union { double z; uint64_t x; } zx[1];
  uint64_t y;
  unsigned char *T = si->S_alg;
#ifdef TRACE_K
  unsigned char *orig_alg_S = alg_S;
#endif
  l = NORM_BITS - (int) ceil (log2 (si->logmax_alg));
  halfI = (double) (si->I >> 1);

  t = si->fijd;
  u = (double*) malloc ((si->degree + 1) * sizeof (double));
  ASSERT(u != NULL);

  /* bucket_region is a multiple of I. Let's find the starting j
   * corresponding to N and the last j.
   */
  j = (N*si->bucket_region) >> si->logI;
  lastj = j + (si->bucket_region >> si->logI);
  for (; j < lastj; j++)
    {
      /* scale by j^(d-k) the coefficients of fij */
      for (k = 0, powj = 1.0; k <= d; k++, powj *= (double) j)
        u[d - k] = t[d - k] * powj;

#ifdef SSE_NORM_INIT
      __v2df u_vec[d];
      init_fpoly_v2df(u_vec, u, d);
      double i0 = 0.0 ;
      unsigned char *S_ptr0 = NULL;
      int cpt = 0;
#endif
      /* now compute norms */
      for (i = -halfI; i < halfI; i += 1.0)
        {
#ifdef TRACE_K
	      if (N * si->bucket_region +  alg_S - orig_alg_S 
		  == TRACE_K)
		{
		  fprintf (stderr, "init_alg_norms_bucket_region: "
			   "rat_S[%d] = %u, -log(prob) = %u\n",
			   TRACE_K, 
			   (unsigned int) *rat_S, 
			   (unsigned int) si->rat_Bound[*rat_S]);
		}
#endif
          if (si->rat_Bound[*rat_S] < 127)
            {
#ifndef SSE_NORM_INIT
              unsigned char n;
              zx->z = fpoly_eval (u, d, i);
              /* 4607182418800017408 = 1023*2^52 */
              y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
              report++;
              n = T[y & mask];
              ASSERT (n > 0);
              *alg_S = n;
#else
              ASSERT(d == 5);
              __v2di mask_vec = { mask, mask };
              __v2di cst_vec = { (uint64_t) 0x3FF0000000000000,
                  (uint64_t) 0x3FF0000000000000 };
              __v2di shift_value = { (uint64_t)(52 - l), (uint64_t)(52 - l) };
              report++;
              if (cpt == 0) {
                  i0 = i;
                  cpt++;
                  S_ptr0 = alg_S;
              } else if (cpt == 1) {
                  cpt--;
                  __v2df i_vec = {i0, i};
                  union { __v2df dble;
                      __v2di intg;
                      struct {uint64_t y0; uint64_t y1; } intpair;
                  } fi_vec;
                  fi_vec.dble = fpoly_eval_v2df_deg5(u_vec, i_vec);
                  fi_vec.intg -= cst_vec;
                  fi_vec.intg = _mm_srl_epi64(fi_vec.intg, shift_value);
                  fi_vec.intg &= mask_vec;
                  *S_ptr0 = T[fi_vec.intpair.y0];
                  *alg_S = T[fi_vec.intpair.y1];
              }
#endif
            }
          else
	    {
	      *alg_S = UCHAR_MAX;
	    }
	  alg_S++;
	  rat_S++;
        }
#ifdef SSE_NORM_INIT
      if (cpt == 1) { // odd number of computations
          zx->z = fpoly_eval (u, d, i0);
          y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
          *S_ptr0 = T[y & mask];
      }
#endif
    }

  free(u);
  return report;
}

typedef struct {
    fbprime_t p;
    fbprime_t r;        // in [ 0, p [
    int next_position;  // start of the sieve for next bucket_region
    unsigned char logp;
} small_typical_prime_data_t;

/* Small primes or powers of small primes p^k with projective root.
   These hit at 
     i*v == j*u (mod p^k) 
   for some u,v in Z, but gcd(v, p^k) > 1.
   We may assume gcd(u,p)==1, or we divide the entire equation by p.
   We store g = gcd(v, p^k), q = p^k / g, and U = u * (v/g)^(-1) (mod q).

   Then we have
     i*v == j*u (mod p^k)  <==>  i == (j/g)*U (mod q)
   with g|j. 
   
   In other words, we can sieve this bad prime (powers) much like a 
   normal prime (power) q with root U, except that after sieving a line 
   we don't advance by one line, but by g lines.
   The case where g = q^k and thus q = 1 can be sieve more efficiently,
   of course, since every entry in each g-th line will be hit, so that
   the sieving should use long word transfers.

   Just like for normal primes, the next_position value points at the first
   position to sieve relative to the start of the current sieve region.

   Within a line that starts at index line_start, for array element of 
   index x, we have x - line_start = i+I/2. 
   We skip j=0, as it contains only the single possible relation 
   (i,j) = (1,0). 
   For j=1*g, we want i=U (mod q), so x - line_start == I/2+U (mod q),
   so we initialise 
     next_position = I*g + (I/2 + U) % q
   to get the first array index in line j=g, 
   then within a line sieve next_position + t*q < I, t in N,
   and update 
     next_position = (next_position - line_start + U) % q + line_start + g*I 
   to get the first position to sieve in the next suitable line. */

typedef struct {
    fbprime_t g, q, U;
    int next_position;
    unsigned char logp;
} small_bad_prime_data_t;

typedef struct {
    // primes with non-projective root
    small_typical_prime_data_t *nice_p;
    // primes with projective root
    small_bad_prime_data_t *bad_p;
    int nb_nice_p;
    int nb_bad_p;
} small_sieve_data_t;

void clear_small_sieve(small_sieve_data_t ssd) {
    free(ssd.nice_p);
    ssd.nice_p = NULL;
    if (ssd.nb_bad_p > 0)
      free(ssd.bad_p);
    ssd.bad_p = NULL;
}

void
clone_small_sieve(small_sieve_data_t *r, const small_sieve_data_t *s) {
    memcpy(r, s, sizeof(small_sieve_data_t));
    r->nice_p = malloc(r->nb_nice_p*sizeof(small_typical_prime_data_t));
    r->bad_p = malloc(r->nb_bad_p*sizeof(small_bad_prime_data_t));
    ASSERT(r->nice_p != NULL && r->bad_p != NULL);
    memcpy(r->nice_p, s->nice_p,
            r->nb_nice_p*sizeof(small_typical_prime_data_t));
    memcpy(r->bad_p, s->bad_p,
            r->nb_bad_p*sizeof(small_bad_prime_data_t));
}

/* Assume q is a prime power p^k with k>=1, return p if k > 1, 0 otherwise. */
static fbprime_t
is_prime_power (fbprime_t q)
{
  unsigned int maxk, k;
  uint32_t p;

  for (maxk = 0, p = q; p > 1; p /= 2, maxk ++);
  for (k = maxk; k >= 2; k--)
    {
      p = (fbprime_t) (pow ((double) q, 1.0 / (double) k) + 0.5);
      if (q % p == 0)
        return p;
    }
  return 0;
}

/* Copy those primes in s to r that need to be resieved, i.e., those
   that are not in trialdiv_primes and that are not prime powers */

static void
copy_small_sieve (small_sieve_data_t *r, const small_sieve_data_t *s, 
		  const fbprime_t *trialdiv_primes)
{
  int i, j, td_idx;

  r->nb_nice_p = 0;
  r->nice_p = NULL;
  if (s->nb_nice_p > 0)
    {
      const size_t size = s->nb_nice_p * sizeof (small_typical_prime_data_t);
      r->nice_p = (small_typical_prime_data_t *) malloc (size);
      ASSERT (r->nice_p != NULL);
      
      td_idx = 0;
      j = 0;
      for (i = 0; i < s->nb_nice_p; i++)
	{
	  if (is_prime_power(s->nice_p[i].p))
	    continue;
	  
	  while (trialdiv_primes[td_idx] != FB_END && 
		 trialdiv_primes[td_idx] < s->nice_p[i].p)
	    td_idx++;
	  
	  if (trialdiv_primes[td_idx] == FB_END || 
	      trialdiv_primes[td_idx] !=s->nice_p[i].p)
	    r->nice_p[j++] = s->nice_p[i];
	}
      
      r->nb_nice_p = j;
      if (j == 0)
	{
	  free (r->nice_p);
	  r->nice_p = NULL;
	}
    }
  
  r->nb_bad_p = 0;
  r->bad_p = NULL;
  if (s->nb_bad_p > 0)
    {
      const size_t size = s->nb_bad_p * sizeof (small_bad_prime_data_t);
      r->bad_p = (small_bad_prime_data_t *) malloc (size);
      ASSERT (r->bad_p != NULL);
      
      td_idx = 0;
      j = 0;
      for (i = 0; i < s->nb_bad_p; i++)
	{
	  /* p^k = q*g, g > 1, so k>1 if g is power or if q > 1 */
	  if (s->bad_p[i].q > 1 || is_prime_power(s->bad_p[i].g))
	    continue;
	  
	  while (trialdiv_primes[td_idx] != FB_END && 
		 trialdiv_primes[td_idx] < s->bad_p[i].g)
	    td_idx++;
	  
	  if ((trialdiv_primes[td_idx] == FB_END || 
	       trialdiv_primes[td_idx] !=  s->bad_p[i].g))
	    r->bad_p[j++] = s->bad_p[i];
	}
      
      r->nb_bad_p = j;
      if (j == 0)
	{
	  free (r->bad_p);
	  r->bad_p = NULL;
	}
    }
}

// Prepare sieving of small primes: initialize a small_sieve_data_t
// structure to be used thereafter during sieving each region.
// next_position points at the next position that will be hit by sieving,
// relative to the start of the next bucket region to sieve. It may exceed I 
// and even BUCKET_REGION
void init_small_sieve(small_sieve_data_t *ssd, const factorbase_degn_t *fb,
                      const sieve_info_t *si, const char side, FILE *output)
{
    const factorbase_degn_t *fb_sav = fb;
    int size = 0, n;
    const unsigned int thresh = si->bucket_thresh;
    const int verbose = 0;
    const int do_bad_primes = 1;

    ASSERT (side == 'a' || side == 'r');

    // Count prime ideals of factor base primes p < thresh
    while (fb->p != FB_END && fb->p < thresh) {
        size += fb->nr_roots;
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    fb = fb_sav;

    // allocate space for these. n is an upper bound, since some of the
    // ideals might become special ones.
    ssd->nice_p = (small_typical_prime_data_t *)
      malloc(size * sizeof(small_typical_prime_data_t));
    n = 0;
    ssd->nb_bad_p = 0;
    ssd->bad_p = NULL;
    // Do another pass on fb and badprimes, to fill in the data
    // while we have any regular primes or bad primes < thresh left
    while (fb->p != FB_END && fb->p < thresh) {
      const fbprime_t p = fb->p;
      int nr;
      fbprime_t r;

      for (nr = 0; nr < fb->nr_roots; nr++)
        {
	  r = fb_root_in_qlattice (p, fb->roots[nr], fb->invp, si);
	  /* If this root is somehow interesting (projective in (a,b) or
	     in (i,j) plane), print a message */
	  if (verbose && (fb->roots[nr] >= p || r >= p))
	    fprintf (output, "# init_small_sieve: Side %c, prime " 
		     FBPRIME_FORMAT " root " FBPRIME_FORMAT " -> " 
		     FBPRIME_FORMAT "\n", side, p, fb->roots[nr], r);

	  /* Handle projective roots */
          if (r >= p)
	    {
	      fbprime_t q, g;
	      r -= p;
	      g = gcd_ul (p, r);
	      if (do_bad_primes && g < si->J)
		{
		  // Fill in data for this bad prime
		  /* Shouldn't happen that often so a realloc is ok */
		  ssd->bad_p = (small_bad_prime_data_t *)
		    realloc (ssd->bad_p, (ssd->nb_bad_p + 1) *
			     sizeof (small_bad_prime_data_t));
		  q = p / g;
		  ssd->bad_p[ssd->nb_bad_p].g = g;
		  ssd->bad_p[ssd->nb_bad_p].q = q;
		  if (q == 1)
		    {
		      ASSERT (r == 0);
		      ssd->bad_p[ssd->nb_bad_p].U = 0;
		      ssd->bad_p[ssd->nb_bad_p].next_position = g * si->I;
		    }
		  else
		    {
		      int rc;
		      /* gcd(r / g, q) = 1 here */
		      uint64_t U = r / g;
		      rc = invmod (&U, q);
		      ASSERT_ALWAYS (rc != 0);
		      ssd->bad_p[ssd->nb_bad_p].U = U;
		      ssd->bad_p[ssd->nb_bad_p].next_position = 
			g * si->I + (si->I / 2 + U) % q;
		    }
		  ssd->bad_p[ssd->nb_bad_p].logp = fb->plog;
		  ssd->nb_bad_p++;
		} 
	      else 
		{
		  if (verbose && !do_bad_primes)
		    fprintf (output, "# init_small_sieve: not adding bad "
			     "prime " FBPRIME_FORMAT " to small sieve because "
			     "do_bad_primes = 0\n", g);
		  else if (verbose && g >= si->J)
		    fprintf (output, "# init_small_sieve: not adding bad "
			     "prime g = " FBPRIME_FORMAT " to small sieve "
			     " because  g >= si->J = %d\n", g, si->J);
		}
	    }
	  else
	    {
	      // Fill in data for this nice prime
	      ASSERT (n < size);
	      ssd->nice_p[n].p = p;
	      ssd->nice_p[n].r = r;
	      ssd->nice_p[n].logp = fb->plog;
	      ssd->nice_p[n].next_position = (si->I >> 1)%p;
	      /* For powers of 2, we sieve only odd lines and 
	         next_position needs to point at line j=1. We assume
	         that in this case (si->I/2) % p == 0 */
	      if (UNLIKELY(p % 2 == 0))
	        {
	          ASSERT (ssd->nice_p[n].next_position == 0);
	          ssd->nice_p[n].next_position = r + si->I;
                }
	      n++;
	    }
        }
      fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    ssd->nb_nice_p = n;
    if (verbose)
      fprintf (output, 
	       "# init_small_sieve: side %c has %d nice and %d bad primes\n",
	       side, ssd->nb_nice_p, ssd->nb_bad_p);
}


// Sieve small primes (up to p < bucket_thresh) of the factor base fb in the
// next sieve region S.
// Information about where we are is in ssd.
void sieve_small_bucket_region(unsigned char *S, const int bucket_nr,
			       small_sieve_data_t ssd, sieve_info_t *si,
			       const unsigned char side)
{
    const uint32_t I = si->I;
    const fbprime_t pattern2_size = 2 * sizeof(unsigned long);
    unsigned long j;
    int n;
    const int test_divisibility = 0; /* very slow, but nice for debugging */
    const unsigned long nj = si->bucket_region >> si->logI; /* Nr. of lines 
                                                        per bucket region */

    ASSERT ((nj & 1) == 0);

    /* Handle powers of 2 up to 2 * sizeof(long) separately. 
       TODO: use SSE2 */
    /* First collect updates for powers of two in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (j = 0; j < nj; j++)
      {
        unsigned long pattern[2];

        /* Prepare the pattern */
        pattern[0] = pattern[1] = 0UL;

        /* This loop assumes that entries in ssd.nice_p are in order of 
           increasing p, or more accurately, that all powers of 2 up to
           2*sizeof(long) appear before any p > 2*sizeof(long) */
        for (n = 0; n < ssd.nb_nice_p && ssd.nice_p[n].p <= pattern2_size; 
             n++)
          if (ssd.nice_p[n].p % 2 == 0) 
            {
              const fbprime_t p = ssd.nice_p[n].p;
              unsigned int i0 = ssd.nice_p[n].next_position;
              /* Sieve only odd lines */
              if (i0 < I)
                {
                  unsigned int i;
                  ASSERT (i0 < p);
                  ASSERT ((nj * bucket_nr + j) % 2 == 1);
                  for (i = i0; i < pattern2_size; i += p)
                    ((unsigned char *)&pattern)[i] += ssd.nice_p[n].logp;
                  i0 = ((i0 + 2 * ssd.nice_p[n].r) & (p - 1)) + 2 * I;
                }
              /* In this loop, next_position gets updated to the first 
                 index to sieve relative to the start of the next line, 
                 but after all lines of this bucket region are processed, 
                 it will point the the first position to sieve relative  
                 to the start of the next bucket region, as required */
              ssd.nice_p[n].next_position = i0 - I;
            }
        
        /* Apply the pattern */
        if (pattern[0] != 0UL || pattern[1] != 0UL)
          {
            unsigned long *S_ptr = (unsigned long *) (S + j * I);
            const unsigned long *end = (unsigned long *)(S + j * I + I);
            
            while (S_ptr < end)
              {
                *(S_ptr) -= pattern[0];
                *(S_ptr + 1) -= pattern[1];
                *(S_ptr + 2) -= pattern[0];
                *(S_ptr + 3) -= pattern[1];
                S_ptr += 4;
              }
          }
      }


    /* Handle 3 */
    /* First collect updates for powers of two in a pattern,
       then apply pattern to sieve line.
       Repeat for each line in bucket region. */
    for (j = 0; j < nj; j++)
      {
        unsigned long pattern[3];

        pattern[0] = pattern[1] = pattern[2] = 0UL;

        for (n = 0; n < ssd.nb_nice_p && ssd.nice_p[n].p <= 3; 
             n++)
          if (ssd.nice_p[n].p == 3) 
            {
              const fbprime_t p = 3;
              unsigned int i0 = ssd.nice_p[n].next_position;
              unsigned int i;
              ASSERT (i0 < p);
              for (i = i0; i < 3 * sizeof(unsigned long); i += p)
                ((unsigned char *)&pattern)[i] += ssd.nice_p[n].logp;
              i0 += ssd.nice_p[n].r;
              if (i0 >= p)
                i0 -= p;
              ssd.nice_p[n].next_position = i0;
            }
        
        if (pattern[0] != 0UL)
          {
            unsigned long *S_ptr = (unsigned long *) (S + j * I);
            const unsigned long *end = (unsigned long *)(S + j * I + I) - 2;
            
            while (S_ptr < end)
              {
                *(S_ptr) -= pattern[0];
                *(S_ptr + 1) -= pattern[1];
                *(S_ptr + 2) -= pattern[2];
                S_ptr += 3;
              }
            
            end += 2;
            if (S_ptr < end)
              *(S_ptr++) -= pattern[0];
            if (S_ptr < end)
              *(S_ptr) -= pattern[1];
          }
      }

    for (n = 0 ; n < ssd.nb_nice_p; ++n) {
        const fbprime_t p = ssd.nice_p[n].p;
        const fbprime_t r = ssd.nice_p[n].r;
        const unsigned char logp = ssd.nice_p[n].logp;
        unsigned char *S_ptr = S;
        fbprime_t twop;
        unsigned int i, i0, linestart = 0;
        /* Always S_ptr = S + linestart. S_ptr is used for the actual array
           updates, linestart keeps track of position relative to start of
           bucket region and is used only for computing i,j-coordinates
           in overflow and divisibility checking, and relation tracing. */

        i0 = ssd.nice_p[n].next_position;

        /* Powers of 2 are treated separately */
        if (UNLIKELY(p % 2 == 0))
          {
            /* Don't sieve powers of 2 again that were pattern-sieved */
            if (p <= pattern2_size)
              continue;
            
            for (j = 0; j < nj; j++)
              {
                if (i0 < I)
                  {
                    ASSERT(i0 < p);
                    ASSERT ((nj * bucket_nr + j) % 2 == 1);
                    for (i = i0; i < I; i += p)
                      {
                        if (test_divisibility)
                          test_divisible_x (p, linestart + i, bucket_nr, si, 
                                            side);
                        sieve_decrease (S_ptr + i, logp, 
                                        linestart + i, bucket_nr, p, si, 2);
                      }
                    i0 = ((i0 + 2 * r) & (p - 1)) + 2 * I;
                  }
                i0 -= I;
                linestart += I;
                S_ptr += I;
              }
            ssd.nice_p[n].next_position = i0;
            continue;
          }

        /* Don't sieve 3 again as it was pattern-sieved */
        if (p == 3)
          continue;
          
        ASSERT(i0 < p);
        for (j = 0; j < nj; j++)
          {
            twop = p;
            i = i0;
            if ((((nj & bucket_nr) ^ j) & 1) == 0) /* (nj*bucket_nr+j)%2 */
              {
                /* for j even, we sieve only odd i. */
                twop += p;
                i += (i0 & 1) ? 0 : p;
              }
            for ( ; i < I; i += twop)
              {
                if (test_divisibility)
                  test_divisible_x (p, linestart + i, bucket_nr, si, side);
                sieve_decrease (S_ptr + i, logp, linestart + i, bucket_nr, 
                                p, si, 3);
              }
            i0 += r;
            if (i0 >= p)
              i0 -= p;
            S_ptr += I;
            linestart += I;
          }
        ssd.nice_p[n].next_position = i0;
    }

    /* Sieve the bad primes. We have p^k | fij(i,j) for i,j such that
       i * g == j * U (mod p^k) where g = p^l and gcd(U, p) = 1. 
       This hits only for g|j, then j = j' * g, and i == j' * U (mod p^(k-l)).
       In every g-th line, we sieve the entries with i == (j/g)*U (mod q).
       In ssd we have stored g, q = p^(k-l), U, and next_position so that
       S + next_position is the next sieve entry that needs to be sieved.
       So if S + next_position is in the current bucket region, we
       update all  S + next_position + n*q  where  next_position + n*q < I,
       then set next_position = ((next_position % I) + U) % q) + I * g.  */
    for (n = 0; n < ssd.nb_bad_p; ++n) {
      if (!test_divisibility && ssd.bad_p[n].q == 1)
        {
	  /* q = 1, therefore U = 0, and we sieve all entries in lines
	     with g|j, beginning with the line starting at S[next_position] */
          unsigned long logps;
          unsigned int i0 = ssd.bad_p[n].next_position;
	  ASSERT (ssd.bad_p[n].U == 0);
          ASSERT (i0 % I == 0);
          ASSERT (I % (4 * sizeof (unsigned long)) == 0);
	  for (j = 0; j < sizeof (unsigned long); j++)
	    ((unsigned char *)&logps)[j] = ssd.bad_p[n].logp;
          while (i0 < (unsigned int) si->bucket_region)
            {
              unsigned long *S_ptr = (unsigned long *) (S + i0);
              unsigned long *end = S_ptr + I / sizeof (unsigned long);
              unsigned long logps2 = logps;
              if ((i0 & I) == 0) /* Even j coordinate? */
                {
                  /* Yes, update only odd i-coordinates */
                  /* Use array indexing to avoid endianness issues. */
		  for (j = 0; j < sizeof (unsigned long); j += 2)
		    ((unsigned char *)&logps2)[j] = 0;
                }
              while (S_ptr < end)
                {
                  *(S_ptr) -= logps2;
                  *(S_ptr + 1) -= logps2;
                  *(S_ptr + 2) -= logps2;
                  *(S_ptr + 3) -= logps2;
                  S_ptr += 4;
                }
              i0 += ssd.bad_p[n].g * I;
            }
          ssd.bad_p[n].next_position = i0 - (1U << LOG_BUCKET_REGION);
        }
      else
	{
	  /* q > 1, more general sieving code. */
	  const unsigned int i0 = ssd.bad_p[n].next_position;
	  const fbprime_t g = ssd.bad_p[n].g;
	  const fbprime_t q = ssd.bad_p[n].q;
	  const fbprime_t U = ssd.bad_p[n].U;
	  const unsigned char logp = ssd.bad_p[n].logp;
	  const fbprime_t evenq = (q % 2 == 0) ? q : 2 * q;
	  unsigned int lineoffset = i0 & (I - 1U),
	               linestart = i0 - lineoffset;
	  ASSERT (U < q);
	  while (linestart < (1U << LOG_BUCKET_REGION))
	    {
	      unsigned int i = lineoffset;
	      if ((linestart & I) == 0) /* Is j even? */
		{
		  /* Yes, sieve only odd i values */
		  if (i % 2 == 0) /* Make i odd */
		    i += q;
		  if (i % 2 == 1) /* If not both i,q are even */
		    for ( ; i < I; i += evenq)
		      {
			if (test_divisibility)
			  test_divisible_x (g*q, linestart + i, bucket_nr, si, side);
                        sieve_decrease (S + linestart + i, logp, 
                                        linestart + i, bucket_nr, q, si, 4);
		      }
		}
	      else
		{
		  for ( ; i < I; i += q)
		    {
		      if (test_divisibility)
			test_divisible_x (g*q, linestart + i, bucket_nr, 
			                  si, side);
                      sieve_decrease (S + linestart + i, logp, linestart + i, 
                                      bucket_nr, q, si, 5);
		    }
		}
              
	      linestart += g * I;
	      lineoffset += U;
	      if (lineoffset >= q)
	        lineoffset -= q;
	    }
	  ssd.bad_p[n].next_position = linestart + lineoffset - 
	                               (1U << LOG_BUCKET_REGION);
	}
    }
}

/* Sieve small primes (p < I, p not in trialdiv_primes list) of the factor
   base fb in the next sieve region S, and add primes and the x position
   where they divide and where there's a sieve report to a bucket (rather
   than subtracting the log norm from S, as during sieving).
   Information about where we are is in ssd.
   Primes in trialdiv_primes must be in increasing order. */
void
resieve_small_bucket_region (bucket_primes_t *BP, unsigned char *S,
			     small_sieve_data_t *ssd,
			     const sieve_info_t *si)
{
  const uint32_t I = si->I;
  unsigned char *S_ptr;
  unsigned long j, nj;
  int n;
  const int resieve_very_verbose = 0, resieve_very_verbose_bad = 0;

  nj = (si->bucket_region >> si->logI);
  ASSERT ((nj & 1) == 0);

  for (n = 0 ; n < ssd->nb_nice_p; ++n) 
    {
      const fbprime_t p = ssd->nice_p[n].p;
      fbprime_t r, twop;
      unsigned int i0;
      
      twop = p + p;
      r = ssd->nice_p[n].r;
      i0 = ssd->nice_p[n].next_position;
      S_ptr = S;
      ASSERT(i0 < p);
      for (j = 0; j < nj; j += 2)
        {
          unsigned int i;
          /* for j even, we sieve only odd i */
          for (i = (i0 & 1) ? i0 : i0 + p ; i < I; i += twop)
            {
              if (S_ptr[i] != 255)
                {
                  bucket_prime_t prime;
                  unsigned int x = (j << (si->logI)) + i;
                  if (resieve_very_verbose) {
                      pthread_mutex_lock(&io_mutex);
                    fprintf (stderr, "resieve_small_bucket_region: root "
                             FBPRIME_FORMAT ",%d divides at x = "
                             "%d = %lu * %u + %d\n",
                             p, r, x, j, 1 << si->logI, i);
                      pthread_mutex_unlock(&io_mutex);
                  }
                  prime.p = p;
                  prime.x = x;
                  push_bucket_prime (BP, prime);
                }
            }
          i0 += r;
          if (i0 >= p)
              i0 -= p;
          S_ptr += I;
          /* j odd */
          for (i = i0 ; i < I; i += p)
              if (S_ptr[i] != 255)
                {
                  bucket_prime_t prime;
                  unsigned int x = ((j + 1) << (si->logI)) + i;
                  if (resieve_very_verbose) {
                      pthread_mutex_lock(&io_mutex);
                    fprintf (stderr, "resieve_small_bucket_region: root "
                             FBPRIME_FORMAT ",%d divides at x = "
                             "%d = %lu * %d + %d\n",
                             p, r, x, j + 1, 1 << si->logI, i);
                      pthread_mutex_unlock(&io_mutex);
                  }
                  prime.p = p;
                  prime.x = x;
                  push_bucket_prime (BP, prime);
                }
          i0 += r;
          if (i0 >= p)
              i0 -= p;
          S_ptr += I;
        }
      ssd->nice_p[n].next_position = i0;
    }


  /* Resieve bad primes */
  if (resieve_very_verbose_bad)
    {
        pthread_mutex_lock(&io_mutex);
      fprintf (stderr, "# %d bad primes to resieve: ", ssd->nb_bad_p);
      for (n = 0; n < ssd->nb_bad_p; ++n)
	fprintf (stderr, "%s" FBPRIME_FORMAT, 
		 (n>0) ? ", " : "", ssd->bad_p[n].g);
      fprintf (stderr, "\n");
      pthread_mutex_unlock(&io_mutex);
    }
  for (n = 0; n < ssd->nb_bad_p; ++n) 
    {
      const fbprime_t g = ssd->bad_p[n].g;

      /* Test every p-th line, starting at S[next_position] */
      unsigned int i, i0 = ssd->bad_p[n].next_position;
      ASSERT (n == 0 || ssd->bad_p[n - 1].g <= ssd->bad_p[n].g);
      ASSERT (i0 % I == 0); /* make sure next_position points at start
                               of line */
      if (resieve_very_verbose_bad) {
          pthread_mutex_lock(&io_mutex);
          fprintf (stderr, "# resieving bad prime " FBPRIME_FORMAT
                  ", i0 = %u\n", g, i0);
          pthread_mutex_unlock(&io_mutex);
      }
      while (i0 < (unsigned int) si->bucket_region)
        {
          unsigned char *S_ptr = S + i0;
          if ((i0 >> si->logI) % 2 == 0) /* Even j coordinate? */
            {
              /* Yes, test only odd i-coordinates */
              for (i = 1; i < I; i += 2)
                {
                  if (S_ptr[i] != 255)
                    {
                      bucket_prime_t prime;
                      const unsigned int x = i0 + i;
                      if (resieve_very_verbose_bad) {
                          pthread_mutex_lock(&io_mutex);
                          fprintf (stderr, "resieve_small_bucket_region even j: root "
                                  FBPRIME_FORMAT ",inf divides at x = %u\n",
                                  g, x);
                          pthread_mutex_unlock(&io_mutex);
                      }
                      prime.p = g;
                      prime.x = x;
                      push_bucket_prime (BP, prime);
                    }
                }
            }
          else
            {
              /* No, test all i-coordinates */
              for (i = 0; i < I; i++)
                {
                  if (S_ptr[i] != 255)
                    {
                      bucket_prime_t prime;
                      const unsigned int x = i0 + i;
                      if (resieve_very_verbose_bad) {
                          pthread_mutex_lock(&io_mutex);
                          fprintf (stderr, "resieve_small_bucket_region odd j: root "
                                  FBPRIME_FORMAT ",inf divides at x = %u\n",
                                  g, x);
                          pthread_mutex_unlock(&io_mutex);
                      }
                      prime.p = g;
                      prime.x = x;
                      push_bucket_prime (BP, prime);
                    }
                }
            }
          i0 += g * I;
        }
      ssd->bad_p[n].next_position = i0 - si->bucket_region;
      if (resieve_very_verbose_bad) {
          pthread_mutex_lock(&io_mutex);
          fprintf (stderr, "# resieving: new i0 = %u, bucket_region = %d, "
                  "new next_position = %d\n",
                  i0, si->bucket_region, ssd->bad_p[n].next_position);
          pthread_mutex_unlock(&io_mutex);
      }
    }
}

typedef struct {
    uint64_t *fac;
    int n;
} factor_list_t;

#define FL_MAX_SIZE 200

void factor_list_init(factor_list_t *fl) {
    fl->fac = (uint64_t *) malloc (FL_MAX_SIZE * sizeof(uint64_t));
    ASSERT_ALWAYS(fl->fac != NULL);
    fl->n = 0;
}

void factor_list_clear(factor_list_t *fl) {
    free(fl->fac);
}

// print a comma-separated list of factors.
// assumes there is at least one factor
void factor_list_fprint(FILE *f, factor_list_t fl) {
    int i;
    for (i = 0; i < fl.n-1; ++i)
        fprintf(f, "%" PRIx64 ",", fl.fac[i]);
    fprintf(f, "%" PRIx64, fl.fac[fl.n-1]);
}


/* The entries in BP must be sorted in order of increasing x */
static void
divide_primes_from_bucket (factor_list_t *fl, mpz_t norm, const int x,
                           bucket_primes_t *BP, const unsigned long fbb)
{
  bucket_prime_t prime;
  while (!bucket_primes_is_end (BP)) {
      prime = get_next_bucket_prime (BP);
      if (prime.x > x)
        {
          rewind_primes_by_1 (BP);
          break;
        }
      if (prime.x == x) {
          unsigned long p = prime.p;
          while (p <= fbb && !mpz_divisible_ui_p (norm, p)) {
              /* It may have been a case of incorrectly reconstructing p
                 from bits 1...16, so let's try if a bigger prime works */
              p += BUCKET_P_WRAP;
          }
          if (p > fbb) {
              pthread_mutex_lock(&io_mutex);
              fprintf (stderr,
                       "# Error, p = %lu does not divide at x = %d\n",
                       (unsigned long) prime.p, x);
              pthread_mutex_unlock(&io_mutex);
              continue;
          }
          do {
              fl->fac[fl->n] = p;
              fl->n++;
              ASSERT_ALWAYS(fl->n <= FL_MAX_SIZE);
              mpz_divexact_ui (norm, norm, p);
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}


NOPROFILE_STATIC void
trial_div (factor_list_t *fl, mpz_t norm, int x,
           factorbase_degn_t *fb, bucket_primes_t *primes,
	   trialdiv_divisor_t *trialdiv_data, const unsigned long fbb)
{
    const int trial_div_very_verbose = 0; // (x == 30878);
    fl->n = 0; /* reset factor list */

    if (trial_div_very_verbose) {
        pthread_mutex_lock(&io_mutex);
        gmp_fprintf (stderr, "# trial_div() entry, x = %d, norm = %Zd\n",
                x, norm);
        pthread_mutex_unlock(&io_mutex);
    }

    // handle 2 separately, if it is in fb
    if (fb->p == 2) {
        int bit = mpz_scan1(norm, 0);
        int i;
        for (i = 0; i < bit; ++i) {
            fl->fac[fl->n] = 2;
            fl->n++;
        }
        if (trial_div_very_verbose) {
            pthread_mutex_lock(&io_mutex);
            gmp_fprintf (stderr, "# x = %d, dividing out 2^%d, norm = %Zd\n",
                    x, bit, norm);
            pthread_mutex_unlock(&io_mutex);
        }
        mpz_tdiv_q_2exp(norm, norm, bit);
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket (fl, norm, x, primes, fbb);
    if (trial_div_very_verbose) {
        pthread_mutex_lock(&io_mutex);
        gmp_fprintf (stderr, "# x = %d, after dividing out bucket/resieved norm = %Zd\n",
                x, norm);
        pthread_mutex_unlock(&io_mutex);
    }

    {
      /* Trial divide primes with precomputed tables */
      int nr_factors, i;
      unsigned long factors[32];
      if (trial_div_very_verbose)
      {
          pthread_mutex_lock(&io_mutex);
          fprintf (stderr, "# Trial division by ");
          for (i = 0; trialdiv_data[i].p != 1; i++)
              fprintf (stderr, " %lu", trialdiv_data[i].p);
          fprintf (stderr, "\n");
          pthread_mutex_unlock(&io_mutex);
      }

      nr_factors = trialdiv (factors, norm, trialdiv_data, 32);
      ASSERT (nr_factors <= 32);

      for (i = 0; i < nr_factors; i++)
      {
          if (trial_div_very_verbose) {
              pthread_mutex_lock(&io_mutex);
              fprintf (stderr, " %lu", factors[i]);
              pthread_mutex_unlock(&io_mutex);
          }
          fl->fac[fl->n++] = factors[i];
      }
      if (trial_div_very_verbose) {
          pthread_mutex_lock(&io_mutex);
          gmp_fprintf (stderr, "\n# After trialdiv(): norm = %Zd\n", norm);
          pthread_mutex_unlock(&io_mutex);
      }
    }
}

/* Return 0 if the leftover norm n cannot yield a relation:
   (a) if n > 2^mfb
   (b) if L < n < B^2
   (c) if L^2 < n < B^3

   FIXME: need to check L^k < n < B^(k+1) too.
*/
static int
check_leftover_norm (mpz_t n, size_t lpb, mpz_t BB, mpz_t BBB, size_t mfb)
{
  size_t s = mpz_sizeinbase (n, 2);

  if (s > mfb)
    return 0;
  if (lpb < s && mpz_cmp (n, BB) < 0)
    return 0;
  if (2 * lpb < s && mpz_cmp (n, BBB) < 0)
    return 0;
  if (lpb < s && mpz_probab_prime_p (n, 1))
    return 0;
  return 1;
}

#ifdef UNSIEVE_NOT_COPRIME
static void
unsieve_one_prime (unsigned char *line_start, const unsigned int p, 
                   const unsigned int y, const sieve_info_t *si)
{
  unsigned int x, np = p; /* if 2|y, np=2p, else np=p */

  x = (si->I / 2U) % p;
  if (y % 2U == 0U)
    {
      np += p;
      if (x % 2U == 0U)
        x += p;
    }
  for ( ; x < si->I; x += np)
    line_start[x] = 255;
}


/* Set locations where gcd(i,j) != 1 to 255*/
void 
unsieve_not_coprime (unsigned char *S, const int N, const sieve_info_t *si)
{
  unsigned int y; /* Line coordiante within bucket region */
  for (y = 0U + (N == 0U ? 1U : 0U); 
       y < 1U << (LOG_BUCKET_REGION - si->logI); y++)
    {
      unsigned int c = y + (N << (LOG_BUCKET_REGION - si->logI));
      unsigned int p;
      unsigned char *line_start = S + (y << si->logI);

      p = si->lpf[c]; /* set p to largest prime factor of c */
      if (p > 3U)
        {
          ASSERT (c % p == 0U);
          unsieve_one_prime (line_start, p, y, si);
          do {c /= p;} while (c % p == 0U);
        }
      
      while (c % 2U == 0U) 
        c >>= 1;
      
      if (c % 3U == 0U)
        {
          const unsigned int I_ul  = si->I / sizeof (unsigned long);
          unsigned int i, x;
          unsigned long s[3] = {0UL, 0UL, 0UL};
          unsigned long * restrict ul_line_start = (unsigned long *) line_start;
          
          /* Sieve only the small pattern */
          for (x = (si->I / 2U) % 3U; x < 3U * sizeof(unsigned long); x += 3U)
            ((unsigned char *) s)[x] = 255;
          
          /* Then apply pattern to array */
          for (i = 0U; i < I_ul - 2U; i += 3U)
            {
              ul_line_start[i] |= s[0];
              ul_line_start[i + 1] |= s[1];
              ul_line_start[i + 2] |= s[2];
            }
          if (i < I_ul - 1U)
            line_start[i] |= s[0];
          if (i < I_ul)
            line_start[i + 1] |= s[1];
          
          do {c /= 3U;} while (c % 3U == 0U);
        }

      for (p = 5U; p * p <= c; p += 2U)
        if (c % p == 0U)
          {
            unsieve_one_prime (line_start, p, y, si);
            do {c /= p;} while (c % p == 0U);
          }

      /* Now c == 1 or c is a prime */
      if (c != 1U)
        {
          ASSERT(c > 3U && si->lpf[c] == c);
          unsieve_one_prime (line_start, c, y, si);
        }
    }
}
#endif /* ifdef UNSIEVE_NOT_COPRIME */


/* Adds the number of sieve reports to *survivors,
   number of survivors with coprime a, b to *coprimes */

NOPROFILE_STATIC int
factor_survivors (const unsigned char *rat_S, unsigned char *alg_S, int N, 
		  bucket_array_t *rat_BA,
                  bucket_array_t *alg_BA, factorbase_degn_t *fb_rat,
                  factorbase_degn_t *fb_alg, cado_poly cpoly, sieve_info_t *si,
                  unsigned long *survivors, unsigned long *coprimes,
		  small_sieve_data_t *srsd_alg, small_sieve_data_t *srsd_rat,
		  unsigned long *report_sizes_a, unsigned long *report_sizes_r,
		  FILE *output)
{
    int x, xul;
    int64_t a;
    uint64_t b;
    int cpt = 0;
    int surv = 0, copr = 0;
    mpz_t alg_norm, rat_norm, BBalg, BBrat, BBBalg, BBBrat, BLPrat;
    factor_list_t alg_factors, rat_factors;
    mpz_array_t *f_r = NULL, *f_a = NULL;    /* large prime factors */
    uint32_array_t *m_r = NULL, *m_a = NULL; /* corresponding multiplicities */
    bucket_primes_t rat_primes, alg_primes;

    f_r = alloc_mpz_array (8);
    f_a = alloc_mpz_array (8);
    m_r = alloc_uint32_array (8);
    m_a = alloc_uint32_array (8);

    factor_list_init(&alg_factors);
    factor_list_init(&rat_factors);
    mpz_init (alg_norm);
    mpz_init (rat_norm);
    mpz_init (BBalg);
    mpz_init (BBrat);
    mpz_init (BBBalg);
    mpz_init (BBBrat);
    mpz_init (BLPrat);
    mpz_ui_pow_ui (BBalg, cpoly->alim, 2);
    mpz_ui_pow_ui (BBrat, cpoly->rlim, 2);
    mpz_mul_ui (BBBalg, BBalg, cpoly->alim);
    mpz_mul_ui (BBBrat, BBrat, cpoly->rlim);
    mpz_set_ui (BLPrat, cpoly->alim);
    mpz_mul_2exp (BLPrat, BLPrat, cpoly->lpba); /* fb bound * lp bound */

#ifdef UNSIEVE_NOT_COPRIME
    unsieve_not_coprime (alg_S, N, si);
#endif
    
    for (x = 0; x < si->bucket_region; ++x)
      {
        unsigned int X;
        unsigned int i, j;

        if (si->alg_Bound[alg_S[x]] + si->rat_Bound[rat_S[x]] >= 127)
          {
            alg_S[x] = 255;
            continue;
          }
        surv++;

        X = x + (N << LOG_BUCKET_REGION);
        i = abs ((int) (X & (si->I - 1)) - si->I / 2);
        j = X >> si->logI;
#ifndef UNSIEVE_NOT_COPRIME
        if (bin_gcd_safe (i, j) != 1)
          {
            alg_S[x] = 255;
            continue;
          }
#endif
      }

    /* Copy those bucket entries that belong to sieving survivors and
       store them with the complete prime */
    /* FIXME: choose a sensible size here */
    alg_primes = init_bucket_primes (si->bucket_region);
    { 
        int i;
        for (i = 0; i < si->nb_threads; ++i) 
            purge_bucket (&alg_primes, alg_BA[i], N, alg_S);
    }

    /* Resieve small primes for this bucket region and store them 
       together with the primes recovered from the bucket updates */
    resieve_small_bucket_region (&alg_primes, alg_S, srsd_alg, si);
    /* Sort the entries to avoid O(n^2) complexity when looking for
       primes during trial division */
    bucket_sortbucket (&alg_primes);

    /* Do the same thing for rational side */
    rat_primes = init_bucket_primes (si->bucket_region);
    {
        int i;
        for (i = 0; i < si->nb_threads; ++i)
            purge_bucket (&rat_primes, rat_BA[i], N, alg_S);
    }
    resieve_small_bucket_region (&rat_primes, alg_S, srsd_rat, si);
    bucket_sortbucket (&rat_primes);

    /* Scan array one long word at a time. If any byte is <255, i.e. if
       the long word is != 0xFFFF...FF, examine the bytes */
    for (xul = 0; xul < si->bucket_region; xul += sizeof (unsigned long))
      if (*(unsigned long *)(alg_S + xul) != (unsigned long)(-1L))
        for (x = xul; x < xul + (int) sizeof (unsigned long); ++x)
          {
            unsigned int i, j;

            if (alg_S[x] == 255)
              continue;

            // Compute algebraic and rational norms.
            xToAB (&a, &b, x + N*si->bucket_region, si);

            /* since a,b both even were not sieved, either a or b should be odd */
            // ASSERT((a | b) & 1);
            if (UNLIKELY(((a | b) & 1) == 0))
              {
                pthread_mutex_lock(&io_mutex);
                fprintf (stderr, "# Error: a and b both even for N = %d, x = %d,\n"
                                 "i = %d, j = %d, a = %ld, b = %lu\n",
                         N, x, ((x + N*si->bucket_region) & (si->I - 1))
                           - (si->I >> 1),
                         (x + N*si->bucket_region) >> si->logI,
                         (long) a, (unsigned long) b);
                pthread_mutex_unlock(&io_mutex);
                continue;
              }

            /* Since the q-lattice is exactly those (a, b) with
               a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
            if (b == 0 || (b >= si->q && b % si->q == 0))
              continue;

            copr++;

            /* For hunting missed relations */
#if 0
            if (a == -6537753 && b == 1264)
              fprintf (stderr, "# Have relation %ld,%lu at bucket nr %d, "
                       "x = %d, K = %lu\n", 
                       a, b, N, x, (unsigned long) N * si->bucket_region + x);
#endif

            // Trial divide rational norm
            eval_fij (rat_norm, (const mpz_t *) cpoly->g, 1, a, b);
            if (si->ratq)
                mpz_divexact_ui (rat_norm, rat_norm, si->q);
            trial_div (&rat_factors, rat_norm, x, fb_rat,
                       &rat_primes, si->trialdiv_data_rat, cpoly->rlim);

            if (!check_leftover_norm (rat_norm, cpoly->lpbr, BBrat, BBBrat,
                                      cpoly->mfbr))
              continue;

            // Trial divide algebraic norm
            eval_fij (alg_norm, (const mpz_t *) cpoly->f, cpoly->degree, a, b);
            if (!si->ratq)
                mpz_divexact_ui (alg_norm, alg_norm, si->q);
            trial_div (&alg_factors, alg_norm, x, fb_alg,
                       &alg_primes, si->trialdiv_data_alg, cpoly->alim);

            if (!check_leftover_norm (alg_norm, cpoly->lpba, BBalg, BBBalg,
                                      cpoly->mfba))
              continue;

            if (mpz_cmp (rat_norm, BLPrat) > 0)
              {
                /* rat_norm might not be smooth, factor it first */
                if (factor_leftover_norm (rat_norm, cpoly->lpbr, f_r, m_r,
                                          si->strategy) == 0)
                  continue;

                if (factor_leftover_norm (alg_norm, cpoly->lpba, f_a, m_a,
                                          si->strategy) == 0)
                  continue;
              }
            else
              {
                /* rat_norm is definitely smooth, factor it last */
                if (factor_leftover_norm (alg_norm, cpoly->lpba, f_a, m_a,
                                          si->strategy) == 0)
                  continue;

                if (factor_leftover_norm (rat_norm, cpoly->lpbr, f_r, m_r,
                                          si->strategy) == 0)
                  continue;
              }

#ifdef UNSIEVE_NOT_COPRIME
            ASSERT (bin_gcd_safe (a, b) == 1);
#endif

            pthread_mutex_lock(&io_mutex);
            if (!si->bench)
                fprintf (output, "%" PRId64 ",%" PRIu64 ":", a, b);
            int first_factor = 0;
            if (rat_factors.n != 0) {
                if (!si->bench)
                    factor_list_fprint (output, rat_factors);
                first_factor = 1;
            }
            for (i = 0; i < f_r->length; ++i)
              for (j = 0; j < m_r->data[i]; j++)
                  if (!first_factor) {
                      if (!si->bench)
                          gmp_fprintf (output, "%Zx", f_r->data[i]);
                      first_factor = 1;
                  } else {
                      if (!si->bench)
                          gmp_fprintf (output, ",%Zx", f_r->data[i]);
                  }
            if (si->ratq) {
                if (!si->bench)
                    fprintf (output, ",%" PRIx64 "", si->q);
            }
            if (!si->bench)
                fprintf (output, ":");
            first_factor = 0;
            if (alg_factors.n != 0) {
                if (!si->bench)
                    factor_list_fprint (output, alg_factors);
                first_factor = 1;
            }
            for (i = 0; i < f_a->length; ++i)
              for (j = 0; j < m_a->data[i]; j++)
                  if (!first_factor) {
                      if (!si->bench)
                          gmp_fprintf (output, "%Zx", f_a->data[i]);
                      first_factor = 1;
                  } else {
                      if (!si->bench)
                          gmp_fprintf (output, ",%Zx", f_a->data[i]);
                  }
            /* print special q */
            if (!si->ratq) {
                if (!si->bench)
                    fprintf (output, ",%" PRIx64 "", si->q);
            }
            if (!si->bench) {
                fprintf (output, "\n");
                fflush (output);
            }
            pthread_mutex_unlock(&io_mutex);
            cpt++;
	    /* Build histogram of lucky S[x] values */
            report_sizes_a[alg_S[x]]++;
	    report_sizes_r[rat_S[x]]++;
          }

    survivors[0] += surv;
    coprimes[0] += copr;
    clear_bucket_primes (&rat_primes);
    clear_bucket_primes (&alg_primes);
    mpz_clear (BBalg);
    mpz_clear (BBrat);
    mpz_clear (BBBalg);
    mpz_clear (BBBrat);
    mpz_clear (BLPrat);
    mpz_clear(alg_norm);
    mpz_clear(rat_norm);
    factor_list_clear(&alg_factors);
    factor_list_clear(&rat_factors);
    clear_mpz_array (f_r);
    clear_mpz_array (f_a);
    clear_uint32_array (m_r);
    clear_uint32_array (m_a);

    return cpt;
}


/****************************************************************************/

// Conversions between different representations for sieve locations:
//   N          is the number of the bucket region, N in [0, nr_buckets[
//   x          is the index in the sieving array. x in [0,I*J[
//   (i,j)      is the coordinates in the q-lattice. i in [-I/2,I/2[
//                                                   j in [0,J[
//   (a,b)      is the original coordinates. a is signed, b is unsigned.

static void 
NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, 
       const sieve_info_t * si)
{
    unsigned int X = x + N * si->bucket_region;
    *i = (X % (si->I)) - (si->I >> 1);
    *j = X / si->I;
}

#if 0
static void
IJTox(int *x, unsigned int *N, int i, int j, sieve_info_t * si)
{
    *x = i + (si->I)*j + (si->I>>1);
}
#endif

static void
IJToAB(int64_t *a, uint64_t *b, const int i, const unsigned int j, 
       const sieve_info_t * si)
{
    int64_t s, t;
    s = (int64_t)i * (int64_t) si->a0 + (int64_t)j * (int64_t)si->a1;
    t = (int64_t)i * (int64_t) si->b0 + (int64_t)j * (int64_t)si->b1;
    if (t >= 0)
      {
        *a = s;
        *b = t;
      }
    else
      {
        *a = -s;
        *b = -t;
      }
}

/* Warning: b might be negative, in which case we return (-a,-b) */
static void
xToAB(int64_t *a, uint64_t *b, const unsigned int x, const sieve_info_t *si)
{
    int i, j;
    int64_t c;
    uint32_t I = si->I;

    i = (x & (I - 1)) - (I >> 1);
    j = x >> si->logI;
    *a = (int64_t) i * (int64_t) si->a0 + (int64_t) j * (int64_t) si->a1;
    c =  (int64_t) i * (int64_t) si->b0 + (int64_t) j * (int64_t) si->b1;
    if (c >= 0)
      *b = c;
    else
      {
        *a = -*a;
        *b = -c;
      }
}

void test_divisible_x (const fbprime_t p, const unsigned long x, 
                       const int n, const sieve_info_t *si, const char side)
{
  const unsigned long X = x + n * si->bucket_region;
  const long i = (long) (X & ((unsigned long) si->I - 1UL)) - (long) si->I / 2;
  const unsigned long j = X >> si->logI;
  mpz_t v;

  mpz_init (v);
  if (side == 'a')
    eval_fij (v, (const mpz_t *) si->fij, si->degree, i, j);
  else if (side == 'r')
    eval_fij (v, (const mpz_t *) si->gij, 1, i, j);
  else
    abort();
  
  if (!mpz_divisible_ui_p (v, (unsigned long) p))
   gmp_fprintf (stderr, "# test_divisible_x (" FBPRIME_FORMAT 
                ", %lu, %d, %c): i = %ld, j = %lu, v = %Zd\n", 
                p, x, n, side, i, j, v);
  
  mpz_clear (v);
}


/*********************** norm computation ************************************/

/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
   Put in fijd[] a double-precision approximation of fij[]/q.
*/
static void
fij_from_f (mpz_t *fij, const sieve_info_t *si, mpz_t *f, const int d)
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;

  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  g = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (g[k]);
  mpz_init (f0);

  /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
     F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
     We use that formula recursively. */

  mpz_set_ui (g[0], 1); /* g = 1 */

  for (k = d - 1; k >= 0; k--)
    {
      /* invariant: we have already translated coefficients of degree > k,
         in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
         g[0..d - (k+1)]:
         f[k]   <- f[k] + a1*f[k+1]
         ...
         f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
         ...
         f[d] <- a0*f[d] */
      mpz_swap (f0, fij[k]); /* save the new constant coefficient */
      mpz_mul_si (fij[k], fij[k + 1], si->a1);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], si->a0);
          mpz_addmul_si (fij[l], fij[l + 1], si->a1);
        }
      mpz_mul_si (fij[d], fij[d], si->a0);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], si->b0);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], si->b1);
          mpz_addmul_si (g[l], g[l-1], si->b0);
        }
      mpz_mul_si (g[0], g[0], si->b1);

      /* now g has degree d-k, and we add f0*g */
      for (l = k; l <= d; l++)
        mpz_addmul (fij[l], g[l - k], f0);
    }

  mpz_clear (f0);
  for (k = 0; k <= d; k++)
    mpz_clear (g[k]);
  free (g);
}

/* return max |g(x)| for 0 <= x <= s,
   where g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double *g, const unsigned int d, double s)
{
  unsigned int k, l, sign_change, new_sign_change;
  double **dg;    /* derivatives of g */
  double a, va, b, vb;
  double *roots, gmax;

  dg = (double**) malloc (d * sizeof (double*));
  dg[0] = g;
  for (k = 1; k < d; k++) /* dg[k] is the k-th derivative, thus has
                             degree d-k, i.e., d-k+1 coefficients */
    dg[k] = (double*) malloc ((d - k + 1) * sizeof (double));
  roots = (double*) malloc (d * sizeof (double));
  for (k = 1; k < d; k++)
    for (l = 0; l <= d - k; l++)
      dg[k][l] = (l + 1) * dg[k - 1][l + 1];
  /* now dg[d-1][0]+x*dg[d-1][1] is the (d-1)-th derivative: it can have at
     most one sign change, iff dg[d-1][0] and dg[d-1][0]+dg[d-1][1] have
     different signs */
  if (dg[d-1][0] * (dg[d-1][0] + dg[d-1][1]) < 0)
    {
      sign_change = 1;
      roots[0] = - dg[d-1][0] / dg[d-1][1]; /* root of (d-1)-th derivative */
    }
  else
    sign_change = 0;
  roots[sign_change] = s; /* end of interval */
  for (k = d - 1; k-- > 1;)
    {
      /* invariant: sign_change is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_change-1], and roots[sign_change] = s. */
      a = 0.0;
      va = dg[k][0]; /* value of dg[k] at x=0 */
      new_sign_change = 0;
      for (l = 0; l <= sign_change; l++)
        {
          b = roots[l]; /* root of dg[k+1], or end of interval */
          vb = fpoly_eval (dg[k], d - k, b);
          if (va * vb < 0) /* root in interval */
            roots[new_sign_change++] = fpoly_dichotomy (dg[k], d - k,
                                                        a, b, va, 20);
          a = b;
          va = vb;
        }
      roots[new_sign_change] = s; /* end of interval */
      sign_change = new_sign_change;
    }
  /* now all extrema of g are 0, roots[0], ..., roots[sign_change] = s */
  gmax = fabs (g[0]);
  for (k = 0; k <= sign_change; k++)
    {
      va = fabs (fpoly_eval (g, d, roots[k]));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  for (k = 1; k < d; k++)
    free (dg[k]);
  free (dg);
  return gmax;
}

/* returns the maximal value of log2 |F(a,b)/q| for
   a = a0 * i + a1 * j, b = b0 * i + b1 * j and q >= q0,
   -I/2 <= i <= I/2, 0 <= j <= I/2*min(s*B/|a1|,B/|b1|)
   where B >= sqrt(2*q/s/sqrt(3)) for all special-q in the current range
   (s is the skewness, and B = si.B, see lattice.tex).

   Since |a0| <= s*B and |b0| <= B, then
   |a0 * i + a1 * j| <= s*B*I and |b0 * i + b1 * j| <= B*I,
   thus it suffices to compute M = max |F(x,y)| in the rectangle
   -s <= x <= s, 0 <= y <= 1, and to multiply M by (B*I)^deg(F).

   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(s,y) for 0 <= y <= 1
   (b) either on F(x,1) for -s <= x <= s
   (c) either on F(-s,y) for 0 <= y <= 1
   (d) or on F(x,0) for -s <= x <= s, but this maximum is f[d]*s^d,
       and is attained in (a) or (c).
*/
static double
get_maxnorm (cado_poly cpoly, sieve_info_t *si, uint64_t q0)
{
  unsigned int d = cpoly->degree, k;
  double *fd; /* double-precision coefficients of f */
  double norm, max_norm, pows, tmp;

  fd = (double*) malloc ((d + 1) * sizeof (double));
  for (k = 0; k <= d; k++)
    fd[k] = mpz_get_d (cpoly->f[k]);

  /* (b1) determine the maximum of |f(x)| for 0 <= x <= s */
  max_norm = get_maxnorm_aux (fd, d, cpoly->skew);

  /* (b2) determine the maximum of |f(-x)| for 0 <= x <= s */
  norm = get_maxnorm_aux (fd, d, -cpoly->skew);
  if (norm > max_norm)
    max_norm = norm;

  for (pows = 1.0, k = 0; k <= d; k++)
    {
      fd[k] *= pows;
      pows *= cpoly->skew;
    }
  /* swap coefficients; if d is odd, we need to go up to k = floor(d/2) */
  for (k = 0; k <= d / 2; k++)
    {
      tmp = fd[k];
      fd[k] = fd[d - k];
      fd[d - k] = tmp;
    }

  /* (a) determine the maximum of |g(y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, 1.0);
  if (norm > max_norm)
    max_norm = norm;

  /* (c) determine the maximum of |g(-y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, -1.0);
  if (norm > max_norm)
    max_norm = norm;

  free (fd);

  /* multiply by (B*I)^d and divide by q0 if sieving on alg side */
  tmp = max_norm * pow (si->B * (double) si->I, (double) d);
  if (!si->ratq)
      tmp /= (double) q0;
  return log2(tmp);
}

/* v <- |f(i,j)|, where f is of degree d */
NOPROFILE_STATIC void
eval_fij (mpz_t v, const mpz_t *f, const unsigned int d, const long i,
	  const unsigned long j)
{
  unsigned int k;
  mpz_t jpow;

  mpz_init_set_ui (jpow, 1);
  mpz_set (v, f[d]);
  for (k = d; k-- > 0;)
    {
      mpz_mul_si (v, v, i);
      mpz_mul_ui (jpow, jpow, j);
      mpz_addmul (v, f[k], jpow);
    }
  mpz_abs (v, v); /* avoids problems with negative norms */
  mpz_clear (jpow);
}


/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

void
SkewGauss (sieve_info_t *si, double skewness)
{
  double a[2], b[2], q, maxab0, maxab1;

  a[0] = (double) si->q;
  ASSERT_ALWAYS(a[0] < 9007199254740992.0); /* si.q should be less than 2^53
                                               so that a[0] is exact */
  b[0] = 0.0;
  a[1] = (double) si->rho;
  b[1] = skewness;
  ASSERT(b[1] != 0);
  while (1)
    {
      /* reduce vector (a[0], b[0]) with respect to (a[1], b[1]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[1] * a[1] + b[1] * b[1]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[0] -= q * a[1];
      b[0] -= q * b[1];

      /* reduce vector (a[1], b[1]) with respect to (a[0], b[0]) */
      q = (a[0] * a[1] + b[0] * b[1]) / (a[0] * a[0] + b[0] * b[0]);
      q = rint (q);
      if (q == 0.0)
        break;
      a[1] -= q * a[0];
      b[1] -= q * b[0];
    }
  ASSERT(fits_int32_t(a[0]));
  ASSERT(fits_int32_t(b[0] / skewness));
  ASSERT(fits_int32_t(a[1]));
  ASSERT(fits_int32_t(b[1] / skewness));
  /* now b[0], b[1] should be of the form i*skewness, but this might not be
     exact due to rounding errors, thus we round them to the nearest integer */
  maxab0 = fabs (a[0]) > fabs (b[0]) ? fabs (a[0]) : fabs (b[0]);
  maxab1 = fabs (a[1]) > fabs (b[1]) ? fabs (a[1]) : fabs (b[1]);
  if (maxab0 <= maxab1)
    {
      si->a0 = (int32_t) a[0];
      si->b0 = (int32_t) rint (b[0] / skewness);
      si->a1 = (int32_t) a[1];
      si->b1 = (int32_t) rint (b[1] / skewness);
    }
  else /* swap (a0,b0) and (a1,b1) */
    {
      si->a1 = (int32_t) a[0];
      si->b1 = (int32_t) rint (b[0] / skewness);
      si->a0 = (int32_t) a[1];
      si->b0 = (int32_t) rint (b[1] / skewness);
      maxab1 = maxab0;
    }

  /* make sure J does not exceed I/2 */
  if (maxab1 >= si->B)
    si->J = (uint32_t) (si->B * skewness / maxab1 * (double) (si->I >> 1));
  else
    si->J = si->I >> 1;

  /* Make sure the bucket region size divides the sieve region size, 
     partly covered bucket regions may lead to problems when 
     reconstructing p from half-empty buckets. */
  /* Compute number of i-lines per bucket region, must be integer */
  ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->logI);
  uint32_t i = 1U << (LOG_BUCKET_REGION - si->logI);
  i *= si->nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */
  si->J = ((si->J - 1U) / i + 1U) * i; /* Round up to multiple of i */
}

/* return max(|a0|,|a1|)/min(|a0|,|a1|) */
double
skewness (sieve_info_t *si)
{
  double a0 = fabs ((double) si->a0);
  double a1 = fabs ((double) si->a1);

  return (a0 > a1) ? a0 / a1 : a1 / a0;
}

/************************ cofactorization ********************************/

/* FIXME: the value of 20 seems large. Normally, a few Miller-Rabin passes
   should be enough. See also http://www.trnicely.net/misc/mpzspsp.html */
#define NMILLER_RABIN 1 /* in the worst case, what can happen is that a
                           composite number is declared as prime, thus
                           a relation might be missed, but this will not
                           affect correctness */
#define IS_PROBAB_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
#define NFACTORS        8 /* maximal number of large primes */

/* This function was contributed by Jerome Milan (and bugs were introduced
   by Paul Zimmermann :-).
   Input: n - the number to be factored (leftover norm). Must be composite!
          l - (large) prime bit size bound is L=2^l
   Assumes n > 0.
   Return value:
          0 if n has a prime factor larger than 2^l
          1 if all prime factors of n are < 2^l
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/
int
factor_leftover_norm (mpz_t n, unsigned int l,
                      mpz_array_t* const factors, uint32_array_t* const multis,
		      facul_strategy_t *strategy)
{
  uint32_t i, nr_factors;
  unsigned long ul_factors[16];
  int facul_code;

  factors->length = 0;
  multis->length = 0;

  /* factoring programs do not like 1 */
  if (mpz_cmp_ui (n, 1) == 0)
    return 1;

  /* If n < L, we know that n is prime, since all primes < B have been
     removed, and L < B^2 in general, where B is the factor base bound,
     thus we only need a primality test when n > L. */
  if (BITSIZE(n) <= l)
    {
      append_mpz_to_array (factors, n);
      append_uint32_to_array (multis, 1);
      return 1;
    }
/* Input is required to be composite!
  else if (IS_PROBAB_PRIME(n))
    {
      return 0;
    }
*/

  /* use the facul library */
  // gmp_printf ("facul: %Zd\n", n);
  facul_code = facul (ul_factors, n, strategy);

  if (facul_code == FACUL_NOT_SMOOTH)
    return 0;

  ASSERT (facul_code == 0 || mpz_cmp_ui (n, ul_factors[0]) != 0);

  if (facul_code > 0)
    {
      nr_factors = facul_code;
      for (i = 0; i < nr_factors; i++)
	{
	  unsigned long r;
	  mpz_t t;
	  if (ul_factors[i] > (1UL << l)) /* Larger than large prime bound? */
	    return 0;
	  r = mpz_tdiv_q_ui (n, n, ul_factors[i]);
	  ASSERT_ALWAYS (r == 0UL);
	  mpz_init (t);
	  mpz_set_ui (t, ul_factors[i]);
	  append_mpz_to_array (factors, t);
	  mpz_clear (t);
	  append_uint32_to_array (multis, 1); /* FIXME, deal with repeated
						 factors correctly */
	}

      if (mpz_cmp_ui (n, 1UL) == 0)
	return 1;
      unsigned int s = BITSIZE(n);
      if (s <= l)
        {
          append_mpz_to_array (factors, n);
          append_uint32_to_array (multis, 1);
          return 1;
        }
      /* If we still have more than two primes (or something non-smooth),
         bail out */
      if (s > 2*l)
        return 0;
      /* We always abort below, so let's skip the prp test
      if (IS_PROBAB_PRIME(n))
        return 0; */
    }
  /* When sieving for 3 large primes, here are so many left over, non-smooth
     numbers here that factoring them all takes a long time, for few
     additional relations */
  return 0;
}


typedef struct {
    bucket_array_t *alg_BA;
    factorbase_degn_t *fb_alg;
    bucket_array_t *rat_BA;
    factorbase_degn_t *fb_rat;
    small_sieve_data_t *ssd_rat;
    small_sieve_data_t *ssd_alg;
    cado_poly_ptr cpoly;
    sieve_info_t *si;
    int id;
    FILE *output;
    int verbose;
} process_bucket_region_arg_t;

typedef struct {
    int reports;
    unsigned long survivors0;
    unsigned long survivors1;
    unsigned long survivors2;
    double tn_rat;
    double tn_alg;
    double ttsm;
    double ttf;
    unsigned long report_sizes_a[256];
    unsigned long report_sizes_r[256];
} process_bucket_region_report_t;

/* id gives the number of the thread: it is supposed to deal with the set
 * of bucket_regions corresponding to that number
 */
process_bucket_region_report_t * 
process_regions_one_thread(process_bucket_region_arg_t *arg)
{
    bucket_array_t *alg_BA      = arg->alg_BA;
    bucket_array_t *rat_BA      = arg->rat_BA;
    small_sieve_data_t *ssd_rat = arg->ssd_rat;
    small_sieve_data_t *ssd_alg = arg->ssd_alg;
    factorbase_degn_t *fb_rat   = arg->fb_rat;
    factorbase_degn_t *fb_alg   = arg->fb_alg;
    cado_poly_ptr cpoly         = arg->cpoly;
    sieve_info_t *si            = arg->si;
    int id                      = arg->id;
    FILE *output                = arg->output;
    int i;

    unsigned long survivors0, survivors1, survivors2;
    int reports = 0;
    unsigned long report_sizes_a[256], report_sizes_r[256];
    double tn_rat = 0.0;
    double tn_alg = 0.0;
    double ttsm = 0.0;
    double ttf = 0.0;

    survivors0 = survivors1 = survivors2 = 0;

    /* make copies of small sieve data: the "next_position" field is
     * specific to each thread.
     */
    small_sieve_data_t lssd_alg, lssd_rat;
    clone_small_sieve (&lssd_alg, ssd_alg);
    clone_small_sieve (&lssd_rat, ssd_rat);

    /* find the first sieve location for each small prime */
    small_sieve_data_t *ssd;
    int side;
    for (side = 0; side < 2; side++) {
        if (side == 0) 
            ssd = &lssd_alg;
        else
            ssd = &lssd_rat;
        int n;
	/* Number of lines to skip to get the first bucket for this thread */
        const unsigned long nj = (si->bucket_region >> si->logI) * 
	  (si->nb_buckets / si->nb_threads) * id;
        ASSERT(si->nb_buckets % si->nb_threads == 0);

        // nice primes
        for (n = 0; n < ssd->nb_nice_p; ++n) {
            unsigned long i0;

            if (ssd->nice_p[n].p % 2 == 0 && nj % 2 == 1)
	      {
		/* Make sure that next_position points to a location
		   where i and j are not both even */
		i0 = ssd->nice_p[n].next_position & (si->I - 1);
		ASSERT (i0 < ssd->nice_p[n].p);
		i0 += (nj + 1) * ssd->nice_p[n].r;
		i0 = i0 & (ssd->nice_p[n].p - 1);
		ssd->nice_p[n].next_position = i0 + 
		  (ssd->nice_p[n].next_position & (~(si->I - 1)));
	      }
	    else
	      {
		/* We want to add nj*r to the offset *relative to the 
		   start of the line*, but next_position may be larger 
		   than I, so we treat the multiple-of-I and mod-I parts
		   separately */
		i0 = ssd->nice_p[n].next_position & (si->I - 1);
		ASSERT (i0 < ssd->nice_p[n].p);
		i0 += nj*ssd->nice_p[n].r;
		i0 = i0 % ssd->nice_p[n].p;
		ssd->nice_p[n].next_position = i0 + 
		  (ssd->nice_p[n].next_position & (~(si->I - 1)));
	      }
        }
	
        // bad primes
        for (n = 0; n < ssd->nb_bad_p; ++n) 
	  {
            /* First line to sieve is the smallest j with g|j and j >= nj,
	       however, if nj == 0 we don't sieve j==0 since it contains
	       only one possible relation (i,j) = (1,0). */
	    unsigned int ng, x, j;
	    ng = iceildiv(nj, ssd->bad_p[n].g);
	    if (ng == 0)
	      ng++;
	    x = (si->I / 2 + ng * ssd->bad_p[n].U) % ssd->bad_p[n].q;
	    j = ng * ssd->bad_p[n].g;
	    ssd->bad_p[n].next_position = (j - nj) * si->I + x;
	  }
    }

    /* Yet another copy: used in factor_survivors for resieving small
     * primes */
    small_sieve_data_t lsrsd_alg, lsrsd_rat;
    copy_small_sieve (&lsrsd_alg, &lssd_alg, si->trialdiv_primes_alg);
    copy_small_sieve (&lsrsd_rat, &lssd_rat, si->trialdiv_primes_rat);

    /* local sieve region */
    unsigned char *alg_S, *rat_S;
    alg_S = (unsigned char *)malloc(si->bucket_region*sizeof(unsigned char));
    ASSERT_ALWAYS (alg_S != NULL);
    rat_S = (unsigned char *)malloc(si->bucket_region*sizeof(unsigned char));
    ASSERT_ALWAYS (rat_S != NULL);

    for (i = 0; i < 256; i++) 
      {
	report_sizes_a[i] = 0;
	report_sizes_r[i] = 0;
      }

    /* loop over appropriate set of sieve regions */
    i = (si->nb_buckets / si->nb_threads) * id;
    for (; i < (si->nb_buckets / si->nb_threads) * (id+1); ++i) 
      {
	int j;
        /* Init rational norms */
        tn_rat -= seconds ();
        init_rat_norms_bucket_region(rat_S, i, cpoly, si);
        tn_rat += seconds ();
        /* Apply rational buckets */
        ttsm -= seconds();
        for (j = 0; j < si->nb_threads; ++j) 
            apply_one_bucket(rat_S, rat_BA[j], i, si);
        ttsm += seconds();
        /* Sieve small rational primes */
        sieve_small_bucket_region(rat_S, i, lssd_rat, si, 'r');
	
        /* Init algebraic norms */
        tn_alg -= seconds ();
        survivors0 += init_alg_norms_bucket_region(alg_S, rat_S, i, 
						   cpoly, si);
        tn_alg += seconds ();
        /* Apply algebraic buckets */
        ttsm -= seconds();
        for (j = 0; j < si->nb_threads; ++j)
            apply_one_bucket(alg_S, alg_BA[j], i, si);
        ttsm += seconds();
        /* Sieve small algebraic primes */
        sieve_small_bucket_region(alg_S, i, lssd_alg, si, 'a');

        /* Factor survivors */
        ttf -= seconds ();
        reports += factor_survivors (rat_S, alg_S, i, rat_BA, alg_BA, 
		fb_rat, fb_alg, cpoly, si, &survivors1,
                &survivors2, &lsrsd_alg, &lsrsd_rat,
                report_sizes_a, report_sizes_r,
                output);
        ttf += seconds ();
      }

    /* clear */
    free(alg_S);
    free(rat_S);
    clear_small_sieve(lssd_alg);
    clear_small_sieve(lssd_rat);
    clear_small_sieve(lsrsd_alg);
    clear_small_sieve(lsrsd_rat);

    process_bucket_region_report_t *rep;
    rep = (process_bucket_region_report_t *)
        malloc(sizeof(process_bucket_region_report_t));
    rep->reports = reports;
    rep->survivors0 = survivors0;
    rep->survivors1 = survivors1;
    rep->survivors2 = survivors2;
    rep->tn_rat = tn_rat;
    rep->tn_alg = tn_alg;
    rep->ttsm = ttsm;
    rep->ttf = ttf;
    for (i = 0; i < 256; i++) {
        rep->report_sizes_a[i] = report_sizes_a[i];
        rep->report_sizes_r[i] = report_sizes_r[i];
    }
    return rep;
}

process_bucket_region_report_t
process_bucket_region_mt(bucket_array_t *alg_BA, bucket_array_t *rat_BA,
        small_sieve_data_t *ssd_rat, small_sieve_data_t *ssd_alg, 
        factorbase_degn_t *fb_rat, factorbase_degn_t *fb_alg,
        cado_poly_ptr cpoly, sieve_info_t *si, FILE *output, int verbose) {
    int i;
    process_bucket_region_report_t report;
    report.reports = 0;
    report.survivors0 = 0;
    report.survivors1 = 0;
    report.survivors2 = 0;
    report.tn_rat = 0.0;
    report.tn_alg = 0.0;
    report.ttsm = 0.0;
    report.ttf = 0.0;
    for (i = 0; i < 256; i++) {
        report.report_sizes_a[i] = 0;
        report.report_sizes_r[i] = 0;
    }

    pthread_t *th;
    th = malloc(si->nb_threads*sizeof(pthread_t));
    ASSERT_ALWAYS(th != NULL);
    process_bucket_region_arg_t *my_arg;
    my_arg = malloc(si->nb_threads*sizeof(process_bucket_region_arg_t));
    ASSERT_ALWAYS(my_arg != NULL);
    for (i = 0; i < si->nb_threads; ++i) {
        int ret;
        my_arg[i].alg_BA = alg_BA;
        my_arg[i].fb_alg = fb_alg;
        my_arg[i].rat_BA = rat_BA;
        my_arg[i].fb_rat = fb_rat;
        my_arg[i].ssd_rat = ssd_rat;
        my_arg[i].ssd_alg = ssd_alg;
        my_arg[i].cpoly = cpoly;
        my_arg[i].si = si;
        my_arg[i].id = i;
        my_arg[i].output = output;
        my_arg[i].verbose = verbose;
        ret = pthread_create(&(th[i]), NULL, 
	      (void * (*)(void *)) &process_regions_one_thread,
                (void *)(&(my_arg[i])));
        ASSERT_ALWAYS(!ret);
    }
    for (i = 0; i < si->nb_threads; ++i) {
        int ret;
        process_bucket_region_report_t *rep;
        ret = pthread_join(th[i], (void **)(&rep));
        ASSERT_ALWAYS(!ret);
        report.reports += rep->reports;
        report.survivors0 += rep->survivors0;
        report.survivors1 += rep->survivors1;
        report.survivors2 += rep->survivors2;
        report.tn_rat += rep->tn_rat;
        report.tn_alg += rep->tn_alg;
        report.ttsm += rep->ttsm;
        report.ttf += rep->ttf;
        int j;
        for (j = 0; j < 256; j++) {
            report.report_sizes_a[j] += rep->report_sizes_a[j];
            report.report_sizes_r[j] += rep->report_sizes_r[j];
        }
        free(rep);
    }

    if (verbose)
    {
        fprintf (output, "# %lu survivors after rational sieve,", report.survivors0);
        fprintf (output, " %lu survivors after algebraic sieve, ", report.survivors1);
        fprintf (output, "coprime: %lu\n", report.survivors2);
    }
    fprintf (output, "# %d relation(s) for (%" PRIu64 ",%" PRIu64
            ")\n", report.reports, si->q, si->rho);
    free(my_arg);
    free(th);
    return report;
}


/*************************** main program ************************************/

static void
usage (const char *argv0, const char * missing)
{
  fprintf (stderr, "Usage: %s [-I I] -poly xxx.poly -fb xxx.roots -q0 q0 [-q1 q1] [-rho rho]\n",
           argv0);
  fprintf (stderr, "          -I i            sieving region has side 2^i [default %u]\n", DEFAULT_I);
  fprintf (stderr, "          -poly xxx.poly  use polynomial xxx.poly\n");
  fprintf (stderr, "          -fb xxx.roots   use factor base xxx.roots\n");
  fprintf (stderr, "          -q0 nnn         left bound of special-q range\n");
  fprintf (stderr, "          -q1 nnn         right bound of special-q range\n");
  fprintf (stderr, "          -rho r          sieve only algebraic root r mod q0\n");
  fprintf (stderr, "          -tdthresh nnn   trial-divide primes p/r <= nnn (r=number of roots)\n");
  fprintf (stderr, "          -bkthresh nnn   bucket-sieve primes p >= nnn\n");
  fprintf (stderr, "          -rlim     nnn   rational factor base bound nnn\n");
  fprintf (stderr, "          -alim     nnn   algebraic factor base bound nnn\n");
  fprintf (stderr, "          -lpbr     nnn   rational large prime bound 2^nnn\n");
  fprintf (stderr, "          -lpba     nnn   algebraic large prime bound 2^nnn\n");
  fprintf (stderr, "          -mfbr     nnn   rational cofactor bound 2^nnn\n");
  fprintf (stderr, "          -mfba     nnn   algebraic cofactor bound 2^nnn\n");
  fprintf (stderr, "          -rlambda  nnn   rational lambda value is nnn\n");
  fprintf (stderr, "          -alambda  nnn   algebraic lambda value is nnn\n");
  fprintf (stderr, "          -S        xxx   skewness value is xxx\n");
  fprintf (stderr, "          -v              be verbose (print some sieving statistics)\n");
  fprintf (stderr, "          -out filename   write relations to filename instead of stdout\n");
  fprintf (stderr, "          -mt nnn   use nnn threads\n");
  fprintf (stderr, "          -ratq           use rational special-q\n");
  fprintf (stderr, "          The following are for benchs:\n");
  fprintf (stderr, "          -bench          activate bench mode\n");
  fprintf (stderr, "          -skfact   xxx   skip factor, default=1.01\n");
  if (missing) {
      fprintf(stderr, "\nError: missing parameter %s\n", missing);
  }
  exit (EXIT_FAILURE);
}

int
main (int argc0, char *argv0[])
{
    sieve_info_t si;
    const char *fbfilename = NULL;
    const char *polyfilename = NULL;
    cado_poly cpoly;
    double t0, tfb, tq, tn_rat, tn_alg, tts, ttsm, ttf;
    uint64_t q0 = 0, q1 = 0, rho = 0;
    uint64_t *roots;
    unsigned long nroots, tot_reports = 0;
    factorbase_degn_t * fb_alg, * fb_rat;
    int td_thresh = 1024; /* cost threshold trialdiv/resieving */
    int bucket_thresh = 0;
    int rpow_lim = 0, apow_lim = 0;
    int I = DEFAULT_I, i;
    int verbose = 0;
    int ratq = 0;
    unsigned long sq = 0;
    double totJ = 0.0;
    unsigned long report_sizes_a[256], report_sizes_r[256];
    /* following command-line values override those in the polynomial file */
    FILE *output;
    const char *outputname = NULL;
    int argc = argc0;
    char **argv = argv0;
    double max_full = 0.;
    int nb_threads = 1;
    int bench = 0;
    double skip_factor = 1.01;  /* next_q = q*skip_factor in bench mode */
    long bench_tot_rep = 0;
    double bench_tot_time = 0.0;

    param_list pl;
    param_list_init(pl);
    cado_poly_init (cpoly);
    param_list_configure_knob(pl, "-v", &verbose);
    param_list_configure_knob(pl, "-ratq", &ratq);
    param_list_configure_knob(pl, "-bench", &bench);
    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0[0],NULL);
    }

    polyfilename = param_list_lookup_string(pl, "poly");
    if (polyfilename) param_list_read_file(pl, polyfilename);
    fbfilename = param_list_lookup_string(pl, "fb");

    if (!cado_poly_set_plist (cpoly, pl)) {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
    }

    param_list_parse_int(pl, "mt", &nb_threads);
    param_list_parse_int(pl, "I", &I);
    param_list_parse_uint64(pl, "q0", &q0);
    param_list_parse_uint64(pl, "q1", &q1);
    param_list_parse_uint64(pl, "rho", &rho);
    param_list_parse_int(pl, "tdthresh", &td_thresh);
    param_list_parse_int(pl, "bkthresh", &bucket_thresh);
    param_list_parse_int(pl, "rpowlim", &rpow_lim);
    param_list_parse_int(pl, "apowlim", &apow_lim);
    param_list_parse_double(pl, "S", &cpoly->skew);
    param_list_parse_double(pl, "skfact", &skip_factor);
    int ok = 1;
    ok = ok && param_list_parse_ulong(pl, "rlim", &cpoly->rlim);
    ok = ok && param_list_parse_ulong(pl, "alim", &cpoly->alim);
    ok = ok && param_list_parse_int(pl, "lpbr", &cpoly->lpbr);
    ok = ok && param_list_parse_int(pl, "lpba", &cpoly->lpba);
    ok = ok && param_list_parse_int(pl, "mfbr", &cpoly->mfbr);
    ok = ok && param_list_parse_int(pl, "mfba", &cpoly->mfba);
    ok = ok && param_list_parse_double(pl, "rlambda", &cpoly->rlambda);
    ok = ok && param_list_parse_double(pl, "alambda", &cpoly->alambda);

    if (!ok) {
        fprintf(stderr, "Some parameters are missing among *lim lpb* mfb* *lambda\n");
        usage(argv0[0],NULL);
    }

    outputname = param_list_lookup_string(pl, "out");
    /* Init output file */
    if (outputname == NULL)
      output = stdout;
    else
      {
	output = gzip_open (outputname, "w");
	if (output == NULL)
	  {
	    fprintf (stderr, "Could not open %s for writing\n", outputname);
	    exit (EXIT_FAILURE);
	  }
      }

    /* Print command line to output */
    fprintf (output, "# %s.r%s", argv0[0], CADO_REV);
    for (i = 1; i < argc0; i++)
      fprintf (output, " %s", argv0[i]);
    fprintf (output, "\n");

    if (fbfilename == NULL) usage(argv0[0], "fb");
    if (q0 == 0) usage(argv0[0], "q0");

    /* if -rho is given, we sieve only for q0, thus -q1 is not allowed */
    if (rho != 0 && q1 != 0)
      {
        fprintf (stderr, "Error, -q1 and -rho are mutually exclusive\n");
        exit (EXIT_FAILURE);
      }

    /* if -q1 is not given, sieve only for q0 */
    if (q1 == 0)
      q1 = q0 + 1;

    /* check that q1 fits into an unsigned long */
    if (q1 > (uint64_t) ULONG_MAX)
      {
        fprintf (stderr, "Error, q1=%" PRIu64 " exceeds ULONG_MAX\n", q1);
        exit (EXIT_FAILURE);
      }

    fprintf (output, "# Sieving parameters: rlim=%lu alim=%lu lpbr=%d lpba=%d\n",
             cpoly->rlim, cpoly->alim, cpoly->lpbr, cpoly->lpba);
    fprintf (output, "#                     mfbr=%d mfba=%d rlambda=%1.1f alambda=%1.1f\n",
             cpoly->mfbr, cpoly->mfba, cpoly->rlambda, cpoly->alambda);
    fprintf (output, "#                     skewness=%1.1f\n",
             cpoly->skew);

    if (cpoly->skew <= 0.0)
      {
        fprintf (stderr, "Error, please provide a positive skewness\n");
        exit (EXIT_FAILURE);
      }

    /* check that nb_threads (-mt nnn) is positive */
    if (nb_threads <= 0)
      {
        fprintf (stderr, "Error, please provide a positive number of threads\n");
        exit (EXIT_FAILURE);
      }

    si.bench=bench;

    /* this does not depend on the special-q */
    si.ratq = ratq;
    sieve_info_init (&si, cpoly, I, q0, bucket_thresh, output, nb_threads);

    factorbase_degn_t **fb_alg_mt;
    fbprime_t *alg_pmax;
    fb_alg_mt = (factorbase_degn_t **)malloc(si.nb_threads*
            sizeof(factorbase_degn_t *));
    alg_pmax = (fbprime_t *)malloc(si.nb_threads*sizeof(fbprime_t));

    /* Read algebraic factor base */
    {
      fbprime_t *leading_div;
      tfb = seconds ();
      leading_div = factor_small (cpoly->f[cpoly->degree], cpoly->alim);
      fb_alg = fb_read_addproj (fbfilename, si.scale_alg * LOG_SCALE, 0,
				leading_div);
      ASSERT_ALWAYS(fb_alg != NULL);
      tfb = seconds () - tfb;
      fprintf (output, 
               "# Reading algebraic factor base of %zuMb took %1.1fs\n", 
               fb_size (fb_alg) >> 20, tfb);
      free (leading_div);
      
      factorbase_degn_t *fb = fb_alg;
      while (fb->p != FB_END && fb->p < (fbprime_t) si.bucket_thresh)
          fb = fb_next (fb); 
      ASSERT (fb->p != FB_END);
      for (i = 0; i < si.nb_threads; ++i) {
          fb_alg_mt[i] = fb;
          fb = fb_next(fb);
          ASSERT (fb->p != FB_END);
          alg_pmax[i] = FBPRIME_MAX;
      }
    }

    /* Prepare rational factor base */
    tfb = seconds ();
    if (rpow_lim >= si.bucket_thresh)
      {
        rpow_lim = si.bucket_thresh - 1;
        printf ("# rpowthresh reduced to %d\n", rpow_lim);
      }
    fb_rat = fb_make_linear ((const mpz_t *) cpoly->g, (fbprime_t) cpoly->rlim,
                             rpow_lim, si.scale_rat * LOG_SCALE, 
                             verbose, 1, output);
    tfb = seconds () - tfb;
    fprintf (output, "# Creating rational factor base of %zuMb took %1.1fs\n",
             fb_size (fb_rat) >> 20, tfb);
    
    factorbase_degn_t **fb_rat_mt;
    fbprime_t *rat_pmax;
    fb_rat_mt = (factorbase_degn_t **)malloc(si.nb_threads*
            sizeof(factorbase_degn_t *));
    rat_pmax = (fbprime_t *)malloc(si.nb_threads*sizeof(fbprime_t));

    factorbase_degn_t *fb = fb_rat;
    while (fb->p != FB_END && fb->p < (fbprime_t) si.bucket_thresh)
        fb = fb_next (fb); 
    ASSERT (fb->p != FB_END);
    for (i = 0; i < si.nb_threads; ++i) {
        fb_rat_mt[i] = fb;
        fb = fb_next(fb);
        ASSERT (fb->p != FB_END);
        rat_pmax[i] = FBPRIME_MAX;
    }

    /* counting factor base elements */
    {
        unsigned long na = 0, nr = 0;
        factorbase_degn_t *fb;
        fb = fb_alg;
        while (fb->p != FB_END) {
            na += fb->nr_roots; 
            fb = fb_next (fb);
        }
        fb = fb_rat;
        while (fb->p != FB_END) {
            nr++; 
            fb = fb_next (fb);
        }
        fprintf (output, "# Number of primes in alg factor base = %lu\n", na);
        fprintf (output, "# Number of primes in rat factor base = %lu\n", nr);
    }

    init_rat_norms (&si);
    init_alg_norms (&si);

    /* Init refactoring stuff */
    int skip2;
    si.trialdiv_primes_rat = fb_extract_bycost (fb_rat, si.bucket_thresh,
                                                td_thresh);
    si.trialdiv_primes_alg = fb_extract_bycost (fb_alg, si.bucket_thresh,
                                                td_thresh);
    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */
    for (i = 0; si.trialdiv_primes_alg[i] != FB_END; i++);
    skip2 = (i > 0 && si.trialdiv_primes_alg[0] == 2) ? 1 : 0;
    si.trialdiv_data_alg = trialdiv_init (si.trialdiv_primes_alg + skip2,
                                          i - skip2);
    for (i = 0; si.trialdiv_primes_rat[i] != FB_END; i++);
    skip2 = (i > 0 && si.trialdiv_primes_rat[0] == 2) ? 1 : 0;
    si.trialdiv_data_rat = trialdiv_init (si.trialdiv_primes_rat + skip2,
                                          i - skip2);
    si.strategy = facul_make_strategy (15, MIN(cpoly->rlim, cpoly->alim),
                                       1UL << MIN(cpoly->lpbr, cpoly->lpba));

    for (i = 0; i < 256; i++)
      {
        report_sizes_a[i] = 0;
        report_sizes_r[i] = 0;
      }

    /* special q (and root rho) */
    roots = (uint64_t *) malloc (cpoly->degree * sizeof (uint64_t));
    q0 --; /* so that nextprime gives q0 if q0 is prime */
    nroots = 0;
    tot_reports = 0;

    tn_rat = tn_alg = tts = ttsm = ttf = 0.0;
    t0 = seconds ();
    fprintf (output, "#\n");
    while (q0 < q1)
      {
        double t_bench = seconds();
        int rep_bench = 0;
        while (nroots == 0) /* go to next prime and generate roots */
          {
            q0 = uint64_nextprime (q0);
            if (q0 >= q1)
              goto end;  // breaks two whiles.
            si.q = q0;
            if (si.ratq)
                nroots = poly_roots_uint64 (roots, cpoly->g, 1, q0);
            else
                nroots = poly_roots_uint64 (roots, cpoly->f, cpoly->degree, q0);
            if (nroots > 0)
              {
                fprintf (output, "### q=%" PRIu64 ": root%s", q0,
                         (nroots == 1) ? "" : "s");
                for (i = 1; i <= (int) nroots; i++)
                  fprintf (output, " %" PRIu64, roots[nroots-i]);
                fprintf (output, "\n");
              }
          }
        tq = seconds ();

        /* computes a0, b0, a1, b1 from q, rho, and the skewness */
        si.rho = roots[--nroots];
        if (rho != 0 && si.rho != rho) /* if -rho, wait for wanted root */
          continue;
        SkewGauss (&si, cpoly->skew);
        /* FIXME: maybe we can discard some special q's if a1/a0 is too large,
           see http://www.mersenneforum.org/showthread.php?p=130478 */
        /* if (skewness (&si) > 32.0)
           continue; */

        fprintf (output, "# Sieving q=%" PRIu64 "; rho=%" PRIu64
                 "; a0=%d; b0=%d; a1=%d; b1=%d\n",
                 si.q, si.rho, si.a0, si.b0, si.a1, si.b1);
        sq ++;
        /* precompute the skewed polynomials of f(x) and g(x) */
        fij_from_f (si.fij, &si, cpoly->f, cpoly->degree);
        fij_from_f (si.gij, &si, cpoly->g, 1);

        /* checks the value of J, updates floating-point fij(x) */
        sieve_info_update (&si, verbose, output);
        totJ += (double) si.J;

        /* Allocate alg buckets */
        ttsm -= seconds();
        bucket_array_t alg_BA[si.nb_threads];
        for (i = 0; i < si.nb_threads; ++i)
            alg_BA[i] = init_bucket_array(si.nb_buckets, si.bucket_limit / si.nb_threads);

        /* Allocate rational buckets */
        bucket_array_t rat_BA[si.nb_threads];
        for (i = 0; i < si.nb_threads; ++i)
            rat_BA[i] = init_bucket_array(si.nb_buckets, si.bucket_limit / si.nb_threads);


        /* Fill in rat and alg buckets */
        fill_in_buckets_mt(alg_BA, fb_alg_mt, alg_pmax,
                rat_BA, fb_rat_mt, rat_pmax,
                &si);
        for (i = 0; i < si.nb_threads; ++i) {
            if (buckets_max_full (alg_BA[i]) > max_full)
                max_full = buckets_max_full (alg_BA[i]);
            if (buckets_max_full (rat_BA[i]) > max_full)
                max_full = buckets_max_full (rat_BA[i]);
        }
        if (max_full >= 1.0) 
            fprintf(stderr, "maxfull=%f\n", max_full);
        ASSERT (max_full < 1.);

        ttsm += seconds();

        /* Initialize data for sieving small primes */
        small_sieve_data_t ssd_alg, ssd_rat;
        init_small_sieve(&ssd_rat, fb_rat, &si, 'r', output);
        init_small_sieve(&ssd_alg, fb_alg, &si, 'a', output);

        /* Process bucket regions in parallel */
        {
            process_bucket_region_report_t rep;
            rep = process_bucket_region_mt(alg_BA, rat_BA,
                &ssd_rat, &ssd_alg, fb_rat, fb_alg, cpoly, &si, output, verbose);
            tot_reports += rep.reports;
            tn_rat      += rep.tn_rat;
            tn_alg      += rep.tn_alg;
            ttsm        += rep.ttsm;
            ttf         += rep.ttf;
            rep_bench = rep.reports;
            for (i = 0; i < 256; ++i) {
                report_sizes_a[i] += rep.report_sizes_a[i];
                report_sizes_r[i] += rep.report_sizes_r[i];
            }
        }
	
        clear_small_sieve(ssd_rat);
        clear_small_sieve(ssd_alg);
        for (i = 0; i < si.nb_threads; ++i) {
            clear_bucket_array(alg_BA[i]);
            clear_bucket_array(rat_BA[i]);
        }
        if (bench) {
            uint64_t newq0 = (uint64_t) (skip_factor*((double) q0));
            uint64_t savq0 = q0;
            // print some estimates for special-q's between q0 and the next
            int nb_q = 1;
            do {
                q0 = uint64_nextprime (q0);
                nb_q ++;
            } while (q0 < newq0);
            q0 = newq0;
            nroots=0;
            t_bench = seconds() - t_bench;
            fprintf(output,
              "# Stats for q=%" PRIu64 ": %d reports in %1.1f s\n",
              savq0, rep_bench, t0);
            fprintf(output,
              "# Estimates for next %d q's: %d reports in %1.0f s, %1.2f s/r\n",
              nb_q, nb_q*rep_bench, t0*nb_q, t0/((double)rep_bench));
            bench_tot_time += t0*nb_q;
            bench_tot_rep += nb_q*rep_bench;
            fprintf(output, "# Cumulative (estimated): %lu reports in %1.0f s, %1.2f s/r\n",
                    bench_tot_rep, bench_tot_time,
		    (double) bench_tot_time / (double) bench_tot_rep);
        }
      } // end of loop over special q ideals.

 end:
    t0 = seconds () - t0;
    fprintf (output, "# Average J=%1.0f for %lu special-q's, max bucket fill %f\n",
             totJ / (double) sq, sq, max_full);
    tts = t0 - (tn_rat + tn_alg + ttf);
    if (verbose)
      facul_print_stats (output);
    if (verbose)
      {
        fprintf (output, "# Histogram of sieve report values that led to "
                 "relations:\n# Algebraic side: ");
        for (i = 0; i < 256; i++)
          if (report_sizes_a[i] > 0)
            fprintf (output, "%d:%lu ", i, report_sizes_a[i]);
        fprintf (output, "\n# Rational side: ");
        for (i = 0; i < 256; i++)
          if (report_sizes_r[i] > 0)
            fprintf (output, "%d:%lu ", i, report_sizes_r[i]);
        fprintf (output, "\n");
      }
    if (si.nb_threads > 1) 
        fprintf (output, "# Total cpu time %1.1fs [precise timings available only for mono-thread]\n", t0);
    else
        fprintf (output, "# Total time %1.1fs [norm %1.2f+%1.1f, sieving %1.1f"
            " (%1.1f + %1.1f),"
             " factor %1.1f]\n", t0, tn_rat, tn_alg, tts, ttsm, tts-ttsm, ttf);
    fprintf (output, "# Total %lu reports [%1.3fs/r, %1.1fr/sq]\n",
             tot_reports, t0 / (double) tot_reports,
             (double) tot_reports / (double) sq);
    if (bench) {
        fprintf(output, "# Total (estimated): %lu reports in %1.1f s\n",
                bench_tot_rep, bench_tot_time);
    }

    facul_clear_strategy (si.strategy);
    si.strategy = NULL;
    trialdiv_clear (si.trialdiv_data_alg);
    trialdiv_clear (si.trialdiv_data_rat);
    free (si.trialdiv_primes_alg);
    free (si.trialdiv_primes_rat);
    free (fb_alg);
    free (fb_alg_mt);
    free (fb_rat);
    free (fb_rat_mt);
    free (alg_pmax);
    free (rat_pmax);
    sieve_info_clear (&si);
    cado_poly_clear (cpoly);
    free (roots);

    if (outputname != NULL)
      gzip_close (output, "");

    param_list_clear(pl);

    return 0;
}
