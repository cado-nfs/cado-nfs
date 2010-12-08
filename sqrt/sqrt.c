#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <inttypes.h>
#include <math.h> /* for log */
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "utils.h"
#include "purgedfile.h"

#define DEBUG 0
#define MAPLE 0

static int verbose = 0;

/* Although the functions in plain_poly are not readily available in the
 * publicized interface of utils.h, it's ok to use them if we explicitly
 * include the corresponding header.
 */
#include "plain_poly.h"

/********** RATSQRT **********/

/* Returns memory usage, in KB 
 * This is the VmSize field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
static long Memusage() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmSize: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

/* same as above, for resident memory (column RES of top) */
static long Memusage2() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmRSS: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

/* Returns peak memory usage, in KB 
 * This is the VmPeak field in the status file of /proc/pid/ dir
 * This is highly non portable.
 * Return -1 in case of failure.
 */
static long PeakMemusage() {
  pid_t pid = getpid();

  char str[1024];
  char *truc;
  snprintf(str, 1024, "/proc/%d/status", pid);

  FILE *file;
  file = fopen(str, "r");
  if (file == NULL)
    return -1;

  long mem;
  for(;;) {
    truc = fgets(str, 1023, file);
    if (truc == NULL) {
      fclose(file);
      return -1;
    }
    int ret = sscanf(str, "VmPeak: %ld", &mem);
    if (ret == 1) {
      fclose(file);
      return mem;
    }
  }
}

static void
my_mpz_mul (mpz_t a, mpz_t b, mpz_t c)
{
  int large, st = 0;

  large = mpz_size (b) + mpz_size (c) >= 5000000;
  if (large)
    {
      fprintf (stderr, "[multiplying %zu*%zu limbs: ",
               mpz_size (b), mpz_size (c));
      fflush (stderr);
      st = cputime ();
    }
  mpz_mul (a, b, c);
  mpz_realloc2 (c, 0);
  if (large)
    {
      fprintf (stderr, "%dms]\n", cputime () - st);
      fflush (stderr);
    }
}

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 + 
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
static mpz_t*
accumulate_fast (mpz_t *prd, mpz_t a, unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  my_mpz_mul (prd[0], prd[0], a);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (mpz_t*) realloc (prd, *lprd * sizeof (mpz_t));
          mpz_init_set_ui (prd[i + 1], 1);
        }
      my_mpz_mul (prd[i + 1], prd[i + 1], prd[i]);
      mpz_set_ui (prd[i], 1);
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
static void
accumulate_fast_end (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    my_mpz_mul (prd[0], prd[0], prd[i]);
}

static size_t
stats (mpz_t *prd, unsigned long lprd)
{
  unsigned long i;
  size_t s = 0;

  for (i = 0; i < lprd; i++)
    s += mpz_size (prd[i]);
  return s;
}

int 
calculateSqrtRat (const char *prefix, int numdep, cado_poly pol, mpz_t Np)
{
  char depname[200];
  char ratname[200];
  snprintf(depname, sizeof(depname), "%s.%03d", prefix, numdep);
  snprintf(ratname, sizeof(ratname), "%s.rat.%03d", prefix, numdep);
  FILE *depfile = NULL;
  FILE *ratfile;
  //int sign;
  long a, b;
  int ret;
  unsigned long ab_pairs = 0, line_number, freerels = 0;
  mpz_t v, *prd;
  unsigned long lprd; /* number of elements in prd[] */
  unsigned long nprd; /* number of accumulated products in prd[] */
  unsigned long res, peakres = 0;

  if (pol->degreeg != 1)
    {
      fprintf (stderr, "Error, calculateSqrtRat called with non-linear polynomial\n");
      exit (EXIT_FAILURE);
    }

#ifdef __MPIR_VERSION
  fprintf (stderr, "Using MPIR %s\n", mpir_version);
#else
  fprintf (stderr, "Using GMP %s\n", gmp_version);
#endif

  mpz_init (v);

  lprd = 1;
  nprd = 0;
  prd = (mpz_t*) malloc (lprd * sizeof (mpz_t));
  mpz_init_set_ui (prd[0], 1);

  depfile = fopen (depname, "r");
  /*if(!depfile) {
    fprintf (stderr, "Error: file %s not exist\n", depname);
    exit(1);
  }*/
  ASSERT_ALWAYS(depfile != NULL);
    
  line_number = 2;
  for (;;)
    {
      ret = fscanf (depfile, "%ld %ld\n", &a, &b);

      if (ret != 2)
        {
          fprintf (stderr, "Invalid line %lu\n", line_number);
          break;
        }

      ab_pairs ++;
      line_number ++;

      if (ab_pairs % 1000000 == 0)
        {
          res = Memusage2 ();
          if (res > peakres)
            peakres = res;
            fprintf (stderr, "%lu pairs: size %zuMb, %dms, VIRT %luM (peak %luM), RES %luM (peak %luM)\n",
                       ab_pairs, stats (prd, lprd) >> 17, cputime (),
                       Memusage () >> 10, PeakMemusage () >> 10,
                       res >> 10, peakres >> 10);
        }

        if (b == 0)
          freerels ++;

        /* accumulate g1*a+g0*b */
        mpz_mul_si (v, pol->g[1], a);
        mpz_addmul_si (v, pol->g[0], b);

        prd = accumulate_fast (prd, v, &lprd, nprd++);
          
        if (feof (depfile))
          break;
      }
    fprintf (stderr, "%lu (a,b) pairs\n", line_number);

    fclose (depfile);

  fprintf (stderr, "Read %lu (a,b) pairs, including %lu free\n", ab_pairs,
           freerels);

  accumulate_fast_end (prd, lprd);

  /* we must divide by g1^ab_pairs: if the number of (a,b) pairs is odd, we
     multiply by g1, and divide by g1^(ab_pairs+1) */
  if (ab_pairs & 1)
    mpz_mul (prd[0], prd[0], pol->g[1]);

  fprintf (stderr, "Size of product = %zu bits\n", mpz_sizeinbase (prd[0], 2));

  if (mpz_sgn (prd[0]) < 0)
    {
      fprintf (stderr, "Error, product is negative: try another dependency\n");
      exit (1);
    }

  fprintf (stderr, "Starting square root at %dms\n", cputime ());

  /* since we know we have a square, take the square root */
  mpz_sqrtrem (prd[0], v, prd[0]);
  
  fprintf (stderr, "Computed square root at %dms\n", cputime ());

  if (mpz_cmp_ui (v, 0) != 0)
    {
      unsigned long p = 2, e;
      mpz_t pp;

      mpz_init (pp);
      fprintf (stderr, "Error, square root remainder is not zero\n");
      /* reconstruct the initial value of prd[0] to debug */
      mpz_mul (prd[0], prd[0], prd[0]);
      mpz_add (prd[0], prd[0], v);
      while (mpz_cmp_ui (prd[0], 1) > 0)
        {
          e = 0;
          if (verbose)
            printf ("Removing p=%lu:", p);
          mpz_set_ui (pp, p);
          e = mpz_remove (prd[0], prd[0], pp);
          if (verbose)
            printf (" exponent=%lu, remaining %zu bits\n", e,
                    mpz_sizeinbase (prd[0], 2));
          if ((e % 2) != 0)
            {
              fprintf (stderr, "Prime %lu appears to odd power %lu\n", p, e);
              if (verbose)
                break;
            }
          p = getprime (p);
        }
      mpz_clear (pp);
      p = getprime (0);
      exit (1);
    }

  mpz_mod (prd[0], prd[0], Np);

  fprintf (stderr, "Reduced mod n at %dms\n", cputime ());

  /* now divide by g1^(ab_pairs/2) if ab_pairs is even, and g1^((ab_pairs+1)/2)
     if ab_pairs is odd */
  
  mpz_powm_ui (v, pol->g[1], (ab_pairs + 1) / 2, Np);
  fprintf (stderr, "Computed g1^(nab/2) mod n at %dms\n", cputime ());

  mpz_invert (v, v, Np);
  mpz_mul (prd[0], prd[0], v);
  mpz_mod (prd[0], prd[0], Np);
  
  ratfile = fopen (ratname, "w");
  gmp_fprintf (ratfile, "%Zd\n", prd[0]);
  fclose(ratfile);

  gmp_fprintf (stderr, "rational square root is %Zd\n", prd[0]);

  fprintf (stderr, "Rational square root time at %dms\n", cputime ());

  mpz_clear (prd[0]);

  mpz_clear (v);
  return 0;
}



/********** ALGSQRT **********/
static void
polymodF_from_ab(polymodF_t tmp, long a, unsigned long b) {
  tmp->v = 0;
  tmp->p->deg = (b != 0) ? 1 : 0;
  mpz_set_si (tmp->p->coeff[1], - (long) b);
  mpz_set_si (tmp->p->coeff[0], a);
}

/* Reduce the coefficients of R in [-m/2, m/2], which are assumed in [0, m[ */
static void
poly_mod_center (poly_t R, const mpz_t m)
{
  int i;
  mpz_t m_over_2;

  mpz_init (m_over_2);
  mpz_div_2exp (m_over_2, m, 2);
  for (i=0; i <= R->deg; i++)
    {
      ASSERT_ALWAYS(mpz_cmp_ui (R->coeff[i], 0) >= 0);
      ASSERT_ALWAYS(mpz_cmp (R->coeff[i], m) < 0);
      if (mpz_cmp (R->coeff[i], m_over_2) > 0)
        mpz_sub (R->coeff[i], R->coeff[i], m);
    }
  mpz_clear (m_over_2);
}

#if 0
/* Check whether the coefficients of R (that are given modulo m) are in
   fact genuine integers. We assume that x mod m is a genuine integer if
   x or |x-m| is less than m/10^6, i.e., the bit size of x or |x-m| is
   less than that of m minus 20.
   Assumes the coefficients x satisfy 0 <= x < m.
*/
static int
poly_integer_reconstruction (poly_t R, const mpz_t m)
{
  int i;
  size_t sizem = mpz_sizeinbase (m, 2), sizer;

  for (i=0; i <= R->deg; ++i)
    {
      sizer = mpz_sizeinbase (R->coeff[i], 2);
      if (sizer + 20 > sizem)
        {
          mpz_sub (R->coeff[i], R->coeff[i], m);
          sizer = mpz_sizeinbase (R->coeff[i], 2);
          if (sizer + 20 > sizem)
            return 0;
        }
    }
  return 1;
}
#endif

// compute res := sqrt(a) in Fp[x]/f(x)
static void
TonelliShanks (poly_t res, const poly_t a, const poly_t F, unsigned long p)
{
  int d = F->deg;
  mpz_t q;
  poly_t delta;  // a non quadratic residue
  poly_t auxpol;
  mpz_t aux;
  mpz_t t;
  int s;
  mpz_t myp;

  mpz_init_set_ui(myp, p);

  mpz_init(aux);
  mpz_init(q);
  poly_alloc(auxpol, d);
  mpz_ui_pow_ui(q, p, (unsigned long)d);

  // compute aux = (q-1)/2
  // and (s,t) s.t.  q-1 = 2^s*t
  mpz_sub_ui(aux, q, 1);
  mpz_divexact_ui(aux, aux, 2);
  mpz_init_set(t, aux);
  s = 1;
  while (mpz_divisible_2exp_p(t, 1)) {
    s++;
    mpz_divexact_ui(t, t, 2);
  }
  // find a non quadratic residue delta
  {
    poly_alloc(delta, d);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    do {
      int i;
      // pick a random delta
      for (i = 0; i < d; ++i)
    mpz_urandomm(delta->coeff[i], state, myp);
      cleandeg(delta, d-1);
      // raise it to power (q-1)/2
      poly_power_mod_f_mod_ui(auxpol, delta, F, aux, p);
    } while ((auxpol->deg != 0) || (mpz_cmp_ui(auxpol->coeff[0], p-1)!= 0));
    gmp_randclear (state);
  }

  // follow the description of Crandall-Pomerance, page 94
  {
    poly_t A, D;
    mpz_t m;
    int i;
    poly_alloc(A, d);
    poly_alloc(D, d);
    mpz_init_set_ui(m, 0);
    poly_power_mod_f_mod_ui(A, a, F, t, p);
    poly_power_mod_f_mod_ui(D, delta, F, t, p);
    for (i = 0; i <= s-1; ++i) {
      poly_power_mod_f_mod_ui(auxpol, D, F, m, p);
      poly_mul_mod_f_mod_mpz(auxpol, auxpol, A, F, myp, NULL);
      mpz_ui_pow_ui(aux, 2, (s-1-i));
      poly_power_mod_f_mod_ui(auxpol, auxpol, F, aux, p);
      if ((auxpol->deg == 0) && (mpz_cmp_ui(auxpol->coeff[0], p-1)== 0))
    mpz_add_ui(m, m, 1UL<<i);
    }
    mpz_add_ui(t, t, 1);
    mpz_divexact_ui(t, t, 2);
    poly_power_mod_f_mod_ui(res, a, F, t, p);
    mpz_divexact_ui(m, m, 2);
    poly_power_mod_f_mod_ui(auxpol, D, F, m, p);

    poly_mul_mod_f_mod_mpz(res, res, auxpol, F, myp, NULL);
    poly_free(D);
    poly_free(A);
    mpz_clear(m);
  }

  poly_free(auxpol);
  poly_free(delta);
  mpz_clear(q);
  mpz_clear(aux);
  mpz_clear(myp);
  mpz_clear (t);
}

// res <- Sqrt(AA) mod F, using p-adic lifting, at prime p.
void
polymodF_sqrt (polymodF_t res, polymodF_t AA, poly_t F, unsigned long p)
{
  poly_t A, *P;
  int v;
  int d = F->deg;
  int k, lk, target_k, logk, K[32];
  size_t target_size; /* target bit size for Hensel lifting */

  /* The size of the coefficients of the square root of A should be about half
     the size of the coefficients of A. Here is an heuristic argument: let
     K = Q[x]/(f(x)) where f(x) is the algebraic polynomial. The square root
     r(x) might be considered as a random element of K: it is smooth, not far
     from an integer, but except that has no relationship with the coefficients
     of f(x). When we square r(x), we obtain a polynomial with coefficients
     twice as large, before reduction by f(x). The reduction modulo f(x)
     produces A(x), however that reduction should not decrease the size of
     the coefficients. */
  target_size = poly_sizeinbase (AA->p, AA->p->deg, 2);
  target_size = target_size / 2;
  target_size += target_size / 10;
  fprintf (stderr, "target_size=%lu\n", (unsigned long int) target_size);

  poly_alloc(A, d-1);
  // Clean up the mess with denominator: if it is an odd power of fd,
  // then multiply num and denom by fd to make it even.
  if (((AA->v)&1) == 0) {
    v = AA->v / 2;
    poly_copy(A, AA->p);
  } else {
    v = (1+AA->v) / 2;
    poly_mul_mpz(A, AA->p, F->coeff[d]);
  }

  // Now, we just have to take the square root of A (without denom) and
  // divide by fd^v.

  // Variables for the lifted values
  poly_t invsqrtA;
  // variables for A and F modulo pk
  poly_t a;
  poly_alloc(invsqrtA, d-1);
  poly_alloc(a, d-1);
  // variable for the current pk
  mpz_t pk, invpk;
  mpz_init (pk);
  mpz_init (invpk);

  /* Jason Papadopoulos's trick: since we will lift the square root of A to at
     most target_size bits, we can reduce A accordingly */
  double st = seconds ();
  target_k = (int) ((double) target_size * log ((double) 2) / log((double) p));
  mpz_ui_pow_ui (pk, p, target_k);
  while (mpz_sizeinbase (pk, 2) <= target_size)
    {
      mpz_mul_ui (pk, pk, p);
      target_k ++;
    }
  poly_reduce_mod_mpz (A, A, pk);
  for (k = target_k, logk = 0; k > 1; k = (k + 1) / 2, logk ++)
    K[logk] = k;
  K[logk] = 1;
  fprintf (stderr, "Reducing A mod p^%d took %2.2lf\n", target_k,
           seconds () - st);

  // Initialize things modulo p:
  mpz_set_ui (pk, p);
  k = 1; /* invariant: pk = p^k */
  lk = 0; /* k = 2^lk */
  st = seconds ();
  P = poly_base_modp_init (A, p, K, logk);
  fprintf (stderr, "poly_base_modp_init took %2.2lf\n", seconds () - st);

  poly_copy (a, P[0]);

  // First compute the inverse square root modulo p
  {
    mpz_t q, aux;
    mpz_init(q);
    mpz_init(aux);
    mpz_ui_pow_ui(q, p, (unsigned long)d);

#if 0
    // compute (q-2)(q+1)/4   (assume q == 3 mod 4, here !!!!!)
    // where q = p^d, the size of the finite field extension.
    // since we work mod q-1, this gives (3*q-5)/4
    mpz_mul_ui(aux, q, 3);
    mpz_sub_ui(aux, aux, 5);
    mpz_divexact_ui(aux, aux, 4);               // aux := (3q-5)/4
    poly_power_mod_f_mod_ui(invsqrtA, a, F, aux, p);
#else
    TonelliShanks(invsqrtA, a, F, p);
    mpz_sub_ui(aux, q, 2);
    poly_power_mod_f_mod_ui(invsqrtA, invsqrtA, F, aux, p);
#endif

    mpz_clear(aux);
    mpz_clear(q);
  }

  // Now, the lift begins
  // When entering the loop, invsqrtA contains the inverse square root
  // of A computed modulo pk.

  poly_t tmp, tmp2;
  poly_alloc(tmp, 2*d-1);
  poly_alloc(tmp2, 2*d-1);
  do {
    double st;

    if (mpz_sizeinbase (pk, 2) > target_size)
      {
        fprintf (stderr, "Failed to reconstruct an integer polynomial\n");
        printf ("Failed\n");
        exit (1);
      }

    /* invariant: invsqrtA = 1/sqrt(A) bmod p^k */

    st = seconds ();
    poly_base_modp_lift (a, P, ++lk, pk);
    st = seconds () - st;

    /* invariant: k = K[logk] */
    ASSERT_ALWAYS(k == K[logk]);

    mpz_set (invpk, pk);
    mpz_mul (pk, pk, pk);   // double the current precision
    logk --;
    if (K[logk] & 1)
      {
        mpz_div_ui (pk, pk, p);
        k --;
      }
    k = K[logk];
    barrett_init (invpk, pk); /* FIXME: we could lift 1/p^k also */
    fprintf (stderr, "start lifting mod p^%d (%lu bits) at %2.2lf\n",
             k, (unsigned long int) mpz_sizeinbase (pk, 2), seconds ());
#ifdef VERBOSE
    fprintf (stderr, "   poly_base_modp_lift took %2.2lf\n", st);
#endif

    // now, do the Newton operation x <- 1/2(3*x-a*x^3)
    st = seconds ();
    poly_sqr_mod_f_mod_mpz(tmp, invsqrtA, F, pk, invpk); /* tmp = invsqrtA^2 */
#ifdef VERBOSE
    fprintf (stderr, "   poly_sqr_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif

    /* Faster version which computes x <- x + x/2*(1-a*x^2).
       However I don't see how to use the fact that the coefficients
       if 1-a*x^2 are divisible by p^(k/2). */
    st = seconds ();
    poly_mul_mod_f_mod_mpz (tmp, tmp, a, F, pk, invpk); /* tmp=a*invsqrtA^2 */
#ifdef VERBOSE
    fprintf (stderr, "   poly_mul_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif
    poly_sub_ui (tmp, 1); /* a*invsqrtA^2-1 */
    poly_div_2_mod_mpz (tmp, tmp, pk); /* (a*invsqrtA^2-1)/2 */
    st = seconds ();
    poly_mul_mod_f_mod_mpz (tmp, tmp, invsqrtA, F, pk, invpk);
#ifdef VERBOSE
    fprintf (stderr, "   poly_mul_mod_f_mod_mpz took %2.2lf\n", seconds () - st);
#endif
    /* tmp = invsqrtA/2 * (a*invsqrtA^2-1) */
    poly_sub_mod_mpz (invsqrtA, invsqrtA, tmp, pk);

  } while (k < target_k);

  /* multiply by a to get an approximation of the square root */
  poly_mul_mod_f_mod_mpz (tmp, invsqrtA, a, F, pk, invpk);
  poly_mod_center (tmp, pk);

  poly_base_modp_clear (P);

  poly_copy(res->p, tmp);
  res->v = v;

  mpz_clear (pk);
  mpz_clear (invpk);
  poly_free(tmp);
  poly_free(tmp2);
  poly_free (A);
  poly_free (invsqrtA);
  poly_free (a);

  size_t sqrt_size = poly_sizeinbase (res->p, F->deg - 1, 2);
  fprintf (stderr, "maximal sqrt bit-size = %zu (%.0f%% of target size)\n",
          sqrt_size, 100.0 * (double) sqrt_size / target_size);
}

static unsigned long
FindSuitableModP (poly_t F, mpz_t N)
{
  unsigned long p = 2;
  int dF = F->deg;

  plain_poly_t fp;

  plain_poly_init (fp, dF);
  while (1)
    {
    int d;

    p = getprime (p);
    if (mpz_gcd_ui(NULL, N, p) != 1)
      continue;
    if (! plain_poly_fits (dF, p))
      {
        fprintf (stderr, "You are in trouble. Please contact the CADO support team at cado-nfs-commits@lists.gforge.inria.fr.\n");
        exit (1);
      }
    d = plain_poly_set_mod (fp, F->coeff, dF, p);
    if (d != dF)
      continue;
    if (plain_poly_is_irreducible (fp, p))
      break;
    }
  plain_poly_clear(fp);
  getprime (0);

  return p;
}

// Products are computed modulo the polynomial F.
polymodF_t*
accumulate_fast_F (polymodF_t *prd, const polymodF_t a, const poly_t F,
                 unsigned long *lprd, unsigned long nprd)
{
  unsigned long i;

  polymodF_mul (prd[0], prd[0], a, F);
  nprd ++;

  for (i = 0; nprd % THRESHOLD == 0; i++, nprd /= THRESHOLD)
    {
      /* need to access prd[i + 1], thus i+2 entries */
      if (i + 2 > *lprd)
        {
          lprd[0] ++;
          prd = (polymodF_t*) realloc (prd, *lprd * sizeof (polymodF_t));
      poly_alloc(prd[i+1]->p, F->deg);
          mpz_set_ui(prd[i + 1]->p->coeff[0], 1);
      prd[i+1]->p->deg = 0;
      prd[i+1]->v = 0;
        }
      polymodF_mul (prd[i+1], prd[i+1], prd[i], F);
      mpz_set_ui(prd[i]->p->coeff[0], 1);
      prd[i]->p->deg = 0;
      prd[i]->v = 0;
    }

  return prd;
}

/* prd[0] <- prd[0] * prd[1] * ... * prd[lprd-1] */
void
accumulate_fast_F_end (polymodF_t *prd, const poly_t F, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    polymodF_mul (prd[0], prd[0], prd[i], F);
}

/* side=0: consider the polynomial f
   side=1: consider the polynomial g
*/
int
calculateSqrtAlg (const char *prefix, int numdep, cado_poly pol, int side, 
        mpz_t Np)
{
  char depname[200];
  char algname[200];
  FILE *depfile = NULL;
  FILE *algfile;
  poly_t F;
  polymodF_t prd, tmp;
  long a;
  unsigned long b;
  unsigned long p;
  double t0 = seconds ();
  mpz_t algsqrt, aux;
  int i, deg;
  mpz_t *f;
  int nab = 0, nfree = 0;

  ASSERT_ALWAYS(side == 0 || side == 1);

  sprintf (depname, "%s.%03d", prefix, numdep);
  if (side == 0)
    sprintf (algname, "%s.alg.%03d", prefix, numdep);
  else
    sprintf (algname, "%s.rat.%03d", prefix, numdep);
  depfile = fopen (depname, "r");
  ASSERT_ALWAYS(depfile != NULL);

  deg = (side == 0) ? pol->degree : pol->degreeg;
  f = (side == 0) ? pol->f : pol->g;

  ASSERT_ALWAYS(deg > 1);
  
  // Init F to be the corresponding polynomial
  poly_alloc (F, deg);
  for (i = deg; i >= 0; --i)
    poly_setcoeff (F, i, f[i]);
  
  // Init prd to 1.
  poly_alloc (prd->p, deg);
  mpz_set_ui (prd->p->coeff[0], 1);
  prd->p->deg = 0;
  prd->v = 0;
  
  // Allocate tmp
  poly_alloc (tmp->p, 1);
  
  // Accumulate product
  #if 0
    // Naive version, without subproduct tree
    while(fscanf(depfile, "%ld %lu", &a, &b) != EOF){
      if(!(nab % 100000))
        fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n",nab,seconds());
      if((a == 0) && (b == 0))
        break;
      polymodF_from_ab (tmp, a, b);
      polymodF_mul (prd, prd, tmp, F);
      nab++;
    }
  #else
    // With a subproduct tree
    {
      polymodF_t *prd_tab;
      unsigned long lprd = 1; /* number of elements in prd_tab[] */
      unsigned long nprd = 0; /* number of accumulated products in prd_tab[] */
      prd_tab = (polymodF_t*) malloc (lprd * sizeof (polymodF_t));
      poly_alloc (prd_tab[0]->p, F->deg);
      mpz_set_ui (prd_tab[0]->p->coeff[0], 1);
      prd_tab[0]->p->deg = 0;
      prd_tab[0]->v = 0;
      while(fscanf(depfile, "%ld %lu", &a, &b) != EOF){
        if(!(nab % 100000))
    fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab, seconds ());
        if((a == 0) && (b == 0))
    break;
        polymodF_from_ab(tmp, a, b);
        prd_tab = accumulate_fast_F (prd_tab, tmp, F, &lprd, nprd++);
        nab++;
        if(b == 0)
      nfree++;
      }
      fprintf (stderr, "# Read %d including %d free relations\n", nab, nfree);
      ASSERT_ALWAYS ((nab & 1) == 0);
      ASSERT_ALWAYS ((nfree & 1) == 0);
      /* nfree being even is forced by a specific character column added
       * by character.c. The correspond assert should not fail.
       *
       * nab being even is kind of a mystery. There is a character column
       * that gives the sign of the rational side. It could well be that
       * with the parameters we usually use, it is negative for all
       * relations; that would force the overall number of relations to
       * be even. Another possibility is when f_d contains a big prime
       * number that does not occur anywhere else, so that the power of
       * this prime is controlled by the usual characters, and since f_d
       * is always present... 
       *
       * But: I wouldn't be surprised if the assert(even(nab)) fails.
       * Then, a patch would be:
       *    - remove the assert (!)
       *    - in the numerator of the big product, eliminate powers of
       *       f_d that divides all coefficients.
       *    - this should finally give an even power of f_d in the
       *      denominator, and the algorithm can continue.
       */
      accumulate_fast_F_end (prd_tab, F, lprd);
      fclose(depfile);
  
      poly_copy(prd->p, prd_tab[0]->p);
      prd->v = prd_tab[0]->v;
      for (i = 0; i < (long)lprd; ++i)
        poly_free(prd_tab[i]->p);
      free(prd_tab);
    }
  #endif
  
    fprintf(stderr, "Finished accumulating the product at %2.2lf\n", seconds());
    fprintf(stderr, "nab = %d, nfree = %d, v = %d\n", nab, nfree, prd->v);
    fprintf (stderr, "maximal polynomial bit-size = %lu\n",
             (unsigned long) poly_sizeinbase (prd->p, deg - 1, 2));
  
    p = FindSuitableModP(F, Np);
    fprintf(stderr, "Using p=%lu for lifting\n", p);
  
    double tm = seconds();
    polymodF_sqrt (prd, prd, F, p);
    fprintf (stderr, "Square root lifted in %2.2lf\n", seconds()-tm);
  
    mpz_init(algsqrt);
    mpz_init(aux);
    poly_eval_mod_mpz(algsqrt, prd->p, pol->m, Np);
    mpz_invert(aux, F->coeff[F->deg], Np);  // 1/fd mod n
    mpz_powm_ui(aux, aux, prd->v, Np);      // 1/fd^v mod n
    mpz_mul(algsqrt, algsqrt, aux);
    mpz_mod(algsqrt, algsqrt, Np);
  
    algfile = fopen (algname, "w");
    gmp_fprintf (algfile, "%Zd\n", algsqrt);
    fclose(algfile);

    gmp_fprintf(stderr, "algebraic square root is: %Zd\n", algsqrt);
    fprintf (stderr, "Algebraic square root time is %2.2lf\n", seconds() - t0);
    mpz_clear(aux);
    mpz_clear(algsqrt);
    poly_free(prd->p);
    poly_free(tmp->p);
    poly_free(F);
    return 0;
}

/*
 * Try to factor input using trial division up to bound B.
 * Found factors are printed (one per line).
 * Returns 1 if input is completely factored, otherwise, returns
 * remaining factor.
 */
unsigned long 
trialdivide_print(unsigned long N, unsigned long B)
{
    ASSERT(N != 0);
    if (N == 1) return 1;
    unsigned long p;
    for (p = 2; p <= B; p = getprime (p)) {
        while ((N%p) == 0) {
            N /= p;
            printf("%ld\n", p);
            if (N == 1) {
                getprime(0);  // free table
                return N;
            }
        }
    }
    getprime(0);
    return N;
}

void print_nonsmall(mpz_t zx)
{
    if (mpz_probab_prime_p(zx, 10))
        gmp_printf("%Zd\n", zx);
    else {
        int pp = mpz_perfect_power_p(zx);
        if (pp) {
            pp = mpz_sizeinbase(zx, 2);
            mpz_t roo;
            mpz_init(roo);
            while (!mpz_root(roo, zx, pp))
                pp--;
            int i;
            for (i = 0; i < pp; ++i)
                gmp_printf("%Zd\n", roo);
            mpz_clear(roo);
        } else
            gmp_printf("%Zd\n", zx);
    }
}

void print_factor(mpz_t N) 
{
    unsigned long xx = mpz_get_ui(N);
    if (mpz_cmp_ui(N, xx) == 0) {
        xx = trialdivide_print(xx, 1000000);
        if (xx != 1) {
            mpz_t zx;
            mpz_init(zx);
            mpz_set_ui(zx, xx);
            print_nonsmall(zx);
            mpz_clear(zx);
        }
    } else 
        print_nonsmall(N);
}


/********** GCD **********/
int
calculateGcd(const char *prefix, int numdep, mpz_t Np)
{
    char ratname[200];
    char algname[200];
    snprintf(ratname, sizeof(ratname), "%s.rat.%03d", prefix, numdep);
    snprintf(algname, sizeof(algname), "%s.alg.%03d", prefix, numdep);
    FILE *ratfile = NULL;
    FILE *algfile = NULL;
    int found = 0;
    mpz_t ratsqrt, algsqrt, g1, g2;

    mpz_init(ratsqrt);
    mpz_init(algsqrt);
    mpz_init(g1);
    mpz_init(g2);
  
    ratfile = fopen(ratname, "r");
    algfile = fopen(algname, "r");
    ASSERT_ALWAYS(ratfile != NULL);
    ASSERT_ALWAYS(algfile != NULL);
    gmp_fscanf(ratfile, "%Zd", ratsqrt);
    gmp_fscanf(algfile, "%Zd", algsqrt);
    fclose(ratfile);
    fclose(algfile);

    // reduce mod Np
    mpz_mod(ratsqrt, ratsqrt, Np);
    mpz_mod(algsqrt, algsqrt, Np);

    // First check that the squares agree
    mpz_mul(g1, ratsqrt, ratsqrt);
    mpz_mod(g1, g1, Np);

    mpz_mul(g2, algsqrt, algsqrt);
    mpz_mod(g2, g2, Np);

    if (mpz_cmp(g1, g2)!=0) {
      fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
      exit(1);
      //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
    }

    mpz_sub(g1, ratsqrt, algsqrt);
    mpz_gcd(g1, g1, Np);
    if (mpz_cmp(g1,Np)) {
      if (mpz_cmp_ui(g1,1)) {
        found = 1;
        print_factor(g1);
      }
    }

    mpz_add(g2, ratsqrt, algsqrt);
    mpz_gcd(g2, g2, Np);
    if (mpz_cmp(g2,Np)) {
      if (mpz_cmp_ui(g2,1)) {
        found = 1;
        print_factor(g2);
      }
    }
    mpz_clear(g1);
    mpz_clear(g2);

    if (!found)
      printf("Failed\n");
  
    mpz_clear(ratsqrt);
    mpz_clear(algsqrt);
  
    return 0;
}

void create_dependencies(const char * prefix, const char * indexname, const char * purgedname, const char * kername)
{
    FILE * ix = fopen(indexname, "r");
    int small_nrows, small_ncols;
    int ret;

    /* small_ncols isn't used here: we don't care. */
    ret = fscanf(ix, "%d %d", &small_nrows, &small_ncols);
    ASSERT(ret == 2);

    FILE * ker;
    size_t ker_stride;
    /* Check that kername has consistent size */
    {
        ker = fopen(kername, "r");
        if (ker == NULL) { perror(kername); exit(errno); }
        struct stat sbuf[1];
        ret = fstat(fileno(ker), sbuf);
        if (ret < 0) { perror(kername); exit(errno); }
        ASSERT_ALWAYS(sbuf->st_size % small_nrows == 0);
        unsigned int ndepbytes = sbuf->st_size / small_nrows;
        fprintf(stderr, "%s contains %u dependencies (including padding)\n",
                kername, 8 * ndepbytes);
        ker_stride = ndepbytes - sizeof(uint64_t);
        if (ker_stride)
            fprintf(stderr, "Considering only the first 64 dependencies\n");
    }

    /* We also initiate the purgedfile_stream structure, it gives us the
     * number of (a,b) pairs -- and we need to allocate a structure of
     * that size early on
     */
    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, purgedname);
    
    uint64_t * abs = malloc(ps->nrows * sizeof(uint64_t));
    memset(abs, 0, ps->nrows * sizeof(uint64_t));

    for(int i = 0 ; i < small_nrows ; i++) {
        uint64_t v;
        ret = fread(&v, sizeof(uint64_t), 1, ker);
        if (ker_stride) fseek(ker, ker_stride, SEEK_CUR);

        /* read the index row */
        int nc;
        ret = fscanf(ix, "%d", &nc); ASSERT_ALWAYS(ret == 1);
        for(int k = 0 ; k < nc ; k++) {
            unsigned int col;
            ret = fscanf(ix, "%x", &col); ASSERT_ALWAYS(ret == 1);
            ASSERT_ALWAYS(col < (unsigned int) ps->nrows);
            abs[col] ^= v;
        }
    }
    fclose(ix);
    fclose(ker);

    unsigned int nonzero_deps = 0;
    uint64_t sanity = 0;
    for(int i = 0 ; i < ps->nrows ; i++) {
        sanity |= abs[i];
    }
    uint64_t dep_masks[64]={0,};
    char * dep_names[64];
    FILE * dep_files[64];
    unsigned int dep_counts[64]={0,};

    for(int i = 0 ; i < 64 ; i++) {
        uint64_t m = UINT64_C(1) << i;
        if (sanity & m)
            dep_masks[nonzero_deps++] = m;
    }
    fprintf(stderr, "Total: %u non-zero dependencies\n", nonzero_deps);
    for(unsigned int i = 0 ; i < nonzero_deps ; i++) {
        ret = asprintf(&dep_names[i], "%s.%03d", prefix, i);
        ASSERT_ALWAYS(ret >= 0);
        dep_files[i] = fopen(dep_names[i], "w");
        ASSERT_ALWAYS(dep_files[i] != NULL);
    }
    ps->parse_only_ab = 1;
    for(int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
        ASSERT_ALWAYS(i < ps->nrows);
        for(unsigned int j = 0 ; j < nonzero_deps ; j++) {
            if (abs[i] & dep_masks[j]) {
                fprintf(dep_files[j], "%"PRId64" %"PRIu64"\n", ps->a, ps->b);
                dep_counts[j]++;
            }
        }
        if (purgedfile_stream_disp_progress_now_p(ps)) {
            fprintf(stderr, "read (a,b) pair # %d / %d at %.1f -- %.1f MB/s -- %.1f pairs / s\n",
                    ps->rrows, ps->nrows, ps->dt, ps->mb_s, ps->rows_s);
        }
    }
    purgedfile_stream_trigger_disp_progress(ps);
    fprintf(stderr, "read (a,b) pair # %d / %d at %.1f -- %.1f MB/s -- %.1f pairs / s\n",
            ps->rrows, ps->nrows, ps->dt, ps->mb_s, ps->rows_s);

    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);

    fprintf(stderr, "Written %u dependencies files\n", nonzero_deps);
    for(unsigned int i = 0 ; i < nonzero_deps ; i++) {
        fprintf(stderr, "%s : %u (a,b) pairs\n", dep_names[i], dep_counts[i]);
        fclose(dep_files[i]);
        free(dep_names[i]);
    }
}


void usage(const char * me)
{
    fprintf(stderr, "Usage: %s [-ab || -rat || -alg || -gcd] -poly polyname -prefix prefix -dep numdep", me);
    fprintf(stderr, " -purged purgedname -index indexname -ker kername\n");
    fprintf(stderr, "or %s (-rat || -alg || -gcd) -poly polyname -prefix prefix -dep numdep\n\n", me);
    fprintf(stderr, "(a,b) pairs of dependency relation 'numdep' will be r/w in file 'prefix.numdep',");
    fprintf(stderr, " rational sqrt in 'prefix.rat.numdep' ...\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    cado_poly pol;
    int numdep = -1, ret, i;

    char * me = *argv;
    /* print the command line */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list pl;
    param_list_init(pl);

    int opt_ab = 0;
    int opt_rat = 0;
    int opt_alg = 0;
    int opt_gcd = 0;
    param_list_configure_knob(pl, "ab", &opt_ab);
    param_list_configure_knob(pl, "rat", &opt_rat);
    param_list_configure_knob(pl, "alg", &opt_alg);
    param_list_configure_knob(pl, "gcd", &opt_gcd);
    param_list_configure_knob(pl, "-v", &verbose);
    argc--,argv++;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (strcmp(*argv, "--help") == 0) {
            usage(me);
            exit(0);
        } else {
            fprintf(stderr, "unexpected argument: %s\n", *argv);
            usage(me);
            exit(1);
        }
    }
    const char * tmp;
    ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
    cado_poly_init(pol);
    ret = cado_poly_read(pol, tmp);
    ASSERT (ret);

    param_list_parse_int (pl, "dep", &numdep);
    const char * purgedname = param_list_lookup_string(pl, "purged");
    const char * indexname = param_list_lookup_string(pl, "index");
    const char * kername = param_list_lookup_string(pl, "ker");
    const char * prefix = param_list_lookup_string(pl, "prefix");
    if (param_list_warn_unused(pl))
        exit(1);

    /* if no options then -ab -rat -alg -gcd */
    if (!(opt_ab || opt_rat || opt_alg || opt_gcd))
        opt_ab = opt_rat = opt_alg = opt_gcd = 1;

    ASSERT_ALWAYS(prefix);

    /*
     * In the case where the number N to factor has a prime factor that
     * divides the leading coefficient of f or g, the reduction modulo N
     * will fail. Let's compute N', the factor of N that is coprime to
     * those leading coefficients.
     */
    mpz_t Np;
    mpz_init(Np);
    {
        mpz_t gg;
        mpz_init(gg);
        mpz_set(Np, pol->n);
        do {
            mpz_gcd(gg, Np, pol->g[pol->degreeg]);
            if (mpz_cmp_ui(gg, 1) != 0) {
                gmp_fprintf(stderr, "Warning: found the following factor of N as a factor of g: %Zd\n", gg);
                print_factor(gg);
                mpz_divexact(Np, Np, gg);
            }
        } while (mpz_cmp_ui(gg, 1) != 0);
        do {
            mpz_gcd(gg, Np, pol->f[pol->degree]);
            if (mpz_cmp_ui(gg, 1) != 0) {
                gmp_fprintf(stderr, "Warning: found the following factor of N as a factor of f: %Zd\n", gg);
                print_factor(gg);
                mpz_divexact(Np, Np, gg);
            }
        } while (mpz_cmp_ui(gg, 1) != 0);
        mpz_clear(gg);
        /* Trial divide Np, to avoid bug if a stupid input is given */
        {
            unsigned long p;
            for (p = 2; p <= 1000000; p = getprime (p)) {
                while (mpz_tdiv_ui(Np, p) == 0) {
                    printf("%lu\n", p);
                    mpz_divexact_ui(Np, Np, p);
                }
            }
            getprime(0);
        }
        if (mpz_cmp(pol->n, Np) != 0) 
            gmp_fprintf(stderr, "Now factoring N' = %Zd\n", Np);
        if (mpz_cmp_ui(Np, 1) == 0) {
            gmp_fprintf(stderr, "Hey N' is 1! Stopping\n");
            cado_poly_clear (pol);
            param_list_clear(pl);
            mpz_clear(Np);
            return 0;
        }
        if (mpz_probab_prime_p(Np, 10) || mpz_perfect_power_p(Np)) {
            gmp_fprintf(stderr, "Hey N' is (power of) prime! Stopping\n");
            print_factor(Np);
            cado_poly_clear (pol);
            param_list_clear(pl);
            mpz_clear(Np);
            return 0;
        }
    }

    if (opt_ab) {
        /* Computing (a,b) pairs is now done in batch for 64 dependencies
         * together -- should be enough for our purposes, even if we do
         * have more dependencies !
         */
        create_dependencies(prefix, indexname, purgedname, kername);
        ASSERT_ALWAYS(indexname != NULL);
        ASSERT_ALWAYS(purgedname != NULL);
        ASSERT_ALWAYS(kername != NULL);
    }

    if (opt_rat) {
        ASSERT_ALWAYS(numdep != -1);
        if (pol->degreeg == 1)
            calculateSqrtRat (prefix, numdep, pol, Np);
        else
            calculateSqrtAlg (prefix, numdep, pol, 1, Np);
    }

    if (opt_alg) {
        ASSERT_ALWAYS(numdep != -1);
        calculateSqrtAlg (prefix, numdep, pol, 0, Np);
    }

    if (opt_gcd) {
        ASSERT_ALWAYS(numdep != -1);
        calculateGcd (prefix, numdep, Np);
    }
    
    cado_poly_clear (pol);
    param_list_clear(pl);
    mpz_clear(Np);
    return 0;
}
  
