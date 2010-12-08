/*
 * Program: algsqrt
 * Authors: P. Gaudry, P. Zimmermann
 * Purpose: computing the squareroots and finishing the factorization
 *
 Possible easy improvements:
 * implement a subquadratic reduction mod F in poly_mod_f_mod_mpz.
   However the main cost now comes from the barrett_mod() calls, i.e.,
   from the reduction of coefficients, not from that of the polynomial.
 * instead of collecting all a-b*x in P(x), computing 1/sqrt(P(x)) mod p^k,
   then multiplying by P(x) to get sqrt(P(x)) mod p^k, avoid the product
   by P(x), and try to recognize a polynomial with rational coefficients
   (using rational reconstruction). If this works, and if the denominator is
   the same, this will trade a cost of (2d-1)M(n) for M(n)log(n) + (d-1)M(n).
 Harder improvements:
 * instead of collecting all a-b*x in the numerator, collect half of them
   in the numerator, and half of them in the denominator (cf paper of Nguyen
   at ANTS3, and thesis of Elkenbracht-Huizing). Some experiments suggest
   that only 1/3rd of the precision is enough to recover sqrt(P/Q), which
   is then a polynomial with *rational* coefficients.
 * cache FFT transforms wherever possible
 * use Karp-Markstein trick to incorporate in the last Hensel iteration for
   1/sqrt(x) the value of x to get an approximation of sqrt(x)
 * use the trick that multiplies p(x) by f'(x)^2 [see Crandall & Pomerance,
   Prime Numbers, A Computational Perspective, Sections 6.2.4, 6.2.5, 6.2.6]
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <math.h> /* for log */
#include "utils.h"

/* Although the functions in plain_poly are not readily available in the
 * publicized interface of utils.h, it's ok to use them if we explicitly
 * include the corresponding header.
 */
#include "plain_poly.h"

//#include "poly.h"

// #define VERBOSE

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
    mpz_divexact_ui(aux, aux, 4);		        // aux := (3q-5)/4
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
    mpz_mul (pk, pk, pk);	// double the current precision
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

unsigned long FindSuitableModP(poly_t F) {
  unsigned long p = 2;
  int dF = F->deg;

  plain_poly_t fp;
  plain_poly_init(fp, dF);
  while (1) {
    int d;
    // select prime congruent to 3 mod 4 to have easy sqrt.
  //  do {
      p = getprime(p);
  //  } while ((p%4)!= 3);
    d = plain_poly_set_mod (fp, F->coeff, dF, p);
    if (d!=dF)
      continue;
    if (plain_poly_is_irreducible(fp, p))
      break;
  }
  plain_poly_clear(fp);
  getprime (0);

  return p;
}

#define THRESHOLD 2 /* must be >= 2 */

/* accumulate up to THRESHOLD products in prd[0], 2^i*THRESHOLD in prd[i].
   nprd is the number of already accumulated values: if nprd = n0 +
   n1 * THRESHOLD + n2 * THRESHOLD^2 + ..., then prd[0] has n0 entries,
   prd[1] has n1*THRESHOLD entries, and so on.
*/
// Products are computed modulo the polynomial F.
polymodF_t*
accumulate_fast (polymodF_t *prd, const polymodF_t a, const poly_t F,
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
accumulate_fast_end (polymodF_t *prd, const poly_t F, unsigned long lprd)
{
  unsigned long i;

  for (i = 1; i < lprd; i++)
    polymodF_mul (prd[0], prd[0], prd[i], F);
}

int
main(int argc, char **argv)
{
  FILE * depfile;
  cado_poly pol;
  poly_t F;
  polymodF_t prd, tmp;
  long a;
  unsigned long b;
  int ret, i;
  unsigned long p;
  double t0 = seconds ();

  if (argc != 4) {
    fprintf(stderr, "usage: %s algdepfile ratdepfile polyfile\n", argv[0]);
    exit(1);
  }

  /* print the command line */
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  depfile = fopen(argv[1], "r");
  ASSERT_ALWAYS(depfile != NULL);
  cado_poly_init(pol);
  ret = cado_poly_read(pol, argv[3]);
  ASSERT_ALWAYS(ret == 1);

  // Init F to be the algebraic polynomial
  poly_alloc(F, pol->degree);
  for (i = pol->degree; i >= 0; --i)
    poly_setcoeff(F, i, pol->f[i]);

  // Init prd to 1.
  poly_alloc(prd->p, pol->degree);
  mpz_set_ui(prd->p->coeff[0], 1);
  prd->p->deg = 0;
  prd->v = 0;

  // Allocate tmp
  poly_alloc(tmp->p, 1);

  // Accumulate product
  int nab = 0, nfree = 0;
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
    poly_alloc(prd_tab[0]->p, F->deg);
    mpz_set_ui(prd_tab[0]->p->coeff[0], 1);
    prd_tab[0]->p->deg = 0;
    prd_tab[0]->v = 0;
    while(fscanf(depfile, "%ld %lu", &a, &b) != EOF){
      if(!(nab % 100000))
	fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n", nab, seconds ());
      if((a == 0) && (b == 0))
	break;
      polymodF_from_ab(tmp, a, b);
      prd_tab = accumulate_fast (prd_tab, tmp, F, &lprd, nprd++);
      nab++;
      if(b == 0)
	  nfree++;
    }
    fprintf (stderr, "# Read %d including %d free relations\n", nab, nfree);
    ASSERT_ALWAYS ((nab & 1) == 0);
    ASSERT_ALWAYS ((nfree & 1) == 0);
    accumulate_fast_end (prd_tab, F, lprd);
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
           (unsigned long) poly_sizeinbase (prd->p, pol->degree - 1, 2));

  p = FindSuitableModP(F);
  fprintf(stderr, "Using p=%lu for lifting\n", p);

  double tm = seconds();
  polymodF_sqrt (prd, prd, F, p);
  fprintf (stderr, "Square root lifted in %2.2lf\n", seconds()-tm);

  mpz_t algsqrt, aux;
  mpz_init(algsqrt);
  mpz_init(aux);
  poly_eval_mod_mpz(algsqrt, prd->p, pol->m, pol->n);
  mpz_invert(aux, F->coeff[F->deg], pol->n);  // 1/fd mod n
  mpz_powm_ui(aux, aux, prd->v, pol->n);      // 1/fd^v mod n
  mpz_mul(algsqrt, algsqrt, aux);
  mpz_mod(algsqrt, algsqrt, pol->n);

  gmp_fprintf(stderr, "Algebraic square root is: %Zd\n", algsqrt);

  depfile = fopen(argv[2], "r");
  ASSERT_ALWAYS(depfile != NULL);
  gmp_fscanf(depfile, "%Zd", aux);
  fclose(depfile);

  gmp_fprintf(stderr, "Rational square root is: %Zd\n", aux);
  fprintf (stderr, "Total square root time is %2.2lf\n", seconds() - t0);

  int found = 0;
  {
    mpz_t g1, g2;
    mpz_init(g1);
    mpz_init(g2);

    // First check that the squares agree
    mpz_mul(g1, aux, aux);
    mpz_mod(g1, g1, pol->n);
    if(mpz_cmp_ui(pol->g[1], 1) != 0){
	// case g(X)=m1*X+m2 with m1 != 1
	// we should have prod (a+b*m2/m1) = A^2 = R^2/m1^(nab-nfree)
	// and therefore nab should be even
	if(nab & 1){
	    fprintf(stderr, "Sorry, but #(a, b) is odd\n");
	    fprintf(stderr, "Bug: this should be patched! Please report your buggy input\n");
	    printf("Failed\n");
	    return 0;
	}
	mpz_powm_ui(g2, pol->g[1], (nab-nfree)>>1, pol->n);
	mpz_mul(algsqrt, algsqrt, g2);
	mpz_mod(algsqrt, algsqrt, pol->n);
    }
    mpz_mul(g2, algsqrt, algsqrt);
    mpz_mod(g2, g2, pol->n);
    if (mpz_cmp(g1, g2)!=0) {
      fprintf(stderr, "Bug: the squares do not agree modulo n!\n");
      //      gmp_printf("g1:=%Zd;\ng2:=%Zd;\n", g1, g2);
    }
    mpz_sub(g1, aux, algsqrt);
    mpz_gcd(g1, g1, pol->n);

    if (mpz_cmp(g1,pol->n)) {
      if (mpz_cmp_ui(g1,1)) {
        found = 1;
        gmp_printf("%Zd\n", g1);
      }
    }
    mpz_add(g2, aux, algsqrt);
    mpz_gcd(g2, g2, pol->n);
    if (mpz_cmp(g2,pol->n)) {
      if (mpz_cmp_ui(g2,1)) {
        found = 1;
        gmp_printf("%Zd\n", g2);
      }
    }
    mpz_clear(g1);
    mpz_clear(g2);
  }
  if (!found)
    printf("Failed\n");

  cado_poly_clear (pol);
  mpz_clear(aux);
  mpz_clear(algsqrt);
  poly_free(F);
  poly_free(prd->p);
  poly_free(tmp->p);

  return 0;
}

