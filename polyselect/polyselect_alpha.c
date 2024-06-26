#include "cado.h"
#include <gmp.h>
#include <math.h>
#include <pthread.h>

#include "mpz_poly.h"
#include "rootfinder.h"
#include "gmp_aux.h"
#include "polyselect_alpha.h"

/************************** polynomial arithmetic ****************************/

/* h(x) <- h(x + r/p), where the coefficients of h(x + r/p) are known to
   be integers */
static void
poly_shift_divp (mpz_t *h, unsigned int d, unsigned long r, unsigned long p)
{
  unsigned int i, k;
  mpz_t t;

  mpz_init (t);
  for (i = 1; i <= d; i++)
    for (k = d - i; k < d; k++)
      { /* h[k] <- h[k] + r/p h[k+1] */
        ASSERT (mpz_divisible_ui_p (h[k+1], p) != 0);
        mpz_divexact_ui (t, h[k+1], p);
        mpz_addmul_ui (h[k], t, r);
      }
  mpz_clear (t);
}

/********************* computation of alpha **********************************/

/* Auxiliary routine for special_valuation(), see below. It returns the
   average p-valuation of the polynomial f. Works recursively.
   Assumes f is square-free, otherwise it will loop. */
static double
special_val0 (mpz_poly_srcptr f, unsigned long p, gmp_randstate_ptr rstate)
{
  double v;
  mpz_t c,  *h;
  unsigned long *roots, r, r0;
  int i, d = f->deg, nroots;
  mpz_poly g, H;

  mpz_init (c);
  mpz_poly_content (c, f);
  for (v = 0.0; mpz_divisible_ui_p (c, p); v++, mpz_divexact_ui (c, c, p));

  mpz_poly_init(g, d);
  g->deg = d;

  /* g <- f/p^v */
  if (v != 0.0)
    {
      mpz_ui_pow_ui (c, p, (unsigned long) v); /* p^v */
      for (i = 0; i <= d; i++)
        mpz_divexact (g->coeff[i], f->coeff[i], c);
    }
  else
    mpz_poly_set (g, f);

  mpz_poly_init (H, d);
  H->deg = d;
  h = H->coeff;
  /* first compute h(x) = g(p*x) = f(p*x)/p^v */
  mpz_set_ui (c, 1);
  for (i = 0; i <= d; i++)
    {
      mpz_mul (h[i], g->coeff[i], c);
      mpz_mul_ui (c, c, p);
    }
  /* Search for roots of g mod p */
  ASSERT (d > 0);
  roots = (unsigned long*) malloc (d * sizeof (unsigned long));
  FATAL_ERROR_CHECK(roots == NULL, "not enough memory");

  nroots = mpz_poly_roots_ulong (roots, g, p, rstate);
  ASSERT (nroots <= d);
  for (r0 = 0, i = 0; i < nroots; i++)
    {
      r = roots[i];
      mpz_poly_eval_diff_ui (c, g, r);
      if (mpz_divisible_ui_p (c, p) == 0) /* g'(r) <> 0 mod p */
        v += 1.0 / (double) (p - 1);
      else /* hard case */
        {
          /* g(p*x+r) = h(x + r/p), thus we can go from h0(x)=g(p*x+r0)
             to h1(x)=g(p*x+r1) by computing h0(x + (r1-r0)/p).
             Warning: we can have h = f, and thus an infinite loop, when
             the p-valuation of f is d, and f has a single root r/(1-p) of
             multiplicity d.
             Moreover if f(x) = c*p^d*(x-r+b*p)^d, where c is coprime to p,
             then h(x) = f(p*x+r)/p^d = c*p^d*(x+b)^d, and most likely after
             at most p iterations we'll go back to f(x), thus we should avoid
             all cases where f(x) has a root of multiplicity d, but how to
             check that efficiently? And which value to return in such a case?
          */
          ASSERT_ALWAYS (r >= r0); /* the roots are sorted */
          poly_shift_divp (h, d, r - r0, p);
          r0 = r;
          v += special_val0 (H, p, rstate) / (double) p;
        }
    }
  free (roots);
  mpz_poly_clear (H);
  mpz_poly_clear (g);
  mpz_clear (c);

  return v;
}

/* Compute the average valuation of F(a,b) for gcd(a,b)=1, for a prime p
   dividing the discriminant of f, using the following algorithm from
   Guillaume Hanrot (which is some kind of p-adic variant of Uspensky's
   algorithm):

   val(f, p)
     return val0(f, p) * p / (p+1) + val0(f(1/(p*x))*(p*x)^d, p) * 1/(p+1)

   val0(f, p).
     v <- valuation (content(f), p);
     f <- f/p^v

     r <- roots mod p(f, p)

     for r_i in r do
         if f'(r_i) <> 0 mod p then v +=  1/(p-1).
         else
              f2 <- f(p*x + r_i)
              v += val0(f2, p) / p.
         endif
     endfor
     Return v.

A special case when:
(a) p^2 does not divide disc(f),
(b) p does not divide lc(f),
then the average valuation is (p q_p - 1)/(p^2 - 1), where q_p is the number
of roots of f mod p. When q_p=1, we get 1/(p+1).

Note: when p does not divide lc(f), the val0(f(1/(p*x))*(p*x)^d, p) call
always returns 0 in val(f,p).

Assumes p divides disc = disc(f), d is the degree of f.
*/
double
special_valuation (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate)
{
    double v;
    int p_divides_lc;
    int pvaluation_disc = 0;
    double pd = (double) p;
    int d = f->deg;

    if (mpz_divisible_ui_p(disc, p)) {
  mpz_t t;
  pvaluation_disc++;
  mpz_init(t);
  mpz_divexact_ui(t, disc, p);
  if (mpz_divisible_ui_p(t, p))
      pvaluation_disc++;
  mpz_clear(t);
    }

    p_divides_lc = mpz_divisible_ui_p(f->coeff[d], p);

    if (pvaluation_disc == 0) {
  /* easy ! */
  int e;
  e = mpz_poly_roots_ulong (NULL, f, p, rstate);
  if (p_divides_lc) {
      /* Or the discriminant would have valuation 1 at least */
      ASSERT(mpz_divisible_ui_p(f->coeff[d - 1], p) == 0);
      e++;
  }
  return (pd * e) / (pd * pd - 1);
    } else if (pvaluation_disc == 1) {
      /* special case where p^2 does not divide disc */
  int e;
  e = mpz_poly_roots_ulong (NULL, f, p, rstate);
        if (p_divides_lc)
          e ++;
  /* something special here. */
  return (pd * e - 1) / (pd * pd - 1);
    } else {
  v = special_val0(f, p, rstate) * pd;
  if (p_divides_lc) {
      /* compute g(x) = f(1/(px))*(px)^d, i.e., g[i] = f[d-i]*p^i */
      /* IOW, the reciprocal polynomial evaluated at px */
      mpz_poly G;
      mpz_t *g;
      mpz_t t;
      int i;

      mpz_poly_init (G, d);
      G->deg = d;
      g = G->coeff;
      mpz_init_set_ui(t, 1);  /* will contains p^i */
      for (i = 0; i <= d; i++) {
        mpz_mul(g[i], f->coeff[d - i], t);
        mpz_mul_ui(t, t, p);
      }
      v += special_val0(G, p, rstate);
      mpz_poly_clear (G);
      mpz_clear(t);
  }
  v /= pd + 1.0;
  return v;
    }
}

/* Compute the value alpha(F) from Murphy's thesis, page 49:
   alpha(F) = sum(prime p <= B, (1 - q_p*p/(p+1)) log(p)/(p-1))
   where q_p is the number of roots of F mod p, including the number of
   projective roots (i.e., the zeros of the reciprocal polynomial mod p).

   alpha(F) is an estimate of the average logarithm of the part removed
   from sieving, compared to a random integer.

   We want alpha as small as possible, i.e., alpha negative with a large
   absolute value. Typical good values are alpha=-4, -5, ...
*/
double
get_alpha (mpz_poly_srcptr f, unsigned long B)
{
  double alpha, e;
  unsigned long p;
  mpz_t disc;

  /* for F linear, we have q_p = 1 for all p, thus
     alpha(F) = sum(prime p <= B, log(p)/(p^2-1)) ~ 0.569959993064325 */
  if (f->deg == 1)
    return 0.569959993064325;

  /* a gmp_randstate init/clear cycle is about 10 times the cost of a
   * one-word gmp random pick. Given that we assume that B is at least
   * in the hundreds, the random picks alone would be a sufficient
   * observation to conclude that it's relatively harmless to do a
   * random initialization here, and present get_alpha as as something
   * deterministic. (and of course, there's tons of arithmetic in there
   * too.
   */

  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);

  mpz_init (disc);
  mpz_poly_discriminant (disc, f);

  /* special_valuation returns the expected average exponent of p in F(a,b)
     for coprime a, b, i.e., e = q_p*p/(p^2-1), thus the contribution for p
     is (1/(p-1) - e) * log(p) */

  /* prime p=2 */
  e = special_valuation (f, 2, disc, rstate);
  alpha = (1.0 - e) * log (2.0);

  /* FIXME: generate all primes up to B and pass them to get_alpha */
  for (p = 3; p <= B; p += 2)
    if (ulong_isprime (p))
      {
        e = special_valuation (f, p, disc, rstate);
        alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
      }
  gmp_randclear(rstate);
  mpz_clear (disc);
  return alpha;
}

/* affine part of the special valution for polynomial f over p. */
double
special_valuation_affine (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate)
{
   double v;
   int pvaluation_disc = 0;
   double pd = (double) p;

   if (mpz_divisible_ui_p(disc, p)) {
      mpz_t t;
      pvaluation_disc++;
      mpz_init(t);
      mpz_divexact_ui(t, disc, p);
      if (mpz_divisible_ui_p(t, p))
         pvaluation_disc++;
      mpz_clear(t);
   }

   if (pvaluation_disc == 0) {
      /* case 1: root must be simple*/
      int e = 0;
      e = mpz_poly_roots_ulong (NULL, f, p, rstate);

      return (pd * e) / (pd * pd - 1);
   }
   /* else if (pvaluation_disc == 1) { */
   /*     /\* case 2: special case where p^2 does not divide disc *\/ */
   /*     int e = 0; */
   /*     e = mpz_poly_roots_ulong (NULL, f, p); */

   /*     /\* something special here. *\/ */
   /*     return (pd * e - 1) / (pd * pd - 1); */

   /* } */
   else {
      v = special_val0(f, p, rstate) * pd;
      v /= pd + 1.0;
      return v;
   }
}


/*
  Find alpha_projective for a poly f. It uses
  some hacks here which need to be changed in future.
  Until now, since this will only be done several
  times, hence the speed is not critical.

  Note that, the returned alpha is the -val * log(p)
  part in the alpha. Hence, we can just add
  this to our affine part.
*/
double
get_alpha_projective (mpz_poly_srcptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;

   /* about random state: see comment in get_alpha */
   gmp_randstate_t rstate;
   gmp_randinit_default(rstate);

   mpz_init (disc);
   mpz_poly_discriminant (disc, f);

   /* prime p=2 */
   e = special_valuation (f, 2, disc, rstate) - special_valuation_affine (f, 2, disc, rstate);

   /* 1/(p-1) is counted in the affine part */
   alpha =  (- e) * log (2.0);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation(f, p, disc, rstate) - special_valuation_affine (f, p, disc, rstate);
         alpha += (- e) * log ((double) p);
      }

   gmp_randclear(rstate);
   mpz_clear (disc);

   return alpha;
}

/*
  Similar to above, but for affine part.
*/
double
get_alpha_affine (mpz_poly_srcptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;

  /* about random state: see comment in get_alpha */
  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);

   mpz_init (disc);
   mpz_poly_discriminant (disc, f);

   /* prime p=2 */
   e = special_valuation_affine (f, 2, disc, rstate);
   alpha =  (1.0 - e) * log (2.0);

   //printf ("\np: %u, val: %f, alpha: %f\n", 2, e, alpha);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation_affine (f, p, disc, rstate);
         alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
         //printf ("\np: %u, val: %f, alpha: %f\n", p, e, alpha);

      }
   mpz_clear (disc);
   gmp_randclear(rstate);
   return alpha;
}

/*
  Similar to above, but for a given prime p.
*/
double
get_alpha_affine_p (mpz_poly_srcptr f, unsigned long p, gmp_randstate_ptr rstate)
{
   double alpha, e;
   mpz_t disc;

   mpz_init (disc);
   mpz_poly_discriminant (disc, f);

   if (p == 2)
     {
       e = special_valuation_affine (f, 2, disc, rstate);
       alpha =  (1.0 - e) * log (2.0);
     }
   else
     {
       ASSERT (ulong_isprime (p));
       e = special_valuation_affine (f, p, disc, rstate);
       alpha = (1.0 / (double) (p - 1) - e) * log ((double) p);
     }
   mpz_clear (disc);
   return alpha;
}

#if 0
/*
  Contribution from a particular multiple root r of the polynomial f
  over p. Note, r must also be a double root of f mod p.
*/
static double
average_valuation_affine_root (mpz_poly_ptr f, unsigned long p, unsigned long r )
{
   unsigned long v = 0UL;
   int i, j;
   mpz_t c, *fv;
   double val;

   mpz_init (c);

   /* init fv */
   fv = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
   if (fv == NULL) {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit (1);
   }

   for (i = 0; i <= d; i++)
      mpz_init_set (fv[i], f[i]);

   /* remove the p-valuations from fv */
   mpz_poly_content (c, f);
   while (mpz_divisible_ui_p(c, p)) {
      v += 1;
      for (i = 0; i <= d; i ++) {
         mpz_fdiv_q_ui (fv[i], f[i], p);
      }
   }

   /* first translate, then scale */
   for (i = d - 1; i >= 0; i--)
      for (j = i; j < d; j++)
         mpz_addmul_ui (fv[j], fv[j+1], r);
   /* t is p^i */
   mpz_set_ui(c, 1);
   for (i = 0; i <= d; i++) {
      mpz_mul(fv[i], fv[i], c);
      mpz_mul_ui(c, c, p);
   }

   /* now c is disc. */
   discriminant (c, fv, d);
   val = special_valuation_affine (fv, d, p, c);
   val = val / (double) p;

   /* clear */
   for (i = 0; i <= d; i++) {
      mpz_clear (fv[i]);
   }

   /* !!! REMEMBER THIS !!! */
   free (fv);
   mpz_clear(c);
   return val;
}
#endif

int alpha_bound = ALPHA_BOUND; /* default value */

/* affecting a global variable in a central state like this should be
 * avoided. It's rather ridiculous to have this, and reflects bad design.
 * This must go.
 */

pthread_mutex_t alpha_bound_lock = PTHREAD_MUTEX_INITIALIZER;

void
set_alpha_bound (unsigned long bound)
{
    pthread_mutex_lock(&alpha_bound_lock);
    alpha_bound = bound;
    pthread_mutex_unlock(&alpha_bound_lock);
}

unsigned long
get_alpha_bound (void)
{
    unsigned long a;
    pthread_mutex_lock(&alpha_bound_lock);
    a = alpha_bound;
    pthread_mutex_unlock(&alpha_bound_lock);
    return a;
}

/* mu and sigma in function exp_alpha() */
#define MU 0.0

#define SIGMA 0.824

/**
 * Asymptotic estimate of minimum order statistics
 * for 2^K many rotations where 0 <= K <= 149.
 * function expected_alpha_est() in alpha.sage
 */
double expected_alpha (double logK) 
{
  if (logK < 0.999)
    return 0.0;
  return MU - SIGMA * (sqrt(2*logK)-(log(logK)+1.3766)/(2*sqrt(2*logK)));
}

#if 0
/* turns out that this function was implemented several times, with
 * different semantics. Let's keep just one.
 *
 * Note how it is not at all clear that we have logK being the log in
 * base 2 or in base e ...
 */
static double
expected_alpha (double S)
{
  double logS, t;

  if (S <= 1.0)
    return 0.0;

  logS = log (S);
  t = sqrt (2 * logS);
  return -0.824 * (t - (log (logS) + 1.3766) / (2 * t));
}
#endif
