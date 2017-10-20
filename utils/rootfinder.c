#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <gmp.h>

#include "macros.h"
#include "rootfinder.h"
#include "mpz_poly.h"
#include "modul_poly.h"
#include "gmp_aux.h"
#include "getprime.h"
#include "portability.h"

/* Entry point for rootfind routines, for an integer n0 not necessarily prime.
   Since we cannot know in advance an easy bound on the number of
   roots, we allocate them in the function: if r = rp[0] at exit,
   the roots are r[0], r[1], ..., r[k-1] and the return value is k.
   Note: the elements r[j] must be mpz_clear'ed by the caller, and the
   array r also.
*/
unsigned long
mpz_poly_roots_gen (mpz_t **rp, mpz_poly_srcptr F, mpz_srcptr n)
{
  unsigned long k, i, j, d = F->deg;
  mpz_t Q, nn, p;

  ASSERT_ALWAYS (mpz_sgn (n) > 0);

  if (mpz_probab_prime_p (n, 1))
    {
      rp[0] = malloc (d * sizeof (mpz_t));
      for (i = 0; i < d; i++)
        mpz_init (rp[0][i]);
      k = mpz_poly_roots (rp[0], F, n);
      /* free the unused roots */
      for (i = k; i < d; i++)
        mpz_clear (rp[0][i]);
      if (k < d)
        rp[0] = realloc (rp[0], k * sizeof (mpz_t));
      return k;
    }

  rp[0] = malloc (sizeof (mpz_t));
  mpz_init_set_ui (rp[0][0], 0);
  mpz_init_set_ui (Q, 1);
  mpz_init_set (nn, n);
  k = 1;

  /* now n is composite */
  mpz_t v, q, x;
  mpz_init (v);
  mpz_init (q);
  mpz_init (x);
  mpz_t *roots_p = malloc (d * sizeof(mpz_t));
  for (i = 0; i < d; i++)
    mpz_init (roots_p[i]);
  for (mpz_init_set_ui (p, 2); mpz_cmp_ui (nn, 1) > 0; mpz_nextprime (p, p))
    {
      if (mpz_probab_prime_p (nn, 1))
        mpz_set (p, nn);
      if (mpz_divisible_p (nn, p))
        {
          unsigned long kp;
          kp = mpz_poly_roots (roots_p, F, p);
          mpz_divexact (nn, nn, p);
          /* lift roots mod p^j if needed */
          mpz_set (q, p);
          while (mpz_divisible_p (nn, p))
            {
              int ii;
              mpz_mul (q, q, p);
              mpz_divexact (nn, nn, p);
              for (i = ii = 0; i < kp; i++)
                {
                  /* FIXME: replace this naive for-loop */
                  mpz_set (roots_p[ii], roots_p[i]);
                  for (mpz_set_ui (x, 0); mpz_cmp (x, p) < 0;
                       mpz_add_ui (x, x, 1))
                    {
                      mpz_poly_eval (v, F, roots_p[ii]);
                      if (mpz_divisible_p (v, q))
                        break;
                      mpz_add (roots_p[ii], roots_p[ii], p);
                    }
                  /* some roots might disappear, for example x^3+2*x^2+3*x-4
                     has two roots mod 2 (0 and 1) but only one mod 4 (0) */
                  ii += (mpz_cmp (x, p) < 0);
                }
              kp = ii;
            }
          /* do a CRT between r[0][0..k-1] mod Q and roots_p[0..kp-1] */
          rp[0] = realloc (rp[0], k * kp * sizeof (mpz_t));
          mpz_invert (x, Q, q); /* x = 1/Q mod q */
          for (i = 0; i < k; i++)
            for (j = kp; j-- > 0;)
              {
                if (j > 0)
                  mpz_init (rp[0][j*k + i]);
                /* x = rp[0][i] mod Q and x = roots_p[j] mod q,
                   thus x = rp[0][i] + Q * t, where
                   t = (roots_p[j] - rp[0][i])/Q mod q */
                mpz_sub (v, roots_p[j], rp[0][i]);
                mpz_mul (v, v, x);
                mpz_mod (v, v, q);
                mpz_mul (v, Q, v);
                mpz_add (rp[0][j*k + i], rp[0][i], v);
              }
          k *= kp;
          mpz_mul (Q, Q, q);
        }
    }
  for (i = 0; i < d; i++)
    mpz_clear (roots_p[i]);
  free (roots_p);
  mpz_clear (v);
  mpz_clear (p);
  mpz_clear (q);
  mpz_clear (x);
  mpz_clear (Q);
  mpz_clear (nn);
  return k;
}

/* Entry point for rootfind routines, for prime p.
   Assume r is an array of deg(F) entries, which are mpz_init'ed. */
int
mpz_poly_roots (mpz_t *r, mpz_poly_srcptr F, mpz_srcptr p)
{
    int d = F->deg;

    if (mpz_cmp_ui(p, ULONG_MAX) <= 0) {
        /* There's a chance of using one of our layers. */
        unsigned long pp = mpz_get_ui(p);
        unsigned long * rr;
        int i;
        int n;
	
        if (r == NULL)
            return mpz_poly_roots_ulong (NULL, F, pp);

        rr = (unsigned long *) malloc(d * sizeof(unsigned long));
        n = mpz_poly_roots_ulong (rr, F, pp);

        for(i = 0 ; i < n ; i++) {
            /* The assumption is that p fits within an unsigned long
             * anyway. So the roots do as well.
             */
            mpz_set_ui(r[i], rr[i]);
        }
        free(rr);
        return n;
    } else {
      int n;
      n = mpz_poly_roots_mpz (r, F, p);
      return n;
    }
}


/* put in r[0], ..., r[n-1] the roots of F modulo p, where p is prime,
   and the return value n is the number of roots (without multiplicities) */
int
mpz_poly_roots_ulong (unsigned long *r, mpz_poly_srcptr F, unsigned long p)
{
    int n;

    residueul_t * rr;
    modulusul_t pp;
    modul_initmod_ul(pp, p);
    int i;
    int d = F->deg;

    if (r == NULL)
      return modul_poly_roots(NULL, F, pp);

    rr = (residueul_t *) malloc(d * sizeof(residueul_t));
    for(i = 0 ; i < d ; i++) {
      modul_init_noset0(rr[i], pp);
    }
    n = modul_poly_roots(rr, F, pp);
    for(int i = 0 ; i < n ; i++) {
      /* The assumption is that p fits within an unsigned long
       * anyway. So the roots do as well.
       */
      r[i] = modul_get_ul(rr[i], pp);
    }
    for(i = 0 ; i < d ; i++) {
      modul_clear(rr[i], pp);
    }
    free(rr);
    modul_clearmod(pp);

    return n;
}



int
mpz_poly_roots_uint64 (uint64_t * r, mpz_poly_srcptr F, uint64_t p)
{
    /* This is glue around poly_roots_ulong, nothing more. When uint64
       is larger than ulong, we call mpz_poly_roots_mpz as a fallback */
    unsigned long *rr;
    int i, n;
    int d = F->deg;

    if (p > (uint64_t) ULONG_MAX)
      {
        mpz_t pp;
        mpz_init (pp);
        mpz_set_uint64 (pp, p);
        if (r == NULL)
          n = mpz_poly_roots_mpz (NULL, F, pp);
        else
          {
            mpz_t *rr;
            rr = malloc ((d + 1) * sizeof (mpz_t));
            for (i = 0; i <= d; i++)
              mpz_init (rr[i]);
            n = mpz_poly_roots_mpz (NULL, F, pp);
            for (i = 0; i <= d; i++)
              {
                if (i < n)
                  r[i] = mpz_get_uint64 (rr[i]);
                mpz_clear (rr[i]);
              }
            free (rr);
          }
        mpz_clear (pp);
        return n;
      }

    if (r == NULL)
      return mpz_poly_roots_ulong (NULL, F, p);

    if (sizeof (unsigned long) != sizeof (uint64_t))
      rr = (unsigned long *) malloc(d * sizeof(unsigned long));
    else
      rr = (unsigned long *) r;
    n = mpz_poly_roots_ulong (rr, F, p);
    if (sizeof (unsigned long) != sizeof (uint64_t))
      {
        for(i = 0 ; i < n ; i++)
          r[i] = rr[i];
        free (rr);
      }
    return n;
}




typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_coeff_cmp(const mpz_t *a, const mpz_t *b) {
  return mpz_cmp(*a, *b) < 0 ? -1 : 1;
}

/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. Return number of roots
   which should be degree of f. Assumes p is odd, and deg(f) >= 1. */
static int
mpz_poly_cantor_zassenhaus (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p, int depth)
{
  mpz_t a, aux;
  mpz_poly q, h;
  int d = f->deg, dq, n, m;

  mpz_init (a);
  mpz_init (aux);

  /* linear polynomial */
  if (d == 1) {
    mpz_neg (aux, f->coeff[1]);
    mpz_invert (a, aux, p);
    mpz_mul (r[0], a, f->coeff[0]);
    mpz_fdiv_r (r[0], r[0], p);
    n = 1;
    goto clear_a;
  }

  /* if f has degree d, then q,h may have up to degree 2d-2 in the
     powering algorithm */
  mpz_poly_init (q, 2 * d - 2);
  mpz_poly_init (h, 2 * d - 2);

  /* random polynomial by a */
  mpz_set_ui (a, lrand48());
  mpz_fdiv_r (a, a, p);
  for (;;)
  {
    /* q=x+a */
    mpz_set_ui (aux, 1);
    mpz_poly_setcoeff (q, 1, aux);
    mpz_poly_setcoeff (q, 0, a);

    /* h=(x+a)^((p-1)/2) mod (f, p) */
    mpz_sub_ui (aux, p, 1);
    mpz_divexact_ui (aux, aux, 2);
    mpz_poly_pow_mod_f_mod_mpz (h, q, f, aux, p);
    mpz_poly_sub_ui (h, h, 1);

    /* q = gcd(f,h) */
    mpz_poly_gcd_mpz (q, f, h, p);
    dq = q->deg;
    ASSERT (dq >= 0);

    /* recursion-split */
    if (0 < dq && dq < d)
    {
      n = mpz_poly_cantor_zassenhaus (r, q, p, depth+1);
      ASSERT (n == dq);

      mpz_poly_divexact (h, f, q, p);
      m = mpz_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
      ASSERT (m == h->deg);
      n += m;
      break;
    }

    mpz_add_ui (a, a, 1); /* no need to reduce mod p, since it will be done
                             in mpz_poly_pow_mod_f_mod_mpz */
  }

  mpz_poly_clear (q);
  mpz_poly_clear (h);

clear_a:
  mpz_clear (a);
  mpz_clear (aux);
  return n;
}


/* Solve f(x)=0 (mod p), where p is a prime. Return the number of roots.
   Assume d (the degree of f) is at least 1.
 */
int
mpz_poly_roots_mpz (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p)
{
  int nr = 0;
  mpz_poly fp, g, h;
  int d = f->deg;

  ASSERT(d >= 1);

  mpz_poly_init (fp, d);
  mpz_poly_init (g, 2*d-1);
  mpz_poly_init (h, 2*d-1);

  /* If f has small coefficients (like in Joux-Lercier polynomial selection)
     don't make f monic, since it might make the coefficients of fp blow up.
     In that case we only reduce coefficients in [-p+1, p-1], to keep
     negative coefficients small in absolute value. */
  if (mpz_poly_sizeinbase (f, 2) < mpz_sizeinbase (p, 2))
    mpz_poly_mod_mpz_lazy (fp, f, p, NULL);
  else
    mpz_poly_makemonic_mod_mpz (fp, f, p);
  if (fp->deg <= 0)
    goto clear_and_exit;
  /* h=x^p-x (mod mpz_poly_fp) */
  mpz_poly_setcoeff_ui (g, 1, 1);
  mpz_poly_pow_mod_f_mod_mpz (h, g, fp, p, p);
  /* FIXME: instead of computing x^p-x, we could compute x^(p-1) - 1 while
     saving the value of h = x^((p-1)/2). If several roots, gcd(h-1, f)
     might help to split them. */
  mpz_poly_sub (h, h, g);
  /* g = gcd (mpz_poly_fp, h) */
  mpz_poly_gcd_mpz (fp, fp, h, p);
  /* fp contains gcd(x^p-x, f) */
  nr = fp->deg;
  ASSERT (nr >= 0);

  /* If r is NULL, we only return the number of roots. */
  if (r != NULL && nr > 0)
  {
    int n MAYBE_UNUSED = mpz_poly_cantor_zassenhaus (r, fp, p, 0);
    ASSERT (n == nr);
  }

 clear_and_exit:
  mpz_poly_clear(fp);
  mpz_poly_clear(g);
  mpz_poly_clear(h);

  /* Sort the roots */
  if (r && nr)
    qsort(r, nr, sizeof(mpz_t), (sortfunc_t) &mpz_poly_coeff_cmp);

  return nr;
}

