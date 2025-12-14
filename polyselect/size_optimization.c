#include "cado.h" // IWYU pragma: keep

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <gmp.h>

#include "auxiliary.h"
#include "gcd.h"
#include "lll.h"
#include "size_optimization.h"
#include "double_poly.h"
#include "polyselect_norms.h"
#include "macros.h"
#include "mpz_poly.h"


/******************************************************************************/
/************************ Internal functions **********************************/
/******************************************************************************/

/***************** List of mpz_t, to handle list of translations **************/
// I've been fighting for years against these countless ad-hoc lists that
// we have here and there in the code. This one has the at least the good
// thing of not being used outside this compilation unit (and I
// definitely do not want it to be exposed outside).

static inline void
list_mpz_realloc (list_mpz_ptr l, uint64_t newalloc)
{
#ifdef DEBUG_LIST_MPZ
  fprintf (stderr, "debug: %s in %s: l->alloc = %" PRIu64 ", l->len = "
                   "%" PRIu64 ", will be reallocated to newalloc = %" PRIu64
                   "\n", __FILE__, __func__, l->alloc, l->len, newalloc);
#endif
  CHECKED_REALLOC(l->tab, newalloc, mpz_t);
  for (uint64_t i = l->alloc; i < newalloc; i++)
    mpz_init (l->tab[i]);
  l->alloc = newalloc;
}

static inline void
list_mpz_init (list_mpz_ptr l, uint64_t init_alloc)
{
  l->len = 0;
  l->alloc = 0;
  l->tab = NULL;
  list_mpz_realloc (l, init_alloc);
}

static inline void
list_mpz_clear (list_mpz_t l)
{
  for (uint64_t i = 0; i < l->alloc; i++)
    mpz_clear (l->tab[i]);
  free (l->tab);
}

static inline void
list_mpz_append (list_mpz_ptr l, mpz_t e)
{
  if (l->len == l->alloc)
    list_mpz_realloc (l, 2*l->alloc);
  mpz_set (l->tab[l->len], e);
  l->len++;
}

static inline void
list_mpz_append_from_rounded_double (list_mpz_ptr l, double e)
{
  if (l->len == l->alloc)
    list_mpz_realloc (l, 2*l->alloc);
  mpz_set_d (l->tab[l->len], e > 0 ? e + 0.5 : e - 0.5);
  l->len++;
}

static inline void
list_mpz_sort_and_remove_dup (list_mpz_ptr l, const int verbose)
{
  /* Sort list_k by increasing order */
  qsort (l->tab, l->len, sizeof (mpz_t), (int(*)(const void*,const void*)) mpz_cmp);

  /* Remove duplicates */
  uint64_t len = l->len;
  l->len = 1;
  for (unsigned int i = 1; i < len; i++)
    if (mpz_cmp (l->tab[i], l->tab[l->len-1]) != 0)
      list_mpz_append (l, l->tab[i]);
    else if (verbose > 2)
      gmp_fprintf (stderr, "# sopt: Remove duplicate %Zd\n", l->tab[i]);
}

/******************************************************************************/

/* store in Q[0..nb_approx-1] the (at most) nb_approx first best rational
   approximations of q with denominator <= bound.
   nb_approx should be greater than 0
   Return the number of found approximations. */
static inline unsigned int
compute_rational_approximation (double *Q, double q, unsigned int nb_approx,
                                double bound)
{
  unsigned int n = 0;
  double best_e = 2.0;

  for (double den = 1.0; n < nb_approx && den <= bound; den += 1.0)
  {
    double num = floor (den * q + 0.5);
    double e = fabs (q - num / den);
    if (e >= best_e)
      continue;
    best_e = e;
    Q[n++] = num / den;
  }
  return n;
}

/* Compute the skewness that is going to be used for LLL */
static inline void
sopt_get_skewness (mpz_t skew, mpz_poly_srcptr f, mpz_poly_srcptr g)
{
  const int d = f->deg;
  /* skew0 = (|g0/ad|)^(1/d): seems close to optimal experimentally */
  double skew0 = pow (fabs (mpz_get_d (mpz_poly_coeff_const(g, 0)) / mpz_get_d (mpz_poly_coeff_const(f, d))),
                      1.0 / (double) d);
  mpz_set_d (skew, skew0);
}

/* For deg(f) = 6 and deg(g) = 1 */
/* Return in list_k good values of k such that f(x+k) has small coefficients of
   degree 4 and 3. */
static void
sopt_find_translations_deg6 (list_mpz_t list_k, mpz_poly_srcptr f,
                             mpz_poly_srcptr g, const int verbose,
                             int sopt_effort)
{
  ASSERT_ALWAYS (f->deg == 6);
  ASSERT_ALWAYS (g->deg == 1);

  unsigned int i, nb_q3roots;

  /* Let f(x) = a6 x^6 + ... + a0 and g(x) = g1 x + g0. */
  double a6, a5, a4, a3, a2, g1, g0;
  a6 = mpz_get_d (mpz_poly_coeff_const(f, 6));
  a5 = mpz_get_d (mpz_poly_coeff_const(f, 5));
  a4 = mpz_get_d (mpz_poly_coeff_const(f, 4));
  a3 = mpz_get_d (mpz_poly_coeff_const(f, 3));
  a2 = mpz_get_d (mpz_poly_coeff_const(f, 2));
  g1 = mpz_get_d (mpz_poly_coeff_const(g, 1));
  g0 = mpz_get_d (mpz_poly_coeff_const(g, 0));

  double kmax;
  /* after translation by k, we have a2 = Theta(binomial(6,2)*a6*k^4)
     = Theta(15*a6*k^4), thus after reduction by x^2*g, a3 = Theta(15*a6*k^4*g1/g0).
     Since we want a3 << g0, a lower bound for k is |g0^2/(15*a6*g1)|^(1/4). */
  kmax = pow (fabs (g0 * g0 / (15.0 * a6 * g1)), 0.25);

  /* Let f~(x,k,q3) = f(x+k)+q3*x^3*g(x+k).
     Let c_i(k,q3) be the coefficient of x^i in f~. Then
        c_3(k,q3) = 20*a6*k^3 + 10*a5*k^2 + 4*a4*k + q3*g1*k + a3 + q3*g0
        c_4(k,q3) = 15*a6*k^2 + 5*a5*k + q3*g1 + a4

     Let res = Res (c_3(k,q3), c_4(k,q3)), where the resultant is taken to
     eliminate the variable k.
     We compute the roots of res and of the derivative of res.

     #Sage code to compute Res(c3(k), c4(k)):
        R.<x,k,q3,g1,g0,a6,a5,a4,a3,a2,a1,a0> = ZZ[]
        f = a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        g = g1*x + g0
        ff=f(x=x+k)+q3*x^3*g(x=x+k)
        c3=ff.coefficient({x:3})
        c4=ff.coefficient({x:4})
        r = c3.resultant(c4,k)//(25*a6)
        SR(r.coefficient({q3:3})).factor()
        SR(r.coefficient({q3:2})).factor()
        SR(r.coefficient({q3:1})).factor()
        SR(r.coefficient({q3:0})).factor()
  */

  double_poly res;
  double roots_q3[3]; /* roots for res */
  double_poly_init (res, 3);

  /* degree 3: a6*g1^3 */
  res->coeff[3] = a6 * g1 * g1 * g1;

 /* degree 2: 135*a6^2*g0^2 - 45*a5*a6*g0*g1 + 10*a5^2*g1^2 - 15*a4*a6*g1^2 */
  res->coeff[2] = 135.0 * a6 * a6 * g0 * g0 - 45.0 * a5 * a6 * g0 * g1
                  + 10.0 * a5 * a5 * g1 * g1 - 15.0 * a4 * a6 * g1 * g1;

  /* degree 1: 50*a5^3*g0 - 180*a4*a5*a6*g0 + 270*a3*a6^2*g0 - 10*a4*a5^2*g1
     + 48*a4^2*a6*g1 - 45*a3*a5*a6*g1 */
  res->coeff[1] = 50.0 * a5 * a5 * a5 * g0 - 180.0 * a4 * a5 * a6 * g0
                  + 270.0 * a3 * a6 * a6 * g0 - 10.0 * a4 * a5 * a5 * g1
                  + 48.0 * a4 * a4 * a6 * g1 - 45.0 * a3 * a5 * a6 * g1;

  /* degree 0: -20*a4^2*a5^2 + 50*a3*a5^3 + 64*a4^3*a6 - 180*a3*a4*a5*a6
     + 135*a3^2*a6^2 */
  res->coeff[0] = -20.0 * a4 * a4 * a5 * a5 + 50.0 * a3 * a5 * a5 * a5
                  + 64.0 * a4 * a4 * a4 * a6 - 180.0 * a3 * a4 * a5 * a6
                  + 135.0 * a3 * a3 * a6 * a6;
  double_poly_cleandeg(res, 3);

  /* compute real roots of res:
     (1) if a6*a4 < 0, usually we have 2 "small" real roots (almost opposite)
         and one huge real root
     (2) if a6*a4 > 0, usually we only have one huge real root */
  nb_q3roots = double_poly_compute_all_roots_with_bound (roots_q3, res,
						 SOPT_MAX_VALUE_FOR_Q_ROOTS);

  if (verbose > 1)
  {
    fprintf (stderr, "# sopt: q-roots of Res(c3,c4) or Res(c3,c4)' = {");
    for (i = 0; i < nb_q3roots; i++)
      fprintf (stderr, " %f%c", roots_q3[i], (i+1==nb_q3roots)?' ':',');
    fprintf (stderr, "}\n");
  }

  double_poly C;
  double_poly_init (C, 4);
  int a = 1, b = 1;
  /* for an average of 2 roots per polynomial, using sopt_effort > 0 will increase
     the size optimization time by a ratio (sopt_effort + 1) on average */
  for (i = 0; i < nb_q3roots + 2 * sopt_effort; i++)
  {
    /* (1) for i < nb_q3roots, we deal with the i-th root nb_q3roots[i], and generate
           up to SOPT_NB_RAT_APPROX_OF_Q_ROOTS rational approximations of it;
       (2) for i >= nb_q3roots, we consider fractions a/b, gcd(a,b) = 1, |a/b| <= 1 */

    unsigned int t;
    double q3_rat_approx[SOPT_NB_RAT_APPROX_OF_Q_ROOTS];

    if (i < nb_q3roots)
      t = compute_rational_approximation (q3_rat_approx, roots_q3[i],
                                          SOPT_NB_RAT_APPROX_OF_Q_ROOTS,
                                          SOPT_MAX_DEN_IN_RAT_APPROX_OF_Q_ROOTS);
    else
      t = SOPT_NB_RAT_APPROX_OF_Q_ROOTS;

    for (unsigned int j = 0; j < t; j++)
    {
      double q3_rat = q3_rat_approx[j];
      double double_roots_k[5];
      unsigned int nb_k_roots;

      double_poly_realloc(C, 5);

      if (i >= nb_q3roots)
        {
          q3_rat = (double) a / (double) b;
          if (a > 0)
            a = -a;
          else /* generate next Farey fraction:
                  1/1 -> -1/1 -> 1/2 -> -1/2 -> 1/3 -> -1/3 -> 2/3 ->
                  -2/3 -> 1/4 -> ... */
            {
              a = (-a) + 1; /* now a > 0 */
              while (gcd_uint64 (a, b) != 1)
                a ++;
              if (a > b)
                {
                  a = 1;
                  b = b + 1;
                }
            }
        }

      if (verbose > 2)
        {
          if (i < nb_q3roots)
            fprintf (stderr, "# sopt: process rational approximation %f of %f\n",
                     q3_rat, roots_q3[i]);
          else
            fprintf (stderr, "# sopt: process rational approximation %f\n", q3_rat);
        }

      /* find roots k of c4(k, q3=q3_rat) */
      C->coeff[2] = 15.0 * a6;
      C->coeff[1] = 5.0 * a5;
      C->coeff[0] = g1 * q3_rat + a4;
      double_poly_cleandeg(C, 2);

      nb_k_roots = double_poly_compute_all_roots_with_bound (double_roots_k, C, kmax);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose > 2)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c4)\n",
                       list_k->tab[list_k->len-1]);
      }

      /* find roots k of c3(k, q3=q3_rat) */
      C->coeff[3] = 20.0 * a6;
      C->coeff[2] = 10.0 * a5;
      C->coeff[1] = g1 * q3_rat + 4.0 * a4;
      C->coeff[0] = q3_rat * g0 + a3;
      double_poly_cleandeg(C, 3);

      nb_k_roots = double_poly_compute_all_roots_with_bound (double_roots_k, C, kmax);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose > 2)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c3)\n",
                       list_k->tab[list_k->len-1]);
      }

      /* Let f~~(x,k,q3,q2) f(x+k) + (q3*x^3 + q2*x^2)*g(x,k)
         Add roots of Res(c3(k, q3=q3_rat, q2), c2(k, q3_q3rat, q2)) and its
         derivative, where the resultant is taken to eliminate the variable q2.

         #Sage code to compute Res(c3(k,q3,q2), c2(k,q3,q2)):
            R.<a6,k,a5,g1,q,q2,g0,a3,a4,a2,x> = ZZ[]
            f = a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2
            g = g1*x+g0
            ff = f(x=x+k)+q*x^3*g(x=x+k)+q2*x^2*g(x=x+k)
            c3 = ff.coefficient({x:3})
            c2 = ff.coefficient({x:2})
            c3.resultant(c2,q2)
      */
      C->coeff[4] = -5.0 * a6 * g1;
      C->coeff[3] = -20.0 * a6 * g0;
      C->coeff[2] = -g1 * g1 * q3_rat - 10.0 * a5 * g0 + 2.0 * g1 * a4;
      C->coeff[1] = 2.0 * g1 * a3 - 2.0 * g1 * q3_rat * g0 - 4.0 * g0 * a4;
      C->coeff[0] = -q3_rat * g0 * g0 - g0 * a3 + g1 * a2;
      double_poly_cleandeg(C, 4);

      nb_k_roots = double_poly_compute_all_roots_with_bound (double_roots_k, C, kmax);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose > 2)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of "
                               "Res(c3,c2))\n", list_k->tab[list_k->len-1]);
      }

      double_poly_derivative (C, C);
      nb_k_roots = double_poly_compute_all_roots_with_bound (double_roots_k, C, kmax);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose > 2)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of derivative "
                               "of Res(c3,c2)\n", list_k->tab[list_k->len-1]);
      }

    }
// TODO use the following code.
#if 0
    double_poly_init (C, 4);
    int n, i, j, t, m, ret = 0;
    //double_poly R, C;
#define MAX_ROOTSK 40
    //double roots_k[MAX_ROOTSK];
    //double *r, roots_q[5];
    //R->deg = 3;
    //r = R->coeff;
      while (t) {
      roots_q[i] = Q[--t];


      /* add roots of Res(c3(q3=roots_q[i](k),c2(k)) (eliminate k) */
      double m0 = -g0;
      double roots_q2[7];
      int ii, nb_roots_q2 = 0;

      C->deg = 4;
      C->coeff[4] = -125*a6*a6*a6*g1*g1*g1*g1;
      C->coeff[3] = (3000*a6*a6*a6*m0*g1*g1*g1 + 500*a5*a6*a6*g1*g1*g1*g1)*roots_q[i] - 160000*a6*a6*a6*a6*m0*m0*m0 - 80000*a5*a6*a6*a6*m0*m0*g1 - 5000*a5*a5*a6*a6*m0*g1*g1 - 20000*a4*a6*a6*a6*m0*g1*g1 - 1000*a4*a5*a6*a6*g1*g1*g1 - 3500*a3*a6*a6*a6*g1*g1*g1;
      C->coeff[2] = -225*a6*a6*g1*g1*g1*g1*g1*roots_q[i]*roots_q[i]*roots_q[i] + (15750*a6*a6*a6*m0*m0*g1*g1 + 5250*a5*a6*a6*m0*g1*g1*g1 - 500*a5*a5*a6*g1*g1*g1*g1 + 2250*a4*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i] + (20000*a5*a5*a6*a6*m0*m0*g1 - 48000*a4*a6*a6*a6*m0*m0*g1 + 10000*a5*a5*a5*a6*m0*g1*g1 - 28000*a4*a5*a6*a6*m0*g1*g1 + 18000*a3*a6*a6*a6*m0*g1*g1 + 2000*a4*a5*a5*a6*g1*g1*g1 - 7200*a4*a4*a6*a6*g1*g1*g1 + 5000*a3*a5*a6*a6*g1*g1*g1 - 4000*a2*a6*a6*a6*g1*g1*g1)*roots_q[i] - 50000*a5*a5*a5*a5*a6*m0*m0 + 240000*a4*a5*a5*a6*a6*m0*m0 - 192000*a4*a4*a6*a6*a6*m0*m0 - 240000*a3*a5*a6*a6*a6*m0*m0 + 480000*a2*a6*a6*a6*a6*m0*m0 - 20000*a4*a5*a5*a5*a6*m0*g1 + 80000*a4*a4*a5*a6*a6*m0*g1 + 10000*a3*a5*a5*a6*a6*m0*g1 - 216000*a3*a4*a6*a6*a6*m0*g1 + 160000*a2*a5*a6*a6*a6*m0*g1 - 2000*a4*a4*a5*a5*a6*g1*g1 + 7200*a4*a4*a4*a6*a6*g1*g1 - 1000*a3*a4*a5*a6*a6*g1*g1 + 5000*a2*a5*a5*a6*a6*g1*g1 - 33750*a3*a3*a6*a6*a6*g1*g1 + 20000*a2*a4*a6*a6*a6*g1*g1;
      C->coeff[1] = (-18000*a6*a6*a6*m0*m0*m0*g1 - 9000*a5*a6*a6*m0*m0*g1*g1 - 3600*a4*a6*a6*m0*g1*g1*g1 - 900*a3*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-15000*a5*a5*a6*a6*m0*m0*m0 + 36000*a4*a6*a6*a6*m0*m0*m0 - 20000*a5*a5*a5*a6*m0*m0*g1 + 63000*a4*a5*a6*a6*m0*m0*g1 - 67500*a3*a6*a6*a6*m0*m0*g1 - 8000*a4*a5*a5*a6*m0*g1*g1 + 25200*a4*a4*a6*a6*m0*g1*g1 - 7500*a3*a5*a6*a6*m0*g1*g1 - 30000*a2*a6*a6*a6*m0*g1*g1 - 2000*a3*a5*a5*a6*g1*g1*g1 + 6300*a3*a4*a6*a6*g1*g1*g1 - 5000*a2*a5*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i] + (8000*a4*a4*a5*a5*a6*m0*g1 - 20000*a3*a5*a5*a5*a6*m0*g1 - 28800*a4*a4*a4*a6*a6*m0*g1 + 84000*a3*a4*a5*a6*a6*m0*g1 - 20000*a2*a5*a5*a6*a6*m0*g1 - 81000*a3*a3*a6*a6*a6*m0*g1 + 48000*a2*a4*a6*a6*a6*m0*g1 + 2000*a3*a4*a5*a5*a6*g1*g1 - 10000*a2*a5*a5*a5*a6*g1*g1 - 7200*a3*a4*a4*a6*a6*g1*g1 + 4500*a3*a3*a5*a6*a6*g1*g1 + 32000*a2*a4*a5*a6*a6*g1*g1 - 36000*a2*a3*a6*a6*a6*g1*g1)*roots_q[i] + 16000*a4*a4*a4*a5*a5*a6*m0- 60000*a3*a4*a5*a5*a5*a6*m0+ 100000*a2*a5*a5*a5*a5*a6*m0- 57600*a4*a4*a4*a4*a6*a6*m0+ 240000*a3*a4*a4*a5*a6*a6*m0+ 15000*a3*a3*a5*a5*a6*a6*m0- 480000*a2*a4*a5*a5*a6*a6*m0- 324000*a3*a3*a4*a6*a6*a6*m0+ 384000*a2*a4*a4*a6*a6*a6*m0+ 480000*a2*a3*a5*a6*a6*a6*m0- 480000*a2*a2*a6*a6*a6*a6*m0+ 4000*a3*a4*a4*a5*a5*a6*g1 - 20000*a3*a3*a5*a5*a5*a6*g1 + 20000*a2*a4*a5*a5*a5*a6*g1 - 14400*a3*a4*a4*a4*a6*a6*g1 + 81000*a3*a3*a4*a5*a6*a6*g1 - 80000*a2*a4*a4*a5*a6*a6*g1 - 10000*a2*a3*a5*a5*a6*a6*g1 - 121500*a3*a3*a3*a6*a6*a6*g1 + 216000*a2*a3*a4*a6*a6*a6*g1 - 80000*a2*a2*a5*a6*a6*a6*g1;
      C->coeff[0] = (3375*a6*a6*a6*m0*m0*m0*m0 + 2250*a5*a6*a6*m0*m0*m0*g1 + 1350*a4*a6*a6*m0*m0*g1*g1 + 675*a3*a6*a6*m0*g1*g1*g1 + 225*a2*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i]*roots_q[i] + (5000*a5*a5*a5*a6*m0*m0*m0 - 18000*a4*a5*a6*a6*m0*m0*m0 + 27000*a3*a6*a6*a6*m0*m0*m0 + 3000*a4*a5*a5*a6*m0*m0*g1 - 10800*a4*a4*a6*a6*m0*m0*g1 + 4500*a3*a5*a6*a6*m0*m0*g1 + 18000*a2*a6*a6*a6*m0*m0*g1 + 1500*a3*a5*a5*a6*m0*g1*g1 - 5400*a3*a4*a6*a6*m0*g1*g1 + 6000*a2*a5*a6*a6*m0*g1*g1 + 500*a2*a5*a5*a6*g1*g1*g1 - 675*a3*a3*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-6000*a4*a4*a5*a5*a6*m0*m0 + 15000*a3*a5*a5*a5*a6*m0*m0 + 21600*a4*a4*a4*a6*a6*m0*m0 - 63000*a3*a4*a5*a6*a6*m0*m0 + 15000*a2*a5*a5*a6*a6*m0*m0 + 60750*a3*a3*a6*a6*a6*m0*m0 - 36000*a2*a4*a6*a6*a6*m0*m0 - 3000*a3*a4*a5*a5*a6*m0*g1 + 15000*a2*a5*a5*a5*a6*m0*g1 + 10800*a3*a4*a4*a6*a6*m0*g1 - 6750*a3*a3*a5*a6*a6*m0*g1 - 48000*a2*a4*a5*a6*a6*m0*g1 + 54000*a2*a3*a6*a6*a6*m0*g1 - 1500*a3*a3*a5*a5*a6*g1*g1 + 3000*a2*a4*a5*a5*a6*g1*g1 + 4050*a3*a3*a4*a6*a6*g1*g1 - 7200*a2*a4*a4*a6*a6*g1*g1 - 3000*a2*a3*a5*a6*a6*g1*g1 + 12000*a2*a2*a6*a6*a6*g1*g1)*roots_q[i]*roots_q[i] + 6000*a3*a3*a4*a4*a5*a5*a6 - 16000*a2*a4*a4*a4*a5*a5*a6 - 20000*a3*a3*a3*a5*a5*a5*a6 + 60000*a2*a3*a4*a5*a5*a5*a6 - 50000*a2*a2*a5*a5*a5*a5*a6 - 21600*a3*a3*a4*a4*a4*a6*a6 + 57600*a2*a4*a4*a4*a4*a6*a6 + 81000*a3*a3*a3*a4*a5*a6*a6 - 240000*a2*a3*a4*a4*a5*a6*a6 - 15000*a2*a3*a3*a5*a5*a6*a6 + 240000*a2*a2*a4*a5*a5*a6*a6 - 91125*a3*a3*a3*a3*a6*a6*a6 + 324000*a2*a3*a3*a4*a6*a6*a6 - 192000*a2*a2*a4*a4*a6*a6*a6 - 240000*a2*a2*a3*a5*a6*a6*a6 + 160000*a2*a2*a2*a6*a6*a6*a6;
      double_poly_cleandeg(C, 4);

      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2, C,
                                                               SOPT_MAX_VALUE_FOR_Q_ROOTS);
      ASSERT_ALWAYS(nb_roots_q2 <= 4);

      /* roots of the derivative of Res(c3,c2) */
      C->coeff[0] = C->coeff[1];
      C->coeff[1] = C->coeff[2] * 2.0;
      C->coeff[2] = C->coeff[3] * 3.0;
      C->coeff[3] = C->coeff[4] * 4.0;
      double_poly_cleandeg(C, 3);
      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2
                                                  + nb_roots_q2, C, SOPT_MAX_VALUE_FOR_Q_ROOTS);
      ASSERT_ALWAYS(nb_roots_q2 <= 7);

      for (ii = 0; ii < nb_roots_q2; ii++)
      {
        /* Compute roots of c2(q2=roots_q2[ii]) */
        C->coeff[4] = 15.0 * a6;
        C->coeff[3] = 10.0 * a5;
        C->coeff[2] = 6.0 * a4;
        C->coeff[1] = g1 * roots_q2[ii] + 3.0 * a3;
        C->coeff[0] = a2 + g0 * roots_q2[ii];
        double_poly_cleandeg(C, 4);

        m += double_poly_compute_all_roots (roots_k + m, C);
        ASSERT_ALWAYS(nb_roots_q2 <= 12+4*(ii+1));
      }
#endif
      }
  double_poly_clear (C);
  double_poly_clear (res);
}

/* For deg(f) = 5 and deg(g) = 1 */
/* Return in list_k good values of k such that f(x+k) has small coefficients of
   degree 3 and 2. */
static void
sopt_find_translations_deg5 (list_mpz_t list_k, mpz_poly_srcptr f,
                             mpz_poly_srcptr g, const int verbose,
                             int sopt_effort)
{
  ASSERT_ALWAYS (f->deg == 5);
  ASSERT_ALWAYS (g->deg == 1);

  unsigned int i, nb_q2roots;

  /* Let f(x) = a5 x^5 + ... + a0 and g(x) = g1 x + g0. */
  double a5, a4, a3, a2, g1, g0;
  a5 = mpz_get_d (mpz_poly_coeff_const(f, 5));
  a4 = mpz_get_d (mpz_poly_coeff_const(f, 4));
  a3 = mpz_get_d (mpz_poly_coeff_const(f, 3));
  a2 = mpz_get_d (mpz_poly_coeff_const(f, 2));
  g1 = mpz_get_d (mpz_poly_coeff_const(g, 1));
  g0 = mpz_get_d (mpz_poly_coeff_const(g, 0));

  /* Let f~(x,k,q2) = f(x+k)+q2*x^2*g(x+k).
     Let c_i(k,q2) be the coefficient of x^i in f~. Then
        c_2(k,q2) = 10*a5*k^3 + 6*a4*k^2 + 3*a3*k + q2*g1*k + a2 + q2*g0
        c_3(k,q2) = 10*a5*k^2 + 4*a4*k + q2*g1 + a3

     Let res = Res (c_2(k,q3), c_3(k,q3)), where the resultant is taken to
     eliminate the variable k.
     We compute the roots of res and of the derivative of res.
     Degree 5 is the only case where res is of degree 2 in q2 (instead of 3)

     #Sage code to compute Res(c2(k), c3(k)):
        R.<x,k,q2,g1,g0,a5,a4,a3,a2,a1,a0> = ZZ[]
        f = a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        g = g1*x + g0
        ff=f(x=x+k)+q2*x^2*g(x=x+k)
        c2=ff.coefficient({x:2})
        c3=ff.coefficient({x:3})
        r = c2.resultant(c3,k)//(40*a5)
        SR(r.coefficient({q2:2})).factor()
        SR(r.coefficient({q2:1})).factor()
        SR(r.coefficient({q2:0})).factor()
  */

  double_poly res;
  double roots_q2[3]; /* 3 = 2 roots for res + 1 roots for the derivative */
  double_poly_init (res, 2);

  /* degree 2: (5*a5*g0 - a4*g1)^2 */
  res->coeff[2] =  (5.0*a5*g0 - a4*g1)*(5.0*a5*g0 - a4*g1);

  /* degree 1: 8*a4^3*g0 - 30*a3*a4*a5*g0 + 50*a2*a5^2*g0 - 2*a3*a4^2*g1
               + 10*a3^2*a5*g1 - 10*a2*a4*a5*g1 */
  res->coeff[1] = 8.0*a4*a4*a4*g0 - 30.0*a3*a4*a5*g0 + 50.0*a2*a5*a5*g0
                  - 2.0*a3*a4*a4*g1 + 10.0*a3*a3*a5*g1 - 10.0*a2*a4*a5*g1;

  /* degree 0: -3*a3^2*a4^2 + 8*a2*a4^3 + 10*a3^3*a5 - 30*a2*a3*a4*a5
               + 25*a2^2*a5^2 */
  res->coeff[0] = -3.0*a3*a3*a4*a4 + 8.0*a2*a4*a4*a4 + 10.0*a3*a3*a3*a5
                  - 30.0*a2*a3*a4*a5 + 25.0*a2*a2*a5*a5;
  double_poly_cleandeg(res, 2);

  /* compute roots of res */
  nb_q2roots = double_poly_compute_all_roots_with_bound (roots_q2, res,
						   SOPT_MAX_VALUE_FOR_Q_ROOTS);
  if (verbose > 1)
  {
    fprintf (stderr, "# sopt: q-roots of Res(c2,c3) or Res(c2,c3)' = {");
    for (i = 0; i < nb_q2roots; i++)
      fprintf (stderr, " %f%c", roots_q2[i], (i+1==nb_q2roots)?' ':',');
    fprintf (stderr, "}\n");
  }

  double_poly C;
  double_poly_init (C, 2);
  int a = 1, b = 1;
  /* for an average of 1 root per polynomial, using sopt_effort > 0 will increase
     the size optimization time by a ratio (sopt_effort + 1) on average */
  for (i = 0; i < nb_q2roots + sopt_effort; i++)
  {
    unsigned int t;
    double q2_rat_approx[SOPT_NB_RAT_APPROX_OF_Q_ROOTS];

    if (i < nb_q2roots)
    t = compute_rational_approximation (q2_rat_approx, roots_q2[i],
                                        SOPT_NB_RAT_APPROX_OF_Q_ROOTS,
                                        SOPT_MAX_DEN_IN_RAT_APPROX_OF_Q_ROOTS);
    else
      t = SOPT_NB_RAT_APPROX_OF_Q_ROOTS;

    for (unsigned int j = 0; j < t; j++)
    {
      double q2_rat = q2_rat_approx[j];
      double double_roots_k[5];
      unsigned int nb_k_roots;

      if (i >= nb_q2roots)
        {
          q2_rat = (double) a / (double) b;
          if (a > 0)
            a = -a;
          else
            {
              a = (-a) + 1; /* now a > 0 */
              while (gcd_uint64 (a, b) != 1)
                a ++;
              if (a > b)
                {
                  a = 1;
                  b = b + 1;
                }
            }
        }

      if (verbose > 2)
        fprintf (stderr, "# sopt: process rational approximation %f of %f\n",
                         q2_rat, roots_q2[i]);

      C->coeff[2] = 10.0 * a5;
      C->coeff[1] = 4.0 * a4;
      C->coeff[0] = g1 * q2_rat + a3;
      double_poly_cleandeg(C, 2);

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose > 2)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c3)\n",
                       list_k->tab[list_k->len-1]);
      }
    }
  }
  double_poly_clear (C);
  double_poly_clear (res);
}

/* Construct the LLL matrix from the polynomial pair. */
static inline void
LLL_set_matrix_from_polys (mat_Z *m, mpz_poly_srcptr ft, mpz_poly_srcptr gt,
                           mpz_srcptr skew, mpz_ptr tmp)
{
  for (int k = 1; k <= m->NumRows; k++)
    for (int l = 1; l <= m->NumCols; l++)
      mpz_set_ui (m->coeff[k][l], 0);
  mpz_set_ui (tmp, 1);
  for (int k = 0; k < m->NumCols; k++)
  {
    if (k > 0)
      mpz_mul (tmp, tmp, skew); /* skew^k */
    if (k <= m->NumRows - 2)
      mpz_set (m->coeff[k+2][k+1], tmp);
    mpz_mul (m->coeff[1][k+1], tmp, mpz_poly_coeff_const(ft, k));
  }
  for (int k = 0; k <= m->NumRows - 2; k++)
  {
    /* m.coeff[k+2][k+1] is already skew^k */
    mpz_mul (m->coeff[k+2][k+2], m->coeff[k+2][k+1], skew); /* skew^(k+1) */
    mpz_mul (m->coeff[k+2][k+1], m->coeff[k+2][k+1], mpz_poly_coeff_const(gt, 0));
    mpz_mul (m->coeff[k+2][k+2], m->coeff[k+2][k+2], mpz_poly_coeff_const(gt, 1));
  }
}

/* Print only the two coefficients of higher degree of a polynomial */
static inline void
mpz_poly_fprintf_short (FILE *out, mpz_poly_srcptr f)
{
  int d = f->deg;
  ASSERT_ALWAYS (d >= 1);
  gmp_fprintf (out, "%Zd*x^%d + %Zd*x^%d%s\n", mpz_poly_coeff_const(f, d), d,
                    mpz_poly_coeff_const(f, d-1), d-1, (d > 1) ? " + ...":"");
}

/* Print poly: use mpz_poly_fprintf or mpz_poly_fprintf_short depending on
   verbose */
static inline void
mpz_poly_fprintf_verbose (FILE *out, mpz_poly_srcptr f, int verbose)
{
  if (verbose > 2)
    mpz_poly_fprintf (out, f);
  else
    mpz_poly_fprintf_short (out, f);
}

/******************************************************************************/
/******************************************************************************/
/********************** End Internal functions ********************************/
/******************************************************************************/

/* Use rotation and translation to find a polynomial with smaller norm.
   Use local descent algorithm to find a local minimum.
   Maximum degree for the rotation is given by deg_rotation (deg_rotation < 0
   means no rotation).
   Translations are used only if use_translation is non-zero.
   Will look for polynomials of the form
   [f_raw + (k[deg_rotation]*x^deg_rotation + ... + k[1]*x + k[0])*g_raw](x+kt)
   and g(x+kt)
   Optimized polynomials are returned in f_opt and g_opt, and the function
   returns the skew lognorm of f_opt.

   To use only translation: use_translation = 1 and deg_rotation = -1
   To use only rotation   : use_translation = 0 and deg_rotation >= 0
   To use both            : use_translation = 1 and deg_rotation >= 0
   TODO: _mp version like old optimize_aux_mp in auxiliary.c
   XXX: Could we replace lognorm computation by norm computation to save
        computation of the log. Compute lognorm only for printing.
*/
double
sopt_local_descent (mpz_poly_ptr f_opt, mpz_poly_ptr g_opt,
                    mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw,
                    int use_translation, int deg_rotation,
                    unsigned int max_iter, int verbose)
{
  mpz_t kt, k[SOPT_MAX_DEGREE_ROTATION+1]; /* current offset */
  mpz_t tmp;
  int changedt, changed[SOPT_MAX_DEGREE_ROTATION+1];
  double logmu_opt, logmu;
  unsigned int iter = 0;
  mpz_poly ftmp;

  ASSERT_ALWAYS(deg_rotation <= SOPT_MAX_DEGREE_ROTATION);

  logmu_opt = L2_skew_lognorm ((mpz_poly_ptr) f_raw);
  mpz_init (tmp);
  mpz_poly_init (ftmp, f_raw->deg);

  if (verbose > 2)
  {
    fprintf (stderr, "# sopt:       starting local descent with polynomials:\n"
                     "# sopt:       f_raw = ");
    mpz_poly_fprintf_verbose (stderr, f_raw, verbose);
    fprintf (stderr, "# sopt:       g_raw = ");
    mpz_poly_fprintf_verbose (stderr, g_raw, verbose);
    fprintf (stderr, "# sopt:       with lognorm = %f\n", logmu_opt);
  }

  /* initialize k[i] to the smallest power of two that increases the lognorm
     by SOPT_LOCAL_DESCENT_GUARD with respect to the initial polynomial, to
     avoid being stuck in a very flat region */
  for (int i = 0; i <= deg_rotation; i++)
  {
    mpz_init_set_ui (k[i], 1);
    while (1)
    {
      mpz_poly_rotation (ftmp, f_raw, g_raw, k[i], i);
      logmu = L2_skew_lognorm (ftmp);
      if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
      {
        mpz_neg (tmp, k[i]);
        mpz_poly_rotation (ftmp, f_raw, g_raw, tmp, i);
        logmu = L2_skew_lognorm (ftmp);
        if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
          break;
      }
    mpz_mul_2exp (k[i], k[i], 1);
    }
  }

  /* initialize kt likewise */
  mpz_init_set_ui (kt, 1);
  while (1)
  {
    mpz_poly_translation (ftmp, f_raw, kt);
    logmu = L2_skew_lognorm (ftmp);
    if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
    {
      mpz_neg (tmp, kt);
      mpz_poly_translation (ftmp, f_raw, tmp);
      logmu = L2_skew_lognorm (ftmp);
      if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
        break;
    }
    mpz_mul_2exp (kt, kt, 1);
  }

  mpz_poly_set (f_opt, f_raw);
  mpz_poly_set (g_opt, g_raw);

  /* abort in cases where the descent procedure takes too long to converge;
   * 300 iterations seems largely enough in most cases */
  while (iter < max_iter)
  {
    iter++;
    changedt = 0;
    if (deg_rotation >= 0)
      memset (changed, 0, (deg_rotation+1) * sizeof (int));

    if (use_translation != 0)
    {
      /* first try translation by kt */
      mpz_poly_translation (ftmp, f_opt, kt); /* f(x+kt) */
      logmu = L2_skew_lognorm (ftmp);
      if (logmu < logmu_opt)
      {
        changedt = 1;
        logmu_opt = logmu;
        mpz_poly_swap (f_opt, ftmp);
        mpz_poly_translation (g_opt, g_opt, kt);
      }
      else
      {
        mpz_neg (tmp, kt);
        mpz_poly_translation (ftmp, f_opt, tmp); /* f(x-kt) */
        logmu = L2_skew_lognorm (ftmp);
        if (logmu < logmu_opt)
        {
          changedt = 1;
          logmu_opt = logmu;
          mpz_poly_swap (f_opt, ftmp);
          mpz_poly_translation (g_opt, g_opt, tmp);
        }
      }
    }

    for (int i = 0; i <= deg_rotation; i++)
    {
      mpz_poly_rotation (ftmp, f_opt, g_opt, k[i], i); /* f + k[i]*x^i*g */
      logmu = L2_skew_lognorm (ftmp);
      if (logmu < logmu_opt)
      {
        changed[i] = 1;
        logmu_opt = logmu;
        mpz_poly_swap (f_opt, ftmp);
      }
      else
      {
        mpz_neg (tmp, k[i]);
        mpz_poly_rotation (ftmp, f_opt, g_opt, tmp, i); /* f - k[i]*x^i*g */
        logmu = L2_skew_lognorm (ftmp);
        if (logmu < logmu_opt)
        {
          changed[i] = 1;
          logmu_opt = logmu;
          mpz_poly_swap (f_opt, ftmp);
        }
      }
    }

    int is_finished = 1;
    if (use_translation)
      is_finished = (changedt == 0 && mpz_cmp_ui (kt, 1) == 0);
    for (int i = 0; (i <= deg_rotation) && is_finished; i++)
      is_finished = (changed[i] == 0 && mpz_cmp_ui (k[i], 1) == 0);

    if (is_finished)
      break;

    if (changedt == 1)
      mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
    else if (mpz_cmp_ui (kt, 1) > 0)
      mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */

    for (int i = 0; i <= deg_rotation; i++)
    {
      if (changed[i] == 1)
        mpz_mul_2exp (k[i], k[i], 1);       /* k[i] <- 2*k[i] */
      else if (mpz_cmp_ui (k[i], 1) > 0)
        mpz_div_2exp (k[i], k[i], 1);       /* k[i] <- k[i]/2 */
    }
  }

  if (verbose > 2)
  {
    fprintf (stderr, "# sopt:       end of local descent after %u iterations "
                     "with polynomials:\n" "# sopt:       f_opt = ", iter);
    mpz_poly_fprintf_verbose (stderr, f_opt, verbose);
    fprintf (stderr, "# sopt:       g_opt = ");
    mpz_poly_fprintf_verbose (stderr, g_opt, verbose);
    fprintf (stderr, "# sopt:       with lognorm = %f\n", logmu_opt);
  }

  for (int i = 0; i <= deg_rotation; i++)
    mpz_clear (k[i]);
  mpz_clear (kt);
  mpz_clear (tmp);
  mpz_poly_clear (ftmp);

  return logmu_opt;
}

/* for given polynomials, skewness and translation k,
   returns the best norm found (using LLL reduction only)
   and the corresponding best polynomials in fopt, gopt */
static double
best_norm (mpz_poly_ptr fopt, mpz_poly_ptr gopt,
           mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw, mpz_t skew, mpz_t k,
           const int max_rot)
{
  const int d = f_raw->deg;
  mat_Z m;
  mpz_t tmp, a, b;
  double min_norm = DBL_MAX, norm;
  mpz_poly ft, gt;

  mpz_poly_init (ft, d);
  mpz_poly_init (gt, 1);

  LLL_init (&m, max_rot + 2, d + 1);
  mpz_init (tmp);
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 1);

  ASSERT_ALWAYS(m.NumCols == d+1);

  mpz_poly_translation (ft, f_raw, k);
  mpz_poly_translation (gt, g_raw, k);
  LLL_set_matrix_from_polys (&m, ft, gt, skew, tmp);
  LLL (tmp, m, NULL, a, b);
  for (int l = 1; l <= m.NumRows; l++)
    if (mpz_sgn (m.coeff[l][d+1]) != 0) /* degree is d */
      {
        /* compute the square norm in tmp */
        mpz_set_ui (tmp, 0);
        for (int j = 1; j <= m.NumCols; j++)
          mpz_addmul (tmp, m.coeff[l][j], m.coeff[l][j]);
        norm = mpz_get_d (tmp);
        if (norm < min_norm)
          {
            mpz_set_ui (tmp, 1);
            for (int j = 0; j <= d; j++)
              {
                /* invariant: tmp = skew^i */
                mpz_divexact (mpz_poly_coeff(ft, j), m.coeff[l][j+1], tmp);
                mpz_mul (tmp, tmp, skew);
              }
            min_norm = norm;
            mpz_poly_set (fopt, ft);
            mpz_poly_set (gopt, gt);
          }
      }

  mpz_clear (tmp);
  mpz_clear (a);
  mpz_clear (b);
  LLL_clear (&m);
  mpz_poly_clear (ft);
  mpz_poly_clear (gt);

  return min_norm;
}

/* Iterates calls to best_norm() and minimization of degree-(d-2) coefficient
   until a local minimum is found.
   Returns in k the best k-value, and in fopt, gopt the corresponding
   polynomials.
   The input min_norm is the best norm found so far.
*/
static double
best_norm2 (mpz_poly_ptr fopt, mpz_poly_ptr gopt,
            mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw, mpz_t skew, mpz_t k,
            const int max_rot, double min_norm)
{
  double norm, roots[3];
  double_poly eq;
  int d = f_raw->deg, nroots;
  mpz_t best_k;
  mpz_poly ft, gt;

  mpz_poly_init (ft, d);
  mpz_poly_init (gt, 1);

  mpz_init_set (best_k, k);
  norm = best_norm (ft, gt, f_raw, g_raw, skew, k, max_rot);
  if (norm >= min_norm)
    goto clear_ft_gt;
  min_norm = norm;
  mpz_set (best_k, k);
  mpz_poly_set (fopt, ft);
  mpz_poly_set (gopt, gt);

  /* solve d*(d-1)/2*f[d]*k^2 + (d-1)*f[d-1]*k + f[d-2] = 0 */
  double_poly_init (eq, 2);
  eq->coeff[2] = (double) (d * (d-1)) / 2.0 * mpz_get_d (mpz_poly_coeff_const(ft, d));
  eq->coeff[1] = (double) (d-1) * mpz_get_d (mpz_poly_coeff(ft, d-1));
  eq->coeff[0] = mpz_get_d (mpz_poly_coeff(ft, d-2));
  double_poly_cleandeg(eq, 2);
  nroots = double_poly_compute_all_roots (roots, eq);
  for (int l = 0; l < nroots; l++)
    {
      /* mpz_set_d rounds towards zero, thus we add or subtract 0.5 to get
         rounding to nearest */
      mpz_set_d (k, roots[l] >= 0 ? roots[l] + 0.5 : roots[l] - 0.5);
      norm = best_norm2 (ft, gt, f_raw, g_raw, skew, k, max_rot, min_norm);
      if (norm < min_norm)
        {
          min_norm = norm;
          mpz_set (best_k, k);
          mpz_poly_set (fopt, ft);
          mpz_poly_set (gopt, gt);
        }
    }
  double_poly_clear (eq);

#if 0
  /* solve d*(d-1)*(d-2)/6*f[d]*k^3 + (d-1)*(d-2)/2*f[d-1]*k^2 +
     (d-2)*f[d-2]*k + f[d-3] = 0 */
  double_poly_init (eq, 3);
  eq->coeff[3] = (double) (d * (d-1) * (d-2)) / 6.0 * mpz_get_d (ft->coeff[d]);
  eq->coeff[2] = (double) ((d-1) * (d-2)) / 2.0 * mpz_get_d (ft->coeff[d-1]);
  eq->coeff[1] = (double) (d-2) * mpz_get_d (ft->coeff[d-2]);
  eq->coeff[0] = mpz_get_d (ft->coeff[d-3]);
  double_poly_cleandeg (eq, 3);
  nroots = double_poly_compute_all_roots (roots, eq);
  for (int l = 0; l < nroots; l++)
    {
      /* mpz_set_d rounds towards zero, thus we add or subtract 0.5 to get
         rounding to nearest */
      mpz_set_d (k, roots[l] >= 0 ? roots[l] + 0.5 : roots[l] - 0.5);
      norm = best_norm2 (ft, gt, f_raw, g_raw, skew, k, max_rot, min_norm);
      if (norm < min_norm)
        {
          min_norm = norm;
          mpz_set (best_k, k);
          mpz_poly_set (fopt, ft);
          mpz_poly_set (gopt, gt);
        }
    }
  double_poly_clear (eq);
#endif

  mpz_set (k, best_k);
 clear_ft_gt:
  mpz_clear (best_k);
  mpz_poly_clear (ft);
  mpz_poly_clear (gt);

  return min_norm;
}

#define SOPT_INIT_SIZE_ALLOCATED_TRANSLATIONS 1024

/* Size optimize the polynomial pair (f_raw, g_raw) with rotations and
   translations.
   Assume that deg(g) = 1.
   Return the size-optimized pair (f_opt, g_opt) and the skew lognorm of f_opt.
   The sopt_effort parameter is used to increase the number of translations
   that are considered. sopt_effort = 0 means that we consider only rotations.
   Consider rotation up to x^max_rot*g(x), usually max_rot = d-2 or d-3.
   found by sopt_find_translation_deg* and the translation k = 0.
   TODO: LLL on gram matrix to be faster
   TODO: precompute skew^i for i in [0..d]
   TODO: return the n better poly not only the best one

   Note: f_opt==f_raw and g_opt==g_raw are used by callers, so should be supported here.

*/
static double
size_optimization_aux (mpz_poly_ptr f_opt, mpz_poly_ptr g_opt,
                       mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw,
                       const unsigned int sopt_effort, const int verbose,
                       const int max_rot)
{
  ASSERT_ALWAYS (f_raw->deg >= 2);
  ASSERT_ALWAYS (g_raw->deg == 1);
  const int d = f_raw->deg;
  double best_lognorm =
      L2_skew_lognorm ((mpz_poly_ptr) f_raw);
  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);
  best_lognorm += expected_rotation_gain (f_raw, g_raw);

  if (verbose > 1)
  {
    fprintf (stderr, "\n# sopt: start size-optimization with polynomials:\n"
                     "# sopt: f_raw = ");
    mpz_poly_fprintf_verbose (stderr, f_raw, verbose);
    fprintf (stderr, "# sopt: g_raw = ");
    mpz_poly_fprintf_verbose (stderr, g_raw, verbose);
    fprintf (stderr, "# sopt: with lognorm = %f\n", best_lognorm);
  }

  /****************************** init *******************************/
  mpz_t tmp, tmp2, skew, a, b;
  mpz_poly ft, gt, flll, fld, gld, fbest, gbest;
  list_mpz_t list_k;
  mat_Z m, U;

  mpz_init (tmp);
  mpz_init (tmp2);
  mpz_init (skew);

  /* 1/4 < delta = a/b <= 1: the closer delta is from 1, the better the
     reduction is. We take delta=1, since in fixed dimension the algorithm is
     still polynomial (see "The optimal LLL algorithm is still polynomial in
     fixed dimension" by Ali Akhavi, Theoretical Computer Science, 2003). */
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 1);

  mpz_poly_init (ft, d);
  mpz_poly_init (gt, 1);
  mpz_poly_init (flll, d);
  mpz_poly_init (fld, d);
  mpz_poly_init (gld, 1);
  mpz_poly_init (fbest, d);
  mpz_poly_init (gbest, 1);

  list_mpz_init (list_k, SOPT_INIT_SIZE_ALLOCATED_TRANSLATIONS);

  /*************** generate list of translations to try *****************/
  if (d == 6 || d == 5)
  {
    if (d == 6)
      sopt_find_translations_deg6 (list_k, f_raw, g_raw, verbose, sopt_effort);
    else if (d == 5)
      sopt_find_translations_deg5 (list_k, f_raw, g_raw, verbose, sopt_effort);
    if (verbose > 1)
      fprintf (stderr, "# sopt: %" PRIu64 " values for the translations were "
                       "computed and added in list_k\n", list_k->len);

    /* Add 0 to the list of translations to try */
    mpz_set_ui (tmp, 0);
    list_mpz_append (list_k, tmp);
    if (verbose > 1)
      fprintf (stderr, "# sopt: k = 0 was added list_k\n");

    /* Sort list_k by increasing order and remove duplicates */
    list_mpz_sort_and_remove_dup (list_k, verbose);
    if (verbose > 1)
      fprintf (stderr, "# sopt: It remains %" PRIu64 " values after sorting"
                       " and removing duplicates\n", list_k->len);
  }
  else
  {
    /* Start with list_k = {0} */
    mpz_set_ui (tmp, 0);
    list_mpz_append (list_k, tmp);
    if (verbose > 1)
      fprintf (stderr, "# sopt: Start with list_k = { 0 }\n");
  }

  /****************************** Main loop *******************************/
  /* For each value of k, call LLL on translated polynomial, and for each
     polynomial returned by LLL, call size_optimize_local_descent */

  /* Compute the skewness that is going to be used in LLL */
  sopt_get_skewness (skew, f_raw, g_raw);

  /* Init matrix for LLL. We consider rotation up to x^max_rot*g(x), so m has
     max_rot+2 row vectors (max_rot+1 for the rotation + 1 for the polynomial f)
     of d+1 coefficients each. */
  LLL_init (&m, max_rot + 2, d+1);
  /* The transformation matrix is a square matrix of size d (because m has d
     row vectors). */
  LLL_init (&U, d, d);

  if (verbose > 1)
    fprintf (stderr, "# sopt: Start processing all k in list_k of length "
                     "%" PRIu64 "\n", list_k->len);

  mpz_poly_set (fbest, f_raw);
  mpz_poly_set (gbest, g_raw);
  mpz_t ki;
  mpz_init (ki);
  list_mpz_t list_k_opt;
    {
      list_mpz_init (list_k_opt, list_k->len);
      for (unsigned int i = 0; i < list_k->len; i++)
        {
      mpz_set (ki, list_k->tab[i]);

      best_norm2 (ft, gt, f_raw, g_raw, skew, ki, max_rot, DBL_MAX);

      /* check if this translation 'ki' was already used */
      int new = 1;
      for (unsigned int k = 0; k < list_k_opt->len; k++)
        if (mpz_cmp (ki, list_k_opt->tab[k]) == 0)
          {
            new = 0;
            break;
          }
      if (new == 0)
        continue;

      list_mpz_append (list_k_opt, ki);

      double lognorm = sopt_local_descent (fld, gld, ft, gt, 1, d-2,
                                           SOPT_DEFAULT_MAX_STEPS, verbose);
      lognorm += expected_rotation_gain (fld, gld);
#ifdef POLYSELECT_CLASSICAL
      /* we only keep polynomial pairs with Res(f,g) = +/-N */
      if (lognorm < best_lognorm &&
          mpz_cmpabs (fld->coeff[d], f_raw->coeff[d]) == 0)
#else
      if (lognorm < best_lognorm)
#endif
        {
          if (verbose > 1)
            {
              gmp_fprintf (stderr, "# sopt:       better lognorm %.2f (previ"
                           "ous was %.2f) for skew = %Zd and k = %Zd\n",
                           lognorm, best_lognorm, skew, ki);
            }
          best_lognorm = lognorm;
          mpz_poly_swap (fbest, fld);
          mpz_poly_swap (gbest, gld);
        }
        }
      list_mpz_clear (list_k_opt);
  }
  mpz_clear (ki);

  mpz_poly_set (f_opt, fbest);
  mpz_poly_set (g_opt, gbest);

  /************************** Clear everything ***************************/
  gmp_randclear(rstate);
  LLL_clear (&m);
  LLL_clear (&U);
  mpz_clear (skew);
  list_mpz_clear (list_k);
  mpz_poly_clear (flll);
  mpz_poly_clear (fld);
  mpz_poly_clear (gld);
  mpz_poly_clear (fbest);
  mpz_poly_clear (gbest);
  mpz_poly_clear (ft);
  mpz_poly_clear (gt);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (tmp);
  mpz_clear (tmp2);

  if (verbose > 1)
  {
    fprintf (stderr, "# sopt: end of size-optimization, best polynomials are\n");
    fprintf (stderr, "# sopt: f_opt = ");
    mpz_poly_fprintf_verbose (stderr, f_opt, verbose);
    fprintf (stderr, "# sopt: g_opt = ");
    mpz_poly_fprintf_verbose (stderr, g_opt, verbose);
    fprintf (stderr, "# sopt: with lognorm = %f\n\n", best_lognorm);
  }

  return best_lognorm;
}

double
size_optimization (mpz_poly_ptr f_opt, mpz_poly_ptr g_opt,
                   mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw,
                   const unsigned int sopt_effort, const int verbose)
{
  int d = f_raw->deg;
  /* rotation up to degree d-2 */
  return size_optimization_aux (f_opt, g_opt, f_raw, g_raw,
                                sopt_effort, verbose, d - 2);
}
