/* usp.c - isolation of real roots of a polynomial using Descartes' rule

Copyright (C) 1998, 2002, 2010 Paul Zimmermann

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

Changes:

- May 2010: modified to enable the use in a library too,
            changed to GNU coding style, got rid of global variables.
- November 2002: remove x factors so that ldegree = 0.

A) to use as a stand-alone program:

compile with -DMAIN

usp < file where file contains (one number per line)
(1) the polynomial degree n
(2) the integer coefficients, from degree 0 to n

Example: the following file represents x^4-5*x^3+3*x^2-x+1.
% usp << EOF
4
1
-1
3
-5
1
EOF
initial interval is -0.16e2..0.16e2
1: 0..4
2: 4..8
2 real root(s)

B) to use within another program: compile without -DMAIN, the main function
   is numberOfRealRoots (mpz_t *p, int n, double T, int verbose):
   - the input polynomial is p[0]+p[1]*x+...+p[n]*x^n
   - n is the degree, and p[n] should not be zero
   - T is either 0.0, or a bound on the absolute value of the real roots
   - if verbose is non-zero, the isolating intervals of the roots are
     printed on stdout
*/

/* define MAIN if you want to compile as a stand-alone program */
/* #define MAIN */

#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "usp.h"
#include "mpz_poly.h"
#include "double_poly.h"
#include "macros.h"

/*
 *
#define DEBUG
#define USP_VERBOSE
 */

void
usp_root_interval_init (usp_root_interval_ptr R) {
  mpz_init (R->a);
  mpz_init (R->b);
}

void
usp_root_interval_clear (usp_root_interval_ptr R) {
  mpz_clear (R->a);
  mpz_clear (R->b);
}

/* returns an approximation of log2|a| */
static double
mpz_ln2 (mpz_srcptr a)
{
   int l;
   double r;

   l = mpz_sizeinbase (a, 2);
   if (l <= 1024) /* a fits in a double */
     r = log (fabs (mpz_get_d (a))) / log (2.0);
   else
     {
       mpz_t b;

       mpz_init (b);
       mpz_tdiv_q_2exp (b, a, l - 900);
       r = mpz_get_d (b);
       r = (double) l - 900.0 + log (fabs (r)) / log (2.0);
       mpz_clear (b);
     }
   return r;
}

/* divides in-place the input polynomial by 2^k*x-a */
static void
divide (mpz_poly_ptr p, mpz_srcptr a, int k)
{
  int i;
  mpz_t u, v, w;

  mpz_init (u);
  mpz_init (v);
  mpz_init (w);
  mpz_set (u, mpz_poly_lc(p));
  for (i = p->deg-1; i >= 0; i--)
    {
      mpz_ptr pi = mpz_poly_coeff(p, i);
      mpz_tdiv_q_2exp (w, u, k); /* p[i] <- u/2^k */
      mpz_mul (v, a, w);
      mpz_add (u, pi, v);
      mpz_set (pi, w);
    }
  mpz_poly_cleandeg(p, p->deg - 1); /* reduces degree by 1 */
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (w);
}

void usp_root_interval_set (usp_root_interval_ptr R, mpz_srcptr a, int ka, mpz_srcptr b, int kb)
{
#ifdef USP_VERBOSE
    gmp_printf ("isolated root in [%Zd/2^%d, %Zd/2^%d]\n", a, ka, b, kb);
#endif
    mpz_set (R->a, a);
    mpz_set (R->b, b);
    R->ka = ka;
    R->kb = kb;
}

void usp_root_interval_set_ui (usp_root_interval_ptr R, unsigned long a, int ka, unsigned long b, int kb)
{
    mpz_set_ui(R->a, a);
    mpz_set_ui(R->b, b);
    R->ka = ka;
    R->kb = kb;
}

/* returns the sign of p(a/2^k) */
static int
signValue (mpz_poly_srcptr p, mpz_srcptr a, int k)
{
  int i, ret;
  mpz_t v, w;

  mpz_init (v);
  mpz_init (w);
  mpz_set (w, mpz_poly_lc(p));
  for (i = p->deg-1; i>=0; i--)
    {
      mpz_mul (w, w, a);
      mpz_mul_2exp (v, mpz_poly_coeff_const(p, i), k * (p->deg-i));
      mpz_add (w, w, v);
    }
  ret = mpz_sgn (w);
  mpz_clear (v);
  mpz_clear (w);
  return ret;
}

int usp_internal_transform(mpz_poly_ptr r, mpz_poly_srcptr p, mpz_srcptr a, mpz_srcptr b, int m)
{
    mpz_t u, v, w;
    mpz_init (w);
    mpz_init (v);
    mpz_init (u);

    /*
     * r[n]   = f[n]
     * r[n-1] = f[n]*b   + f[n-1] / 2^m
     * r[n-2] = f[n]*b^2 + f[n-1] / 2^m * b + f[n-2] / 2^(2m)
     * ...
     */
    int n = p->deg;
    mpz_poly_set_zero(r);
    mpz_poly_setcoeff(r, n, mpz_poly_lc(p));
    for (int i = n-1; i >= 0; i--) {
        mpz_ptr ri = mpz_poly_coeff(r, i);
        mpz_srcptr ri1 = mpz_poly_coeff_const(r, i+1);
        mpz_srcptr pi = mpz_poly_coeff_const(p, i);

        mpz_mul (ri, b, ri1); 
        mpz_mul_2exp (w, pi, (n-i) * m);
        mpz_add (ri, ri, w);
    }


    mpz_sub (v, a, b);
    mpz_set (u, v);
    for (int k = 1; k <= n; k++) {
        mpz_ptr rk = mpz_poly_coeff(r, k);
        for (int i = n-1; i >= k; i--) {
            mpz_ptr ri = mpz_poly_coeff(r, i);
            mpz_srcptr ri1 = mpz_poly_coeff_const(r, i+1);
            mpz_addmul (ri, b, ri1);
        }
        /* and then adjust by u = (a-b)^k */
        mpz_mul (rk, rk, u);
        mpz_mul (u, u, v);
    }

    mpz_clear (v);
    mpz_clear (w);
    mpz_clear (u);
    return mpz_sgn (mpz_poly_coeff_const(r, 0));
}

/* returns number of real roots (isolated) in a/2^m..b/2^m of polynomial
   p[0]+p[1]*x+...+p[n]*x^n
   r[0..n] is an auxiliary array.
   If R is not NULL, puts the roots in R.
*/
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static int usp (mpz_ptr a, mpz_ptr b, int m, int up, int va, int vb, int n, int *nroots,
     mpz_poly_srcptr p, mpz_poly_ptr r, int verbose, usp_root_interval *R)
{
   int lmi, i, k, c, s, smi, last, d;
   mpz_t mi;

#ifdef DEBUG
   printf ("enter usp with "); mpz_poly_fprintf(stdout, p);
   gmp_printf ("looking at interval %Zd/2^%d..%Zd/2^%d\n", a, m, b, m);
   printf ("up=%d va=%d vb=%d\n", up, va, vb);
#endif
   if (va * vb == 2 * (up % 2) - 1)
     up--;
   if (up == 0)
     return 0;
   else if (up == 1)
     {
       if (R) usp_root_interval_set(R[*nroots], a, m, b, m);
       (*nroots)++;
       return 1;
     }
   mpz_init (mi);
   mpz_add (mi, a, b);
   lmi = m;
   if (mpz_fdiv_ui (mi, 2) == 0)
     mpz_tdiv_q_2exp (mi, mi, 1);
   else
     lmi++;
   smi = signValue (p, mi, lmi);
   if (smi == 0)
     { /* rational root at mi */
       int le, ri;

       mpz_poly q;
       mpz_poly_init(q, -1);
       mpz_poly_set(q, p);
       /* we cannot divide in-place, otherwise we will modify the input
          polynomial for the rest of the algorithm */
       while (smi == 0)
         {
#ifdef DEBUG
           printf ("rational root at ");
           mpz_out_str (stdout, 10, mi);
           printf ("/2^%d\n", lmi);
#endif
           if (R) usp_root_interval_set(R[*nroots], mi, lmi, mi, lmi);
           (*nroots)++;

           divide (q, mi, lmi);
           n--;
#ifdef DEBUG
           printf ("new input polynomial is ");
           mpz_poly_fprintf(stdout, q);
#endif
           smi = signValue (q, mi, lmi);
         }
       if (lmi > m)
         {
           mpz_mul_2exp (a, a, 1);
           mpz_mul_2exp (b, b, 1);
         }
       le = usp (a, mi, lmi, n,
		 signValue (q, a, lmi),
                 signValue (q, mi, lmi),
		 n, nroots, q, r, verbose, R);
       ri = usp (mi, b, lmi, n,
		 signValue (q, mi, lmi),
                 signValue (q, b, lmi),
		 n, nroots, q, r, verbose, R);
       /* restore original a, b */
       if (lmi > m)
         {
           mpz_div_2exp (a, a, 1);
           mpz_div_2exp (b, b, 1);
         }
       mpz_clear (mi);
       mpz_poly_clear(q);
       return 1 + le + ri;
   }



   if (va * smi < 0)
     {
       if (up == 2)
         {
             if (R) usp_root_interval_set(R[*nroots], a, m, mi, lmi);
             (*nroots)++;
             if (R) usp_root_interval_set(R[*nroots], mi, lmi, b, m);
             (*nroots)++;
             mpz_clear (mi);
             return 2;
         }
       else if (vb * smi < 0)
         {
           mpz_t aa;
 
           mpz_init (aa);
           if (lmi > m)
             mpz_mul_2exp (aa, a, 1);
           else
             mpz_set (aa, a);
           c = usp (aa, mi, lmi, up, va, smi, n, nroots, p, r, verbose, R);
           if (c < up)
             {
               if (lmi > m)
                 mpz_mul_2exp (aa, b, 1);
               else
                 mpz_set (aa, b);
               c += usp (mi, aa, lmi, up - c, smi, vb, n, nroots, p, r, verbose, R);
             }
           mpz_clear (mi);
           mpz_clear (aa);
           return c;
         }
       }

   last = usp_internal_transform(r, p, a, b, m);


   d = n-1;
   for (c = k = s = 0; k <= n && c < 2; k++)
     {
       /* invariant: all signs in r[0]..r[n-(d+1)] are equal */
       while (d > k && mpz_sgn (mpz_poly_coeff_const(r, n-d)) * last >= 0)
         d--;
       if (d < k)
         {
           /* d+1 <= k, thus all signs in r[0]..r[n-k] are equal,
              thus only one more sign change is possible */
           c += (mpz_sgn (mpz_poly_coeff_const(r, n-k)) * s < 0);
           k = n;
         }
       else
         {
           for (i = n-1; i >= k; i--)
             mpz_add (mpz_poly_coeff(r, n-i),
                     mpz_poly_coeff(r, n-i),
                     mpz_poly_coeff_const(r, n-i-1));
           i = mpz_cmp_ui (mpz_poly_coeff_const(r, n-k), 0);
           if (s * i < 0)
             {
               c++;
               if (va * vb > 0)
                 c = 2;
             }
           if (i != 0)
             s = i; /* s is the last sign */
         }
       /* when k=n-1 here and c=1, necessarily va * vb < 0, otherwise
          we would have c>=2 already, thus when we exit we cannot have
          c = 2 and k=n+1 */
     }
   if (c == 1) {
       if (R) usp_root_interval_set(R[*nroots], a, m, b, m);
       (*nroots)++;
   } else if (c > 1)
     {
       mpz_t aa;

       mpz_init (aa);
       ASSERT(k <= n);
       if (lmi > m)
         mpz_mul_2exp (aa, a, 1);
       else
         mpz_set (aa, a);
       c = usp (aa, mi, lmi, up, va, smi, n, nroots, p, r, verbose, R);
       if (c < up)
         {
           if (lmi > m)
             mpz_mul_2exp (aa, b, 1);
           else
             mpz_set (aa, b);
           c += usp (mi, aa, lmi, up-c, smi, vb, n, nroots, p, r, verbose, R);
         }
       mpz_clear (aa);
     }
   mpz_clear (mi);
   return c;
}

/* return the number of real roots of the polynomial p[0]+p[1]*x+...+p[n]*x^n
   (where n is orig_n below).
   Assume p[n] is not zero.
   T (if not zero) is a bound on the absolute value of the real roots.
   If verbose is non-zero, print the isolating intervals for the roots.
   If Roots is not NULL, put the isolating intervals in Roots[0..nroots-1].
*/
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
int mpz_poly_number_of_real_roots_extra (mpz_poly_srcptr f, double T, usp_root_interval *Roots)
{
  const int verbose = 0;

#ifdef DEBUG
  printf ("enter numberOfRealRoots with ");
  mpz_poly_fprintf (stdout, f);
#endif

  ASSERT_ALWAYS(mpz_cmp_ui (mpz_poly_lc(f), 0) != 0);

  if (mpz_cmp_ui (mpz_poly_coeff_const(f, 0), 0) == 0) /* root at 0 */
    {
      mpz_poly q;
      mpz_t a;

      mpz_poly_init(q, -1);
      mpz_init (a);

      mpz_poly_set(q, f);
      mpz_set_ui (a, 0);

      if (Roots) usp_root_interval_set(Roots[0], a, 0, a, 0);
      int val = 0;
      while (mpz_cmp_ui (mpz_poly_coeff_const(f, val), 0) == 0) {
          // divide p by X, but we do not want to touch p !!
          divide (q, a, 0);
          val++;
      }

      int nroots = 1 + mpz_poly_number_of_real_roots_extra(q, T, Roots ? (Roots + 1) : NULL);
      
      mpz_clear(a);
      mpz_poly_clear(q);

      return nroots;
    }

  int nroots = 0; /* initialize number of roots found */

  if (T != 0.0)
    T = log (T - 0.5) / log (2.0);
  else {
      /* we use here the near optimal Fujiwara bound
       * 2*max(|p[n-1]/p[n]|, |p[n-2]/p[n]|^{1/2}, ..., |p[0]/(2p[n])|^{1/n})
       * */
      double pn = mpz_ln2 (mpz_poly_lc(f)); /* leading coefficient */
      T = 0.0;
      for (int i = 1; i <= f->deg; i++) {
          mpz_srcptr h = mpz_poly_coeff_const(f, f->deg - i);
          if (mpz_cmp_ui (h, 0) != 0) {
              double x = (mpz_ln2 (h) - pn - (i >= f->deg)) / (double) i;
              if (x > T)
                  T = x;
          }
      }
      T += 1.0; /* factor of 2 in Fujiwara bound */
  }

#ifdef DEBUG
  printf ("root bound is 2^%f\n", T);
#endif

  mpz_t a;
  mpz_t R, R1;
  mpz_poly r;

  mpz_init(a);
  mpz_init (R);
  mpz_init (R1);
  mpz_poly_init(r, -1);

  mpz_set_ui (a, 1);
  mpz_mul_2exp (a, a, 1 + (int) T);
  mpz_set (R, a);
  mpz_neg (R1, R);


  if (verbose)
    {
      mpf_t aa;
      mpf_init2 (aa, 10);
      mpf_set_z (aa, a);
      printf ("initial interval is -");
      mpf_out_str (stdout, 10, 0, aa);
      printf ("..");
      mpf_out_str (stdout, 10, 0, aa);
      printf ("\n");
      mpf_clear (aa);
    }

  usp (R1, R, 0, f->deg, signValue (f, R1, 0), signValue (f, R, 0),
           f->deg, &nroots, f, r, verbose, Roots);

  mpz_clear (a);
  mpz_clear (R);
  mpz_clear (R1);

  mpz_poly_clear(r);


  return nroots;
}

/* refine the root interval r[0] for the polynomial p of degree n,
   and return a double-precision approximation of the corresponding root,
   with a maximal error <= precision.
   Warning: precision is an absolute value, not a relative.
   Use it only if you know what you do!
*/
double
usp_root_interval_refine (usp_root_interval_ptr r, mpz_poly_srcptr P, double precision)
{
  double a, b, c;
  double sa, sb, sc;
  double_poly q;

  /* Note: if precision = 0.0, rootRefine will stop when the bound a and b
     are two adjacent floating-point numbers. */

  a = ldexp (mpz_get_d (r[0].a), -r[0].ka); /* a/2^ka */
  b = ldexp (mpz_get_d (r[0].b), -r[0].kb); /* b/2^kb */
  c = (a + b) * .5;

  ASSERT_ALWAYS (a <= b);

  if (b - a <= precision) /* includes the case a = b */
    return c;

  double_poly_init (q, P->deg);
  double_poly_set_mpz_poly (q, P);
  sa = double_poly_eval (q, a);
  sb = double_poly_eval (q, b);
  /* due to truncation of the initial coefficients, and rounding error in
     evaluation of q, it might be that sa and sb do not have opposite signs */
  if (sa * sb >= 0)
    {
      sa = double_poly_eval_safe (q, a);
      sb = double_poly_eval_safe (q, b);
    }
  if (sa == 0.0)
    {
      c = a;
      goto end_refine;
    }
  if (sb == 0.0)
    {
      c = b;
      goto end_refine;
    }
  ASSERT_ALWAYS(sa * sb < 0);
  while (b - a > precision) {
    /* Warning: with precision == 0. + x86 32 bits + gcc mathematical default
       -mfpmath=387, the computation of c = (a + b) * .5 is done on the
       top of the i387 stack in 80 bits precision (extended precision).
       The next comparison c == a is done in this case between a memory
       double (64 bits) and this 80 bits value, and might fail forever,
       unless we convert c to binary64. */
#if defined(__i386)
    { volatile double ms = (a + b) * 0.5; c = ms; }
#else
    c = (a + b) * .5;
#endif
    if (c == a || c == b) break; /* avoids infinite loops if precision = 0
                                    or precision < ulp(a) */

    /* Note: in principle we should also use double_poly_eval_safe here,
     because due to rounding errors double_poly_eval() might return a value
     with the wrong sign, and thus we might search for a root in the wrong
     half-interval. However this should happen rarely, thus for efficiency
     reasons we keep double_poly_eval() here (rootRefine is critical in the
     norm initialization in las and in the skewness computation for degree 6).
     Another solution would be to translate the input polynomial at x=a before
     the loop, which should reduce the cancellations when evaluating p(c). */
    sc = double_poly_eval (q, c);
    if (sa * sc < 0.) b = c; else a = c;
  }
 end_refine:
  double_poly_clear (q);
  return c;
}
#undef MAX_LOOPS

#ifdef MAIN
int main(int argc, char const * argv[])
{
  int i, s, n, nroots, verbose = 1;
  double T = 0.0;
  char c;
  mpz_t *p;

  scanf ("%d\n", &n); /* degree of polynomial */
  p = (mpz_t*) malloc ((n+1) * sizeof (mpz_t));
  for (i = 0; i <= n; i++)
    {
      mpz_init (p[i]);
      do
        c = getchar ();
      while (!(isdigit (c) || c=='+' || c=='-'));
      if (c=='+')
        s = 1;
      else if (c=='-')
        s = -1;
    else if (isdigit (c))
      {
        s = 1;
        ungetc (c, stdin);
      }
      mpz_inp_str (p[i], stdin, 0);
      if (s < 0)
        mpz_neg (p[i], p[i]);
    }
#ifdef DEBUG
  printf ("input polynomial is ");
  printPol (p, n);
#endif
  if (argc >= 4){
      /* ./a.out a ka b kb precision to test rootRefine */
      int a = atoi(argv[1]), ka = atoi(argv[2]);
      int b = atoi(argv[3]), kb = atoi(argv[4]);
      double precision = atof(argv[5]);
      usp_root_interval r;
      usp_root_interval_init(r);
      mpz_init_set_ui(r->a, a);
      r->ka = ka;
      mpz_init_set_ui(r->b, b);
      r->kb = kb;
      printf("rf=%lf\n", usp_root_interval_refine(r, p, precision));
      usp_root_interval_clear(r);
  }
  else{
      if (argc >= 2){
	  T = atof (argv[1]);
	  nroots = mpz_poly_number_of_real_roots (p, T, NULL);
	  printf ("%d real root(s)\n", nroots);
      }
  }

  for (i = 0; i <= n; i++)
    mpz_clear (p[i]);
  free (p);

  return 0;
}
#endif
