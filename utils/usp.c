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

#include "cado.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "usp.h"
#include "double_poly.h"
#include "portability.h"
#include "macros.h"

/* #define DEBUG */

void
root_struct_init (root_struct *R) {
  mpz_init (R->a);
  mpz_init (R->b);
}

void
root_struct_clear (root_struct *R) {
  mpz_clear (R->a);
  mpz_clear (R->b);
}

static int
sign (mpz_t a)
{
  int i;

  i = mpz_cmp_ui (a, 0);
  if (i > 0)
    return 1;
  else if (i < 0)
    return -1;
  else
    return 0;
}

static double
ln2 (mpz_t a)
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
divide (mpz_t a, int k, int n, mpz_t *p)
{
  int i;
  mpz_t u, v, w;

  mpz_init (u);
  mpz_init (v);
  mpz_init (w);
  mpz_set (u, p[n]);
  for (i = n-1; i >= 0; i--)
    {
      mpz_tdiv_q_2exp (w, u, k); /* p[i] <- u/2^k */
      mpz_mul (v, a, w);
      mpz_add (u, p[i], v);
      mpz_set (p[i], w);
    }
  n--; /* reduces degree by 1 */
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (w);
}

/* isolating interval is [a/2^ka, b/2^kb] */
static void
printInt (mpz_t a, int ka, mpz_t b, int kb, int *nroots, root_struct *R,
          int verbose)
{
  if (verbose)
    gmp_printf ("isolated root in [%Zd/2^%d, %Zd/2^%d]\n", a, ka, b, kb);
  if (R != NULL)
    {
      mpz_set (R[*nroots].a, a);
      R[*nroots].ka = ka;
      mpz_set (R[*nroots].b, b);
      R[*nroots].kb = kb;
    }
  *nroots += 1;
}

/* returns the sign of p(a/2^k) */
static int
signValue (mpz_t a, int k, int n, mpz_t *p)
{
  int i, ret;
  mpz_t v, w;

  mpz_init (v);
  mpz_init (w);
  mpz_set (w, p[n]);
  for (i = n-1; i>=0; i--)
    {
      mpz_mul (w, w, a);
      mpz_mul_2exp (v, p[i], k * (n-i));
      mpz_add (w, w, v);
    }
  ret = sign (w);
  mpz_clear (v);
  mpz_clear (w);
  return ret;
}

#ifdef DEBUG
static void
printQ (mpz_t a, int k)
{
  mpz_t w;

  mpz_init (w);
  mpz_out_str (stdout, 10, a);
  if (k > 0)
    {
      /* printf ("/%d",1<<k); may overflow */
      mpz_set_ui (w, 1);
      mpz_mul_2exp (w, w, k);
      printf ("/");
      mpz_out_str (stdout, 10, w);
    }
  mpz_clear (w);
}

static void
printPol (mpz_t *p, int n)
{
  int i;

  for (i = 0; i <= n; i++)
    {
      if (i > 0 && mpz_cmp_ui (p[i], 0) >= 0)
        printf ("+");
      printQ (p[i], 0);
      printf ("*x^%d", i);
    }
  printf ("\n");
}
#endif

/* returns number of real roots (isolated) in a/2^m..b/2^m of polynomial
   p[0]+p[1]*x+...+p[n]*x^n
   r[0..n] is an auxiliary array.
   If R is not NULL, puts the roots in R.
*/
static int
usp (mpz_t a, mpz_t b, int m, int up, int va, int vb, int n, int *nroots,
     mpz_t *p, mpz_t *r, int verbose, root_struct *R)
{
   int lmi, i, k, c, s, smi, last, d;
   mpz_t mi, u, v, w;

#ifdef DEBUG
   gmp_printf ("looking at interval %Zd/2^%d..%Zd/2^%d\n", a, m, b, m);
   printf ("up=%d va=%d vb=%d\n", up, va, vb);
#endif
   if (va * vb == 2 * (up % 2) - 1)
     up--;
   if (up == 0)
     return 0;
   else if (up == 1)
     {
       printInt (a, m, b, m, nroots, R, verbose);
       return 1;
     }
   mpz_init (mi);
   mpz_add (mi, a, b);
   lmi = m;
   if (mpz_fdiv_ui (mi, 2) == 0)
     mpz_tdiv_q_2exp (mi, mi, 1);
   else
     lmi++;
   smi = signValue (mi, lmi, n, p);
   if (smi == 0)
     { /* rational root at mi */
       int le, ri, i, n0 = n;
       mpz_t *q;
       /* we cannot divide in-place, otherwise we will modify the input
          polynomial for the rest of the algorithm */
       q = malloc ((n + 1) * sizeof (mpz_t));
       for (i = 0; i <= n0; i++)
         mpz_init_set (q[i], p[i]);
       while (smi == 0)
         {
#ifdef DEBUG
           printf ("rational root at ");
           mpz_out_str (stdout, 10, mi);
           printf ("/2^%d\n", lmi);
#endif
           printInt (mi, lmi, mi, lmi, nroots, R, verbose);
           divide (mi, lmi, n, q);
           n --;
#ifdef DEBUG
           printf ("new input polynomial is ");
           printPol (q, n);
#endif
           smi = signValue (mi, lmi, n, q);
         }
       if (lmi > m)
         {
           mpz_mul_2exp (a, a, 1);
           mpz_mul_2exp (b, b, 1);
         }
       le = usp (a, mi, lmi, n, signValue (a, lmi, n, q),
                 signValue (mi, lmi, n, q), n, nroots, q, r, verbose, R);
       ri = usp (mi, b, lmi, n, signValue (mi, lmi, n, q),
                 signValue (b, lmi, n, q), n, nroots, q, r, verbose, R);
       mpz_clear (mi);
       for (i = 0; i <= n0; i++)
         mpz_clear (q[i]);
       free (q);
       return 1 + le + ri;
   }
   if (va * smi < 0)
     {
       if (up == 2)
         {
           printInt (a, m, mi, lmi, nroots, R, verbose);
           printInt (mi, lmi, b, m, nroots, R, verbose);
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
   mpz_set (r[n], p[n]);
   mpz_init (w);
   for (i = n-1; i >= 0; i--)
     {
       mpz_mul (r[i], b, r[i+1]); 
       mpz_mul_2exp (w, p[i], (n-i) * m);
       mpz_add (r[i], r[i], w);
     }
   mpz_init (v);
   mpz_sub (v, a, b);
   mpz_init (u);
   mpz_set (u, v);
   for (k = 1; k < n; k++)
     {
       for (i = n-1; i >= k; i--)
         {
           mpz_mul (w, b, r[i+1]);
           mpz_add (r[i], r[i], w);
         }
       mpz_mul (r[k], r[k], u);
       mpz_mul (u, u, v);
     }
   mpz_clear (v);
   mpz_clear (w);
   mpz_mul (r[n], r[n], u);
   mpz_clear (u);
   last = sign (r[0]);
   d = n-1;
   for (c = k = s = 0; k <= n && c < 2; k++)
     {
       /* invariant: all signs in r[0]..r[n-(d+1)] are equal */
       while (d > k && sign (r[n-d]) * last >= 0)
         d--;
       if (d < k)
         {
           /* d+1 <= k, thus all signs in r[0]..r[n-k] are equal,
              thus only one more sign change is possible */
           c += (sign (r[n-k]) * s < 0);
           k = n;
         }
       else
         {
           for (i = n-1; i >= k; i--)
             mpz_add (r[n-i], r[n-i], r[n-i-1]);
           i = mpz_cmp_ui (r[n-k], 0);
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
   if (c == 1)
     printInt (a, m, b, m, nroots, R, verbose);
   else if (c > 1)
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
   Assume p[n] is not zero.
   T (if not zero) is a bound on the absolute value of the real roots.
   If verbose is non-zero, print the isolating intervals for the roots.
   If Roots is not NULL, put the isolating intervals in Roots[0..nroots-1].
*/
int
numberOfRealRoots (mpz_t *p, const int orig_n, double T, int verbose, root_struct *Roots)
{
  int i, nroots, n = orig_n;
  mpz_t a, R, R1, *r;
  double C, pn, x;
  mpf_t aa;

  mpz_init (a);
  r = (mpz_t*) malloc ((n+1) * sizeof (mpz_t));
  for (i = 0; i <= n; i++)
    mpz_init (r[i]);
  ASSERT_ALWAYS(mpz_cmp_ui (p[n], 0) != 0);
  nroots = 0; /* initialize number of roots found */
  if (mpz_cmp_ui (p[0], 0) == 0) /* root at 0 */
    {
      mpz_set_ui (a, 0);
      printInt (a, 0, a, 0, &nroots, Roots, verbose);
      while (mpz_cmp_ui (p[0], 0) == 0)
        {
          divide (a, 0, n, p);
          n--;
        }
    }
  if (T != 0.0)
    T = log (T - 0.5) / log (2.0);
  else
    {
      pn = ln2 (p[n]); /* leading coefficient */
      C = (mpz_cmp_ui (p[n-1], 0) == 0) ? 0.0
        : ln2 (p[n-1]) - pn - log (1.0 * n) / log (2.0);
      T = 0.0;
      for (i = 1; i <= n; i++)
        {
          if (mpz_cmp_ui (p[n-i], 0))
            {
              /* for i=1, we get ln2 (p[n-1]) - pn, which is always
                 larger than C = ln2 (p[n-1]) - pn - log2 (n) */
              x = (ln2 (p[n-i]) - pn) / (double) i;
              if (x > T)
                T = x;
            }
        }
      T += 1.0;
      T = T + log (1 + exp ((C - T) / log (2.0))) / log (2.0);
    }
  i = 1 + (int) T;
  mpz_set_ui (a, 1);
  mpz_mul_2exp (a, a, i);
  mpz_init (R);
  mpz_set (R, a);

  if (verbose)
    {
      mpf_init2 (aa, 10);
      mpf_set_z (aa, a);
      printf ("initial interval is -");
      mpf_out_str (stdout, 10, 0, aa);
      printf ("..");
      mpf_out_str (stdout, 10, 0, aa);
      printf ("\n");
      mpf_clear (aa);
    }

  mpz_init (R1);
  mpz_neg (R1, R);
  i = usp (R1, R, 0, n, signValue (R1, 0, n, p), signValue (R, 0, n, p),
           n, &nroots, p, r, verbose, Roots);

  mpz_clear (a);
  mpz_clear (R);
  mpz_clear (R1);
  for (i = 0; i <= orig_n; i++)
    mpz_clear (r[i]);
  free (r);

  return nroots;
}

/* refine the root interval r[0] for the polynomial p of degree n,
   and return a double-precision approximation of the corresponding root */
double
rootRefine (root_struct *r, mpz_t *p, int n)
{
  volatile double a, b, c; /* we need 'volatile' here to force c to be cast
                              to double precision in c = (a + b) * 0.5 below,
                              otherwise if c is computed in extended precision,
                              the comparison c == a || c == b might fail
                              forever and we might enter an infinite loop */
  double sa, sb, sc;
  double_poly_t q;
  mpz_poly_t P;
  unsigned long count = 0;

  P->coeff = p;
  P->deg = n;
  double_poly_init (q, n);
  double_poly_set_mpz_poly (q, P);
  a = ldexp (mpz_get_d (r[0].a), -r[0].ka); /* a/2^ka */
  b = ldexp (mpz_get_d (r[0].b), -r[0].kb); /* b/2^kb */
  sa = double_poly_eval (q, a);
  sb = double_poly_eval (q, b);
  ASSERT_ALWAYS (sa * sb < 0.0);
#define MAX_LOOPS (1024 + 1074) /* difference between the maximal and minimal
                                   exponents of the double format */
  while (count++ < MAX_LOOPS)
    {
      c = (a + b) * 0.5;
      if (c == a || c == b)
        break;
      sc = double_poly_eval (q, c);
      if (sa * sc < 0.0)
        b = c;
      else
        a = c;
    }
  if (count >= MAX_LOOPS)
    {
      mpz_poly_fprintf (stdout, P);
      printf ("a=%.16e\n", ldexp (mpz_get_d (r[0].a), -r[0].ka));
      printf ("b=%.16e\n", ldexp (mpz_get_d (r[0].b), -r[0].kb));
      printf ("sa=%.16e sb=%.16e\n", sa, sb);
      exit (EXIT_FAILURE);
    }
  double_poly_clear (q);
  return c;
}

#ifdef MAIN
int
main (int argc, char *argv[])
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
  if (argc >= 2)
    T = atof (argv[1]);
  nroots = numberOfRealRoots (p, n, T, verbose);
  printf ("%d real root(s)\n", nroots);

  for (i = 0; i <= n; i++)
    mpz_clear (p[i]);
  free (p);

  return 0;
}
#endif
