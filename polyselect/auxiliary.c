/* Auxiliary routines for polynomial selection

Copyright 2008, 2009, 2010, 2013 Emmanuel Thome, Paul Zimmermann

This file is part of CADO-NFS.

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
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h> /* for DBL_MAX */
#include <math.h>
#include "gmp.h"
#include "macros.h" /* for ASSERT_ALWAYS */
#include "portability.h"
#include "utils.h"
#include "auxiliary.h"
#include "rootsieve.h"
#include "murphyE.h"

/* define OPTIMIZE_MP to perform computations in multiple-precision */
//#define OPTIMIZE_MP

//#define DEBUG_OPTIMIZE_AUX

/* for the rotation, we try (j*x+k) for |k| <= 2^MAX_k */
int MAX_k = 16;

/************************* norm and skewness *********************************/

/* Same as L2_lognorm, but takes 'double' instead of 'mpz_t' as coefficients.
   Returns 1/2*log(int(int(F(r*cos(t)*s,r*sin(t))^2*r/s^d, r=0..1), t=0..2*Pi))
   (circular method). Cf Remark 3.2 in Kleinjung paper, Math. of Comp., 2006.

   Maple code for degree 2:
   f:=x->a2*x^2+a1*x+a0: d:=degree(f(x),x):
   int(int((f(s*cos(t)/sin(t))*(r*sin(t))^d)^2*r/s^d, r=0..1), t=0..2*Pi);
 */
static double
L2_lognorm_d (double_poly_srcptr p, double s)
{
  double n;
  double *a = p->coeff;
  unsigned int d = p->deg;

  ASSERT_ALWAYS(1 <= d && d <= 7);

  if (d == 1)
  {
    double a1, a0;
    a1 = a[1] * s;
    a0 = a[0];
    /* use circular integral (Sage code):
       var('a1,a0,x,y,r,s,t')
       f = a1*x+a0
       F = expand(f(x=x/y)*y)
       F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
       v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
       (s*v).expand().collect(pi)
    */
    n = a0 * a0 + a1 * a1;
    n = n * 0.785398163397448310; /* Pi/4 */
    return 0.5 * log (n / s);
  }
  else if (d == 2)
  {
    double a2, a1, a0;
    a2 = a[2] * s * s;
    a1 = a[1] * s;
    a0 = a[0];
    /* use circular integral (Sage code):
       var('a2,a1,a0,x,y,r,s,t')
       f = a2*x^2+a1*x+a0
       F = expand(f(x=x/y)*y^2)
       F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
       v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
       (s^2*v).expand().collect(pi)
    */
    n = 3.0 * (a2 * a2 + a0 * a0) + 2.0 * a0 * a2 + a1 * a1;
    n = n * 0.130899693899574704; /* Pi/24 */
    return 0.5 * log(n / (s * s));
  }
  else if (d == 3)
    {
      double a3, a2, a1, a0;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('a3,a2,a1,a0,x,y,r,s,t')
         f = a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^3)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^3*v).expand().collect(pi)
      */
      n = 5.0 * (a3 * a3 + a0 * a0) + 2.0 * (a3 * a1 + a0 * a2)
        + a1 * a1 + a2 * a2;
      n = n * 0.049087385212340519352; /* Pi/64 */
      return 0.5 * log(n / (s * s * s));
    }
  else if (d == 4)
    {
      double a4, a3, a2, a1, a0;

      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('a4,a3,a2,a1,a0,x,r,s,t')
         f = a4*x^4+a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^4)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^4*v).expand().collect(pi)
      */
      n = 35.0 * (a4 * a4 + a0 * a0) + 10.0 * (a4 * a2 + a2 * a0)
        + 5.0 * (a3 * a3 + a1 * a1) + 6.0 * (a4 * a0 + a3 * a1)
        + 3.0 * a2 * a2;
      n = n * 0.0049087385212340519352; /* Pi/640 */
      return 0.5 * log(n / (s * s * s * s));
    }
  else if (d == 5)
    {
      double a5, a4, a3, a2, a1, a0, s2 = s * s, s3 = s2 * s, s4 = s3 * s,
        s5 = s4 * s;

      /*
        f := a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
        F := expand(y^5*subs(x=x/y,f));
        int(int(subs(x=x*s,F)^2/s^5, x=-1..1), y=-1..1);
       */
      a0 = a[0];
      a1 = a[1] * s;
      a2 = a[2] * s2;
      a3 = a[3] * s3;
      a4 = a[4] * s4;
      a5 = a[5] * s5;
      /* use circular integral (Sage code):
         var('a5,a4,a3,a2,a1,a0,x,r,s,t')
         f = a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
         F = expand(f(x=x/y)*y^5)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^5*v).expand().collect(pi)
      */
      n = 6.0 * (a3 * a1 + a1 * a5 + a4 * a2 + a0 * a4)
        + 14.0 * (a0 * a2 + a3 * a5) + 63.0 * (a0 * a0 + a5 * a5)
        + 7.0 * (a4 * a4 + a1 * a1) + 3.0 * (a3 * a3 + a2 * a2);
      n = n * 0.0020453077171808549730; /* Pi/1536 */
      return 0.5 * log(n / s5);
    }
  else if (d == 6)
    {
      double a6, a5, a4, a3, a2, a1, a0;

      a6 = a[6] * s * s * s * s * s * s;
      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=6; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^d*v).expand().collect(pi)
      */
      n = 231.0 * (a6 * a6 + a0 * a0) + 42.0 * (a6 * a4 + a2 * a0)
        + 21.0 * (a5 * a5 + a1 * a1) + 7.0 * (a4 * a4 + a2 * a2)
        + 14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
        + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 5.0 * a3 * a3;
      n = n * 0.00043828022511018320850; /* Pi/7168 */
      return 0.5 * log(n / (s * s * s * s * s * s));
    }
  else /* d == 7 */
    {
      double a7, a6, a5, a4, a3, a2, a1, a0;

      a7 = a[7] * s * s * s * s * s * s * s;
      a6 = a[6] * s * s * s * s * s * s;
      a5 = a[5] * s * s * s * s * s;
      a4 = a[4] * s * s * s * s;
      a3 = a[3] * s * s * s;
      a2 = a[2] * s * s;
      a1 = a[1] * s;
      a0 = a[0];
      /* use circular integral (Sage code):
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=7; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         (s^d*v).expand().collect(pi)
      */
      n = 429.0*(a0*a0+a7*a7) + 33.0*(a1*a1+a6*a6) + 66.0*(a0*a2+a5*a7)
        + 9*(a2*a2+a5*a5) + 18*(a1*a3+a0*a4+a4*a6+a3*a7) + 5*(a3*a3+a4*a4)
        + 10*(a2*a4+a1*a5+a3*a5+a0*a6+a2*a6+a1*a7);
      n = n * 0.000191747598485705154; /* Pi/16384 */
      return 0.5 * log(n / (s * s * s * s * s * s * s));
    }
}

/* Returns the logarithm of the L2-norm as defined by Kleinjung, i.e.,
   log(1/2 sqrt(int(int((F(sx,y)/s^(d/2))^2, x=-1..1), y=-1..1))).
   Since we only want to compare norms, we don't consider the log(1/2) term,
   and compute only 1/2 log(int(int(...))) [here the 1/2 factor is important,
   since it is added to the alpha root property term].

   Circular method: integrate over the unit circle.
*/
double
L2_lognorm (mpz_poly_ptr f, double s)
{
  double res;
  double_poly_t a;
  double a_coeffs[MAXDEGREE];

  a->coeff = a_coeffs;
  double_poly_set_mpz_poly (a, f);
  res = L2_lognorm_d (a, s);
  return res;
}

#ifdef OPTIMIZE_MP
/* The name is misleading. It returns the before-log part in
   1/2 log(int(int(...)))
*/
void
L2_lognorm_mp (mpz_poly_ptr f, mpz_t s, int method, mpz_t norm)
{
  mpz_t n, tmp, tmp1, tmpsum;
  unsigned long i;
  unsigned int d = f->deg;

  mpz_init_set_ui (n, 1);
  mpz_init_set_ui (tmp, 1);
  mpz_init_set_ui (tmp1, 1);
  mpz_init_set_ui (tmpsum, 1);

  if (d != 6)
    {
      fprintf (stderr, "not yet implemented for degree %u\n", d);
      exit (1);
    }
  else {
    mpz_t a[d+1];
    mpz_init_set_ui (a[0], 1);
    for (i=1; i<=d; i++) {
      mpz_init (a[i]);
      mpz_mul (a[i], a[i-1], s);
    }
    for (i=0; i<=d; i++)
      mpz_mul (a[i], a[i], f->coeff[i]);

    // n = 231.0 * (a6 * a6 + a0 * a0)
    mpz_mul (tmp, a[6], a[6]);
    mpz_mul (tmp1, a[0], a[0]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul_ui(n, tmpsum, 231);
    //   + 42.0 * (a6 * a4 + a2 * a0)
    mpz_mul (tmp, a[6], a[4]);
    mpz_mul (tmp1, a[2], a[0]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 42);
    //   + 21.0 * (a5 * a5 + a1 * a1)
    mpz_mul (tmp, a[5], a[5]);
    mpz_mul (tmp1, a[1], a[1]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 21);
    // + 7.0 * (a4 * a4 + a2 * a2)
    mpz_mul (tmp, a[4], a[4]);
    mpz_mul (tmp1, a[2], a[2]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_addmul_ui (n, tmpsum, 7);
    // +  14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
    mpz_mul (tmp, a[6], a[2]);
    mpz_mul (tmp1, a[5], a[3]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul (tmp, a[4], a[0]);
    mpz_mul (tmp1, a[3], a[1]);
    mpz_add (tmp, tmp, tmp1);
    mpz_add (tmpsum, tmp, tmpsum);
    mpz_addmul_ui (n, tmpsum, 14);
    // + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2)
    mpz_mul (tmp, a[6], a[0]);
    mpz_mul (tmp1, a[5], a[1]);
    mpz_add (tmpsum, tmp, tmp1);
    mpz_mul (tmp, a[4], a[2]);
    mpz_add (tmpsum, tmpsum, tmp);
    mpz_addmul_ui (n, tmpsum, 10);
    //  + 5.0 * a3 * a3;
    mpz_mul (tmp, a[3], a[3]);
    mpz_addmul_ui (n, tmp, 5);

    for (i=0; i<=d; i++)
      mpz_clear (a[i]);
  }

  mpz_set (norm, n);
  mpz_clear (tmp);
  mpz_clear (tmp1);
  mpz_clear (tmpsum);
  mpz_clear (n);
}
#endif

/* Use derivative test, with ellipse regions */
double
L2_skewness (mpz_poly_ptr f, int prec)
{
  double_poly_t ff, df;
  double s = 0.0, a = 0.0, b = 0.0, c, nc, *fd, *dfd,
    s1, s2, s3, s4, s5, s6, s7;
  unsigned int d = f->deg;

  ASSERT_ALWAYS(1 <= d && d <= 7);

  double_poly_init (ff, d);
  double_poly_init (df, d);

  /* convert once for all to double's to avoid expensive mpz_get_d() */
  double_poly_set_mpz_poly (ff, f);
  fd = ff->coeff;
  dfd = df->coeff;
  if (d == 7)
    {
      /* Sage code:
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=7; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (16384*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^8).expand().collect(s)
         3003*a7^2*s^14 + 165*a6^2*s^12 + 330*a5*a7*s^12 + 27*a5^2*s^10
         + 54*a4*a6*s^10 + 54*a3*a7*s^10 + 5*a4^2*s^8 + 10*a3*a5*s^8
         + 10*a2*a6*s^8 + 10*a1*a7*s^8 - 5*a3^2*s^6 - 10*a2*a4*s^6
         - 10*a1*a5*s^6 - 10*a0*a6*s^6 - 27*a2^2*s^4 - 54*a1*a3*s^4
         - 54*a0*a4*s^4 - 165*a1^2*s^2 - 330*a0*a2*s^2 - 3003*a0^2
      */
      dfd[7] = 3003.0 * fd[7] * fd[7];
      dfd[6] = 165.0 * (fd[6] * fd[6] + 2.0 * fd[5] * fd[7]);
      dfd[5] = 27.0 * (fd[5]*fd[5] + 2.0*fd[4]*fd[6] + 2.0*fd[3]*fd[7]);
      dfd[4] = 5.0*(fd[4]*fd[4]+2.0*fd[3]*fd[5]+2.0*fd[2]*fd[6]+2.0*fd[1]*fd[7]);
      dfd[3] = 5.0*(fd[3]*fd[3]+2.0*fd[2]*fd[4]+2.0*fd[1]*fd[5]+2.0*fd[0]*fd[6]);
      dfd[2] = 27.0 * (fd[2]*fd[2] + 2.0*fd[1]*fd[3] + 2.0*fd[0]*fd[4]);
      dfd[1] = 165.0 * (fd[1]*fd[1] + 2.0*fd[0]*fd[2]);
      dfd[0] = 3003 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[7] + dfd[6] + dfd[5] + dfd[4] - dfd[3] - dfd[2] - dfd[1]
        - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s3 * s3; /* s^12 */
          s7 = s6 * s1; /* s^14 */
          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s3 * s3; /* s^12 */
          s7 = s6 * s1; /* s^14 */
          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;
          s2 = s1 * s1;
          s3 = s2 * s1;
          s4 = s2 * s2;
          s5 = s4 * s1;
          s6 = s3 * s3;
          s7 = s6 * s1;

          nc = dfd[7] * s7 + dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[3] * s3 - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 6)
    {
      /* Sage code:
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=6; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (7168*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^7/14).expand().collect(s)
         99*a6^2 * s^12 +
         6*(2*a4*a6 + a5^2) * s^10 +
         (2*a2*a6 + 2*a3*a5 + a4^2) * s^8 -
         (2*a0*a4 + 2*a1*a3 + a2^2) * s^4 -
         6*(2*a0*a2 + a1^2) * s^2 -
         99*a0^2
      */
      dfd[6] = 99.0 * fd[6] * fd[6];
      dfd[5] = 6.0 * ( 2.0 * fd[4] * fd[6] + fd[5] * fd[5] );
      dfd[4] = 2.0 * ( fd[2] * fd[6] + fd[3] * fd[5] ) + fd[4] * fd[4];
      dfd[2] = 2.0 * ( fd[0] * fd[4] + fd[1] * fd[3] ) + fd[2] * fd[2];
      dfd[1] = 6.0 * ( 2.0 * fd[0] * fd[2] + fd[1] * fd[1] );
      dfd[0] = 99.0 * fd[0] * fd[0] ;
      s = 1.0;
      nc = dfd[6] + dfd[5] + dfd[4] - dfd[2] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s5 * s1; /* s^12 */
          nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          s6 = s5 * s1; /* s^12 */
          nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;
          s2 = s1 * s1;
          s4 = s2 * s2;
          s5 = s4 * s1;
          s6 = s5 * s1;

          nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 5)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=5; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (1536*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^6/3).expand().collect(s)
      */
      dfd[5] = 105.0 * fd[5] * fd[5];
      dfd[4] = 7.0 * (2.0 * fd[3] * fd[5] + fd[4] * fd[4]);
      dfd[3] = 2.0 * (fd[1] * fd[5] + fd[2] * fd[4]) + fd[3] * fd[3];
      dfd[2] = 2.0 * (fd[0] * fd[4] + fd[1] * fd[3]) + fd[2] * fd[2];
      dfd[1] = 7.0 * (2.0 * fd[0] * fd[2] + fd[1] * fd[1]);
      dfd[0] = 105.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[5] + dfd[4] + dfd[3] - dfd[2] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          s5 = s4 * s1; /* s^10 */
          nc = dfd[5] * s5 + dfd[4] * s4 + dfd[3] * s3
            - dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 4)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=4; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (640*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^5/10).expand().collect(s)
      */
      dfd[4] = 14.0 * fd[4] * fd[4];
      dfd[3] = 2.0 * fd[2] * fd[4] + fd[3] * fd[3];
      dfd[1] = 2.0 * fd[0] * fd[2] + fd[1] * fd[1];
      dfd[0] = 14.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[4] + dfd[3] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          s4 = s2 * s2; /* s^8 */
          nc = dfd[4] * s4 + dfd[3] * s3 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 3)
    {
      /* Sage code:
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=3; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         var('r,s,t')
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (64*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^4).expand().collect(s)
      */
      dfd[3] = 15.0 * fd[3] * fd[3];
      dfd[2] = 2.0 * fd[1] * fd[3] + fd[2] * fd[2];
      dfd[1] = 2.0 * fd[0] * fd[2] + fd[1] * fd[1];
      dfd[0] = 15.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[3] + dfd[2] - dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
        }
      while (nc < 0);

      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          s3 = s2 * s1; /* s^6 */
          nc = dfd[3] * s3 + dfd[2] * s2 - dfd[1] * s1 - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else if (d == 2)
    {
      /* Sage code:
         var('r,s,t,y')
         R.<x> = PolynomialRing(ZZ)
         S.<a> = InfinitePolynomialRing(R)
         d=2; f = SR(sum(a[i]*x^i for i in range(d+1)))
         F = expand(f(x=x/y)*y^d)
         F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
         v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
         v = (24*v/pi).expand()
         dv = v.diff(s)
         dv = (dv*s^3).expand().collect(s)
      */
      dfd[2] = 6.0 * fd[2] * fd[2];
      dfd[1] = 0.0;
      dfd[0] = 6.0 * fd[0] * fd[0];
      s = 1.0;
      nc = dfd[2] + dfd[1] - dfd[0];
      /* first isolate the minimum in an interval [s, 2s] by dichotomy */
      while (nc > 0)
        {
          s = 0.5 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
        }
      do
        {
          s = 2.0 * s;
          s1 = s * s;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
        }
      while (nc < 0);
      /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
      a = (s == 2.0) ? 1.0 : 0.5 * s;
      b = s;
      /* use dichotomy to refine the root */
      while (prec--)
        {
          c = (a + b) * 0.5;
          s1 = c * c;   /* s^2 */
          s2 = s1 * s1; /* s^4 */
          nc = dfd[2] * s2 + dfd[1] - dfd[0];
          if (nc > 0)
            b = c;
          else
            a = c;
        }
    }
  else /* d == 1 */
    a = b = (fd[0] / fd[1] >= 0) ? fd[0] / fd[1] : -fd[0] / fd[1];

  s = (a + b) * 0.5;

  double_poly_clear (ff);
  double_poly_clear (df);

  return s;
}


#ifdef OPTIMIZE_MP

/* Use derivative test, with ellipse regions */
void
L2_skewness_derivative_mp (mpz_poly_ptr f, int prec, int method, mpz_t skewness)
{
  mpz_t s, s1, s2, s3, s4, s5, s6, a, b, c, nc;
  mpz_init (s);
  mpz_init (s1);
  mpz_init (s2);
  mpz_init (s3);
  mpz_init (s4);
  mpz_init (s5);
  mpz_init (s6);
  mpz_init (a);
  mpz_init (b);
  mpz_init (c);
  mpz_init (nc);

  int i;

  if (d == 6) {
    mpz_t df[d+1];
    mpz_t tmp;
    mpz_init_set_ui (tmp, 1);
    for (i=0; i<=d; i++)
      mpz_init (df[i]);

    // dfd[6] = 99.0 * fd[6] * fd[6];
    mpz_mul (df[6], f[6], f[6]);
    mpz_mul_ui(df[6], df[6], 99);
    // dfd[5] = 6.0 * ( 2.0 * fd[4] * fd[6] + fd[5] * fd[5] );
    mpz_mul (df[5], f[5], f[5]);
    mpz_mul (tmp, f[4], f[6]);
    mpz_add (df[5], df[5], tmp);
    mpz_add (df[5], df[5], tmp);
    mpz_mul_ui (df[5], df[5], 6);
    // dfd[4] = 2.0 * ( fd[2] * fd[6] + fd[3] * fd[5] ) + fd[4] * fd[4];
    mpz_mul (df[4], f[4], f[4]);
    mpz_mul (tmp, f[3], f[5]);
    mpz_add (df[4], df[4], tmp);
    mpz_add (df[4], df[4], tmp);
    mpz_mul (tmp, f[2], f[6]);
    mpz_add (df[4], df[4], tmp);
    mpz_add (df[4], df[4], tmp);
    //  dfd[2] = 2.0 * ( fd[0] * fd[4] + fd[1] * fd[3] ) + fd[2] * fd[2];
    mpz_mul (df[2], f[2], f[2]);
    mpz_mul (tmp, f[1], f[3]);
    mpz_add (df[2], df[2], tmp);
    mpz_add (df[2], df[2], tmp);
    mpz_mul (tmp, f[0], f[4]);
    mpz_add (df[2], df[2], tmp);
    mpz_add (df[2], df[2], tmp);
    // dfd[1] = 6.0 * ( 2.0 * fd[0] * fd[2] + fd[1] * fd[1] );
    mpz_mul (df[1], f[1], f[1]);
    mpz_mul (tmp, f[0], f[2]);
    mpz_add (df[1], df[1], tmp);
    mpz_add (df[1], df[1], tmp);
    mpz_mul_ui (df[1], df[1], 6);
    // dfd[0] = 99.0 * fd[0] * fd[0] ;
    mpz_mul (df[0], f[0], f[0]);
    mpz_mul_ui(df[0], df[0], 99);
    /*
    gmp_fprintf (stderr, "df[6]: %Zd\n", df[6]);
    gmp_fprintf (stderr, "df[5]: %Zd\n", df[5]);
    gmp_fprintf (stderr, "df[4]: %Zd\n", df[4]);
    gmp_fprintf (stderr, "df[2]: %Zd\n", df[2]);
    gmp_fprintf (stderr, "df[1]: %Zd\n", df[1]);
    gmp_fprintf (stderr, "df[0]: %Zd\n", df[0]);
    */

    mpz_set_si (nc, -1);
    mpz_set_ui (s, 1);

    /* first isolate the minimum in an interval [s, 2s] by dichotomy */
    while ( mpz_cmp_ui(nc, 0) < 0 ) {

      mpz_add (s, s, s); /* s = 2.0 * s */
      mpz_mul (s1, s, s); /* s^2 */
      mpz_mul (s2, s1, s1); /* s^4 */
      mpz_mul (s4, s2, s2); /* s^8 */
      mpz_mul (s5, s4, s1); /* s^10 */
      mpz_mul (s6, s5, s1); /* s^12 */

      /* nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
         - dfd[2] * s2 - dfd[1] * s1 - dfd[0];         */
      mpz_mul (s6, s6, df[6]);
      mpz_mul (s5, s5, df[5]);
      mpz_mul (s4, s4, df[4]);
      mpz_mul (s2, s2, df[2]);
      mpz_mul (s1, s1, df[1]);
      mpz_add (nc, s6, s5);
      mpz_add (nc, nc, s4);
      mpz_sub (nc, nc, s2);
      mpz_sub (nc, nc, s1);
      mpz_sub (nc, nc, df[0]);

    }

    /* now dv(s/2) < 0 < dv(s) thus the minimum is in [s/2, s] */
    mpz_cdiv_q_2exp (a, s, 1);
    mpz_set (b, s);
    /* use dichotomy to refine the root */
    while (prec--)
    {
      mpz_add (tmp, a, b);
      mpz_cdiv_q_2exp (c, tmp, 1);
      mpz_mul (s1, c, c); //s1 = c * c;
      mpz_mul (s2, s1, s1); //s2 = s1 * s1;
      mpz_mul (s4, s2, s2); //s4 = s2 * s2;
      mpz_mul (s5, s4, s1); //s5 = s4 * s1;
      mpz_mul (s6, s5, s1); //s6 = s5 * s1;

      /* nc = dfd[6] * s6 + dfd[5] * s5 + dfd[4] * s4
         - dfd[2] * s2 - dfd[1] * s1 - dfd[0]; */
      mpz_mul (s6, s6, df[6]);
      mpz_mul (s5, s5, df[5]);
      mpz_mul (s4, s4, df[4]);
      mpz_mul (s2, s2, df[2]);
      mpz_mul (s1, s1, df[1]);
      mpz_add (nc, s6, s5);
      mpz_add (nc, nc, s4);
      mpz_sub (nc, nc, s2);
      mpz_sub (nc, nc, s1);
      mpz_sub (nc, nc, df[0]);

      if (mpz_cmp_ui (nc, 0) > 0)
        mpz_set (b, c);
      else
        mpz_set (a, c);
    }

    mpz_clear (tmp);
    for (i=0; i<=d; i++)
      mpz_clear (df[i]);

  } // end
  else  {
    fprintf (stderr, "L2_skewness_derivative_mp not yet implemented for degree %d\n", d);
    exit (1);
  }

  mpz_add (s, a, b);
  mpz_cdiv_q_2exp (skewness, s, 1);

  mpz_clear (s);
  mpz_clear (s1);
  mpz_clear (s2);
  mpz_clear (s3);
  mpz_clear (s4);
  mpz_clear (s5);
  mpz_clear (s6);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (c);
  mpz_clear (nc);
}
#endif

/************************** polynomial arithmetic ****************************/

/* g <- content(f) where deg(f)=d */
void
content_poly (mpz_t g, mpz_poly_ptr f)
{
  unsigned int i;
  unsigned int d = f->deg;

  ASSERT(d >= 1);

  mpz_gcd (g, f->coeff[0], f->coeff[1]);
  for (i = 2; i <= d; i++)
    mpz_gcd (g, g, f->coeff[i]);
}

/* v <- f(r), where f is of degree d */
void
eval_poly_ui (mpz_t v, mpz_t *f, unsigned int d, unsigned long r)
{
  int i;

  mpz_set (v, f[d]);
  for (i = d - 1; i >= 0; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_add (v, v, f[i]);
    }
}

/* v <- f'(r), where f is of degree d */
void
eval_poly_diff_ui (mpz_t v, mpz_t *f, unsigned int d, unsigned long r)
{
  int i;

  mpz_mul_ui (v, f[d], d);
  for (i = d - 1; i >= 1; i--)
    {
      mpz_mul_ui (v, v, r);
      mpz_addmul_ui (v, f[i], i); /* v <- v + i*f[i] */
    }
}

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
   average p-valuation of the polynomial f. Works recursively. */
double
special_val0 (mpz_poly_ptr f, unsigned long p)
{
  double v;
  mpz_t c,  *h;
  int v0;
  unsigned long *roots, r, r0;
  int i, d = f->deg, nroots;
  mpz_poly_t g, H;

  mpz_init (c);
  content_poly (c, f);
  for (v = 0.0; mpz_divisible_ui_p (c, p); v++, mpz_divexact_ui (c, c, p));
  v0 = (int) v;

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
    mpz_poly_copy (g, f);

  mpz_poly_init (H, d);
  H->deg = d;
  h = H->coeff;
  /* first compute h(x) = g(px) */
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

  nroots = mpz_poly_roots_ulong (roots, g, p);
  ASSERT (nroots <= d);
  for (r0 = 0, i = 0; i < nroots; i++)
    {
      r = roots[i];
      eval_poly_diff_ui (c, g->coeff, d, r);
      if (mpz_divisible_ui_p (c, p) == 0) /* g'(r) <> 0 mod p */
  v += 1.0 / (double) (p - 1);
      else /* hard case */
  {
    /* g(px+r) = h(x + r/p), thus we can go from h0(x)=g(px+r0)
       to h1(x)=g(px+r1) by computing h0(x + (r1-r0)/p).
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
    if (v0 != d) /* this should catch all the cases where f(x) has a
        root of multiplicity d, but also more cases.
        In those cases we avoid an infinite loop, but the
        result is probably wrong. */
      v += special_val0 (H, p) / (double) p;
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
special_valuation (mpz_poly_ptr f, unsigned long p, mpz_t disc)
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
  e = mpz_poly_roots_ulong (NULL, f, p);
  if (p_divides_lc) {
      /* Or the discriminant would have valuation 1 at least */
      ASSERT(mpz_divisible_ui_p(f->coeff[d - 1], p) == 0);
      e++;
  }
  return (pd * e) / (pd * pd - 1);
    } else if (pvaluation_disc == 1) {
      /* special case where p^2 does not divide disc */
  int e;
  e = mpz_poly_roots_ulong (NULL, f, p);
        if (p_divides_lc)
          e ++;
  /* something special here. */
  return (pd * e - 1) / (pd * pd - 1);
    } else {
  v = special_val0(f, p) * pd;
  if (p_divides_lc) {
      /* compute g(x) = f(1/(px))*(px)^d, i.e., g[i] = f[d-i]*p^i */
      /* IOW, the reciprocal polynomial evaluated at px */
      mpz_poly_t G;
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
      v += special_val0(G, p);
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
get_alpha (mpz_poly_ptr f, unsigned long B)
{
  double alpha, e;
  unsigned long p;
  mpz_t disc;
  unsigned int d = f->deg;

  /* for F linear, we have q_p = 1 for all p, thus
     alpha(F) = sum(prime p <= B, log(p)/(p^2-1)) ~ 0.569959993064325 */
  if (d == 1)
    return 0.569959993064325;

  mpz_init (disc);
  discriminant (disc, f->coeff, d);

  /* special_valuation returns the expected average exponent of p in F(a,b)
     for coprime a, b, i.e., e = q_p*p/(p^2-1), thus the contribution for p
     is (1/(p-1) - e) * log(p) */

  /* prime p=2 */
  e = special_valuation (f, 2, disc);
  alpha = (1.0 - e) * log (2.0);

  /* FIXME: generate all primes up to B and pass them to get_alpha */
  for (p = 3; p <= B; p += 2)
    if (ulong_isprime (p))
      {
        e = special_valuation (f, p, disc);
        alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
      }
  mpz_clear (disc);
  return alpha;
}

/* affine part of the special valution for polynomial f over p. */
double
special_valuation_affine (mpz_poly_ptr f, unsigned long p, mpz_t disc )
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
      e = mpz_poly_roots_ulong (NULL, f, p);

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
      v = special_val0(f, p) * pd;
      v /= pd + 1.0;
      return v;
   }
}


/*
  Find biased alpha_projective for a poly f. It uses
  some hacks here which need to be changed in future.
  Until now, since this will only be done several
  times, hence the speed is not critical.

  Note that, the returned alpha is the  -val * log(p)
  biased part in the alpha. Hence, we can just add
  this to our affine part.
*/
double
get_biased_alpha_projective (mpz_poly_ptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;
   int d = f->deg;

   mpz_init (disc);
   discriminant (disc, f->coeff, d);

   /* prime p=2 */
   e = special_valuation (f, 2, disc) - special_valuation_affine (f, 2, disc);

   /* 1/(p-1) is counted in the affine part */
   alpha =  (- e) * log (2.0);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation(f, p, disc) - special_valuation_affine (f, p, disc);
         alpha += (- e) * log ((double) p);
      }

   mpz_clear (disc);

   return alpha;
}

#if 0
/*
  Similar to above, but for affine part.
*/
double
get_biased_alpha_affine (mpz_poly_ptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;

   mpz_init (disc);
   discriminant (disc, f, d);

   /* prime p=2 */
   e = special_valuation_affine (f, 2, disc);
   alpha =  (1.0 - e) * log (2.0);

   //printf ("\np: %u, val: %f, alpha: %f\n", 2, e, alpha);

   /* FIXME: generate all primes up to B and pass them to get_alpha */
   for (p = 3; p <= B; p += 2)
      if (ulong_isprime (p)) {
         e = special_valuation_affine (f, p, disc);
         alpha += (1.0 / (double) (p - 1) - e) * log ((double) p);
         //printf ("\np: %u, val: %f, alpha: %f\n", p, e, alpha);

      }
   mpz_clear (disc);
   return alpha;
}


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
      fprintf (stderr, "Error, cannot allocate memory in average_valuation_affine_root.\n");
      exit (1);
   }

   for (i = 0; i <= d; i++)
      mpz_init_set (fv[i], f[i]);

   /* remove the p-valuations from fv */
   content_poly (c, f);
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


/**************************** rotation ***************************************/

/* D <- discriminant (f+k*g), which has degree d */
static void
discriminant_k (mpz_t *D, mpz_poly_ptr f, mpz_t m, mpz_t b)
{
  int i, j, k;
  uint32_t **M, pivot;
  int d = f->deg;

  ASSERT_ALWAYS(d <= 9);

  /* we first put in D[i] the value of disc(f + i*g) for 0 <= i <= d,
     thus if disc(f + k*g) = a[d]*k^d + ... + a[0], then
          D[0] = a[0]
          D[1] = a[0] + a[1] + ... + a[d]
          ...
          D[d] = a[0] + a[1]*d + ... + a[d]*d^d */

  discriminant (D[0], f->coeff, d); /* a[0] */
  for (i = 1; i <= d; i++)
    {
      /* add b*x - m */
      mpz_add (f->coeff[1], f->coeff[1], b);
      mpz_sub (f->coeff[0], f->coeff[0], m);
      discriminant (D[i], f->coeff, d);
    }

  /* initialize matrix coefficients */
  M = (uint32_t**) malloc ((d + 1) * sizeof(uint32_t*));
  for (i = 0; i <= d; i++)
    {
      M[i] = (uint32_t*) malloc ((d + 1) * sizeof(uint32_t));
      M[i][0] = 1;
      for (j = 1; j <= d; j++)
        M[i][j] = i * M[i][j-1];
    }

  for (j = 0; j < d; j++)
    {
      /* invariant: D[i] = M[i][0] * a[0] + ... + M[i][d] * a[d]
         with M[i][k] = 0 for k < j and k < i */
      for (i = j + 1; i <= d; i++)
        {
          /* eliminate M[i][j] */
          pivot = M[i][j] / M[j][j];
          mpz_submul_ui (D[i], D[j], pivot);
          for (k = j; k <= d; k++)
            M[i][k] -= pivot * M[j][k];
        }
    }

  /* now we have an upper triangular matrix */
  for (j = d; j > 0; j--)
    {
      for (k = j + 1; k <= d; k++)
        mpz_submul_ui (D[j], D[k], M[j][k]);
      ASSERT_ALWAYS(mpz_divisible_ui_p (D[j], M[j][j]));
      mpz_divexact_ui (D[j], D[j], M[j][j]);
    }

  /* restore the original f[] */
  mpz_submul_ui (f->coeff[1], b, d);
  mpz_addmul_ui (f->coeff[0], m, d);

  for (i = 0; i <= d; i++)
    free (M[i]);
  free (M);
}

/* replace f + k0 * x^t * (b*x - m) by f + k * x^t * (b*x - m), and return k */
long
rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t)
{
  mpz_addmul_si (f[t + 1], b, k - k0);
  mpz_submul_si (f[t], m, k - k0);
  return k;
}

/* add k*x^t*g to f */
static void
rotate_auxg_si (mpz_t *f, mpz_t *g, long k, unsigned int t)
{
  mpz_addmul_si (f[t + 1], g[1], k);
  mpz_addmul_si (f[t], g[0], k);
}

/* replace f by f + k * x^t * (b*x + g0) */
void
rotate_auxg_z (mpz_t *f, mpz_t b, mpz_t g0, mpz_t k, unsigned int t)
{
  mpz_addmul (f[t + 1], b, k);
  mpz_addmul (f[t], g0, k);
}

unsigned long
rotate_area (long K0, long K1, long J0, long J1)
{
  return (unsigned long) (J1 - J0 + 1) * (unsigned long) (K1 - K0 + 1);
}

/* Use Emmanuel Thome's idea: assuming alpha(f + k*g) admits a Gaussian
   distribution with mean 'm' and standard deviation 's', then the probability
   that one polynomial has a value of alpha >= A is
   1/2*(1 - erf((A-m)/s/sqrt(2))), thus the probability that K polynomials
   have all their alpha values >= A is [1/2*(1 - erf((A-m)/s/sqrt(2)))]^K.
   For 'm' and 's' fixed, and a given K, we define the value of A for which
   this probability is 1/2 to be the expected value of A for K polynomials.

   We assume lognorm(f + k*g) + E(alpha(f + k*g)) is first decreasing, then
   increasing, then the optimal K corresponds to the minimum of that function.
*/
void
rotate_bounds (mpz_poly_ptr f, mpz_t b, mpz_t m, long *K0, long *K1,
               long *J0, long *J1, int verbose)
{
  /* exp_alpha[i] corresponds to K=2^i polynomials for m=0, s=1
     f := x -> 1/2*(1 - erf(x/sqrt(2)));
     seq(fsolve(f(x)^(2^k) = 1/2, x=-log[2](1.0+k)), k=0..30);
   */
  double exp_alpha[] = {0.0 /* K=1 */, -0.5449521356 /* 2 */,
                        -0.9981488825 /* 4 */, -1.385198061 /* 8 */,
                        -1.723526050 /* 16 */, -2.025111894 /* 32 */,
                        -2.298313131 /* 64 */, -2.549067054 /* 128 */,
                        -2.781676726 /* 256 */, -2.999326227 /* 512 */,
                        -3.204420841 /* 1024 */, -3.398814100 /* 2048 */,
                        -3.583961388 /* 4096 */, -3.761025425 /* 8192 */,
                        -3.930949902 /* 16384 */, -4.094511733 /* 32768 */,
                        -4.252358774 /* 65536 */, -4.405037486 /* 131072 */,
                        -4.553013560 /* 262144 */, -4.696687518 /* 524288 */,
                        -4.836406692 /* 2^20 */, -4.972474538 /* 2^21 */,
                        -5.105157963 /* 2^22 */, -5.234693169 /* 2^23 */,
                        -5.361290351 /* 2^24 */, -5.485137511 /* 2^25 */,
                        -5.606403590 /* 2^26 */, -5.725241052 /* 2^27 */,
                        -5.841788041 /* 2^27 */, -5.956170181 /* 2^29 */,
                        DBL_MAX};
  int i;
  long k0 = 0, j0 = 0;
  double lognorm, alpha, E0, E, best_E;
  double skewness = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
  long jmax = (long) ((double) (1L << MAX_k) / skewness);
  unsigned long max_area = 1UL << MAX_k;

#define MARGIN 0.12 /* we allow a small error margin in the expected lognorm
                       + alpha values, to get a larger search range */

  E0 = L2_lognorm (f, skewness);

  *K0 = -2;
  *J0 = -16;
  *J1 = 16;

  /* look for negative k: -2, -4, -8, ... */
  best_E = E0;
  for (i = 1; rotate_area (*K0, -*K0, *J0, *J1) < max_area; i++, *K0 *= 2)
    {
      k0 = rotate_aux (f->coeff, b, m, k0, *K0, 0);
      lognorm = L2_lognorm (f, L2_skewness (f, SKEWNESS_DEFAULT_PREC));
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
    }
  /* go back to k=0 */
  k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);
  *K1 = -*K0;

  /* now try negative j: -1, -3, -7, ... */
  for (i++; exp_alpha[i] != DBL_MAX && rotate_area (*K0, *K1, *J0, -*J0) < max_area; i++, *J0 = 2 * *J0 - 1)
    {
      j0 = rotate_aux (f->coeff, b, m, j0, *J0, 1);
      lognorm = L2_lognorm (f, L2_skewness (f, SKEWNESS_DEFAULT_PREC));
      alpha = exp_alpha[i];
      E = lognorm + alpha;
      if (E < best_E + MARGIN)
        {
          if (E < best_E)
            best_E = E;
        }
      else
        break;
      if (1 - 2 * *J0 > jmax)
        break;
    }
  *J1 = -*J0;

  if (verbose > 0)
    fprintf (stderr, "# Rotate bounds: %ld <= j <= %ld, %ld <= k <= %ld\n",
             *J0, *J1, *K0, *K1);

  /* rotate back to j=k=0 */
  rotate_aux (f->coeff, b, m, k0, 0, 0);
  rotate_aux (f->coeff, b, m, j0, 0, 1);
}

/* Return the smallest value of lognorm + alpha(f + (j*x+k)*(b*x-m)) for
   j and k small enough such that the norm does not increase too much, and
   modify f[] accordingly.
   The parameter "multi" means that several polynomials are wanted. If
   multi=0 or 1, then only 1 polynomial is returned (classical behavior).
   Otherwise, multi polynomials are stored in jmin and kmin (that
   must be initialized arrays with at least multi elements). This option
   might be useful for Coppersmith variant (a.k.a. MNFS).
   In the multi case, the smallest of the returned values of lognorm + alpha
   is returned (and f[] accordingly).
   Warning: the caller is responsible to update the skewness if needed.

   FIXME: it would be better to pass the polynomial g to rotate, instead of
   b and m.
   */
double
rotate (mpz_poly_ptr f, unsigned long alim, mpz_t m, mpz_t b,
        long *jmin, long *kmin, int multi, int verbose)
{
  mpz_array_t *D;
  long K0, K1, J0, J1, k0, k, i, j, j0, bestk;
  double *A, alpha, lognorm, best_alpha = DBL_MAX, best_lognorm = DBL_MAX;
  double corr = 0.0;
  double alpha0;
  unsigned long p;
  double *best_E = NULL; /* set to NULL to avoid warning... */
  double time_alpha = 0.0, time_norm = 0.0;
  int d = f->deg;
  unsigned long pp;
  double one_over_pm1, logp, average_alpha = 0.0;

  /* allocate best_E, to store the best (lognorm+alpha) in multi mode */
  if (multi > 1)
    {
      best_E = (double *) malloc (multi * sizeof (double));
      for (i = 0; i < multi; ++i)
        best_E[i] = DBL_MAX;
    }

  /* allocate D(k) = disc(f + (j*x+k)*g, x) */
  D = alloc_mpz_array (d + 1);

  /* compute range for k */
  rotate_bounds (f, b, m, &K0, &K1, &J0, &J1, verbose);
  ASSERT_ALWAYS(K0 <= 0 && 0 <= K1);

  /* allocate sieving zone for computing alpha */
  A = (double*) malloc ((K1 + 1 - K0) * sizeof (double));
  j0 = k0 = 0; /* the current coefficients f[] correspond to f+(j*x+k)*g */

  *jmin = *kmin = 0;

  alpha0 = get_alpha (f, alim); /* value of alpha without rotation */

  ASSERT_ALWAYS(J0 < 0 && 0 < J1);
  for (j = 0;;)
    {
      /* we consider j=0, 1, ..., J1, then J0, J0+1, ..., -1 */
      j0 = rotate_aux (f->coeff, b, m, j0, j, 1);
      /* go back to k=0 for the discriminant */
      k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);
      /* D(k) = disc(f + (j*x+k)*g, x) (j is now fixed) */
      discriminant_k (D->data, f, m, b);

      for (k = K0; k <= K1; k++)
  A[k - K0] = 0.0; /* A[k - K0] will store the value alpha(f + k*g) */

  for (p = 2; p <= alim; p += 1 + (p & 1))
    if (ulong_isprime (p))
      {
        int i;
        /* We skip primes which divide all coefficients of f, since then
           f mod p is zero. This can only happen when p divides N, which is
           silly in GNFS, but maybe the user is stupid. */
        for (i = d; i >= 0 && mpz_divisible_ui_p (f->coeff[i], p); i--);
        if (i < 0)
          continue;

  if (k0 != 0)
    k0 = rotate_aux (f->coeff, b, m, k0, 0, 0);

  time_alpha -= seconds ();

  one_over_pm1 = 1.0 / (double) (p - 1);
  logp = log ((double) p);
  for (pp = p; pp <= alim; pp *= p)
    {
      /* Murphy (page 48) defines cont_p(F) = q_p*p/(p^2-1)
         = q_p*p/(p+1)*(1/p+1/p^2+...)
         The contribution for p^k is thus q_p*p/(p+1)/p^k. */
      alpha = logp / (double) (p + 1) * (double) p / (double) pp;
      /* the following is the average contribution for a prime not
         dividing the discriminant, cf alpha.sage, function alpha_p.
         We take it into account only for p, not for p^2, p^3, ... */
      if (p == pp)
        average_alpha += logp * one_over_pm1;
      /* we do not consider the projective roots here, since the
         corresponding correction will be considered separately for each
         row below */
      /* + alpha_p_projective (f, d, (D->data)[0], p); */
      update_table (f->coeff, d, m, b, A, K0, K1, pp, alpha);
    }

  time_alpha += seconds ();
      } /* end of loop on primes p */

  /* determine the best alpha in each row */
  bestk = K0;
  for (k = K0 + 1; k <= K1; k++)
    if (A[k - K0] < A[bestk - K0])
      bestk = k;

  /* Correction to apply to the current row (takes into account the projective
     roots). FIXME: we are lazy here, we should only consider the contribution
     from the projective roots. */
  k0 = rotate_aux (f->coeff, b, m, k0, bestk, 0);
  corr = get_alpha (f, alim) - A[bestk - K0];

  if (verbose > 1)
    fprintf (stderr, "# best alpha for j=%ld: k=%ld with %f\n",
             j, bestk, A[bestk - K0] + corr);

  /* now finds the best lognorm+alpha */
  time_norm -= seconds ();
  for (k = K0; k <= K1; k++)
    {
      alpha = A[k - K0] + corr;
      if (alpha < best_alpha + 2.0)
        {
          /* check lognorm + alpha < best_lognorm + best_alpha */

          /* translate from k0 to k */
          k0 = rotate_aux (f->coeff, b, m, k0, k, 0);
          lognorm = L2_lognorm (f, L2_skewness (f, SKEWNESS_DEFAULT_PREC));
          if (multi <= 1) {
              if (lognorm + alpha < best_lognorm + best_alpha) {
                  best_lognorm = lognorm;
                  best_alpha = alpha;
                  *kmin = k;
                  *jmin = j;
              }
          } else { /* multi mode */
              /* Rem: best_lognorm + best_alpha is the worse of the
                 preselected */
              double newE = lognorm + alpha;
              if (newE < best_E[multi-1]) {
                  int ii;
                  /* find position; assume list of preselected is sorted */
                  for (ii = 0; ii < multi; ++ii) {
                      if (best_E[ii] > newE)
                          break;
                  }
                  ASSERT_ALWAYS(ii < multi);
                  /* insert */
                  for (i = multi - 1; i > ii; --i) {
                      kmin[i] = kmin[i-1];
                      jmin[i] = jmin[i-1];
                      best_E[i] = best_E[i-1];
                  }
                  kmin[ii] = k;
                  jmin[ii] = j;
                  best_E[ii] = newE;
              }
          }
        }
    }
  time_norm += seconds ();

  j++;
  if (j > J1)
    j = J0;
  else if (j == 0)
    break;

    } /* end of loop on j */

  /* we now have f + (j0*x+k0)*(bx-m) and we want f + (jmin*x+kmin)*(bx-m),
     thus we have to add ((jmin-j0)*x+(kmin-k0)*(bx-m) */
  /* if you are in multi-mode, we use the best polynomial */
  rotate_aux (f->coeff, b, m, k0, *kmin, 0);
  rotate_aux (f->coeff, b, m, j0, *jmin, 1);

  if ((verbose > 0) && (multi <= 1))
    {
      fprintf (stderr, "# Rotate by ");
      if (*jmin != 0)
        {
          if (*jmin == -1)
            fprintf (stderr, "-");
          else if (*jmin != 1)
            fprintf (stderr, "%ld*", *jmin);
          fprintf (stderr, "x");
          if (*kmin >= 0)
            fprintf (stderr, "+");
        }
      fprintf (stderr, "%ld: alpha improved from %1.2f to %1.2f (alpha %1.2fs, norm %1.2fs)\n", *kmin, alpha0, best_alpha, time_alpha, time_norm);
    }

  if (verbose && (multi > 1)) {
      fprintf(stderr, "Found the following polynomials  (j, k, E):\n");
      for (i = 0; i < multi; ++i) {
          fprintf(stderr, "  %ld\t%ld\t%1.2f\n", jmin[i], kmin[i], best_E[i]);
      }
  }

  free (A);

  clear_mpz_array (D);

  {
      double ret_val = best_lognorm + best_alpha;
      if (multi>1) {
          ret_val = best_E[0]; /* we return the smallest */
          free(best_E);
      }
      return ret_val;
  }
}


/*
   Print f, g only.
   Note: it's a backend for print_cadopoly().
*/
void
print_cadopoly_fg (FILE *fp, mpz_t *f, int df, mpz_t *g, int dg, mpz_t n )
{
   int i;

   /* n */
   gmp_fprintf (fp, "\nn: %Zd\n", n);

   /* Y[i] */
   for (i = dg; i >= 0; i--)
     gmp_fprintf (fp, "Y%d: %Zd\n", i, g[i]);

   /* c[i] */
   for (i = df; i >= 0; i--)
     gmp_fprintf (fp, "c%d: %Zd\n", i, f[i]);
}


/*
   Print f, g only, lognorm, skew, alpha, MurphyE.
   Note:  it's a backend for print_cadopoly_extra().
*/
double
print_cadopoly (FILE *fp, cado_poly p)
{
   unsigned int nroots = 0;
   double alpha, alpha_proj, logmu, e;
   mpz_poly_t F, G;

   F->coeff = p->alg->coeff;
   F->deg = p->alg->deg;
   G->coeff = p->rat->coeff;
   G->deg = p->rat->deg;

   /* print f, g only*/
   print_cadopoly_fg (fp, F->coeff, F->deg, G->coeff, G->deg, p->n);

#ifdef DEBUG
   fprintf (fp, "# ");
   fprint_polynomial (fp, F->coeff, F->deg);
   fprintf (fp, "# ");
   fprint_polynomial (fp, G->coeff, G->deg);
#endif

   fprintf (fp, "skew: %1.3f\n", p->skew);

   if (G->deg > 1)
   {
    logmu = L2_lognorm (G, p->skew);
    alpha = get_alpha (G, ALPHA_BOUND);
    alpha_proj = get_biased_alpha_projective (G, ALPHA_BOUND);
    nroots = numberOfRealRoots (G->coeff, G->deg, 0, 0, NULL);
    fprintf (fp, "# lognorm: %1.2f, alpha: %1.2f (proj: %1.2f), E: %1.2f, "
                 "nr: %u\n", logmu, alpha, alpha_proj, logmu + alpha, nroots);
   }

   logmu = L2_lognorm (F, p->skew);
   alpha = get_alpha (F, ALPHA_BOUND);
   alpha_proj = get_biased_alpha_projective (F, ALPHA_BOUND);
   nroots = numberOfRealRoots (F->coeff, F->deg, 0, 0, NULL);
   fprintf (fp, "# lognorm: %1.2f, alpha: %1.2f (proj: %1.2f), E: %1.2f, "
                "nr: %u\n", logmu, alpha, alpha_proj, logmu + alpha, nroots);

   e = MurphyE (p, bound_f, bound_g, area, MURPHY_K);
   fprintf (fp, "# MurphyE(Bf=%.2e,Bg=%.2e,area=%.2e)=%.2e\n",
        bound_f, bound_g, area, e);

   return e;
}


/*
   Print f, g, lognorm, skew, alpha, MurphyE, REV, time duration.
*/
void
print_cadopoly_extra (FILE *fp, cado_poly cpoly, int argc, char *argv[], double st)
{
   int i;

   print_cadopoly (fp, cpoly);
   /* extra info */
   fprintf (fp, "# generated by %s: %s", CADO_REV, argv[0]);
   for (i = 1; i < argc; i++)
      fprintf (fp, " %s", argv[i]);
   fprintf (fp, " in %.2fs\n", (seconds () - st));
}


/*
  Call print_cadopoly, given f, g and return MurphyE.
*/
double
print_poly_fg (mpz_poly_ptr f, mpz_t *g, mpz_t N, int mode)
{
   double e;
   int i;
   int d = f->deg;

   cado_poly cpoly;
   cado_poly_init(cpoly);
   for (i = 0; i < (d + 1); i++)
      mpz_set(cpoly->alg->coeff[i], f->coeff[i]);
   for (i = 0; i < 2; i++)
      mpz_set(cpoly->rat->coeff[i], g[i]);
   mpz_set(cpoly->n, N);
   cpoly->skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
   cpoly->alg->deg = d;
   cpoly->rat->deg = 1;

   if (mode == 1)
     {
       e = print_cadopoly (stdout, cpoly);
       fflush(stdout);
     }
   else
     e = MurphyE (cpoly, bound_f, bound_g, area, MURPHY_K);

   cado_poly_clear (cpoly);
   return e;
}

/* f <- f(x+k), g <- g(x+k) */
void
do_translate_z (mpz_poly_ptr f, mpz_t *g, mpz_t k)
{
  int i, j;
  int d = f->deg;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul (f->coeff[j], f->coeff[j+1], k);
  mpz_addmul (g[0], g[1], k);
}

/* Use rotation and translation to find a polynomial with smaller norm
   (local minimum). Modify f and g accordingly.
   If use_rotation is non zero, use also rotation.
*/
void
optimize_aux (mpz_poly_ptr f, mpz_t *g, int verbose, int use_rotation)
{
  mpz_t kt, k2, k1, k, l; /* current offset */
  mpz_t ktot, khitot, lamtot, mutot;
  /* compute total translation and rotation, such that current polynomials are
     [f+(khitot*x^2+lamtot*x+mutot)*g](x+k) and g(x+k) */
  mpz_t tmp;
  int changed, changedt, changed2, changed1;
  double logmu00, logmu0, logmu, skew;
  int prec = SKEWNESS_DEFAULT_PREC;
  int count = 0;
  int d = f->deg;
  mpz_poly_t G;

  G->coeff = g;
  G->deg = 1;
  skew = L2_skewness (f, prec);
  logmu00 = logmu0 = L2_lognorm (f, skew);
  mpz_init_set_ui (k, 1);
  mpz_init_set_ui (k2, 1);
  mpz_init_set_ui (k1, 1);
  mpz_init_set_ui (kt, 1);
  mpz_init (l);
  mpz_init (ktot);
  mpz_init (khitot);
  mpz_init (lamtot);
  mpz_init (mutot);
  mpz_init (tmp);
  while (1)
    {
      changed = changedt = changed2 = changed1 = 0;

      /* first try translation by kt */
      do_translate_z (f, g, kt); /* f(x+kt) */
      mpz_add (ktot, ktot, kt);
      skew = L2_skewness (f, prec);
      logmu = L2_lognorm (f, skew);
      if (logmu < logmu0)
        {
#ifdef DEBUG_OPTIMIZE_AUX
          gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t ktot=%Zd\n", logmu, logmu0, ktot);
#endif
          changedt = 1;
          logmu0 = logmu;
        }
      else
        {
          mpz_mul_si (l, kt, -2); /* l = -2*kt */
          do_translate_z (f, g, l); /* f(x-kt) */
          mpz_add (ktot, ktot, l);
          skew = L2_skewness (f, prec);
          logmu = L2_lognorm (f, skew);
          if (logmu < logmu0)
          {
#ifdef DEBUG_OPTIMIZE_AUX
            gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t ktot=%Zd\n", logmu, logmu0, ktot);
#endif
              changedt = 1;
              logmu0 = logmu;
            }
          else
            {
              do_translate_z (f, g, kt); /* original f */
              mpz_add (ktot, ktot, kt);
            }
        }

      /* then do rotation by k2*x^2*g if d >= 6 */
      if (d >= 6 && use_rotation)
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k2, 2);
          /* k2*(x-ktot)^2 = k2*x^2 - 2*k2*ktot*x + k2*ktot^2 */
          mpz_mul (tmp, ktot, k2);
          mpz_submul_ui (lamtot, tmp, 2);
          mpz_addmul (mutot, tmp, ktot);
          mpz_add (khitot, khitot, k2);
          skew = L2_skewness (f, prec);
          logmu = L2_lognorm (f, skew);
          if (logmu < logmu0)
            {
#ifdef DEBUG_OPTIMIZE_AUX
              gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t khitot=%Zd\n", logmu, logmu0, khitot);
#endif
              changed2 = 1;
              logmu0 = logmu;
            }
          else
            {
              mpz_mul_si (l, k2, -2); /* l = -2*k2 */
              rotate_auxg_z (f->coeff, g[1], g[0], l, 2); /* f - k2*x^2*g */
              mpz_mul (tmp, ktot, l);
              mpz_submul_ui (lamtot, tmp, 2);
              mpz_addmul (mutot, tmp, ktot);
              mpz_add (khitot, khitot, l);
              skew = L2_skewness (f, prec);
              logmu = L2_lognorm (f, skew);
              if (logmu < logmu0)
                {
#ifdef DEBUG_OPTIMIZE_AUX
                  gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t khitot=%Zd\n", logmu, logmu0, khitot);
#endif
                  changed2 = 1;
                  logmu0 = logmu;
                }
              else
                {
                  rotate_auxg_z (f->coeff, g[1], g[0], k2, 2); /* previous f */
                  mpz_mul (tmp, ktot, k2);
                  mpz_submul_ui (lamtot, tmp, 2);
                  mpz_addmul (mutot, tmp, ktot);
                  mpz_add (khitot, khitot, k2);
                }
            }
        }

      if (use_rotation)
        {
          /* then do rotation by k1*x*g */
          rotate_auxg_z (f->coeff, g[1], g[0], k1, 1);
          mpz_submul (mutot, ktot, k1);
          mpz_add (lamtot, lamtot, k1);
          skew = L2_skewness (f, prec);
          logmu = L2_lognorm (f, skew);
          if (logmu < logmu0)
            {
#ifdef DEBUG_OPTIMIZE_AUX
              gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t lamtot=%Zd\n", logmu, logmu0, lamtot);
#endif
              changed1 = 1;
              logmu0 = logmu;
            }
          else
            {
              mpz_mul_si (l, k1, -2); /* l = -2*k1 */
              rotate_auxg_z (f->coeff, g[1], g[0], l, 1); /* f - k1*x*g */
              mpz_submul (mutot, ktot, l);
              mpz_add (lamtot, lamtot, l);
              skew = L2_skewness (f, prec);
              logmu = L2_lognorm (f, skew);
              if (logmu < logmu0)
                {
#ifdef DEBUG_OPTIMIZE_AUX
                  gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t lamtot=%Zd\n", logmu, logmu0, lamtot);
#endif
                  changed1 = 1;
                  logmu0 = logmu;
                }
              else
                {
                  rotate_auxg_z (f->coeff, g[1], g[0], k1, 1); /* previous f */
                  mpz_submul (mutot, ktot, k1);
                  mpz_add (lamtot, lamtot, k1);
                }
            }

          /* then do rotation by k*g */
          rotate_auxg_z (f->coeff, g[1], g[0], k, 0);
          mpz_add (mutot, mutot, k);
          skew = L2_skewness (f, prec);
          logmu = L2_lognorm (f, skew);
          if (logmu < logmu0)
            {
#ifdef DEBUG_OPTIMIZE_AUX
              gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t mutot=%Zd\n", logmu, logmu0, mutot);
#endif
              changed = 1;
              logmu0 = logmu;
            }
          else
            {
              mpz_mul_si (l, k, -2); /* l = -2*k */
              rotate_auxg_z (f->coeff, g[1], g[0], l, 0); /* f - k*g */
              mpz_add (mutot, mutot, l);
              skew = L2_skewness (f, prec);
              logmu = L2_lognorm (f, skew);
              if (logmu < logmu0)
                {
#ifdef DEBUG_OPTIMIZE_AUX
                  gmp_fprintf (stderr, "%.10f (logmu) < %.10f (logmu0),\t mutot=%Zd\n", logmu, logmu0, mutot);
#endif
                  changed = 1;
                  logmu0 = logmu;
                }
              else
                {
                  rotate_auxg_z (f->coeff, g[1], g[0], k, 0); /* previous f */
                  mpz_add (mutot, mutot, k);
                }
            }
        } /* use_rotation */

      if (changedt == 0 && changed == 0 && changed2 == 0 && changed1 == 0 &&
          mpz_cmp_ui (kt, 1) == 0 && mpz_cmp_ui (k, 1) == 0 &&
          mpz_cmp_ui (k2, 1) == 0 && mpz_cmp_ui (k1, 1) == 0)
        break;

      if (changedt == 1)
        mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
      else if (mpz_cmp_ui (kt, 1) > 0)
        mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */

      if (changed == 1)
        mpz_mul_2exp (k, k, 1);       /* k <- 2*k */
      else if (mpz_cmp_ui (k, 1) > 0)
        mpz_div_2exp (k, k, 1);       /* k <- k/2 */

      if (changed2 == 1)
        mpz_mul_2exp (k2, k2, 1);       /* k2 <- 2*k2 */
      else if (mpz_cmp_ui (k2, 1) > 0)
        mpz_div_2exp (k2, k2, 1);       /* k2 <- k2/2 */

      if (changed1 == 1)
        mpz_mul_2exp (k1, k1, 1);       /* k1 <- 2*k1 */
      else if (mpz_cmp_ui (k1, 1) > 0)
        mpz_div_2exp (k1, k1, 1);       /* k1 <- k1/2 */

      if (count++ > 10000) /* avoid an infinite loop due to the random
                              choices when logmu=logmu0 */
        break;
    }

  if (verbose > 0)
    {
      gmp_fprintf (stderr, "# ad=%Zd: optimized lognorm from %.2f to %.2f\n",
                   f->coeff[d], logmu00, logmu0);
      gmp_fprintf (stderr, "# (rotation %Zd*x^2+%Zd*x+%Zd, translation %Zd)\n",
                   khitot, lamtot, mutot, ktot);
      if (verbose > 1)
        {
          fprintf (stderr, "# "); mpz_poly_fprintf (stderr, f);
          fprintf (stderr, "# "); mpz_poly_fprintf (stderr, G);
        }
    }

  mpz_clear (k2);
  mpz_clear (k1);
  mpz_clear (k);
  mpz_clear (kt);
  mpz_clear (l);
  mpz_clear (ktot);
  mpz_clear (khitot);
  mpz_clear (lamtot);
  mpz_clear (mutot);
  mpz_clear (tmp);
}

#ifdef OPTIMIZE_MP
/*
   use mpz
*/
void
optimize_aux_mp (mpz_poly_ptr f, mpz_t *g, int verbose, int use_rotation,
                 int method)
{
  mpz_t kt, k2, k1, k, l; /* current offset */
  mpz_t ktot, khitot, lamtot, mutot; /* compute total translation and rotation,
                                        such that current polynomial are
                                        [f+(khitot*x^2+lamtot*x+mutot)*g](x+k)
                                        and g(x+k) */
  mpz_t tmp, tmp0, logmu00, logmu0, logmu, skew, skew0;
  int changed, changedt, changed2, changed1;
  int prec = SKEWNESS_DEFAULT_PREC;
  int count = 0;

  mpz_init (tmp);
  mpz_init (tmp0);
  mpz_init (skew);
  mpz_init (skew0);
  mpz_init (logmu);
  mpz_init (logmu0);
  mpz_init (logmu00);
  mpz_init (l);
  mpz_init (ktot);
  mpz_init (lamtot);
  mpz_init (khitot);
  mpz_init (mutot);
  mpz_init_set_ui (k, 1);
  mpz_init_set_ui (k2, 1);
  mpz_init_set_ui (k1, 1);
  mpz_init_set_ui (kt, 1);

  // init norm
  L2_skewness_derivative_mp (f, prec, skew0);
  L2_lognorm_mp (f, skew0, logmu0);
  mpz_set (logmu00, logmu0);

  while (1)
  {
    changed = changedt = changed2 = changed1 = 0;

    /* first try translation by kt */
    do_translate_z (f, g, kt); /* f(x+kt) */
    mpz_add (ktot, ktot, kt);
    L2_skewness_derivative_mp (f, prec, skew);
    L2_lognorm_mp (f, skew, logmu);

    mpz_pow_ui (tmp, skew0, d); // s_old^6
    mpz_mul (tmp, tmp, logmu);
    mpz_pow_ui (tmp0, skew, d); // s^6
    mpz_mul (tmp0, tmp0, logmu0);
    if (mpz_cmp (tmp, tmp0) < 0)
    {
      changedt = 1;
      mpz_set (logmu0, logmu);
      mpz_set (skew0, skew);
    }
    else
    {
      mpz_mul_si (l, kt, -2); /* l = -2*kt */
      do_translate_z (f, g, l); /* f(x-kt) */
      mpz_add (ktot, ktot, l);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changedt = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        do_translate_z (f, g, kt); /* original f */
        mpz_add (ktot, ktot, kt);
      }
    }

    /* then do rotation by k2*x^2*g if d >= 6 */
    if (d >= 6 && use_rotation)
    {
      rotate_auxg_z (f->coeff, g[1], g[0], k2, 2);
      mpz_mul (tmp, ktot, k2);
      mpz_submul_ui (lamtot, tmp, 2);
      mpz_addmul (mutot, tmp, ktot);
      mpz_add (khitot, khitot, k2);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed2 = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k2, -2); /* l = -2*k2 */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 2); /* f - k2*x^2*g */
        mpz_mul (tmp, ktot, l);
        mpz_submul_ui (lamtot, tmp, 2);
        mpz_addmul (mutot, tmp, ktot);
        mpz_add (khitot, khitot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed2 = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k2, 2); /* previous f */
          mpz_mul (tmp, ktot, k2);
          mpz_submul_ui (lamtot, tmp, 2);
          mpz_addmul (mutot, tmp, ktot);
          mpz_add (khitot, khitot, k2);
        }
      }
    }

    if (use_rotation)
    {
      /* then do rotation by k1*x*g */
      rotate_auxg_z (f->coeff, g[1], g[0], k1, 1);
      mpz_submul (mutot, ktot, k1);
      mpz_add (lamtot, lamtot, k1);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed1 = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k1, -2); /* l = -2*k1 */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 1); /* f - k1*x*g */
        mpz_submul (mutot, ktot, l);
        mpz_add (lamtot, lamtot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed1 = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k1, 1); /* previous f */
          mpz_submul (mutot, ktot, k1);
          mpz_add (lamtot, lamtot, k1);
        }
      }

      /* then do rotation by k*g */
      rotate_auxg_z (f->coeff, g[1], g[0], k, 0);
      mpz_add (mutot, mutot, k);
      L2_skewness_derivative_mp (f, prec, skew);
      L2_lognorm_mp (f, skew, logmu);

      mpz_pow_ui (tmp, skew0, d); // s_old^6
      mpz_mul (tmp, tmp, logmu);
      mpz_pow_ui (tmp0, skew, d); // s^6
      mpz_mul (tmp0, tmp0, logmu0);
      if (mpz_cmp (tmp, tmp0) < 0)
      {
        changed = 1;
        mpz_set (logmu0, logmu);
        mpz_set (skew0, skew);
      }
      else
      {
        mpz_mul_si (l, k, -2); /* l = -2*k */
        rotate_auxg_z (f->coeff, g[1], g[0], l, 0); /* f - k*g */
        mpz_add (mutot, mutot, l);
        L2_skewness_derivative_mp (f, prec, skew);
        L2_lognorm_mp (f, skew, logmu);

        mpz_pow_ui (tmp, skew0, d); // s_old^6
        mpz_mul (tmp, tmp, logmu);
        mpz_pow_ui (tmp0, skew, d); // s^6
        mpz_mul (tmp0, tmp0, logmu0);
        if (mpz_cmp (tmp, tmp0) < 0)
        {
          changed = 1;
          mpz_set (logmu0, logmu);
          mpz_set (skew0, skew);
        }
        else
        {
          rotate_auxg_z (f->coeff, g[1], g[0], k, 0); /* previous f */
          mpz_add (mutot, mutot, k);
        }
      }
    } /* use_rotation */

    if (changedt == 1)
      mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
    else if (mpz_cmp_ui (kt, 1) > 0)
      mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */
    if (changed == 1)
      mpz_mul_2exp (k, k, 1);       /* k <- 2*k */
    else if (mpz_cmp_ui (k, 1) > 0)
      mpz_div_2exp (k, k, 1);       /* k <- k/2 */
    if (changed2 == 1)
      mpz_mul_2exp (k2, k2, 1);       /* k2 <- 2*k2 */
    else if (mpz_cmp_ui (k2, 1) > 0)
      mpz_div_2exp (k2, k2, 1);       /* k2 <- k2/2 */
    if (changed1 == 1)
      mpz_mul_2exp (k1, k1, 1);       /* k1 <- 2*k1 */
    else if (mpz_cmp_ui (k1, 1) > 0)
      mpz_div_2exp (k1, k1, 1);       /* k1 <- k1/2 */
    if (changedt == 0 && changed == 0 && changed2 == 0 && changed1 == 0 &&
        mpz_cmp_ui (kt, 1) == 0 && mpz_cmp_ui (k, 1) == 0 &&
        mpz_cmp_ui (k2, 1) == 0 && mpz_cmp_ui (k1, 1) == 0)
      break;
    if (count++ > 10000) /* avoid an infinite loop due to the random
                            choices when logmu=logmu0 */
      break;
  }

  if (verbose > 0)
  {
    gmp_fprintf (stderr, "# ad=%Zd: optimized lognorm from %Zd to %Zd\n",
                 f[d], logmu00, logmu0);
    gmp_fprintf (stderr, "# (rotation %Zd*x^2+%Zd*x+%Zd, translation %Zd)\n",
                 khitot, lamtot, mutot, ktot);
    if (verbose > 1)
    {
      fprintf (stderr, "# "); fprint_polynomial (stderr, f->coeff, d);
      fprintf (stderr, "# "); fprint_polynomial (stderr, g, d);
    }
  }
  mpz_clear (tmp);
  mpz_clear (tmp0);
  mpz_clear (logmu);
  mpz_clear (logmu0);
  mpz_clear (logmu00);
  mpz_clear (skew);
  mpz_clear (skew0);
  mpz_clear (k2);
  mpz_clear (k1);
  mpz_clear (k);
  mpz_clear (kt);
  mpz_clear (l);
  mpz_clear (ktot);
  mpz_clear (khitot);
  mpz_clear (lamtot);
  mpz_clear (mutot);
}
#endif

static void
eval_poly (mpz_t b, mpz_poly_ptr h, mpz_t k)
{
  unsigned int n = h->deg;

  mpz_mul (b, h->coeff[n], k);
  mpz_add (b, b, h->coeff[--n]);
  while (n > 0)
    {
      mpz_mul (b, b, k);
      mpz_add (b, b, h->coeff[--n]);
    }
}

/* puts in h[3], ..., h[0] the coefficients (in k) of the degree-3 (in x)
   coefficient of f(x+k) */
static void
fdminus3_translated (mpz_poly_ptr h, mpz_poly_ptr f)
{
  unsigned int d = f->deg;
  ASSERT_ALWAYS(h->deg == 3);
  mpz_mul_ui (h->coeff[3], f->coeff[d], (d * (d-1) * (d-2)) / 6);
  mpz_mul_ui (h->coeff[2], f->coeff[d-1], ((d-1) * (d-2)) / 2);
  mpz_mul_ui (h->coeff[1], f->coeff[d-2], d-2);
  mpz_set (h->coeff[0], f->coeff[d-3]);
}

/* assuming h(a)*h(b) < 0 where h(k) = h[n]*k^n+...+h[0] and a < b,
   refines the root such that a + 1 = b.
   The refined root is in a.
*/
static void
root_refine (mpz_t a, mpz_t b, mpz_poly_ptr h)
{
  mpz_t c, v;
  int sa, i;
  unsigned int n = h->deg;

  mpz_init (c);
  mpz_init (v);
  eval_poly (v, h, a);
  sa = mpz_sgn (v);
  while (1)
    {
      mpz_add (c, a, b);
      mpz_fdiv_q_2exp (c, c, 1);
      if (mpz_cmp (c, a) == 0)
        break;
      eval_poly (v, h, c);
      if (mpz_sgn (v) == sa)
        mpz_swap (a, c);
      else
        mpz_swap (b, c);
    }

  /* now b = a+1, we do one more iteration to round to nearest, i.e.,
     we check the sign of 2^n*h(a+1/2) = h[n]*(2a+1)^n + ... + 2^n*h[0] */

  mpz_mul_2exp (c, a, 1);
  mpz_add_ui (c, c, 1); /* c = 2a+1 */
  mpz_mul (v, h->coeff[n], c);
  mpz_addmul_ui (v, h->coeff[n-1], 2);
  for (i = n-2; i >= 0; i--)
    {
      mpz_mul (v, v, c);
      mpz_addmul_ui (v, h->coeff[i], 1UL << (n - i));
    }
  if (mpz_sgn (v) == sa) /* root is in [a+1/2,b], round to b */
    mpz_add_ui (a, a, 1);

  mpz_clear (c);
  mpz_clear (v);
}

/* put in r[0], r[1], r[2] integer approximations of the real roots of
   h[3]*x^3+...+h[0], and return the number of real roots */
static int
roots3 (mpz_t *r, mpz_poly_ptr H)
{
  int n = 0;
  mpz_t v, k;
  int s1, s2;
  mpz_t *h = H->coeff;

  ASSERT_ALWAYS(H->deg == 3);
  /* if h2^2-3*h1*h3 >= 0, the derivative of h has two real roots */
  mpz_mul (r[2], h[2], h[2]);
  mpz_mul (r[1], h[1], h[3]);
  mpz_mul_ui (r[1], r[1], 3);
  mpz_sub (r[2], r[2], r[1]);
  if (mpz_sgn (r[2]) >= 0) /* 1 or 3 real roots */
    {
      /* the roots of h' are (-h2 +/- sqrt(h2^2-3*h1*h3))/(3*h3) */
      mpz_sqrt (r[2], r[2]);
      mpz_neg (r[1], h[2]);
      mpz_sub (r[1], r[1], r[2]);
      mpz_add (r[2], r[2], h[2]);
      mpz_tdiv_q_ui (r[1], r[1], 3);
      mpz_tdiv_q (r[1], r[1], h[3]);
      mpz_tdiv_q_ui (r[2], r[2], 3);
      mpz_tdiv_q (r[2], r[2], h[3]);
      if (mpz_cmp (r[1], r[2]) > 0)
        mpz_swap (r[1], r[2]);
      /* now r[1] < r[2] are approximations of the two real roots of h' */
    }
  else /* only 1 real root */
    {
      mpz_set_ui (r[1], 0);
      mpz_set_ui (r[2], 0);
    }
  /* now we have four control points -Inf < r[1] < r[2] < +Inf */
  mpz_init (v);
  mpz_init (k);
  eval_poly (v, H, r[1]);
  s1 = mpz_sgn (v);
  if (-mpz_sgn (h[3]) * s1 < 0) /* one root in -Inf..r[1] */
    {
      mpz_set_si (k, -1);
      while (mpz_cmp (r[1], k) <= 0)
        mpz_mul_2exp (k, k, 1);
      while (1)
        {
          eval_poly (v, H, k);
          if (mpz_sgn (v) * s1 < 0)
            break;
          mpz_mul_2exp (k, k, 1);
        }
      root_refine (k, r[1], H);
      mpz_swap (r[n++], k);
    }
  eval_poly (v, H, r[2]);

  s2 = mpz_sgn (v);
  if (s1 * s2 < 0) /* one root in r[1]..r[2] */
    {
      root_refine (r[1], r[2], H);
      mpz_swap (r[n++], r[1]);
    }
  if (s2 * mpz_sgn (h[3]) < 0) /* one root in r[2]..+Inf */
    {
      mpz_set_ui (k, 1);
      while (mpz_cmp (k, r[2]) <= 0)
        mpz_mul_2exp (k, k, 1);
      while (1)
        {
          eval_poly (v, H, k);
          if (mpz_sgn (v) * s2 < 0)
            break;
          mpz_mul_2exp (k, k, 1);
        }
      root_refine (r[2], k, H);
      mpz_swap (r[n++], r[2]);
    }
  mpz_clear (v);
  mpz_clear (k);

  return n;
}

#if 0
/* put in r[0..n-1] the real roots of x^3 + a*x^2 + b*x + c and return n */
static int
cubic_roots (double *r, double a, double b, double c)
{
  double Q = (a * a - 3.0 * b) / 9.0;
  double R = ((2.0 * a * a - 9.0 * b) * a + 27.0 * c) / 54.0;
  double R2 = R * R;
  double Q3 = Q * Q * Q;
  double theta;

  a = a / 3.0;
  if (R2 < Q3)
    {
      double twopioverthree = 2.0943951023931954923;

      theta = acos (R / sqrt (Q3)) / 3.0;
      Q = -2.0 * sqrt (Q);
      r[0] = Q * cos (theta) - a;
      r[1] = Q * cos (theta + twopioverthree) - a;
      r[2] = Q * cos (theta - twopioverthree) - a;
      return 3;
    }
  else
    {
      double A = cbrt (fabs (R) + sqrt (R2 - Q3)), B;
      A = (R >= 0) ? -A : A;
      B = Q / A;
      r[0] = (A + B) - a;
      return 1;
    }
}

/* for degree 6, and a given l, estimate the smallest degree-3 coefficient
   of f + l*x^3*g */
static double
optimize_c3 (mpz_t *f, mpz_t *g, long l)
{
  double c6 = mpz_get_d (f[6]);
  double c5 = mpz_get_d (f[5]);
  double c4 = mpz_get_d (f[4]);
  double c3 = mpz_get_d (f[3]);
  double g1 = mpz_get_d (g[1]);
  double g0 = mpz_get_d (g[0]);
  double roots[3], k, v, vmin = DBL_MAX;
  int n;

  c6 = 20.0 * c6;
  c5 = 10.0 * c5;
  c4 = 4.0 * (c4 + (double) l * g1);
  c3 = c3 + (double) l * g0;
  n = cubic_roots (roots, c5 / c6, c4 / c6, c3 / c6);
  while (n)
    {
      k = round (roots[--n]);
      v = ((c6 * k + c5) * k + c4) * k + c3;
      v = fabs (v);
      if (v < vmin)
        vmin = v;
    }
  return vmin;
}
#endif

/* if use_rotation is non-zero, also use rotation */
void
optimize (mpz_poly_ptr f, mpz_t *g, int verbose, int use_rotation)
{
  int d = f->deg;

  /* Pre-optimize routine: we assume that f[d], f[d-1] and f[d-2] are small,
     and we want to force f[d-3] to be small. */
  if (d == 6)
  {
    mpz_t k, r[3], f_copy[MAXDEGREE], g0_copy, best_k;
    int i, j, n;
    long l, best_l = LONG_MAX;
    mpz_poly_t h;

    mpz_init (k);
    mpz_init (best_k);

#ifdef OPTIMIZE_MP
    mpz_t skew, logmu, best_logmu;
    mpz_init (skew);
    mpz_init_set_ui (logmu, 1);
    mpz_init_set_ui (best_logmu, 1);
    mpz_ui_pow_ui (best_logmu, 10, 100);
    double best_logmud = DBL_MAX;
#else
    double skew, logmu, best_logmu = DBL_MAX;
#endif

    mpz_poly_init (h, 3);
    h->deg = 3;
    for (i = 0; i < 3; i++)
      mpz_init (r[i]);

    /* f[d] and g[1] are not changed below */
    for (i = 0; i < d; i++)
      mpz_init_set (f_copy[i], f->coeff[i]);
    mpz_init_set (g0_copy, g[0]);

    /* We use here an idea suggested by Thorsten Kleinjung, namely
       rotating by l*x^3*g for several values of l, and keeping the
       best value of l.
       With RSA-768, P=10^5, admax=2500, incr=60, lq=3, nq=1000, seed=1, we
       get the following min/av/max lognorms (248 hits):
       LMAX=0: 68.97/71.66/73.12
       LMAX=1: 68.97/71.26/72.80
       LMAX=2: 68.97/71.04/72.47
       LMAX=4: 68.97/70.78/72.44
       LMAX=8: 68.59/70.46/72.33
       LMAX=16: 68.15/70.18/71.91
       LMAX=32: 68.05/69.90/71.28
       LMAX=64: 67.94/69.66/71.04
       LMAX=128: 67.94/69.44/70.86
       LMAX=256: 67.94/69.30/70.51
       LMAX=512: 67.94/69.23/70.36
       LMAX=1024: 67.94/69.20/70.10
    */
#define LMAX 256
    for (l = -LMAX; l <= LMAX; l++) /* we consider f + l*x^(d-3)*g */
#undef LMAX
    {
      for (i = 0; i < d; i++)
        mpz_set (f->coeff[i], f_copy[i]);
      mpz_set (g[0], g0_copy);

      rotate_auxg_si (f->coeff, g, l, 3);
      fdminus3_translated (h, f);
      n = roots3 (r, h);

      for (j = 0; j < n; j++)
      {
        for (i = 0; i < d; i++)
          mpz_set (f->coeff[i], f_copy[i]);
        mpz_set (g[0], g0_copy);
        rotate_auxg_si (f->coeff, g, l, 3);
        mpz_set (k, r[j]);

        do_translate_z (f, g, k);

        /* now reduce coefficients f[0], f[1], f[2] using rotation */
        mpz_ndiv_q (k, f->coeff[0], g[0]);
        mpz_neg (k, k);
        rotate_auxg_z (f->coeff, g[1], g[0], k, 0);
        mpz_ndiv_q (k, f->coeff[1], g[0]);
        mpz_neg (k, k);
        rotate_auxg_z (f->coeff, g[1], g[0], k, 1);
        mpz_ndiv_q (k, f->coeff[2], g[0]);
        mpz_neg (k, k);
        rotate_auxg_z (f->coeff, g[1], g[0], k, 2);

#ifdef OPTIMIZE_MP
        optimize_aux_mp (f, g, verbose, use_rotation, 0);
        L2_skewness_derivative_mp (f, SKEWNESS_DEFAULT_PREC, 0, skew);
        L2_lognorm_mp (f, skew, 0, logmu);

        double logmud = 0.0, skewd = 0.0, skewf = 0.0;
        logmud = mpz_get_d (logmu);
        skewf = mpz_get_d (skew);
        mpz_pow_ui (skew, skew, d); // s^6
        skewd = mpz_get_d (skew);
        logmud = logmud / skewd * 0.00043828022511018320850; /* Pi/7168 */
        logmud = 0.5 * log(logmud);

        if (logmud < best_logmud) {
          best_logmud = logmud;
          best_l = l;
          mpz_set (best_k, r[j]);
        }
#else
        optimize_aux (f, g, verbose, use_rotation);
        skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
        logmu = L2_lognorm (f, skew);

        if (logmu < best_logmu)
          {
            best_logmu = logmu;
            best_l = l;
            mpz_set (best_k, r[j]);
        }
#endif

      } // #roots
    }  // #translations

    ASSERT(best_l != LONG_MAX);

    /* now consider the best (l,j) */
    for (i = 0; i < d; i++)
      mpz_set (f->coeff[i], f_copy[i]);
    mpz_set (g[0], g0_copy);
    rotate_auxg_si (f->coeff, g, best_l, 3);
    mpz_set (k, best_k);
    do_translate_z (f, g, k);

    /* now reduce coefficients f[0], f[1], f[2] using rotation */
    mpz_ndiv_q (k, f->coeff[0], g[0]);
    mpz_neg (k, k);
    rotate_auxg_z (f->coeff, g[1], g[0], k, 0);
    mpz_ndiv_q (k, f->coeff[1], g[0]);
    mpz_neg (k, k);
    rotate_auxg_z (f->coeff, g[1], g[0], k, 1);
    mpz_ndiv_q (k, f->coeff[2], g[0]);
    mpz_neg (k, k);
    rotate_auxg_z (f->coeff, g[1], g[0], k, 2);

    mpz_clear (k);
    mpz_clear (best_k);
    mpz_poly_clear (h);
    for (i = 0; i < 3; i++)
      mpz_clear (r[i]);
    for (i = 0; i < d; i++)
      mpz_clear (f_copy[i]);
    mpz_clear (g0_copy);

#ifdef OPTIMIZE_MP
    mpz_clear (logmu);
    mpz_clear (skew);
    mpz_clear (best_logmu);
#endif

  }

#ifdef OPTIMIZE_MP
  optimize_aux_mp (f, g, verbose, use_rotation, 0);
#else
  optimize_aux (f, g, verbose, use_rotation);
#endif
}
