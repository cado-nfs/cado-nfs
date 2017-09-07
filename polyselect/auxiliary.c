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
#include "murphyE.h"

/* define OPTIMIZE_MP to perform computations in multiple-precision */
//#define OPTIMIZE_MP

//#define DEBUG_OPTIMIZE_AUX

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
    if (isnan (n) || isinf (n))
      n = DBL_MAX;
    return 0.5 * log (n / s);
  }
  else if (d == 2)
  {
    double a2, a1, a0;
    a2 = a[2] * s;
    a1 = a[1];
    a0 = a[0] / s;
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
    if (isnan (n) || isinf (n))
      n = DBL_MAX;
    return 0.5 * log(n);
  }
  else if (d == 3)
    {
      double a3, a2, a1, a0, invs = 1.0 / s;
      a3 = a[3] * s * s;
      a2 = a[2] * s;
      a1 = a[1];
      a0 = a[0] * invs;
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
      if (isnan (n) || isinf (n))
        n = DBL_MAX;
      return 0.5 * log(n * invs);
    }
  else if (d == 4)
    {
      double a4, a3, a2, a1, a0, invs = 1.0 / s;

      a4 = a[4] * s * s;
      a3 = a[3] * s;
      a2 = a[2];
      a1 = a[1] * invs;
      a0 = a[0] * invs * invs;
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
      if (isnan (n) || isinf (n))
        n = DBL_MAX;
      return 0.5 * log(n);
    }
  else if (d == 5)
    {
      double a5, a4, a3, a2, a1, a0, invs = 1.0 / s;

      /*
        f := a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0:
        F := expand(y^5*subs(x=x/y,f));
        int(int(subs(x=x*s,F)^2/s^5, x=-1..1), y=-1..1);
       */
      a0 = a[0] * invs * invs;
      a1 = a[1] * invs;;
      a2 = a[2];
      a3 = a[3] * s;
      a4 = a[4] * s * s;
      a5 = a[5] * s * s * s;
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
      if (isnan (n) || isinf (n))
        n = DBL_MAX;
      return 0.5 * log(n * invs);
    }
  else if (d == 6)
    {
      double a6, a5, a4, a3, a2, a1, a0, invs;

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
      invs = 1.0 / s;
      a6 = a[6] * s * s * s;
      a5 = a[5] * s * s;
      a4 = a[4] * s;
      a3 = a[3];
      a2 = a[2] * invs;
      a1 = a[1] * invs * invs;
      a0 = a[0] * invs * invs * invs;
      n = 231.0 * (a6 * a6 + a0 * a0) + 42.0 * (a6 * a4 + a2 * a0)
        + 21.0 * (a5 * a5 + a1 * a1) + 7.0 * (a4 * a4 + a2 * a2)
        + 14.0 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
        + 10.0 * (a6 * a0 + a5 * a1 + a4 * a2) + 5.0 * a3 * a3;
      n = n * 0.00043828022511018320850; /* Pi/7168 */
      if (isnan (n) || isinf (n))
        n = DBL_MAX;
      return 0.5 * log(n);
    }
  else /* d == 7 */
    {
      double a7, a6, a5, a4, a3, a2, a1, a0;
      double invs = 1.0 / s;

      a7 = a[7] * s * s * s * s;
      a6 = a[6] * s * s * s;
      a5 = a[5] * s * s;
      a4 = a[4] * s;
      a3 = a[3];
      a2 = a[2] * invs;
      a1 = a[1] * invs * invs;
      a0 = a[0] * invs * invs * invs;
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
      if (isnan (n) || isinf (n))
        n = DBL_MAX;
      return 0.5 * log(n * invs);
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
L2_lognorm (mpz_poly_srcptr f, double s)
{
  double res;
  double_poly a;
  double_poly_init(a, f->deg);
  double_poly_set_mpz_poly (a, f);

  res = L2_lognorm_d (a, s);

  double_poly_clear(a);
  return res;
}

#ifdef OPTIMIZE_MP
/* The name is misleading. It returns the before-log part in
   1/2 log(int(int(...)))
*/
void
L2_lognorm_mp (mpz_poly_ptr f, mpz_t s, mpz_t norm)
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

static double
L2_skewness_deg6_approx (mpz_poly_srcptr f MAYBE_UNUSED, double_poly_ptr ff,
                         double_poly_ptr dff, int prec)
{
  double *dfd = dff->coeff;
  double *fd = ff->coeff;
  double s, nc, a, b, c, smin, smax;
  int sign_changes = 0;
  double q[7], logmu, best_logmu = DBL_MAX, best_s = DBL_MAX;

  dfd[6] = 99.0 * fd[6] * fd[6];
  dfd[5] = 6.0 * (2.0 * fd[4] * fd[6] + fd[5] * fd[5]);
  dfd[4] = 2.0 * (fd[2] * fd[6] + fd[3] * fd[5]) + fd[4] * fd[4];
  dfd[2] = -2.0 * (fd[0] * fd[4] + fd[1] * fd[3]) - fd[2] * fd[2];
  dfd[1] = -6.0 * (2.0 * fd[0] * fd[2] + fd[1] * fd[1]);
  dfd[0] = -99.0 * fd[0] * fd[0];
  if (dfd[1] > 0)
    sign_changes ++; /* since dfd[0] < 0 */
  if (dfd[1] * dfd[2] < 0)
    sign_changes ++;
  if (dfd[2] * dfd[4] < 0)
    sign_changes ++; /* since dfd[3] = 0 */
  if (dfd[4] * dfd[5] < 0)
    sign_changes ++;
  if (dfd[5] < 0)
    sign_changes ++; /* since dfd[6] > 0 */
  /* since dfd[6] and dfd[0] have opposite signs, we have an odd number of
     roots on [0,+inf[. Moreover since dfd[3]=0, we can't have 5 positive
     roots, thus we have either 1 or 3. */

  q[6] = dfd[6];
  q[5] = (dfd[5] < 0) ? dfd[5] : 0.0;
  q[4] = (dfd[4] < 0) ? dfd[4] : 0.0;
  q[2] = (dfd[2] < 0) ? dfd[2] : 0.0;
  q[1] = (dfd[1] < 0) ? dfd[1] : 0.0;
  q[0] = dfd[0]; /* always negative */
  s = 1.0;
  while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] < 0)
    s = s + s;
  if (s == 1.0)
    {
      while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] > 0)
        s = s * 0.5;
      s = s + s;
    }
  smax = s;

  q[6] = dfd[6]; /* always positive */
  q[5] = (dfd[5] > 0) ? dfd[5] : 0.0;
  q[4] = (dfd[4] > 0) ? dfd[4] : 0.0;
  q[2] = (dfd[2] > 0) ? dfd[2] : 0.0;
  q[1] = (dfd[1] > 0) ? dfd[1] : 0.0;
  q[0] = dfd[0];
  s = smax;
  while ((((((q[6]*s)+q[5])*s+q[4])*s*s+q[2])*s+q[1])*s+q[0] > 0)
    s = s * 0.5;
  smin = s;

  /* positive roots are in [smin, smax] */

  double v = -1.0, oldv;
  for (double t = 2.0 * smin; t <= smax; t = 2.0 * t)
    {
      /* invariant: q(smin) < 0 */
      oldv = v;
      v = (((((dfd[6]*t)+dfd[5])*t+dfd[4])*t*t+dfd[2])*t+dfd[1])*t+dfd[0];
      if (v == 0.0)
        {
          best_s = t;
          break;
        }
      if (oldv < 0 && v > 0)
        {
          /* the derivative has a root in [smin,b] and is increasing, thus
             the norm has a minimum in [smin,b], we refine it by dichotomy */
          a = smin;
          b = t;
          for (int i = 0; i < prec; i++)
            {
              c = (a + b) * 0.5;

              nc = ((((dfd[6] * c + dfd[5]) * c + dfd[4]) * c * c + dfd[2]) * c
                    + dfd[1]) * c + dfd[0];
              if (nc > 0)
                b = c;
              else
                a = c;
            }
          s = sqrt ((a + b) * 0.5);
          if (sign_changes == 1)
            {
              best_s = s;
              break;
            }
          logmu = L2_lognorm_d (ff, s);
          if (logmu < best_logmu)
            {
              best_logmu = logmu;
              best_s = s;
            }
        }
      smin = t; /* to avoid hitting twice the same root of the derivative */
    }

  return best_s;
}

double
L2_skewness_deg6 (mpz_poly_ptr f MAYBE_UNUSED, double_poly_srcptr ff,
                  double_poly_srcptr dff MAYBE_UNUSED, int prec MAYBE_UNUSED)
{
  double s, logmu, logmu_min = DBL_MAX, s_min = 1.0;
  mpz_poly df;
  usp_root_data Roots[6];
  int i, k;

  mpz_poly_init (df, 6);
  for (i = 0; i < 6; i++)
    usp_root_data_init (Roots + i);

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
#if 1 /* using numberOfRealRoots */
  mpz_mul (df->coeff[6], f->coeff[6], f->coeff[6]);
  mpz_mul_ui (df->coeff[6], df->coeff[6], 99); /* 99*a6^2 */
  mpz_mul (df->coeff[5], f->coeff[4], f->coeff[6]);
  mpz_mul_2exp (df->coeff[5], df->coeff[5], 1);
  mpz_addmul (df->coeff[5], f->coeff[5], f->coeff[5]);
  mpz_mul_ui (df->coeff[5], df->coeff[5], 6); /* 6*(2*a4*a6 + a5^2) */
  mpz_mul (df->coeff[4], f->coeff[2], f->coeff[6]);
  mpz_addmul (df->coeff[4], f->coeff[3], f->coeff[5]);
  mpz_mul_2exp (df->coeff[4], df->coeff[4], 1);
  mpz_addmul (df->coeff[4], f->coeff[4], f->coeff[4]); /*2*a2*a6+2*a3*a5+a4^2*/
  mpz_set_ui (df->coeff[3], 0);
  mpz_mul (df->coeff[2], f->coeff[0], f->coeff[4]);
  mpz_addmul (df->coeff[2], f->coeff[1], f->coeff[3]);
  mpz_mul_2exp (df->coeff[2], df->coeff[2], 1);
  mpz_addmul (df->coeff[2], f->coeff[2], f->coeff[2]);
  mpz_neg (df->coeff[2], df->coeff[2]); /* -(2*a0*a4+2*a1*a3+a2^2) */
  mpz_mul (df->coeff[1], f->coeff[0], f->coeff[2]);
  mpz_mul_2exp (df->coeff[1], df->coeff[1], 1);
  mpz_addmul (df->coeff[1], f->coeff[1], f->coeff[1]);
  mpz_mul_si (df->coeff[1], df->coeff[1], -6); /* -6*(2*a0*a2 + a1^2) */
  mpz_mul (df->coeff[0], f->coeff[0], f->coeff[0]);
  mpz_mul_si (df->coeff[0], df->coeff[0], -99); /* -99*a0^2 */

  df->deg = 6;
  k = numberOfRealRoots (df->coeff, 6, 0.0, 0, Roots);
  int kpos = 0;
  for (i = 0; i < k; i++)
    if (mpz_sgn (Roots[i].b) > 0)
      {
        kpos ++;
        s = rootRefine (Roots + i, df->coeff, 6, ldexp (1.0, -prec));
        s = sqrt (s);
        logmu = L2_lognorm_d (ff, s);
        if (logmu < logmu_min)
          {
            logmu_min = logmu;
            s_min = s;
          }
      }
#else /* using double_poly_compute_roots */
  double *dfd = dff->coeff;
  double *fd = ff->coeff;
  double roots[6];

  dfd[6] = 99.0 * fd[6] * fd[6];
  dfd[5] = 6.0 * ( 2.0 * fd[4] * fd[6] + fd[5] * fd[5] );
  dfd[4] = 2.0 * ( fd[2] * fd[6] + fd[3] * fd[5] ) + fd[4] * fd[4];
  dfd[3] = 0.0;
  dfd[2] = -2.0 * ( fd[0] * fd[4] + fd[1] * fd[3] ) - fd[2] * fd[2];
  dfd[1] = -6.0 * ( 2.0 * fd[0] * fd[2] + fd[1] * fd[1] );
  dfd[0] = -99.0 * fd[0] * fd[0];
  double B = double_poly_bound_roots (dff);
  k = double_poly_compute_roots (roots, dff, B);
  ASSERT_ALWAYS(k > 0);
  for (i = 0; i < k; i++)
    {
      s = sqrt (roots[i]);
      logmu = L2_lognorm_d (ff, s);
      if (logmu < logmu_min)
        {
          logmu_min = logmu;
          s_min = s;
        }
    }
#endif

  mpz_poly_clear (df);
  for (i = 0; i < 6; i++)
    usp_root_data_clear (Roots + i);

  return s_min;
}

/* return the skewness giving the best lognorm sum for two polynomials,
   by using trichotomy between the optimal skewness of both polynomials */
double
L2_combined_skewness2 (mpz_poly_srcptr f, mpz_poly_srcptr g, int prec)
{
  double a, b, c, d, va, vb, vc, vd;

  a = L2_skewness (f, prec);
  b = L2_skewness (g, prec);

  if (b < a)
    {
      c = b;
      b = a;
      a = c;
    }

  ASSERT(a <= b);

  va = L2_lognorm (f, a) + L2_lognorm (g, a);
  vb = L2_lognorm (f, b) + L2_lognorm (g, b);

  while (b - a > ldexp (a, -prec))
    {
      c = (2.0 * a + b) / 3.0;
      vc = L2_lognorm (f, c) + L2_lognorm (g, c);

      d = (a + 2.0 * b) / 3.0;
      vd = L2_lognorm (f, d) + L2_lognorm (g, d);

      if (va < vd && va < vb && vc < vd && vc < vb) /* minimum is in a or c */
        {
          b = d;
          vb = vd;
        }
      else /* the minimum is in d or b */
        {
          a = c;
          va = vc;
        }
    }
  return (a + b) * 0.5;
}

/* Use derivative test, with ellipse regions */
double
L2_skewness (mpz_poly_srcptr f, int prec)
{
  double_poly ff, df;
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
      s = L2_skewness_deg6_approx (f, ff, df, prec);
      goto end;
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
         var('r,s,t,y')
         F = expand(f(x=x/y)*y^d)
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
         var('r,s,t,y')
         F = expand(f(x=x/y)*y^d)
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
         We get dv = 6*a_2^2*s^4 - 6*a_0^2
         thus the optimal skewness is sqrt(|a0|/|a2|).
      */
      a = b = sqrt (fabs (fd[0] / fd[2]));
    }
  else /* d == 1 */
    a = b = fabs (fd[0] / fd[1]);

  s = (a + b) * 0.5;

 end:
  double_poly_clear (ff);
  double_poly_clear (df);

  return s;
}

double L2_skew_lognorm (mpz_poly_srcptr f, int prec)
{
  return L2_lognorm (f, L2_skewness (f, prec));
}

#ifdef OPTIMIZE_MP

/* Use derivative test, with ellipse regions */
void
L2_skewness_derivative_mp (mpz_poly_ptr F, int prec, mpz_t skewness)
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

  int i, d = F->deg;
  mpz_t *f = F->coeff;

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
static double
special_val0 (mpz_poly_srcptr f, unsigned long p)
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
      mpz_poly_eval_diff_ui (c, g, r);
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
special_valuation (mpz_poly_srcptr f, unsigned long p, mpz_t disc)
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
get_alpha (mpz_poly_srcptr f, unsigned long B)
{
  double alpha, e;
  unsigned long p;
  mpz_t disc;

  /* for F linear, we have q_p = 1 for all p, thus
     alpha(F) = sum(prime p <= B, log(p)/(p^2-1)) ~ 0.569959993064325 */
  if (f->deg == 1)
    return 0.569959993064325;

  mpz_init (disc);
  mpz_poly_discriminant (disc, f);

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
special_valuation_affine (mpz_poly_srcptr f, unsigned long p, mpz_t disc )
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
get_biased_alpha_projective (mpz_poly_srcptr f, unsigned long B)
{
   double alpha, e;
   unsigned long p;
   mpz_t disc;

   mpz_init (disc);
   mpz_poly_discriminant (disc, f);

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


/**************************** rotation ***************************************/

/* replace f + k0 * x^t * (b*x + m) by f + k * x^t * (b*x + m), and return k */
long
rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t)
{
  /* Warning: k - k0 might not be representable in a long! */
  unsigned long diff;
  if (k >= k0)
    {
      diff = k - k0; /* k - k0 always fits in an unsigned long */
      mpz_addmul_ui (f[t + 1], b, diff);
      mpz_addmul_ui (f[t], m, diff);
    }
  else
    {
      diff = k0 - k;
      mpz_submul_ui (f[t + 1], b, diff);
      mpz_submul_ui (f[t], m, diff);
    }
  return k;
}

/* replace f by f + k * x^t * (b*x + g0) */
void
rotate_auxg_z (mpz_t *f, const mpz_t b, const mpz_t g0, const mpz_t k, unsigned int t)
{
  mpz_addmul (f[t + 1], b, k);
  mpz_addmul (f[t], g0, k);
}

/* replace f by f - k * x^t * (b*x + g0) */
void
derotate_auxg_z (mpz_t *f, const mpz_t b, const mpz_t g0, const mpz_t k, unsigned int t)
{
  mpz_submul (f[t + 1], b, k);
  mpz_submul (f[t], g0, k);
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
   mpz_poly F, G;

   F->coeff = p->pols[ALG_SIDE]->coeff;
   F->deg = p->pols[ALG_SIDE]->deg;
   G->coeff = p->pols[RAT_SIDE]->coeff;
   G->deg = p->pols[RAT_SIDE]->deg;

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
   fprintf (fp, "# MurphyE(Bf=%.0f,Bg=%.0f,area=%.3e)=%.2e\n",
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
   fprintf (fp, "# generated by %s: %s", cado_revision_string, argv[0]);
   for (i = 1; i < argc; i++)
      fprintf (fp, " %s", argv[i]);
   fprintf (fp, " in %.2fs\n", (seconds () - st));
}


/*
  Call print_cadopoly, given f, g and return MurphyE.
*/
double
print_poly_fg (mpz_poly_srcptr f, mpz_t *g, mpz_t N, int mode)
{
   double e;
   int i;
   int d = f->deg;

   cado_poly cpoly;
   cado_poly_init(cpoly);
   for (i = 0; i < (d + 1); i++)
      mpz_set(cpoly->pols[ALG_SIDE]->coeff[i], f->coeff[i]);
   for (i = 0; i < 2; i++)
      mpz_set(cpoly->pols[RAT_SIDE]->coeff[i], g[i]);
   mpz_set(cpoly->n, N);
   cpoly->skew = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
   cpoly->pols[ALG_SIDE]->deg = d;
   cpoly->pols[RAT_SIDE]->deg = 1;

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
do_translate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k)
{
  int i, j;
  int d = f->deg;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_addmul (f->coeff[j], f->coeff[j+1], k);
  mpz_addmul (g[0], g[1], k);
}

/* f <- f(x-k), g <- g(x-k) */
void
do_detranslate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k)
{
  int i, j;
  int d = f->deg;

  for (i = d - 1; i >= 0; i--)
    for (j = i; j < d; j++)
      mpz_submul (f->coeff[j], f->coeff[j+1], k);
  mpz_submul (g[0], g[1], k);
}


/* If final <> 0, print the real value of E (root-optimized polynomial),
   otherwise print the expected value of E. Return E or exp_E accordingly.
   TODO: adapt for more than 2 polynomials and two algebraic polynomials */
double
cado_poly_fprintf_with_info (FILE *fp, cado_poly_ptr poly, const char *prefix,
                             int final)
{
  unsigned int nrroots;
  double lognorm, alpha, alpha_proj, exp_E;

  nrroots = numberOfRealRoots (poly->pols[ALG_SIDE]->coeff, poly->pols[ALG_SIDE]->deg, 0, 0, NULL);
  if (poly->skew <= 0.0) /* If skew is undefined, compute it. */
    poly->skew = L2_skewness (poly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);
  lognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew);
  alpha = get_alpha (poly->pols[ALG_SIDE], ALPHA_BOUND);
  alpha_proj = get_biased_alpha_projective (poly->pols[ALG_SIDE], ALPHA_BOUND);
  exp_E = (final) ? 0.0 : lognorm
    + expected_rotation_gain (poly->pols[ALG_SIDE], poly->pols[RAT_SIDE]);

  cado_poly_fprintf (stdout, poly, prefix);
  cado_poly_fprintf_info (fp, lognorm, exp_E, alpha, alpha_proj, nrroots,
                          prefix);
  return (final) ? lognorm + alpha : exp_E;
}

/* TODO: adapt for more than 2 polynomials and two algebraic polynomials */
double
cado_poly_fprintf_with_info_and_MurphyE (FILE *fp, cado_poly_ptr poly,
                                         double MurphyE, double bound_f,
                                         double bound_g, double area,
                                         const char *prefix)
{
  double exp_E;
  exp_E = cado_poly_fprintf_with_info (fp, poly, prefix, 1);
  cado_poly_fprintf_MurphyE (fp, MurphyE, bound_f, bound_g, area, prefix);
  return exp_E;
}

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

/* compute largest interval kmin <= k <= kmax such that when we add k*x^i*g(x)
   to f(x), the lognorm does not increase more than margin */
void
expected_growth (rotation_space *r, mpz_poly_srcptr f, mpz_poly_srcptr g, int i,
                 double margin)
{
  double s = L2_skewness (f, SKEWNESS_DEFAULT_PREC);
  double n = L2_lognorm (f, s), n2;
  mpz_t fi, fip1, kmin, kmax, k;

  mpz_init_set (fi, f->coeff[i]);
  mpz_init_set (fip1, f->coeff[i+1]);
  mpz_init (kmin);
  mpz_init (kmax);
  mpz_init (k);

  /* negative side */
  mpz_set_si (kmin, -1);
  for (;;)
    {
      mpz_set (f->coeff[i], fi);
      mpz_set (f->coeff[i+1], fip1);
      rotate_auxg_z (f->coeff, g->coeff[1], g->coeff[0], kmin, i);
      n2 = L2_lognorm (f, s);
      if (n2 > n + margin)
        break;
      mpz_mul_2exp (kmin, kmin, 1);
    }
  /* now kmin < k < kmin/2 */
  mpz_tdiv_q_2exp (kmax, kmin, 1);
  while (1)
    {
      mpz_add (k, kmin, kmax);
      mpz_div_2exp (k, k, 1);
      if (mpz_cmp (k, kmin) == 0 || mpz_cmp (k, kmax) == 0)
        break;
      mpz_set (f->coeff[i], fi);
      mpz_set (f->coeff[i+1], fip1);
      rotate_auxg_z (f->coeff, g->coeff[1], g->coeff[0], k, i);
      n2 = L2_lognorm (f, s);
      if (n2 > n + margin)
        mpz_set (kmin, k);
      else
        mpz_set (kmax, k);
    }
  r->kmin = mpz_get_d (kmax);

  /* positive side */
  mpz_set_ui (kmax, 1);
  for (;;)
    {
      mpz_set (f->coeff[i], fi);
      mpz_set (f->coeff[i+1], fip1);
      rotate_auxg_z (f->coeff, g->coeff[1], g->coeff[0], kmax, i);
      n2 = L2_lognorm (f, s);
      if (n2 > n + margin)
        break;
      mpz_mul_2exp (kmax, kmax, 1);
    }
  /* now kmax < k < kmax/2 */
  mpz_tdiv_q_2exp (kmin, kmax, 1);
  while (1)
    {
      mpz_add (k, kmin, kmax);
      mpz_div_2exp (k, k, 1);
      if (mpz_cmp (k, kmin) == 0 || mpz_cmp (k, kmax) == 0)
        break;
      mpz_set (f->coeff[i], fi);
      mpz_set (f->coeff[i+1], fip1);
      rotate_auxg_z (f->coeff, g->coeff[1], g->coeff[0], k, i);
      n2 = L2_lognorm (f, s);
      if (n2 > n + margin)
        mpz_set (kmax, k);
      else
        mpz_set (kmin, k);
    }
  r->kmax = mpz_get_d (kmin);

  /* reset f[i] and f[i+1] */
  mpz_set (f->coeff[i], fi);
  mpz_set (f->coeff[i+1], fip1);

  mpz_clear (fi);
  mpz_clear (fip1);
  mpz_clear (kmin);
  mpz_clear (kmax);
  mpz_clear (k);
}

/* for a given pair (f,g), tries to estimate the value of alpha one might
   expect from rotation (including the projective alpha) */
double
expected_rotation_gain (mpz_poly_srcptr f, mpz_poly_srcptr g)
{
  double S = 1.0, s, incr = 0.0;
  rotation_space r;
  double proj_alpha = get_biased_alpha_projective (f, ALPHA_BOUND_SMALL);

  for (int i = 0; 2 * i < f->deg; i++)
    {
      expected_growth (&r, f, g, i, NORM_MARGIN);
      s = r.kmax - r.kmin + 1.0;
      S *= s;
      /* assume each non-zero rotation increases on average by NORM_MARGIN/2 */
      if (s >= 2.0)
        incr += NORM_MARGIN / 2.0;
    }
  return proj_alpha + expected_alpha (S) + incr;
}
