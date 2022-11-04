#include "cado.h"
#include <float.h>
#include <math.h>
#include "double_poly.h"
#include "usp.h"

#include "polyselect_norms.h"

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
  k = numberOfRealRoots ((const mpz_t *) df->coeff, 6, 0.0, 0, Roots);
  for (i = 0; i < k; i++)
    if (mpz_sgn (Roots[i].b) > 0)
      {
        s = rootRefine (Roots + i, (const mpz_t *) df->coeff, 6, ldexp (1.0, -prec));
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

