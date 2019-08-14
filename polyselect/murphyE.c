/* Compute Murphy's E-value.

Copyright 2010-2019 Paul Zimmermann and Nicolas David.

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

/* Murphy's E-value is defined on pages 86 and 87 of Murphy's thesis:

   E(f,g) = sum(rho(u_f(theta_i))*rho(u_g(theta_i)), i=1..K)
   
   where theta_i = Pi/K*(i-1/2)

   and u_f(theta_i) = (log(|F(cos(theta_i)*s^(1/2),sin(theta_i)/s^(1/2))|)
                       + alpha_f)/log(B_f)

   where s is the skewness, F(x,y) is the bivariate polynomial associated to f,
   alpha_f is the alpha-value for f, B_f is the smoothness bound associated
   to f (idem for g).

   This bound depends on the smoothness bounds B_f and B_g.

   In his pol51opt code, Kleinjung uses B_f = 1e7 and B_g = 5e6.
   He also scales x=cos(theta_i)*s^(1/2) and y=sin(theta_i)/s^(1/2) by
   a factor sqrt(area), with area = 1e16.
*/

#include "cado.h"
#define PI 3.14159265358979324

#include <math.h>
#include "utils.h"
#include "auxiliary.h"
#include "rho.h"

/* define USE_VARIANT to use the simple variant of E', that consists in using
   MurphyE, with alpha replaced by alpha - std, where std is the standard
   deviation corresponding to the measure whose mean is alpha */
// #define USE_VARIANT

double
MurphyE (cado_poly cpoly, double Bf, double Bg, double area, int K)
{
  double E = 0, x, y, ti;
  double alpha_f, alpha_g, xi, yi, vf, vg;
  double one_over_logBf, one_over_logBg;
  double_poly f, g;

  x = sqrt (area * cpoly->skew);
  y = sqrt (area / cpoly->skew);
  double_poly_init (f, cpoly->pols[ALG_SIDE]->deg);
  double_poly_init (g, cpoly->pols[RAT_SIDE]->deg);
  double_poly_set_mpz_poly (f, cpoly->pols[ALG_SIDE]);
  double_poly_set_mpz_poly (g, cpoly->pols[RAT_SIDE]);
  alpha_f = get_alpha (cpoly->pols[ALG_SIDE], ALPHA_BOUND);
#ifdef USE_VARIANT
  double mu, v;
  mu = dist_alpha (cpoly->pols[ALG_SIDE], ALPHA_BOUND, &v);
  ASSERT_ALWAYS(alpha_f == mu);
  ASSERT_ALWAYS(v >= 0);
  /* patch: replace alpha_f by alpha_f - stddev */
  alpha_f -= sqrt (v);
#endif
  alpha_g = get_alpha (cpoly->pols[RAT_SIDE], ALPHA_BOUND);
  one_over_logBf = 1.0 / log (Bf);
  one_over_logBg = 1.0 / log (Bg);
  for (int i = 0; i < K; i++)
    {
      ti = PI / (double) K * ((double) i + 0.5);
      xi = x * cos (ti);
      yi = y * sin (ti);

      vf = double_poly_eval (f, xi / yi) * pow (yi, f->deg);
      vg = double_poly_eval (g, xi / yi) * pow (yi, g->deg);

      vf = log (fabs (vf)) + alpha_f;
      vg = log (fabs (vg)) + alpha_g;

      vf *= one_over_logBf;
      vg *= one_over_logBg;

      E += dickman_rho (vf) * dickman_rho (vg);
    }
  double_poly_clear (f);
  double_poly_clear (g);

  return E / (double) K;
}

/* Return bessel_I (a, y) with tolerance 'epsilon':
   bessel_I (a, y) = (y/2)^a * sum ((y^2/4)^j/factorial(j)/gamma(a+j+1),
   j, 0, infinity).
*/
static double
bessel_I (double a, double y, double epsilon)
{
  double s, c, z, j = 1.0;

  if (y < epsilon * sqrt(a + 1.0))
    return pow (y / 2, a) / tgamma (a + 1.0); /* order-1 expansion around 0 */
  else if (epsilon * y > fabs (a * a - 0.25))
    return exp (y) / sqrt (2 * PI * y); /* asymptotic expansion */
  else
    {
      s = c = 1.0 / tgamma (a + 1.0);
      z = y * y * 0.25;
      while (fabs (c) > epsilon * fabs (s))
        {
          c *= z / (j * (a + j));
          s = s + c;
          j += 1.0;
        }
      return pow (y / 2, a) * s;
    }
}

/* return the pdf (probability density function) at x of the chi2 noncentral
   distribution with parameters 'k' and 'lam', with tolerance 'epsilon':
   1/2*exp(-(x+lam)/2)*(x/lam)^(k/4-1/2)*bessel_I(k/2-1,sqrt(lam*x)).
   With epsilon=0, we should find the same values as with Sage:
   from scipy.stats import ncx2
   ncx2.pdf(x, k, lam) */
double
ncx2_pdf (double x, double k, double lam, double epsilon)
{
  double I = bessel_I (0.5 * k - 1.0, sqrt (lam * x), epsilon);
  double ret = 0.5 * exp (-0.5 * (x + lam)) * pow (x / lam, 0.25 * k - 0.5) * I;
  return ret;
}





/* return the E_value with a non central chi2 density for \alpha.
   It is an alternative of the MurphyE function. */
double
MurphyE_chi2 (cado_poly cpoly, double Bf, double Bg, double area, int K)
{
#ifdef USE_VARIANT
  /* temporary patch: we compute the classical MurphyE, but with alpha
     replaced by alpha - std(alpha) */
  return MurphyE (cpoly, Bf, Bg, area, K);
#else
  double E = 0.0, x, y, ti, h;
  double alpha_f, alpha_g, xi, yi, vf, vg, cof;
  double one_over_logBf, one_over_logBg;
  double_poly f, g;
  double k, lam, mu, v;
  unsigned long p;

  /* sum([1.0/(p-1)*log(p*1.0) for p in prime_range(ALPHA_BOUND)]) */
  for (cof = 0.0, p = 2; p < ALPHA_BOUND; p += 1 + (p > 2))
    if (ulong_isprime (p))
      cof += log ((double) p) / (double) (p - 1);

  mu = dist_alpha (cpoly->pols[ALG_SIDE], ALPHA_BOUND, &v);
  mu = cof - mu;
  lam = v / 2 - mu;
  k = mu - lam;

  /* For K=200, the time spent in MurphyE_chi2() within polyselect_ropt is only
     about 10%, thus since K=1000 by default, we divide it by 5. */
  K /= 5;

  x = sqrt (area * cpoly->skew);
  y = sqrt (area / cpoly->skew);
  h = (double) K / 100.0;
  double_poly_init (f, cpoly->pols[ALG_SIDE]->deg);
  double_poly_init (g, cpoly->pols[RAT_SIDE]->deg);
  double_poly_set_mpz_poly (f, cpoly->pols[ALG_SIDE]);
  double_poly_set_mpz_poly (g, cpoly->pols[RAT_SIDE]);
  alpha_g = get_alpha (cpoly->pols[RAT_SIDE], ALPHA_BOUND);
  one_over_logBf = 1.0 / log (Bf);
  one_over_logBg = 1.0 / log (Bg);
  /* here j / h fits with the integration parameter with a non central
     chi2 measure */
  double *pdf = malloc (K * sizeof (double));
  ASSERT_ALWAYS(pdf != NULL);
  for (int j = 1; j < K; j++)
    pdf[j] = ncx2_pdf (j / h, k, lam, 0.001);
  for (int i = 0; i < K; i++)
    {
	  ti = PI / (double) K * ((double) i + 0.5);
	  xi = x * cos (ti);
	  yi = y * sin (ti);

	  vf = double_poly_eval (f, xi / yi) * pow (yi, f->deg);
	  vg = double_poly_eval (g, xi / yi) * pow (yi, g->deg);

          vf = log (fabs (vf));

	  vg = log (fabs (vg)) + alpha_g;
	  vg *= one_over_logBg;
          vg = dickman_rho (vg);

          for (int j = 1; j < K; j++)
            {
              double vfj;

              alpha_f = cof - (double) j / h;
              vfj = vf + alpha_f;
              vfj *= one_over_logBf;

              double t = dickman_rho (vfj) * vg * pdf[j];

              E += t;
            }
    }
  free (pdf);
  double_poly_clear (f);
  double_poly_clear (g);

  E = E / (K * h);

  return E;
#endif
}
