/* Compute Murphy's E-value.

Copyright 2010, 2013 Paul Zimmermann

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

#include "cado.h" // IWYU pragma: keep
#define PI 3.14159265358979324

#include <math.h>
#include "auxiliary.h"
#include "rho.h"
#include "murphyE.h"
#include "double_poly.h"
#include "mpz_poly.h"

double
MurphyE (cado_poly_srcptr cpoly, double Bf, double Bg, double area, int K,
         unsigned long B)
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
  alpha_f = get_alpha (cpoly->pols[ALG_SIDE], B);
  alpha_g = get_alpha (cpoly->pols[RAT_SIDE], B);
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
