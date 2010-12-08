/* Compute Murphy's E-value.

Copyright 2010 Paul Zimmermann

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

double
MurphyE (cado_poly cpoly, double Bf, double Bg, double area, int K)
{
  double E = 0, x, y, ti;
  int i, deg_g = 1;
  double *f, *g, alpha_f, alpha_g, xi, yi, vf, vg;
  double one_over_logBf, one_over_logBg;

  x = sqrt (area * cpoly->skew);
  y = sqrt (area / cpoly->skew);
  f = malloc ((cpoly->degree + 1) * sizeof (double));
  g = malloc ((deg_g + 1) * sizeof (double));
  for (i = 0; i <= cpoly->degree; i++)
    f[i] = mpz_get_d (cpoly->f[i]);
  for (i = 0; i <= deg_g; i++)
    g[i] = mpz_get_d (cpoly->g[i]);
  alpha_f = get_alpha (cpoly->f, cpoly->degree, ALPHA_BOUND);
  alpha_g = get_alpha (cpoly->g, deg_g, ALPHA_BOUND);
  one_over_logBf = 1.0 / log (Bf);
  one_over_logBg = 1.0 / log (Bg);
  for (i = 0; i < K; i++)
    {
      ti = PI / (double) K * ((double) i + 0.5);
      xi = x * cos (ti);
      yi = y * sin (ti);

      vf = fpoly_eval (f, cpoly->degree, xi / yi) * pow (yi, cpoly->degree);
      vg = fpoly_eval (g, deg_g, xi / yi) * pow (yi, deg_g);

      vf = log (fabs (vf)) + alpha_f;
      vg = log (fabs (vg)) + alpha_g;

      vf *= one_over_logBf;
      vg *= one_over_logBg;

      E += dickman_rho (vf) * dickman_rho (vg);
    }
  free (f);
  free (g);
  return E / (double) K;
}
