/* Rootsieve algorithm.

Copyright 2010 Paul Zimmermann (using notes from Emmanuel Thome)

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
#include <math.h> /* for log */
#include "utils.h"
#include "modul_poly.h"
#include "auxiliary.h"

// #define KK -16346
// #define PP 7

/* if gpn and pn are not coprime, let g = gcd(gpn, pn):
   we must have g*g'*v = r + j*g*p' with g' = gpn/g and p'=pn/g:
   (a) if g does not divide r, there is no solution
   (b) if g divides r, then we must have:
   g'*v = r' + j*p' with r'=r/g, let v0 = r'/g' mod p',
   then v0, v0+p', ..., v0+(g-1)*p' are the g solutions. */
static void
special_update (double *A, long K0, long K1, const residueul_t gpn,
		modulusul_t Pn, residueul_t r, double alpha_p)
{
  modulusul_t Pp;
  residueul_t rp, gp;
  int ret;
  long k0, k;
  modintul_t gg;
  unsigned long g, rr, pp, pn;

  pn = modul_getmod_ul (Pn);
  modul_gcd (gg, gpn, Pn);
  g = gg[0];
  rr = modul_get_ul (r, Pn);
  if (rr % g == 0)
    {
      pp = pn / g;
      modul_initmod_ul (Pp, pp);
      modul_init (gp, Pp);
      modul_init (rp, Pp);
      modul_set_ul (gp, modul_get_ul (gpn, Pn) / g, Pp);
      modul_set_ul (rp, rr / g, Pp);
      ret = modul_inv (gp, gp, Pp);
      ASSERT_ALWAYS (ret != 0);
      modul_mul (rp, rp, gp, Pp);
      k0 = modul_get_ul (rp, Pp);
      for (k = k0; k <= K1; k += pp)
	{
	  A[k - K0] -= alpha_p;
#if (defined(KK) && defined(PP))
	  if (k == KK && (pn % PP) == 0)
	    fprintf (stderr, "subtract %f to AA[%d] for pp=%lu, now %f\n",
		     alpha_p, KK, pp, A[k - K0]);
#endif
	}
      for (k = k0 - pp; k >= K0; k -= pp)
	{
	  A[k - K0] -= alpha_p;
#if (defined(KK) && defined(PP))
	  if (k == KK && (pn % PP) == 0)
	    fprintf (stderr, "subtract %f to AA[%d] for pp=%lu, now %f\n",
		     alpha_p, KK, pp, A[k - K0]);
#endif
	}
      modul_clear (gp, Pp);
      modul_clear (rp, Pp);
      modul_clearmod (Pp);
    }
}

/* f(x) of degree d is the polynomial f0(x) + j*x*g(x)
   where g(x) = b*x-m.
   A is a table of K1-K0+1 entries, where A[k-K0] corresponds to the alpha
   contribution for f(x) + k*g(x)
   p is the current prime, we consider pn = p^n here */
void
update_table (mpz_t *f, int d, mpz_t m, mpz_t b, double *A, long K0, long K1,
	      unsigned long pn, double alpha_p)
{
  long k, k0;
  modul_poly_t fpn; /* f mod p^n */
  modulusul_t Pn;
  residueul_t l, gpn, bpn, r, v;
  int ret;

  modul_initmod_ul (Pn, pn);
  modul_init (gpn, Pn);
  modul_init (bpn, Pn);
  modul_init (r, Pn);
  modul_init (l, Pn);
  modul_init (v, Pn);

  modul_poly_init (fpn, d);

  /* first reduce f(x) and g(x) mod p^n */
  modul_poly_set_mod_raw (fpn, f, d, Pn);

  modul_set_ul (gpn, mpz_fdiv_ui (m, pn), Pn);
  /* invariant: gpn = -g(l) */
  modul_set_ul (bpn, mpz_fdiv_ui (b, pn), Pn);
  modul_set_ul (l, 0, Pn);

  for (;;)
    {
      modul_poly_eval (r, fpn, l, Pn);
      /* invariant: gpn = -g(l) */
      /* we want r = v*gpn, i.e., v = r/gpn; r is never 0 otherwise f(x) would
	 be divisible by p^n, but gpn = g(l) can be 0 */
      ret = modul_intcmp_ul (gpn, 0);
      if (ret != 0)
	{
	  ret = modul_inv (v, gpn, Pn); /* FIXME: use batch inversion */
	  if (ret == 0) /* gpn and pn are not coprime */
	    special_update (A, K0, K1, gpn, Pn, r, alpha_p);
	  else
	    {
	      modul_mul (v, v, r, Pn);

	      k0 = modul_get_ul (v, Pn);
	      for (k = k0; k <= K1; k += pn)
		{
		  A[k - K0] -= alpha_p;
#if (defined(KK) && defined(PP))
		  if (k == KK && (pn % PP) == 0)
		    fprintf (stderr, "subtract %f to AA[%d] for l=%lu, pn=%lu, now %f\n",
			     alpha_p, KK, modul_get_ul (l, Pn), pn, A[k - K0]);
#endif
		}
	      for (k = k0 - (long) pn; k >= K0; k -= (long) pn)
		{
		  A[k - K0] -= alpha_p;
#if (defined(KK) && defined(PP))
		  if (k == KK && (pn % PP) == 0)
		    fprintf (stderr, "subtract %f to AA[%d] for l=%lu, pn=%lu, now %f\n",
			     alpha_p, KK, modul_get_ul (l, Pn), pn, A[k - K0]);
#endif
		}
	    }
	}

      /* since g(x) = b*x-m, -g(l) = m-b*l */
      modul_sub (gpn, gpn, bpn, Pn);
      modul_add_ul (l, l, 1, Pn);
      if (modul_intcmp_ul (l, 0) == 0)
	break;
    }

  modul_clear (gpn, Pn);
  modul_clear (bpn, Pn);
  modul_clear (r, Pn);
  modul_clear (l, Pn);
  modul_clear (v, Pn);
  modul_clearmod (Pn);
  modul_poly_clear (fpn);
}

#define average_valuation_affine special_val0

/* Computes the contribution at p of the alpha value of f (projective part
   only). The return value is always <= 0. */
double
alpha_p_projective (mpz_t *f, int d, mpz_t disc, unsigned long p)
{
  double e;
  int i;

  if (mpz_divisible_ui_p (disc, p))
    {
      mpz_array_t *q;
      mpz_t c;

      /* compute q = reverse(f)(px): q[0] = f[d], q[1] = f[d-1]*p, ...,
	 q[d] = f[0]*p^d */
	q = alloc_mpz_array (d + 1);
	mpz_init_set_ui (c, 1);
	for (i = 0; i <= d; i++)
	  {
	    mpz_mul ((q->data)[i], f[d-i], c);
	    mpz_mul_ui (c, c, p);
	  }
	e = average_valuation_affine (q->data, d, p) / (double) (p + 1);
	mpz_clear (c);
	clear_mpz_array (q);
	return -e * log ((double) p);
    }
  else
    {
      for (i = d; i >= 0 && mpz_divisible_ui_p (f[i], p); i--);
      ASSERT_ALWAYS (i >= 0);
      ASSERT_ALWAYS (d - i <= 1);
      return - (double) (d - i) * (double) p / (double) (p + 1)
	* log ((double) p) / (double) (p - 1);
    }
}
