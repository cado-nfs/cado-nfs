/* This uses Paul Zimmermann implementation. I kept the header here:

   LLL using exact multiprecision arithmetic.
   Translated from NTL 4.1a <http://www.shoup.net/>
   into GMP <http://www.swox.se/gmp/> 
   by Paul Zimmermann, July 2000.

   Revised April 4, 2002 (bug found by Jens Franke <franke (at) math
   (dot) uni-bonn (dot) de>).

   This program is open-source software distributed under the terms 
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   Usage: lll <mat_size> [a] [b] < file

   mat_size - size of the input matrix
   a, b - positive integer coefficients used for swapping vectors. 
   We should have 1/4 < delta=a/b <= 1. The closest delta is from 1,
   the shortest are the output vectors, but the computation takes longer.
 */

#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "lll.h"

#define swap(x, y) { long _tmp = (x); (x) = (y); (y) = _tmp; }
#define mpz_swap_n(x, y) { mpz_t *_tmp = (x); (x) = (y); (y) = _tmp; }

#define ROW /* each row represents a vector */

void
LLL_init (mat_Z *m, int nRows, int nCols)
{
  m->NumRows = nRows;
  m->NumCols = nCols;
  /* indices in the LLL code start at 1, thus we allocate one more row */
  m->coeff = (mpz_t **) malloc ((m->NumRows + 1) * sizeof(mpz_t*));
  ASSERT_ALWAYS (m->coeff != NULL);
  for (int j = 1; j <= m->NumRows; j++)
  {
    m->coeff[j] = (mpz_t *) malloc ((m->NumCols + 1) * sizeof(mpz_t));
    ASSERT_ALWAYS (m->coeff[j] != NULL);
    for (int i = 1; i <= m->NumCols; i++)
      mpz_init (m->coeff[j][i]);
  }
}

void
LLL_clear (mat_Z *m)
{
  for (int j = 1; j <= m->NumRows; j++)
  {
    for (int i = 1; i <= m->NumCols; i++)
      mpz_clear (m->coeff[j][i]);
    free (m->coeff[j]);
  }
  free (m->coeff);
}

static void
ident (mat_Z X, long n)
{
   long i, j;
   X.NumRows = X.NumCols = n;
  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= n; j++)  
         if (i == j)
	    mpz_set_ui (X.coeff[i][j], 1);
         else  
            mpz_set_ui (X.coeff[i][j], 0);
} 

static void
InnerProduct (mpz_t x, mpz_t *a, mpz_t *b, long n)
{
   long i;

   mpz_mul (x, a[1], b[1]);
   for (i = 2; i <= n; i++)
     mpz_addmul (x, a[i], b[i]);
}

static void
IncrementalGS (mat_Z B, long *P, mpz_t *D, mpz_t **lam, long *s, long k)
{
   long n = B.NumCols;
   mpz_t u, t1;
   long i, j, posj;

   mpz_init (u);
   mpz_init (t1);

   for (j = 1; j <= k-1; j++) {
      posj = P[j];
      if (posj == 0) continue;

      InnerProduct (u, B.coeff[k], B.coeff[j], n);
      for (i = 1; i <= posj-1; i++) {
         mpz_mul (t1, D[i], u);
         // coverity[copy_paste_error]
         mpz_submul (t1, lam[k][i], lam[j][i]);
         mpz_div (t1, t1, D[i-1]);
         mpz_set (u, t1);
      }

      mpz_set (lam[k][posj], u);
   }

   InnerProduct (u, B.coeff[k], B.coeff[k], n);

   for (i = 1; i <= *s; i++) {
      mpz_mul (t1, D[i], u);
      mpz_submul (t1, lam[k][i], lam[k][i]);
      mpz_div (t1, t1, D[i-1]);
      mpz_set (u, t1);
   }

   if (mpz_cmp_ui (u, 0) == 0)
     {
       P[k] = 0;
     }
   else
     {
       (*s)++;
       P[k] = *s;
       mpz_set (D[*s], u);
     }

   mpz_clear (u);
   mpz_clear (t1);
}


static void
BalDiv (mpz_t q, mpz_t a, mpz_t d, mpz_t r)
/*  rounds a/d to nearest integer, breaking ties
    by rounding towards zero.  Assumes d > 0. */
{
   long cmp;

   mpz_fdiv_qr (q, r, a, d);

   mpz_mul_2exp (r, r, 1);

   cmp = mpz_cmp (r, d);
   if (cmp > 0 || (cmp == 0 && mpz_cmp_ui (q, 0) < 0))
     mpz_add_ui (q, q, 1);
}


static void
MulSubN (mpz_t *c, mpz_t *c2, mpz_t x, long n, mpz_t tmp)
/* c = c - x*c2 */
{
   long i;
   signed long int x0;

   x0 = mpz_get_si (x);
   if (mpz_cmp_si (x, x0) == 0 && 
       x0 != ((signed long int) 1 << (mp_bits_per_limb - 1))) {
     if (x0 > 0)
       for (i = 1; i <= n; i++) {
	 mpz_mul_ui (tmp, c2[i], x0);
	 mpz_sub (c[i], c[i], tmp);
       }
     else if (x0 < 0) {
       x0 = -x0;
       for (i = 1; i <= n; i++)
	 mpz_addmul_ui (c[i], c2[i], x0);
     }
   }
   else
     {
       for (i = 1; i <= n; i++)
         {
           mpz_mul (tmp, c2[i], x);
           mpz_sub (c[i], c[i], tmp);
         }
     }
}


static void
reduce ( long k, long l, mat_Z B, long *P, mpz_t *D, 
         mpz_t **lam, mat_Z* U, mpz_t t1, mpz_t r )
/* t1 and r are temporary variables */
{
    long j;

   if (P[l] == 0) return;

   mpz_mul_2exp (t1, lam[k][P[l]], 1);
   mpz_abs (t1, t1);
   if (mpz_cmp (t1, D[P[l]]) <= 0)
     return;

   BalDiv (r, lam[k][P[l]], D[P[l]], t1);
   MulSubN (B.coeff[k], B.coeff[l], r, B.NumCols, t1);

   if (U)
     MulSubN (U->coeff[k], U->coeff[l], r, B.NumRows, t1);

   for (j = 1; j <= l-1; j++)
     if (P[j] != 0)
       mpz_submul (lam[k][P[j]], lam[l][P[j]], r);

   mpz_submul (lam[k][P[l]], D[P[l]], r);
}


static long
SwapTest (mpz_t d0, mpz_t d1, mpz_t d2, mpz_t lam,
          mpz_srcptr a, mpz_srcptr b, mpz_t t1, mpz_t t2)
/* test if a*d1^2 > b*(d0*d2 + lam^2)
   t1 and t2 are temporary variables */
{
    mpz_mul (t1, d0, d2);
    mpz_mul (t2, lam, lam);
    mpz_add (t1, t1, t2);
    mpz_mul (t1, t1, b);

    mpz_mul (t2, d1, d1);
    mpz_mul (t2, t2, a);

    return (mpz_cmp (t2, t1) > 0);
}


static void
MulAddDiv (mpz_t c, mpz_t c1, mpz_t c2, 
           mpz_t x, mpz_t y, mpz_t z, mpz_t t1, mpz_t t2)

/* c = (x*c1 + y*c2)/z
   warning: c and z can be the same variable
   t1 and t2 are temporary variables */
{
    mpz_mul(t1, x, c1);
    mpz_mul(t2, y, c2);
    mpz_add(t1, t1, t2);
    mpz_divexact(c, t1, z);
}


static void
MulSubDiv (mpz_t c, mpz_t c1, mpz_t c2, 
           mpz_t x, mpz_t y, mpz_t z, mpz_t t1)
/* c = (x*c1 - y*c2)/z
   t1 is a temporary variable */
{
    mpz_mul(t1, x, c1);
    mpz_mul(c, y, c2);
    mpz_sub(t1, t1, c);
    mpz_divexact(c, t1, z);
}


static void
RowTransform ( mpz_t c1,
               mpz_t c2,
               mpz_t x,
               mpz_t y,
               mpz_t u,
               mpz_t v )
/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
    mpz_t t1, t2;

    mpz_init(t1);
    mpz_init(t2);

    mpz_mul(t1, x, c1);
    mpz_mul(t2, y, c2);
    mpz_add(t1, t1, t2);

    mpz_mul(t2, u, c1);
    mpz_set(c1, t1);
    mpz_mul(t1, v, c2);
    mpz_add(c2, t1, t2);

    mpz_clear(t1);
    mpz_clear(t2);
}


static void
RowTransformN ( mpz_t *c1,
                mpz_t *c2,
                mpz_t x,
                mpz_t y,
                mpz_t u,
                mpz_t v,
                long n )
/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
   mpz_t t1, t2;
   long i;

   mpz_init(t1);
   mpz_init(t2);

   for (i = 1; i <= n; i++)
     {
       mpz_mul(t1, x, c1[i]);
       mpz_mul(t2, y, c2[i]);
       mpz_add(t1, t1, t2);

       mpz_mul(t2, u, c1[i]);
       mpz_set(c1[i], t1);
       mpz_mul(t1, v, c2[i]);
       mpz_add(c2[i], t1, t2);
     }

   mpz_clear(t1);
   mpz_clear(t2);
}


static void
swapLLL (long k, mat_Z B, long *P, mpz_t *D,
         mpz_t **lam, mat_Z* U, long m)
/* swaps vectors k-1 and k;  assumes P(k-1) != 0 */
{
   long i, j;
   mpz_t t1, t2, t3, e, x, y;

   mpz_init(t1);
   mpz_init(t2);
   mpz_init(t3);
   mpz_init(e);
   mpz_init(x);
   mpz_init(y);

   if (P[k] != 0) {
      mpz_swap_n(B.coeff[k-1], B.coeff[k]);
      if (U)
        mpz_swap_n(U->coeff[k-1], U->coeff[k]);
      
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            mpz_swap(lam[k-1][P[j]], lam[k][P[j]]);

      for (i = k+1; i <= m; i++) {
         MulAddDiv(t1, lam[i][P[k]-1], lam[i][P[k]],
                   lam[k][P[k]-1], D[P[k]-2], D[P[k]-1], t2, t3);
         MulSubDiv(lam[i][P[k]], lam[i][P[k]-1], lam[i][P[k]], 
                   D[P[k]], lam[k][P[k]-1], D[P[k]-1], t2);
         mpz_set(lam[i][P[k]-1], t1);
      }

      MulAddDiv(D[P[k]-1], D[P[k]], lam[k][P[k]-1],
                D[P[k]-2], lam[k][P[k]-1], D[P[k]-1], t2, t3);
   }
   else if (mpz_cmp_ui (lam[k][P[k-1]], 0) != 0) {
      mpz_gcdext(e, x, y, lam[k][P[k-1]], D[P[k-1]]);

      mpz_divexact(t1, lam[k][P[k-1]], e);
      mpz_divexact(t2, D[P[k-1]], e);

      mpz_set(t3, t2);
      mpz_neg(t2, t2);
      RowTransformN(B.coeff[k-1], B.coeff[k], t1, t2, y, x, B.NumCols);
      if (U)
        RowTransformN(U->coeff[k-1], U->coeff[k], t1, t2, y, x, B.NumCols);
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            RowTransform(lam[k-1][P[j]], lam[k][P[j]], t1, t2, y, x);

      mpz_mul(t2, t2, t2);
      mpz_divexact(D[P[k-1]], D[P[k-1]], t2);

      for (i = k+1; i <= m; i++)
         if (P[i] != 0) {
            mpz_divexact(D[P[i]], D[P[i]], t2);
            for (j = i+1; j <= m; j++) {
               mpz_divexact(lam[j][P[i]], lam[j][P[i]], t2);
            }
         }

      for (i = k+1; i <= m; i++) {
         mpz_divexact(lam[i][P[k-1]], lam[i][P[k-1]], t3);
      }

      swap(P[k-1], P[k]);
   }
   else {
      mpz_swap_n(B.coeff[k-1], B.coeff[k]);
      if (U)
        mpz_swap_n(U->coeff[k-1], U->coeff[k]);
   
      for (j = 1; j <= k-2; j++)
         if (P[j] != 0)
            mpz_swap(lam[k-1][P[j]], lam[k][P[j]]);

      swap(P[k-1], P[k]);
   }

   mpz_clear(t1);
   mpz_clear(t2);
   mpz_clear(t3);
   mpz_clear(e);
   mpz_clear(x);
   mpz_clear(y);
}

/* LLL-reduce the matrix B (whose rows represent vectors, with indices
   starting at 1):
 * det (output) is the square of the determinant
 * U (output) is the transformation matrix (NULL if not needed)
 * a, b are parameters (delta = a/b = 3/4 classically, we must have
   1/4 < delta < 1, the closer delta is from 1, the better is the reduction)
 m is the number of vectors (i.e., number of rows)
 n is the number of columns (i.e., length of each vector)
 */
long LLL ( mpz_t det, mat_Z B, mat_Z* U, mpz_srcptr a, mpz_srcptr b )
{
   long m, n, *P, j, s, k, max_k;
   mpz_t *D, **lam, tmp1, tmp2;

   mpz_init (tmp1);
   mpz_init (tmp2);

   m = B.NumRows;
   n = B.NumCols;
   ASSERT_ALWAYS(n >= m);

   P = (long*) malloc((m+1) * sizeof(long));
   memset(P, 0, (m+1) * sizeof(long));

   D = (mpz_t*) malloc((m+1) * sizeof(mpz_t));
   for (j=0; j<=m; j++)
     mpz_init_set_ui(D[j], j==0);

   lam = (mpz_t**) malloc((m+1) * sizeof(mpz_t*));
   for (j = 0; j <= m; j++)
     {
       lam[j] = (mpz_t*) malloc((m + 1) * sizeof(mpz_t));
       for (k = 0; k <= m; k++) mpz_init_set_ui(lam[j][k], 0);
     }

   if (U) ident(*U, m);

   s = 0;

   k = 1;
   max_k = 0;

   while (k <= m) {
      if (k > max_k)
        {
          IncrementalGS (B, P, D, lam, &s, k);
          max_k = k;
        }

      if (k == 1) {
         k++;
         continue;
      }

      reduce (k, k-1, B, P, D, lam, U, tmp1, tmp2);

      if (P[k-1] != 0 && 
          (P[k] == 0 || 
           SwapTest (D[P[k]], D[P[k]-1], D[P[k]-2], lam[k][P[k]-1], a, b, tmp1, tmp2))) {
         swapLLL (k, B, P, D, lam, U, max_k);
         k--;
      }
      else {	
         for (j = k-2; j >= 1; j--) 
            reduce (k, j, B, P, D, lam, U, tmp1, tmp2);
         k++;
      }
   }

   mpz_set (det, D[s]);
   for (j = 0; j <= m; j++)
     mpz_clear (D[j]);
   free (D);
   for (j = 0; j <= m; j++)
     {
       for (k = 0; k <= m; k++)
         mpz_clear (lam[j][k]);
       free (lam[j]);
     }
   free (lam);

   mpz_clear(tmp1);
   mpz_clear(tmp2);

   free(P);
   return s;
}


