#include "cado.h" // IWYU pragma: keep

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <iostream>
#include <sstream>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_poly_parallel.hpp"
#include "mpz_polymodF.h"
#include "portability.h" //  IWYU pragma: keep
#include "tests_common.h"

static void mpz_poly_setcoeffs_ui_var(mpz_poly f, int d, ...)
{
    va_list ap;
    va_start(ap, d);
    mpz_poly_realloc(f, d + 1);
    for(int i = 0 ; i <= d ; i++) {
        mpz_poly_setcoeff_ui(f, i, va_arg(ap, int));
    }
    mpz_poly_cleandeg(f, d);
    va_end(ap);
}

static void mpz_poly_setcoeffs_si_var(mpz_poly f, int d, ...)
{
    va_list ap;
    va_start(ap, d);
    mpz_poly_realloc(f, d + 1);
    for(int i = 0 ; i <= d ; i++) {
        mpz_poly_setcoeff_si(f, i, va_arg(ap, int));
    }
    mpz_poly_cleandeg(f, d);
    va_end(ap);
}

/* Return f=g*h, where g has degree r, and h has degree s. */
static int
mpz_poly_mul_basecase (mpz_t *f, mpz_t *g, int r, mpz_t *h, int s) {
  int i, j;
  ASSERT(f != g && f != h);
  for (i = 0; i <= r + s; i++)
    mpz_set_ui (f[i], 0);
  for (i = 0; i <= r; ++i)
    for (j = 0; j <= s; ++j)
      mpz_addmul(f[i+j], g[i], h[j]);
  return r + s;
}

/* check mpz_poly_mul_tc against basecase code */
static void
test_mpz_poly_mul_tc (unsigned long iter)
{
  mpz_poly g, h, f0, f1;
  int r, s;

  mpz_poly_init (g, -1);
  mpz_poly_init (h, -1);
  mpz_poly_init (f0, -1);
  mpz_poly_init (f1, -1);

  while (iter--)
    for (r = 0; r <= MAX_TC_DEGREE; r++)
      for (s = 0; r + s <= MAX_TC_DEGREE; s++)
        {
          mpz_poly_set_randomb (g, r, state, 10,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
          mpz_poly_set_randomb (h, s, state, 10,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
          mpz_poly_mul (f0, g, h);
          mpz_poly_realloc (f1, r + s + 1);
          f1->deg = r + s;
          mpz_poly_mul_basecase (f1->_coeff, g->_coeff, r, h->_coeff, s);
          if (mpz_poly_cmp (f0, f1) != 0)
            {
              printf ("Error, mpz_poly_mul and mpz_poly_mul_basecase differ\n");
              printf ("g="); mpz_poly_fprintf (stdout, g);
              printf ("h="); mpz_poly_fprintf (stdout, h);
              printf ("mpz_poly_mul gives ");
              mpz_poly_fprintf (stdout, f0);
              printf ("mpz_poly_mul_basecase gives ");
              mpz_poly_fprintf (stdout, f1);
              abort ();
            }
        }

  mpz_poly_clear (f0);
  mpz_poly_clear (f1);
  mpz_poly_clear (g);
  mpz_poly_clear (h);
}

/* check mpz_poly_sqr_tc against basecase code */
static void
test_mpz_poly_sqr_tc (unsigned long iter)
{
  mpz_poly g, f0, f1;
  int r;

  mpz_poly_init (g, -1);
  mpz_poly_init (f0, -1);
  mpz_poly_init (f1, -1);

  while (iter--)
    for (r = 0; r <= MAX_TC_DEGREE; r++)
      {
        mpz_poly_set_randomb (g, r, state, 10,
                MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
        mpz_poly_mul(f0, g, g);
        mpz_poly_realloc (f1, r + r + 1);
        f1->deg = r + r;
        mpz_poly_mul_basecase (f1->_coeff, g->_coeff, r, g->_coeff, r);
        if (mpz_poly_cmp (f0, f1) != 0)
          {
            printf ("Error, mpz_poly_sqr and mpz_poly_sqr_basecase differ\n");
            printf ("g="); mpz_poly_fprintf (stdout, g);
            printf ("mpz_poly_mul gives ");
            mpz_poly_fprintf (stdout, f0);
            printf ("mpz_poly_mul_basecase gives ");
            mpz_poly_fprintf (stdout, f1);
            abort ();
          }
      }

  mpz_poly_clear (f0);
  mpz_poly_clear (f1);
  mpz_poly_clear (g);
}

static void
test_mpz_polymodF_mul ()
{
  int d1, d2, d;
  mpz_poly F, T, U;
  mpz_polymodF P1, P2, Q, P1_saved;
  int k = 2 + gmp_urandomm_ui(state, 127), count = 0;
  mpz_t c;

  mpz_poly_init (T, -1);
  mpz_poly_init (U, -1);
  mpz_init (c);
  for (d = 1; d <= 10; d++)
    {
      mpz_poly_init (F, d);
      mpz_poly_init (Q->p, d-1);
      do {
          mpz_poly_set_randomb (F, d, state, k,
                MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_UPPER_BOUND);
      } while (F->deg == -1);
      for (d1 = 1; d1 <= 10; d1++)
        {
          mpz_poly_init (P1->p, d1);
          mpz_poly_init (P1_saved->p, d1);
          mpz_poly_set_randomb (P1->p, d1, state, k,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
          mpz_poly_set (P1_saved->p, P1->p);
          P1->v = 0;
          for (d2 = 1; d2 <= 10; d2++)
            {
              mpz_poly_init (P2->p, d2);
              mpz_poly_set_randomb (P2->p, d2, state, k,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
              P2->v = 0;
              if ((++count % 3) == 0)
                mpz_polymodF_mul (Q, P1, P2, F);
              else if ((count % 3) == 1)
                {
                  mpz_polymodF_mul (P1, P1, P2, F);
                  mpz_poly_set (Q->p, P1->p);
                  Q->v = P1->v;
                  mpz_poly_set (P1->p, P1_saved->p);
                  P1->v = 0;
                }
              else
                {
                  mpz_polymodF_mul (P1, P2, P1, F);
                  mpz_poly_set (Q->p, P1->p);
                  Q->v = P1->v;
                  mpz_poly_set (P1->p, P1_saved->p);
                  P1->v = 0;
                }
              /* check that Q->p = lc(F)^Q->v * P1 * P1 mod F */
              ASSERT_ALWAYS (Q->p->deg < F->deg);
              mpz_poly_mul (T, P1->p, P2->p);
              mpz_pow_ui (c, mpz_poly_lc(F), Q->v);
              mpz_poly_mul_mpz (T, T, c);
              mpz_poly_sub (T, T, Q->p);
              /* T should be a multiple of F */
              while (T->deg >= F->deg)
                {
                  int const oldd = T->deg;
                  if (!mpz_divisible_p (mpz_poly_lc(T), mpz_poly_lc(F)))
                    {
                      printf ("Error in test_mpz_polymodF_mul\n");
                      printf ("F="); mpz_poly_fprintf (stdout, F);
                      printf ("P1="); mpz_poly_fprintf (stdout, P1->p);
                      printf ("P2="); mpz_poly_fprintf (stdout, P2->p);
                      printf ("Q="); mpz_poly_fprintf (stdout, Q->p);
                      exit (1);
                    }
                  mpz_divexact (c, mpz_poly_lc(T), mpz_poly_lc(F));
                  mpz_poly_mul_mpz (U, F, c);
                  /* multiply U by x^(T->deg - F->deg) */
                  mpz_poly_mul_xi (U, U, T->deg - F->deg);
                  mpz_poly_sub (T, T, U);
                  ASSERT_ALWAYS (T->deg < oldd);
                }
              if (T->deg != -1)
                {
                  printf ("count=%d\n", count);
                }
              ASSERT_ALWAYS (T->deg == -1);
              mpz_poly_clear (P2->p);
            }
          mpz_poly_clear (P1->p);
          mpz_poly_clear (P1_saved->p);
        }
      mpz_poly_clear (F);
      mpz_poly_clear (Q->p);
    }
  mpz_poly_clear (T);
  mpz_poly_clear (U);
  mpz_clear (c);
}

#if 0
void 
test_mpz_poly_roots_mpz (unsigned long iter) 
{ 
  mpz_t r[10], f[10], p, res; 
  unsigned long i, n, d; 
  mpz_poly F; 

  for (i = 0; i < 10; i++) 
    { 
      mpz_init (r[i]); 
      mpz_init (f[i]); 
    } 
  mpz_init (p); 
  mpz_init (res); 

  /* -16*x^2 - x - 2 mod 17 */ 
  mpz_set_si (f[2], -16); 
  mpz_set_si (f[1], -1); 
  mpz_set_si (f[0], -2); 
  F->coeff = f; 
  F->deg = 2; 
  mpz_set_ui (p, 17); 
  n = mpz_poly_roots_mpz (r, F, p); 
  ASSERT_ALWAYS(n == 2); 
  ASSERT_ALWAYS(mpz_cmp_ui (r[0], 2) == 0); 
  ASSERT_ALWAYS(mpz_cmp_ui (r[1], 16) == 0); 

  /* 9*x^2 + 6*x + 3 mod 3 */ 
  mpz_set_si (f[2], 9); 
  mpz_set_si (f[1], 6); 
  mpz_set_si (f[0], 3); 
  mpz_set_ui (p, 3); 
  n = mpz_poly_roots_mpz (r, F, p, state); 
  ASSERT_ALWAYS(n == 0); 

  /* try random polynomials */ 
  for (i = 0; i < iter; i++) 
    { 
      d = 1 + gmp_urandomm_ui(state, 7); 
      for (n = 0; n <= d; n++) 
        mpz_set_si (f[n], INT_MIN + gmp_urandomb_ui(state, 32)); 
      mpz_urandomb (p, state, 128); 
      mpz_nextprime (p, p); 
      while (mpz_divisible_p (f[d], p)) 
        mpz_set_si (f[d], INT_MIN + gmp_urandomb_ui(state, 32)); 
      F->coeff = f; 
      F->deg = d; 
      n = mpz_poly_roots_mpz (r, F, p, state); 
      ASSERT_ALWAYS (n <= d); 
      while (n-- > 0) 
        { 
          mpz_poly_eval (res, F, r[n]); 
          ASSERT_ALWAYS (mpz_divisible_p (res, p)); 
        } 
    } 

  for (i = 0; i < 10; i++) 
    { 
      mpz_clear (r[i]); 
      mpz_clear (f[i]); 
    } 
  mpz_clear (p); 
  mpz_clear (res); 
} 
#endif

/* also exercises mpz_poly_mul */
static void
test_mpz_poly_sqr_mod_f_mod_mpz (unsigned long iter)
{
  while (iter--)
    {
      mpz_poly Q, P, f;
      mpz_t m, invm;
      int const k = 2 + gmp_urandomm_ui(state, 127);
      int const d = 1 + gmp_urandomm_ui(state, 7);

      mpz_init (m);
      do mpz_urandomb (m, state, k); while (mpz_tstbit (m, 0) == 0);
      mpz_poly_init (f, d);
      mpz_init (invm);
      while (1)
        {
          mpz_poly_set_randomb (f, d, state, k,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
          mpz_gcd (invm, m, mpz_poly_lc(f));
          if (mpz_cmp_ui (invm, 1) == 0)
            break;
        }
      mpz_poly_init (P, d - 1);
      if (iter)
        mpz_poly_set_randomb (P, d - 1, state, k,
                  MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
      else
        P->deg = -1; /* P=0 */
      mpz_poly_init (Q, d - 1);
      mpz_poly_sqr_mod_f_mod_mpz (Q, P, f, m, nullptr, nullptr);
      if (iter == 0)
        ASSERT_ALWAYS(Q->deg == -1);
      mpz_poly_mul_mod_f_mod_mpz (Q, P, P, f, m, nullptr, nullptr);
      if (iter == 0)
        ASSERT_ALWAYS(Q->deg == -1);
      mpz_poly_clear (f);
      mpz_poly_clear (P);
      mpz_poly_clear (Q);
      mpz_clear (m);
      mpz_clear (invm);
    }
}

/* Also exercises mpz_poly_getcoeff, mpz_poly_setcoeff_int64,
   mpz_poly_setcoeff_si, mpz_poly_cmp, mpz_poly_eval,
   mpz_poly_eval_mod_mpz and mpz_poly_eval_several_mod_mpz */
static void
test_mpz_poly_fprintf()
{
  mpz_poly f, g;
  mpz_t c, v[2], m, invm;
  int res;
  mpz_poly_srcptr F[2];
  mpz_ptr V[2];

  mpz_poly_init (f, 1);
  mpz_poly_init (g, 1);
  F[0] = f;
  F[1] = g;
  V[0] = (mpz_ptr) v[0];
  V[1] = (mpz_ptr) v[1];
  mpz_init (c);
  mpz_init (v[0]);
  mpz_init (v[1]);
  mpz_init_set_ui (m, 11);
  mpz_init (invm);

  f->deg = -1;
  mpz_poly_fprintf (stdout, f);
  mpz_set(c, mpz_poly_coeff_const(f, 0));
  ASSERT_ALWAYS (mpz_cmp_ui (c, 0) == 0);
  mpz_set_ui (c, 17);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 0) == 0);
  mpz_poly_eval_mod_mpz (v[0], f, c, m);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 0) == 0);

  f->deg = 0;
  mpz_poly_setcoeff_ui(f, 0, 17); /* f = 17 */
  mpz_poly_fprintf (stdout, f);
  mpz_set_ui (c, 42);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 17) == 0);
  mpz_poly_eval_mod_mpz (v[0], f, c, m);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 6) == 0);

  mpz_poly_setcoeff_int64 (f, 1, 42); /* f = 42*x+17 */
  mpz_poly_fprintf (stdout, f);
  mpz_set_ui (c, 1);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 59) == 0);
  mpz_set_si (c, -1);
  mpz_poly_eval (v[0], f, c);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], -25) == 0);
  mpz_poly_eval_mod_mpz (v[0], f, c, m);
  ASSERT_ALWAYS (mpz_cmp_ui (v[0], 8) == 0);

  mpz_poly_setcoeff_si (f, 2, -3); /* f = -3*x^2+42*x+17 */
  mpz_poly_fprintf (stdout, f);

  mpz_poly_set (g, f);
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res == 0);
  mpz_add_ui (mpz_poly_lc_w(g), mpz_poly_lc(g), 1); /* g = -2*x^2+42*x+17 */
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res != 0);
  mpz_set_si (c, 3);
  mpz_poly_eval_several_mod_mpz (V, F, 1, c, m);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  mpz_poly_eval_several_mod_mpz (V, F, 2, c, m);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 4) == 0);
  mpz_poly_setcoeff_si (g, g->deg + 1, 1); /* g = x^3-2*x^2+42*x+17 */
  res = mpz_poly_cmp (f, g);
  ASSERT_ALWAYS (res != 0);
  mpz_set_si (c, -3);
  mpz_poly_eval_several_mod_mpz (V, F, 1, c, m);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 7) == 0);
  mpz_poly_eval_several_mod_mpz (V, F, 2, c, m);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 7) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 0) == 0);
  /* test with one zero polynomial */
  g->deg = -1;
  mpz_set_si (c, 3);
  mpz_poly_eval_several_mod_mpz (V, F, 2, c, m);
  ASSERT_ALWAYS (mpz_cmp_si (v[0], 6) == 0);
  ASSERT_ALWAYS (mpz_cmp_si (v[1], 0) == 0);

  mpz_poly_clear (f);
  mpz_poly_clear (g);
  mpz_clear (c);
  mpz_clear (v[0]);
  mpz_clear (v[1]);
  mpz_clear (m);
  mpz_clear (invm);
}

static void
test_mpz_poly_div_2_mod_mpz ()
{
  mpz_poly f;
  mpz_t m;

  mpz_init_set_ui (m, 17);
  mpz_poly_init (f, -1);
  mpz_poly_setcoeff_si (f, 0, 1);
  mpz_poly_setcoeff_si (f, 1, -2);
  mpz_poly_setcoeff_si (f, 2, -3);
  mpz_poly_setcoeff_si (f, 3, 4);
  mpz_poly_div_2_mod_mpz (f, f, m);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(f, 0), 9) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(f, 1), -1) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(f, 2), 7) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(f, 3), 2) == 0);
  mpz_poly_clear (f);
  mpz_clear (m);
}

static void
test_mpz_poly_derivative ()
{
  mpz_poly f, df;

  mpz_poly_init (f, -1);
  mpz_poly_init (df, 1);

  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == -1);

  mpz_poly_setcoeff_si (f, 0, 17); /* f = 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == -1);

  mpz_poly_setcoeff_si (f, 1, 42); /* f = 42*x + 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(df, 0), 42) == 0);

  mpz_poly_setcoeff_si (f, 2, -3); /* f = -3*x^2 + 42*x + 17 */
  mpz_poly_derivative (df, f);
  ASSERT_ALWAYS(df->deg == 1);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(df, 0), 42) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(df, 1), -6) == 0);

  mpz_poly_clear (f);
  mpz_poly_clear (df);
}

/* also exercises mpz_poly_pow_mod_f_mod_mpz */
static void
test_mpz_poly_pow_mod_f_mod_ui ()
{
  mpz_poly Q, P, f;
  mpz_t a, pp;
  unsigned long const p = 4294967291UL;

  mpz_poly_init (Q, -1);
  mpz_poly_init (P, -1);
  mpz_poly_init (f, -1);
  mpz_init (a);
  mpz_init_set_ui (pp, p);
  mpz_poly_setcoeff_si (f, 4, 60);
  mpz_poly_setcoeff_si (f, 3, 165063);
  mpz_poly_setcoeff_int64 (f, 2, (int64_t) 2561596016);
  mpz_poly_setcoeff_int64 (f, 1, (int64_t) -4867193837504);
  mpz_poly_setcoeff_int64 (f, 0, (int64_t) -9292909378109715);
  mpz_poly_setcoeff_si (P, 0, 0);
  mpz_poly_setcoeff_si (P, 1, 1); /* P = x */

  mpz_set_ui (a, 0);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 1) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 1) == 0);

  mpz_set_ui (a, 1);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(mpz_poly_cmp (Q, P) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(mpz_poly_cmp (Q, P) == 0);

  mpz_set_ui (a, 2);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 1) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 1) == 0);

  mpz_set_ui (a, 3);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 1) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 1) == 0);

  mpz_set_ui (a, 4);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 2081229567) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 3524154901) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 1102631344) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 2362229259) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 2081229567) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 3524154901) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 1102631344) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 2362229259) == 0);

  mpz_set_ui (a, 999999);
  mpz_poly_pow_mod_f_mod_ui (Q, P, f, a, p);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 4223801964) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 502704799) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 3358125388) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 1722383279) == 0);
  mpz_poly_pow_mod_f_mod_mpz (Q, P, f, a, pp);
  ASSERT_ALWAYS(Q->deg == 3);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 0), 4223801964) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 1), 502704799) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 2), 3358125388) == 0);
  ASSERT_ALWAYS(mpz_cmp_si (mpz_poly_coeff_const(Q, 3), 1722383279) == 0);

  mpz_clear (a);
  mpz_clear (pp);
  mpz_poly_clear (Q);
  mpz_poly_clear (P);
  mpz_poly_clear (f);
}

/* also test mpz_poly_base_modp_lift and mpz_poly_sizeinbase */
static void
test_mpz_poly_base_modp_init (unsigned long iter)
{
  int p, l, d, i;
  unsigned long k, K[65];
  mpz_poly f, *P, g;
  mpz_t pk;
  size_t s;
  mpz_poly_parallel_info const pinf;

  mpz_poly_init (f, -1);
  mpz_poly_init (g, -1);
  mpz_init (pk);
  while (iter--)
    {
      p = gmp_urandomb_ui(state, 31);
      if (p < 2)
        p = 2;
      k = gmp_urandomb_ui(state, 10);
      if (k < 2)
        k = 2; /* ensures l > 0 */
      for (K[0] = k, l = 0; K[l] > 1; K[l+1] = (K[l] + 1) >> 1, l++);
      d = gmp_urandomm_ui(state, 10);
      int m = (1 + gmp_urandomm_ui(state, 9)) * k;
      if (iter == 0) /* exercise bug found on 32-bit MinGW */
        {
          p = 2048;
          // k = 119;
          d = 1;
          m = 833;
        }
      mpz_poly_set_randomb (f, d, state, m,
              MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
      s = mpz_poly_sizeinbase (f, 2);
      for (i = 0; i <= f->deg; i++)
        ASSERT_ALWAYS(mpz_sizeinbase (mpz_poly_coeff_const(f, i), 2) <= s);
      P = pinf.mpz_poly_base_modp_init(f, p, K, l);
      /* check f = P[0] + p^K[l]*P[1] + p^K[l-1]*P[2] + ... + p^K[1]*P[l] */
      for (i = 1; i <= l; i++)
        {
          mpz_ui_pow_ui (pk, p, K[l + 1 - i]);
          pinf.mpz_poly_base_modp_lift (P[0], P, i, pk);
        }
      ASSERT_ALWAYS(mpz_poly_cmp (f, P[0]) == 0);
      mpz_poly_base_modp_clear (P, l);
    }
  mpz_poly_clear (f);
  mpz_poly_clear (g);
  mpz_clear (pk);
}

static void test_mpz_poly_is_root(unsigned long iter)
{
    mpz_t p, r;
    mpz_poly f, ell;
    mpz_init(p);
    mpz_init(r);
    mpz_poly_init(f, 10);
    mpz_poly_init(ell, 1);

    for( ; iter--; ) {
        mpz_poly_set_randomb(f, 10, state, 100,
              MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
        mpz_urandomb(p, state, 100);
        mpz_rrandomb(r, state, 100);
        mpz_poly_setcoeff_si(ell, 1, 1);
        mpz_neg(mpz_poly_coeff(ell, 0), r);
        mpz_poly_mul(f, f, ell);
        mpz_poly_mod_mpz(f, f, p, nullptr);
        mpz_mod(r, r, p);
    }

    ASSERT_ALWAYS(mpz_poly_is_root(f, r, p));
    mpz_poly_clear(f);
    mpz_poly_clear(ell);
    mpz_clear(r);
    mpz_clear(p);
}

static void test_mpz_poly_factor(unsigned long iter)
{
    mpz_t p;
    mpz_poly_factor_list lf;
    mpz_poly f;

    mpz_init(p);
    mpz_poly_init(f, -1);
    mpz_poly_factor_list_init(lf);

    
    mpz_set_ui(p, 2);

    mpz_poly_set_from_expression(f, "x+2*x^2+2*x^3+2*x^4");
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 1);
    mpz_poly_mod_mpz(f, f, p, nullptr);
    ASSERT_ALWAYS(mpz_poly_cmp(lf->factors[0]->f, f) == 0);

    mpz_poly_set_from_expression(f, "1+2*x^2+2*x^3+2*x^4");
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 0);

    mpz_poly_set_from_expression(f, "x-2*x^2-2*x^3-2*x^4");
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 1);
    mpz_poly_mod_mpz(f, f, p, nullptr);
    ASSERT_ALWAYS(mpz_poly_cmp(lf->factors[0]->f, f) == 0);


    mpz_poly_set_from_expression(f, "x^10 + x^9 + x^8 + 1");
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 3);

    mpz_poly_set_from_expression(f, "x^10 + x^8 + x^7 + x^6 + x^2 + x + 1");
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 1);

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 21, 3, 8, 5, 5, 6, 1, 9, 4, 3, 3, 8, 7, 7, 7, 0, 12, 5, 11, 11, 1, 7, 10);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 2);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 20);

    ASSERT_ALWAYS(mpz_poly_is_irreducible(lf->factors[1]->f, p));

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 25, 0, 0, 0, 0, 5, 1, 1, 3, 7, 9, 7, 4, 11, 8, 2, 10, 8, 11, 0, 5, 12, 12, 1, 11, 2, 3);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 6);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[0]->m == 4);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 2);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 6);
    ASSERT_ALWAYS(lf->factors[5]->f->deg == 9);

    mpz_set_ui(p, 13);
    mpz_poly_setcoeffs_ui_var(f, 21, 0, 7, 4, 3, 6, 11, 3, 6, 12, 2, 5, 4, 7, 6, 4, 6, 12, 8, 0, 3, 3, 6);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 5);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[2]->m == 2);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 5);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 8);

    mpz_set_ui(p, 5);
    mpz_poly_setcoeffs_ui_var(f, 8, 3, 3, 3, 3, 2, 2, 2, 2, 3);
    printf ("f="); mpz_poly_fprintf (stdout, f);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 2);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 4);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 4);

    mpz_set_ui(p, 3);
    mpz_poly_setcoeffs_ui_var(f, 20, 1, 1, 1, 2, 1, 0, 0, 1, 0, 0, 0, 1, 2, 0, 2, 0, 2, 0, 1, 1, 2);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 4);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[1]->m == 3);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 6);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 10);

    mpz_set_ui(p, 3);
    mpz_poly_setcoeffs_ui_var(f, 20, 2, 0, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 0, 2, 2, 1, 1, 2, 0, 2, 2);
    mpz_poly_factor(lf, f, p, state);
    ASSERT_ALWAYS(lf->size == 5);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 1);
    ASSERT_ALWAYS(lf->factors[0]->m == 6);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 3);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 4);
    ASSERT_ALWAYS(lf->factors[4]->f->deg == 4);

    mpz_set_ui(p, 2);
    mpz_poly_setcoeffs_ui_var(f, 12, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1);
    mpz_poly_factor_sqf(lf, f, p);
    /* beware, mpz_poly_factor_sqf output is a bit peculiar */
    ASSERT_ALWAYS(lf->size == 4);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 0);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 9);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 0);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 1);

    /* same, with entries not reduced */
    mpz_set_ui(p, 2);
    mpz_poly_setcoeffs_ui_var(f, 12, -122661, -9, -9, 9, -3, 0, 0, 0, 0, 0, 0, 0, 1);

    mpz_poly_factor_sqf(lf, f, p);
    /* beware, mpz_poly_factor_sqf output is a bit peculiar */
    ASSERT_ALWAYS(lf->size == 4);
    ASSERT_ALWAYS(lf->factors[0]->f->deg == 0);
    ASSERT_ALWAYS(lf->factors[1]->f->deg == 9);
    ASSERT_ALWAYS(lf->factors[2]->f->deg == 0);
    ASSERT_ALWAYS(lf->factors[3]->f->deg == 1);

    mpz_set_ui(p, 3);
    mpz_poly_setcoeffs_ui_var(f, 12, -122661, -9, -9, 9, -3, 0, 0, 0, 0, 0, 0, 0, 1);
    mpz_poly_factor(lf, f, p, state);

    {
        /* Now factor x^(p-1)-1 */
        unsigned long const pp = 7;
        mpz_set_ui(p, pp);
        mpz_poly_set_xi(f, pp - 1);
        mpz_poly_sub_ui(f, f, 1);
        mpz_poly_factor(lf, f, p, state);
    }

    for( ; iter-- ; ) {
        // fprintf(stderr, "%lu ", iter);
        mpz_rrandomb(p, state, 20);
        mpz_nextprime(p, p);
        mpz_poly_set_randomb(f, 10, state, 10,
              MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
        mpz_poly_mod_mpz(f, f, p, nullptr);
        // mpz_poly_fprintf(stderr, f);
        mpz_poly_factor(lf, f, p, state);
        mpz_poly g;
        mpz_poly_init(g, f->deg);
        mpz_poly_set_xi(g, 0);
        for(int i = 0 ; i < lf->size ; i++) {
            mpz_poly_with_m_ptr fx = lf->factors[i];
            mpz_poly_pow_ui_mod_f_mod_mpz(fx->f, fx->f, nullptr, fx->m, p); 
            mpz_poly_mul(g, g, fx->f);
            mpz_poly_mod_mpz(g, g, p, nullptr);
        }
        mpz_poly_makemonic_mod_mpz(f, f, p);
        ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);
        mpz_poly_clear(g);
    }

    mpz_poly_factor_list_clear(lf);
    mpz_clear(p);
    mpz_poly_clear(f);
}

static void test_mpz_poly_factor_padic(unsigned long iter)
{
    cxx_mpz_poly f;
    cxx_mpz p = 1;
    cxx_gmp_randstate rstate;

    for(unsigned long i = 0 ; i < iter ; ) {
        /* pick a prime between 4 and 128 bits */
        unsigned long const pbits = gmp_urandomm_ui(rstate, 124) + 4;
        for(p = 1 ; !mpz_probab_prime_p(p, 2) ; mpz_urandomb(p, rstate, pbits));

        cxx_mpz disc = 0;
        for( ; mpz_cmp_ui(disc, 0) == 0 ; ) {
            /* pick a degree between 2 and 10 */
            int const deg = 2 + (int) gmp_urandomm_ui(rstate, 8);
            mpz_poly_set_randomm_ui(f, deg, state, 100,
                    MPZ_POLY_URANDOM | MPZ_POLY_DEGREE_EXACT | MPZ_POLY_MONIC);
            mpz_poly_mod_mpz(f, f, p, nullptr);

            mpz_poly_discriminant(disc, f);
            mpz_mod(disc, disc, p);
        }

        int const prec = MAX(2, 256 / pbits);
        auto lf = mpz_poly_factor_and_lift_padically(f, p, prec, state);

        cxx_mpz px;
        mpz_pow_ui(px, p, prec);
        auto F = prod(lf, px);

        mpz_poly_sub_mod_mpz(F, F, f, px);
        ASSERT_ALWAYS(F.degree() == -1);

        /* only for sagemath testing. the test above is good enough I
         * think.
        fmt::print("R.<x>=pAdicRing({},{})[];\nA={}\nB=[\n",
                p, prec, f.print_poly("x"));
        for(auto const & fm : lf) {
            fmt::print("\t({})^{},\n", fm.first.print_poly("x"), fm.second);
        }
        fmt::print("\t]\nA == prod(B)\n");
        */
        i++;
    }
}

static void test_mpz_poly_trivialities()
{
    mpz_t p;
    mpz_poly f, g, q, r;
    mpz_poly_factor_list lf;
    int rc;

    mpz_init_set_ui(p, 13);
    mpz_poly_init(f, -1);
    mpz_poly_init(g, -1);
    mpz_poly_init(q, -1);
    mpz_poly_init(r, -1);

    {
        mpz_t a[5];
        mpz_init_set_ui(a[0], 1);
        mpz_init_set_ui(a[1], 4);
        mpz_init_set_ui(a[2], 6);
        mpz_init_set_ui(a[3], 4);
        mpz_init_set_ui(a[4], 1);
        mpz_poly_setcoeffs(f, a, 4);
        mpz_clear(a[4]);
        mpz_clear(a[3]);
        mpz_clear(a[2]);
        mpz_clear(a[1]);
        mpz_clear(a[0]);
    }

    {
        mpz_poly_factor_list_init(lf);
        mpz_poly_factor(lf, f, p, state);
        ASSERT_ALWAYS(lf->factors[0]->m == 4);
        mpz_poly_swap(f, lf->factors[0]->f);
        mpz_poly_factor_list_clear(lf);
    }

    /* we expect to have f == x + 1 */
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 1);
    ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_coeff_const(f, 0), 1) == 0);
    ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_coeff_const(f, 1), 1) == 0);
    mpz_poly_sub_mod_mpz(f, f, g, p);
    mpz_poly_set_zero(g);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* check that (x+1)^5 div x is irreducible mod 13 */
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 1);
    mpz_poly_pow_ui_mod_f_mod_mpz(g, g, nullptr, 5, p);
    mpz_poly_div_xi(f, g, 6);
    ASSERT_ALWAYS(f->deg == -1);
    mpz_poly_div_xi(f, g, 0);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);
    mpz_poly_div_xi(g, g, 1);
    ASSERT_ALWAYS(mpz_poly_is_irreducible(g, p));

    /* check that x+1 - (x+1)^2 = -x^2-x */
    mpz_poly_set_zero(f);
    mpz_poly_add_ui(f, f, 1);
    mpz_poly_set_xi(g, 1);
    mpz_poly_add(f, f, g);
    mpz_poly_mul(g, f, f);
    mpz_poly_sub(f, f, g);
    mpz_poly_set_zero(g);
    mpz_poly_sub_ui(g, g, 1);
    mpz_poly_mul(f, f, g);
    mpz_poly_set_zero(g);
    mpz_poly_set_xi(g, 2);
    mpz_poly_setcoeff_ui(g, 1, 1);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* multiply by zero */
    mpz_poly_set_randomb(f, 10, state, 10,
            MPZ_POLY_SIGNED_COEFFICIENTS | MPZ_POLY_RRANDOM | MPZ_POLY_DEGREE_EXACT);
    mpz_poly_set_zero(g);
    mpz_poly_mul(f, f, g);
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    /* div modulo non-prime N should get factor */
    mpz_set_ui(p, 143);
    mpz_poly_setcoeffs_ui_var(f, 2, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 11);
    rc = mpz_poly_div_r_mod_mpz_clobber(f, g, p);
    ASSERT_ALWAYS(rc == 0);
    mpz_poly_setcoeffs_ui_var(f, 2, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 11);
    /* same test. Of course it doesn't exactly divide, but we do expect
     * failure nonetheless */
    rc = mpz_poly_divexact(f, f, g, p);
    ASSERT_ALWAYS(rc == 0);

    /* test div_qr */
    mpz_poly_setcoeffs_ui_var(f, 4, 1, 1, 1, 1, 3);
    mpz_poly_setcoeffs_ui_var(g, 2, 1, 2, 11);
    rc = mpz_poly_div_qr_mod_mpz(q, r, f, g, p);
    ASSERT_ALWAYS(rc == 0);
    mpz_set_ui(p, 13);
    rc = mpz_poly_div_qr_mod_mpz(q, r, f, g, p);
    mpz_poly_setcoeffs_ui_var(f, 2, 0, 11, 5);
    mpz_poly_setcoeffs_ui_var(g, 1, 1, 3);
    ASSERT_ALWAYS(rc);
    ASSERT_ALWAYS(mpz_poly_cmp(f, q) == 0);
    ASSERT_ALWAYS(mpz_poly_cmp(g, r) == 0);

    /* multiply by p, then reduce mod p */
    mpz_poly_set_randomb(f, 10, state, 10,
            MPZ_POLY_SIGNED_COEFFICIENTS |
            MPZ_POLY_RRANDOM |
            MPZ_POLY_DEGREE_EXACT);
    mpz_poly_mul_mpz (f, f, p);
    mpz_poly_makemonic_mod_mpz(f, f, p);
    ASSERT_ALWAYS(f->deg < 0);

    /* a non-irreducible polynomial */
    mpz_poly_setcoeffs_ui_var(g, 2, 1, 2, 1);
    rc = mpz_poly_is_irreducible(g, p);
    ASSERT_ALWAYS(rc == 0);

    mpz_poly_setcoeffs_ui_var(f, 10, 0, 6, 6, 7, 9, 1, 7, 6, 1, 9, 5);
    mpz_poly_set_from_expression(g, "6*x*(x^6+x+1)+7*(x^3+x^6)+(x^4+x^9)*9+x^5+x^8+5*x^10");
    ASSERT_ALWAYS(mpz_poly_cmp(f, g) == 0);

    mpz_poly_clear(f);
    mpz_poly_clear(g);
    mpz_poly_clear(q);
    mpz_poly_clear(r);
    mpz_clear(p);
}

static void test_mpz_poly_resultant()
{
  mpz_poly f, g;
  mpz_poly_init(f, 10);
  mpz_poly_init(g, 10);
  mpz_t res;
  mpz_init(res);
  mpz_t val;
  mpz_init(val);

  /*f=6+7*x^1+0*x^2+9*x^3+13*x^4+13*x^5+1*x^6+4*x^7+8*x^8+4*x^9+6*x^10*/
  /*g=1+10*x^1+7*x^2+2*x^3+9*x^4+5*x^5+0*x^6+10*x^7+7*x^8+5*x^9+4*x^10*/
  mpz_poly_setcoeffs_ui_var(f, 10, 6, 7, 0, 9, 13, 13, 1, 4, 8, 4, 6);
  mpz_poly_setcoeffs_ui_var(g, 10, 1, 10, 7, 2, 9, 5, 0, 10, 7, 5, 4);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "3787840596130306882", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=12+11*x^1+6*x^2+9*x^3+11*x^4+13*x^5+2*x^6+14*x^7+14*x^8+1*x^9+1*x^10*/
  /*g=0+10*x^1+13*x^2+5*x^3+4*x^4+1*x^5+1*x^6+9*x^7+6*x^8+5*x^9+13*x^10*/
  mpz_poly_setcoeffs_ui_var(f, 10, 12, 11, 6, 9, 11, 13, 2, 14, 14, 1, 1);
  mpz_poly_setcoeffs_ui_var(g, 10, 0, 10, 13, 5, 4, 1, 1, 9, 6, 5, 13);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "52543088043796652195928", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=0+6*x^1+6*x^2+7*x^3+9*x^4+1*x^5+7*x^6+6*x^7+1*x^8+9*x^9+5*x^10*/
  /*g=5+7*x^1+11*x^2+0*x^3+13*x^4+9*x^5+5*x^6+0*x^7+4*x^8+1*x^9+5*x^10*/
  mpz_poly_setcoeffs_ui_var(f, 10, 0, 6, 6, 7, 9, 1, 7, 6, 1, 9, 5);
  mpz_poly_setcoeffs_ui_var(g, 10, 5, 7, 11, 0, 13, 9, 5, 0, 4, 1, 5);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "137271514893118787175", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=-10-11*x^1+13*x^2+6*x^3-13*x^4+3*x^5+5*x^6-13*x^7+11*x^8-6*x^9-11*x^10*/
  /*g=-8-14*x^1-9*x^2+2*x^3-4*x^4+0*x^5-7*x^6-10*x^7-3*x^8-3*x^9-11*x^10*/
  mpz_poly_setcoeffs_si_var(f, 10, -10, -11, 13, 6, -13, 3, 5, -13, 11, -6,
      -11);
  mpz_poly_setcoeffs_si_var
    (g, 10, -8, -14, -9, 2, -4, 0, -7, -10, -3, -3, -11);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "408310242047874808370080", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=-3-10*x^1-12*x^2+1*x^3-3*x^4+5*x^5+0*x^6+0*x^7+10*x^8-12*x^9+14*x^10*/
  /*g=13-11*x^1-8*x^2-13*x^3-14*x^4-9*x^5+10*x^6-5*x^7-3*x^8-11*x^9+12*x^10*/
  mpz_poly_setcoeffs_si_var(f, 10, -3, -10, -12, 1, -3, 5, 0, 0, 10, -12, 14);
  mpz_poly_setcoeffs_si_var(g, 10, 13, -11, -8, -13, -14, -9, 10, -5, -3, -11,
      12);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "-36491329842163368368577782", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=-3-15*x^1-9*x^2+3*x^3-13*x^4-12*x^5-1*x^6-1*x^7-12*x^8-14*x^9+4*x^10*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7-5*x^8-4*x^9*/
  mpz_poly_setcoeffs_si_var(f, 10, -3, -15, -9, 3, -13, -12, -1, -1, -12, -14,
      4);
  mpz_poly_cleandeg(g, 9);
  mpz_poly_setcoeffs_si_var(g, 9, -6, -13, 9, 7, -5, -5, 11, 2, -5, -4);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "3719519175576976543932", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=0*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7-5*x^8-4*x^9*/
  mpz_poly_cleandeg(f, -1);
  mpz_poly_resultant(res, f, g);
  ASSERT_ALWAYS(mpz_cmp_ui(res, 0) == 0);

  /*f=-3-15*x^1-9*x^2+3*x^3-13*x^4-12*x^5-1*x^6-1*x^7-12*x^8-14*x^9+4*x^10*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7*/
  mpz_poly_setcoeffs_si_var(f, 10, -3, -15, -9, 3, -13, -12, -1, -1, -12, -14,
      4);
  mpz_poly_cleandeg(g, 7);
  mpz_poly_setcoeffs_si_var(g, 7, -6, -13, 9, 7, -5, -5, 11, 2);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "-26778351555137831424", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=-3-15*x^1-9*x^2+3*x^3-12*x^4-12*x^5-3*x^6-3*x^7-12*x^8-15*x^9+6*x^10*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7*/
  mpz_poly_setcoeffs_si_var(f, 10, -3, -15, -9, 3, -12, -12, -3, -3, -12, -15,
      6);
  mpz_poly_cleandeg(g, 7);
  mpz_poly_setcoeffs_si_var(g, 7, -6, -13, 9, 7, -5, -5, 11, 2);
  mpz_poly_resultant(res, f, g);
  mpz_set_str(val, "-61519394185549843500", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  mpz_poly_resultant(res, g, f);
  mpz_set_str(val, "-61519394185549843500", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);

  /*f=-3-15*x^1-9*x^2+3*x^3-12*x^4-12*x^5-3*x^6-3*x^7-12*x^8-15*x^9*/
  /*g=-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7*/
  mpz_poly_cleandeg(f, 9);
  mpz_poly_setcoeffs_si_var(f, 9, -3, -15, -9, 3, -12, -12, -3, -3, -12, -15);
  mpz_poly_resultant(res, g, f);
  mpz_set_str(val, "-4936496053264331049", 10);
  ASSERT_ALWAYS(mpz_cmp(res, val) == 0);


  mpz_clear(res);
  mpz_clear(val);
  mpz_poly_clear(f);
  mpz_poly_clear(g);
}

static void
test_mpz_poly_discriminant (unsigned long iter)
{
    mpz_poly f;
    mpz_t D;
    mpz_poly_init(f, -1);
    mpz_init(D);

    mpz_poly_setcoeffs_si_var(f, 3, 4, 3, 0, 1);
    mpz_poly_discriminant(D, f);
    ASSERT_ALWAYS(mpz_cmp_ui (D, 540) == 0);

    mpz_poly_setcoeffs_si_var(f, 3, 0, 0, 0, 2);
    mpz_poly_discriminant(D, f);
    ASSERT_ALWAYS(mpz_cmp_ui (D, 0) == 0);

    while (iter--)
    {
        const int N = 10;
        int d = 1 + (int) gmp_urandomm_ui(state, N-1);
        mpz_poly_set_randomm_ui(f, d, state, 2,
            MPZ_POLY_URANDOM | MPZ_POLY_DEGREE_EXACT);

        mpz_poly_discriminant (D, f);
    }

    mpz_poly_clear(f);
    mpz_clear (D);
}

static void test_mpz_poly_discriminant2(unsigned long iter)
{
    cxx_mpz_poly h, fk;
    unsigned long k;
    cxx_mpz D, E;

    cxx_mpz_poly f = "-37200*x^4-301641500*x^3+11679049396284*x^2-7023696347014750305*x+30546672287719745916994";
    cxx_mpz_poly g = "8498629835017307*x-2606392756281442341909";
    mpz_poly_discriminant_of_linear_combination(h, f, g);

    for(unsigned long i = 0 ; i < iter ; i++) {
        k = gmp_urandomm_ui(state, 1<<16);
        mpz_poly_rotation_ui(fk, f, g, k, 0);

        mpz_poly_discriminant(D, fk);
        mpz_poly_eval_ui(E, h, k);

        ASSERT_ALWAYS(mpz_cmp(D, E) == 0);
    }
}


static void test_mpz_poly_infinity_norm()
{
  mpz_poly f;
  mpz_poly_init(f, -1);

  mpz_t inf_norm;
  mpz_init(inf_norm);

  mpz_poly_infinity_norm(inf_norm, f);
  ASSERT_ALWAYS(mpz_cmp_ui(inf_norm, 0) == 0);

  mpz_poly_setcoeffs_si_var(f, 3, 4, 3, 0, 1);
  mpz_poly_infinity_norm(inf_norm, f);
  ASSERT_ALWAYS(mpz_cmp_ui(inf_norm, 4) == 0);

  mpz_poly_setcoeffs_si_var(f, 5, -4, 3, -5, 1, 6, -8);
  mpz_poly_infinity_norm(inf_norm, f);
  ASSERT_ALWAYS(mpz_cmp_ui(inf_norm, 8) == 0);

  mpz_poly_setcoeffs_si_var(f, 2, -4, 3, -5);
  mpz_poly_infinity_norm(inf_norm, f);
  ASSERT_ALWAYS(mpz_cmp_ui(inf_norm, 5) == 0);

  mpz_clear(inf_norm);
  mpz_poly_clear(f);
}

static void test_mpz_poly_interpolation(unsigned long iter)
{
    {
        std::vector<cxx_mpz> points;
        std::vector<cxx_mpz> evaluations;
        cxx_mpz_poly f = "x+x^2-x^17+3";
        for(int i = 0 ; i < 24 ; i++) { 
            cxx_mpz t(i), e;
            mpz_poly_eval(e, f, t);
            points.push_back(t);
            evaluations.push_back(e);
        }
        cxx_mpz_poly tmp;
        int ok = mpz_poly_interpolate(tmp, points, evaluations);
        ASSERT_ALWAYS(tmp == f);
        ASSERT_ALWAYS(ok);
    }

    for(unsigned long i = 0 ; i < iter ; i++) {
        const int d = 4 + (int) gmp_urandomm_ui(state, 97);
        cxx_mpz_poly f;
        mpz_poly_set_randomb(f, d, state, 5,
                MPZ_POLY_UNSIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_DEGREE_EXACT);

        std::vector<cxx_mpz> points;
        std::vector<cxx_mpz> evaluations;

        const int neval = d + 1;

        for(int i = 0 ; i < neval ; i++) {
            cxx_mpz t(i), e;
            mpz_poly_eval(e, f, t);
            points.push_back(t);
            evaluations.push_back(e);
        }
        cxx_mpz_poly tmp;
        const int ok = mpz_poly_interpolate(tmp, points, evaluations);
        if (!ok || (tmp != f)) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << f << "\n";
            std::cerr << tmp << "\n";
        }
        ASSERT_ALWAYS(tmp == f);
        ASSERT_ALWAYS(ok);
    }
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    unsigned long iter = 500;
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    test_mpz_poly_mul_tc (iter / 5);
    test_mpz_poly_sqr_tc (iter / 5);
    test_mpz_polymodF_mul ();
    /* test_mpz_poly_roots_mpz (iter); */
    test_mpz_poly_sqr_mod_f_mod_mpz (iter);
    test_mpz_poly_fprintf();
    test_mpz_poly_div_2_mod_mpz ();
    test_mpz_poly_derivative ();
    test_mpz_poly_pow_mod_f_mod_ui ();
    test_mpz_poly_base_modp_init (iter / 25);
    test_mpz_poly_is_root(iter);
    test_mpz_poly_factor(2 + iter / 5);
    test_mpz_poly_factor_padic(2 + iter / 5);
    test_mpz_poly_trivialities ();
    test_mpz_poly_resultant();
    test_mpz_poly_discriminant(iter);
    test_mpz_poly_discriminant2(10);
    test_mpz_poly_infinity_norm();
    test_mpz_poly_interpolation(iter / 10);
    tests_common_clear ();
    exit (EXIT_SUCCESS);
}
