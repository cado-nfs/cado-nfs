#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <limits.h>
#include "mpz_poly_bivariate.h"
#include "mpz_poly.h"
#include "macros.h" // ASSERT

// ----- lagrange_poly -----
/*
 * lagrange_poly = p / denom;
 */
typedef struct {
  mpz_poly p;
  mpz_t denom;
} lagrange_poly_struct_t;

typedef lagrange_poly_struct_t lagrange_poly_t[1];
typedef lagrange_poly_struct_t * lagrange_poly_ptr;
typedef const lagrange_poly_struct_t * lagrange_poly_srcptr;

static void lagrange_poly_init(lagrange_poly_ptr p, int deg)
{
  mpz_poly_init(p->p, deg);
  mpz_init(p->denom);
  mpz_set_ui(p->denom, 1);
}

static void lagrange_poly_clear(lagrange_poly_ptr p)
{
  mpz_poly_clear(p->p);
  mpz_clear(p->denom);
}

MAYBE_UNUSED static void lagrange_poly_fprintf(FILE * f, lagrange_poly_srcptr p)
{
  fprintf(f, "(");
  mpz_poly_fprintf(f, p->p);
  gmp_fprintf(f, ") / %Zd\n", p->denom);
}

// r = p + q
static void lagrange_poly_add(lagrange_poly_ptr r, lagrange_poly_srcptr p,
    lagrange_poly_srcptr q)
{
  mpz_t lcm;
  mpz_init(lcm);
  mpz_lcm(lcm, p->denom, q->denom);
  mpz_t lcm_red;
  mpz_init(lcm_red);
  // r = p * lcm_red
  mpz_divexact(lcm_red, lcm, p->denom);
  mpz_poly_mul_mpz(r->p, p->p, lcm_red);
  // r = r + q * lcm_red
  mpz_divexact(lcm_red, lcm, q->denom);
  mpz_poly tmp;
  mpz_poly_init(tmp, q->p->deg);
  mpz_poly_mul_mpz(tmp, q->p, lcm_red);
  mpz_poly_add(r->p, r->p, tmp);
  mpz_set(r->denom, lcm);
  mpz_poly_clear(tmp);
  mpz_clear(lcm_red);
  mpz_clear(lcm);
}

static inline void lagrange_poly_mul_mpz(lagrange_poly_ptr r,
    lagrange_poly_srcptr p, mpz_t a)
{
  mpz_poly_mul_mpz(r->p, p->p, a);
}

static void interpolate(mpz_poly_ptr p, mpz_t * x, mpz_t * y, unsigned int nb)
{
  lagrange_poly_t * l = (lagrange_poly_t * )
    malloc(sizeof(lagrange_poly_t) * nb);
  mpz_poly poly_tmp;
  mpz_poly_init(poly_tmp, 1);
  mpz_poly_setcoeff_int64(poly_tmp, 1, 1);
  mpz_t Z_tmp;
  mpz_init(Z_tmp);
  for (unsigned int i = 0; i < nb; i++) {
    lagrange_poly_init(l[i], 0);
    mpz_poly_setcoeff_int64(l[i]->p, 0, 1);
    mpz_set_ui(l[i]->denom, 1);
    for (unsigned int j = 0; j < nb; j++) {
      if (i != j) {
        mpz_poly_setcoeff(poly_tmp, 0, x[j]);
        mpz_mul_si(poly_tmp->coeff[0], poly_tmp->coeff[0], -1);
        mpz_poly_mul(l[i]->p, l[i]->p, poly_tmp);
        mpz_sub(Z_tmp, x[i], x[j]);
        mpz_mul(l[i]->denom, l[i]->denom, Z_tmp);
      }
    }
  }
  mpz_clear(Z_tmp);
  mpz_poly_clear(poly_tmp);
  lagrange_poly_t interpol;
  lagrange_poly_init(interpol, (int) nb);
  for (unsigned int i = 0; i < nb; i++) {
    lagrange_poly_mul_mpz(l[i], l[i], y[i]);
    lagrange_poly_add(interpol, interpol, l[i]);
  }
  for (int i = 0; i <= interpol->p->deg; i++) {
    mpz_poly_setcoeff(p, i, interpol->p->coeff[i]);
    ASSERT(mpz_divisible_p(p->coeff[i], interpol->denom) != 0);
    /*gmp_printf("%d, %Zd\n", i, p->coeff[i]);*/
    mpz_divexact(p->coeff[i], p->coeff[i], interpol->denom);
  }
  lagrange_poly_clear(interpol);
  for (unsigned int i = 0; i < nb; i++) {
    lagrange_poly_clear(l[i]);
  }
  free(l);
}


// ---- mpz_poly_bivariate -----

#define MPZ_POLY_BIVARIATE_DEG_X 2

void mpz_poly_bivariate_init_y_x(mpz_poly_bivariate_ptr f, int dy, int dx)
{
  f->deg_y = -1;
  f->deg_x = -1;
  if (dy < 0) {
    f->alloc = 0;
    f->coeff = (mpz_poly *) NULL;
  }
  else {
    int i;
    f->alloc = dy + 1;
    f->coeff = (mpz_poly *) malloc ((dy + 1)*sizeof(mpz_poly));
    FATAL_ERROR_CHECK(f->coeff == NULL, "not enough memory");
    for (i = 0; i <= dy; ++i) {
      mpz_poly_init(f->coeff[i], dx);
    }
  }
}

void mpz_poly_bivariate_init(mpz_poly_bivariate_ptr f, int d)
{
  mpz_poly_bivariate_init_y_x(f, d, MPZ_POLY_BIVARIATE_DEG_X);
}

void mpz_poly_bivariate_realloc_x(mpz_poly_bivariate_ptr f, unsigned int nc, int dx)
{
  ASSERT_ALWAYS(nc <= (unsigned int) INT_MAX);
  if (f->alloc < nc) {
    f->coeff = (mpz_poly*) realloc(f->coeff, nc * sizeof(mpz_poly));
    FATAL_ERROR_CHECK(f->coeff == NULL, "not enough memory");
    for (unsigned int i = f->alloc; i < nc; i++) {
      mpz_poly_init(f->coeff[i], dx);
    }
    f->alloc = nc;
  }
}

void mpz_poly_bivariate_realloc(mpz_poly_bivariate_ptr f, unsigned int nc)
{
  mpz_poly_bivariate_realloc_x(f, nc, MPZ_POLY_BIVARIATE_DEG_X);
}

void mpz_poly_bivariate_clear(mpz_poly_bivariate_ptr f) 
{
  for (unsigned i = 0; i < f->alloc; ++i) {
    mpz_poly_clear(f->coeff[i]);
  }
  if (f->coeff != NULL) {
    free(f->coeff);
  }
  f->coeff = NULL; /* to avoid a double-free */
  memset(f, 0, sizeof(mpz_poly_bivariate_t));
  f->deg_y = -1;
  f->deg_x = -1;
  f->alloc = 0; /* to avoid a double-free */
}

void mpz_poly_bivariate_cleandeg(mpz_poly_bivariate_ptr f, int deg_y)
{
  ASSERT(deg_y >= -1);

  f->deg_x = -1;
  mpz_poly_cleandeg(f->coeff[deg_y], f->coeff[deg_y]->deg);
  if (f->deg_x < f->coeff[deg_y]->deg) {
    f->deg_x = f->coeff[deg_y]->deg;
  }
  while ((deg_y >= 0) && f->coeff[deg_y]->deg == -1) {
    deg_y--;
    mpz_poly_cleandeg(f->coeff[deg_y], f->coeff[deg_y]->deg);
  }
  f->deg_y = deg_y;
  for (int i = 0; i < deg_y; i++) {
    mpz_poly_cleandeg(f->coeff[i], f->coeff[i]->deg);
    if (f->deg_x < f->coeff[i]->deg) {
      f->deg_x = f->coeff[i]->deg;
    }
  }
}

void mpz_poly_bivariate_setcoeff(mpz_poly_bivariate_ptr f, int i,
    mpz_poly_srcptr z)
{
  mpz_poly_bivariate_realloc(f, i + 1);
  mpz_poly_set(f->coeff[i], z);
  if (i >= f->deg_y) {
    mpz_poly_bivariate_cleandeg(f, i);
  } else {
    if (z->deg >= f->deg_x) {
      f->deg_x = z->deg;
    } else {
      f->deg_x = -1;
      for (int i = 0; i < f->deg_y; i++) {
        mpz_poly_cleandeg(f->coeff[i], f->coeff[i]->deg);
        if (f->deg_x < f->coeff[i]->deg) {
          f->deg_x = f->coeff[i]->deg;
        }
      }
    }
  }
}

void mpz_poly_bivariate_fprintf(FILE * fp, mpz_poly_bivariate_srcptr f)
{
  if (f->deg_y == -1) {
    fprintf(fp, "0\n");
    return;
  }
  for (int i = 0, printed = 0; i <= f->deg_y; ++i) {
    if (f->coeff[i]->deg == -1) {
      continue;
    }

    if (printed++) {
      fprintf(fp, "+");
    }

    fprintf(fp, "(");
    mpz_poly_fprintf_endl(fp, f->coeff[i], 0);
    fprintf(fp, ")");
    if (i) {
      fprintf(fp, "*y");
    }

    if (i > 1) {
      fprintf(fp, "^%d", i);
    }
  }
  fprintf(fp, "\n");
}

void mpz_poly_bivariate_eval_y(mpz_poly_ptr res,
    mpz_poly_bivariate_srcptr f, mpz_srcptr y)
{
  int i, d;
  d = f->deg_y;
  if (d == -1) {
    res->deg = -1;
    return;
  }
  mpz_poly_set(res, f->coeff[d]);
  for (i = d - 1; i >= 0; --i) {
    mpz_poly_mul_mpz(res, res, y);
    mpz_poly_add(res, res, f->coeff[i]);
  }
}

void mpz_poly_bivariate_eval_x(mpz_poly_ptr res,
    mpz_poly_bivariate_srcptr f, mpz_srcptr x)
{
  int i, d;
  d = f->deg_y;
  if (d == -1) {
    res->deg = -1;
    return;
  }

  mpz_t r;
  mpz_init(r);
  for (i = 0; i <= d; i++) {
    mpz_poly_eval(r, f->coeff[i], x);
    mpz_poly_setcoeff(res, i, r);
  }
  mpz_clear(r);
}

void mpz_poly_bivariate_resultant_y(mpz_poly_ptr resultant,
    mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g)
{
  int nb_eval = MAX(f->deg_x, g->deg_x) * (f->deg_y + g->deg_y) + 1;
  mpz_t * resultants = (mpz_t * ) malloc(sizeof(mpz_t) * nb_eval);
  mpz_t * x = (mpz_t * ) malloc(sizeof(mpz_t) * nb_eval);
  mpz_poly f_eval;
  mpz_poly g_eval;
  mpz_poly_init(f_eval, -1);
  mpz_poly_init(g_eval, -1);
  int point = -nb_eval / 2;

  for (int i = 0; i < nb_eval; i++) {
    mpz_init(resultants[i]);
    mpz_init(x[i]);
    mpz_set_si(x[i], point);
    mpz_poly_bivariate_eval_x(f_eval, f, x[i]);
    mpz_poly_bivariate_eval_x(g_eval, g, x[i]);
    while (f_eval->deg != f->deg_y || g_eval->deg != g->deg_y) {
      point++;
      mpz_set_si(x[i], point);
      mpz_poly_bivariate_eval_x(f_eval, f, x[i]);
      mpz_poly_bivariate_eval_x(g_eval, g, x[i]);
    }
    mpz_poly_resultant(resultants[i], f_eval, g_eval);
    point++;
  }
  interpolate(resultant, x, resultants, nb_eval);

  mpz_poly_clear(f_eval);
  mpz_poly_clear(g_eval);
  for (int i = 0; i < nb_eval; i++) {
    mpz_clear(resultants[i]);
    mpz_clear(x[i]);
  }
  free(resultants);
  free(x);
}

void mpz_poly_bivariate_resultant_x(mpz_poly_ptr resultant,
    mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g)
{
  int nb_eval = MAX(f->deg_y, g->deg_y) * (f->deg_x + g->deg_x) + 1;
  mpz_t * resultants = (mpz_t * ) malloc(sizeof(mpz_t) * nb_eval);
  mpz_t * y = (mpz_t * ) malloc(sizeof(mpz_t) * nb_eval);
  mpz_poly f_eval;
  mpz_poly g_eval;
  mpz_poly_init(f_eval, -1);
  mpz_poly_init(g_eval, -1);
  int point = -nb_eval / 2;

  for (int i = 0; i < nb_eval; i++) {
    mpz_init(resultants[i]);
    mpz_init(y[i]);
    mpz_set_si(y[i], point);
    mpz_poly_bivariate_eval_y(f_eval, f, y[i]);
    mpz_poly_bivariate_eval_y(g_eval, g, y[i]);
    while (f_eval->deg != f->deg_x || g_eval->deg != g->deg_x) {
      point++;
      mpz_set_si(y[i], point);
      mpz_poly_bivariate_eval_y(f_eval, f, y[i]);
      mpz_poly_bivariate_eval_y(g_eval, g, y[i]);
    }
    mpz_poly_resultant(resultants[i], f_eval, g_eval);
    point++;
  }
  interpolate(resultant, y, resultants, nb_eval);

  mpz_poly_clear(f_eval);
  mpz_poly_clear(g_eval);
  for (int i = 0; i < nb_eval; i++) {
    mpz_clear(resultants[i]);
    mpz_clear(y[i]);
  }
  free(resultants);
  free(y);
}
