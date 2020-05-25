#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <string.h>
#include <gmp.h>
#include "int64_poly.h"
#include "macros.h"
#include "sieving_bound.h"
#include "macros.h"
#include "utils_int64.h"

void int64_poly_init(int64_poly_ptr f, int d)
{
  f->deg = -1;
  if (d < 0) {
    f->alloc = 0;
    f->coeff = (int64_t *) NULL;
  } else {
    f->alloc = d + 1;
    f->coeff = (int64_t *) malloc((d + 1) * sizeof(int64_t));
    memset(f->coeff, 0, (d + 1) * sizeof(int64_t));
    FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
  }
}

void int64_poly_realloc(int64_poly_ptr f, int nc)
{
  if (f->alloc < nc)
  {
    f->coeff = (int64_t*) realloc (f->coeff, nc * sizeof(int64_t));
    FATAL_ERROR_CHECK (f->coeff == NULL, "not enough memory");
    memset(f->coeff + f->alloc, 0, (nc - f->alloc) * sizeof(int64_t));
    f->alloc = nc;
  }
}

void int64_poly_copy(int64_poly_ptr g, int64_poly_srcptr f)
{
  g->deg = f->deg;
  for (int i = f->deg; i >= 0; --i) {
    int64_poly_setcoeff(g, i, f->coeff[i]);
  }
}

void int64_poly_swap(int64_poly_ptr f, int64_poly_ptr g)
{
  int i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->deg;
  f->deg = g->deg;
  g->deg = i;
  int64_t * t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

void int64_poly_clear(int64_poly_ptr f)
{
  free(f->coeff);
  memset(f, 0, sizeof(int64_poly_t));
  f->deg = -1;
}

void int64_poly_cleandeg(int64_poly_ptr f, int deg)
{
  ASSERT(deg >= -1);

  while ((deg >= 0) && (f->coeff[deg] == 0))
    deg--;
  f->deg = deg;
}

void int64_poly_set_zero(int64_poly_ptr f)
{
  f->deg = -1;
}

void int64_poly_set_xi(int64_poly_ptr f, int i)
{
  int64_poly_realloc (f, i + 1);
  for(int j = 0 ; j < i + 1 ; j++) {
    f->coeff[j] = (j == i);
  }
  f->deg = i;
}

void int64_poly_set_bxi(int64_poly_ptr f, int i, int64_t b)
{
  int j = 0;
  int64_poly_realloc (f, i + 1);
  for( ; j < i ; j++) {
    f->coeff[j] = 0;
  }
  f->coeff[i] = b;
  f->deg = i;
}

int int64_poly_equal(int64_poly_srcptr a, int64_poly_srcptr b)
{
  int r = (a->deg > b->deg) - (b->deg > a->deg);
  if (r) return 0;
  for(int d = a->deg; d >= 0 ; d--) {
    r = a->coeff[d] - b->coeff[d];
    if (r) return 0;
  }
  return 1;
}

int int64_poly_normalized(int64_poly_srcptr f)
{
  return (f->deg == -1) || f->coeff[f->deg] !=  0;
}

void int64_poly_fprintf(FILE * file, int64_poly_srcptr f)
{
  int i;
  if (f->deg == -1) {
    fprintf(file, "0\n");
    return;
  }
  else if (f->deg == 0) {
    printf("%" PRId64 "\n", f->coeff[0]);
    return;
  }
  fprintf(file, "%" PRId64 "", f->coeff[0]);
  for (i = 1; i <= f->deg; i++) {
    if (f->coeff[i] >= 0) {
      fprintf(file, "+%" PRId64 "*x^%d", f->coeff[i], i);
    } else {
      fprintf(file, "%" PRId64 "*x^%d", f->coeff[i], i);
    }
  }
  fprintf(file, "\n");
}

int64_t int64_poly_max(int64_poly_srcptr f)
{
  int i = 1;
  int64_t max = f->coeff[0];
  for ( ; i < f->deg + 1; i++) {
    max = MAX(max, f->coeff[i]);
  }
  return max;
}

uint64_t int64_poly_infinity_norm(int64_poly_srcptr f)
{
  int i = 1;
  uint64_t in = (uint64_t)ABS(f->coeff[0]);
  for ( ; i < f->deg + 1; i++) {
    in = MAX(in, (uint64_t)ABS(f->coeff[i]));
  }
  return in;
}

void int64_poly_to_mpz_poly(mpz_poly_ptr a, int64_poly_srcptr b)
{
  mpz_t coeff;
  mpz_init(coeff);
  if (b->deg != -1) {
    for (int i = 0; i <= b->deg; i++) {
      mpz_poly_setcoeff_int64(a, i, b->coeff[i]);
    }
  }
  mpz_clear(coeff);
}
