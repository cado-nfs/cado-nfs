#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include "double_poly.h"
#include "mpz_poly.h"
#include "sieving_bound.h"
#include "int64_vector.h"
#include "mat_int64.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))

/*[[255, -50, 47],*/
/*[-83, -8, 416],*/
/*[-147, -285, 128]]*/

/*mpz_poly_setcoeff_si(f_Z, 0, 1);*/
/*mpz_poly_setcoeff_si(f_Z, 1, 112148);*/
/*mpz_poly_setcoeff_si(f_Z, 2, 280355);*/
/*mpz_poly_setcoeff_si(f_Z, 3, -20);*/
/*mpz_poly_setcoeff_si(f_Z, 4, -280370);*/
/*mpz_poly_setcoeff_si(f_Z, 5, -112142);*/
/*mpz_poly_setcoeff_si(f_Z, 6, 1);*/


void resultant()
{
  unsigned int dim = 3;

  mpz_poly a_Z;
  mpz_poly_init(a_Z, dim);

  mpz_poly f_Z;
  mpz_poly_init(f_Z, 6);

  mat_int64_t M;
  mat_int64_init(M, dim, dim);

  M->coeff[1][1] = 255;
  M->coeff[1][2] = -50;
  M->coeff[1][3] = 47;
  M->coeff[2][1] = -83;
  M->coeff[2][2] = -8;
  M->coeff[2][3] = 416;
  M->coeff[3][1] = -147;
  M->coeff[3][2] = -285;
  M->coeff[3][3] = 128;

  mpz_poly_setcoeff_si(f_Z, 0, 56035);
  mpz_poly_setcoeff_si(f_Z, 1, 323734);
  mpz_poly_setcoeff_si(f_Z, 2, -31190);
  mpz_poly_setcoeff_si(f_Z, 3, -1120700);
  mpz_poly_setcoeff_si(f_Z, 4, -809335);
  mpz_poly_setcoeff_si(f_Z, 5, 12476);
  mpz_poly_setcoeff_si(f_Z, 6, 56035);

  double_poly a_d;
  double_poly_init(a_d, dim);

  double_poly f_d;
  double_poly_init(f_d, 6);

  f_d->coeff[0] = (double) 56035;
  f_d->coeff[1] = (double) 323734;
  f_d->coeff[2] = (double) -31190;
  f_d->coeff[3] = (double) -1120700;
  f_d->coeff[4] = (double) -809335;
  f_d->coeff[5] = (double) 12476;
  f_d->coeff[6] = (double) 56035;

  sieving_bound_t H;
  sieving_bound_init(H, dim);
  for (unsigned int i = 0; i < dim; i++) {
    H->h[i] = 128;
  }

  uint64_t size = 1;

  int64_vector_t v;
  int64_vector_init(v, dim);
  for (unsigned int i = 0; i < dim - 1; i++) {
    v->c[i] = -(int64_t)H->h[i];
    size = size * (2 * (uint64_t)H->h[i] - 1);
  }
  v->c[0] = v->c[0] - 1;
  v->c[dim - 1] = 0;
  size = size * ((uint64_t)H->h[dim - 1] - 1);

  mpz_t norm_Z;
  mpz_init(norm_Z);

  double norm_d = 0.0;

  int size_Z = 0;
  int size_d = 0;

  for (uint64_t i = 0; i < size; i++) {

    int64_vector_add_one(v, H);

    int deg = -1;
    for (unsigned int i = 0; i < dim; i++) {
      int64_t tmp = 0;
      for (unsigned int k = 0; k < M->NumCols; k++) {
        tmp = tmp + M->coeff[i + 1][k + 1] * (int64_t) v->c[k];
      }
      mpz_poly_setcoeff_int64(a_Z, i, tmp);
      a_d->coeff[i] = (double) tmp;
      if (tmp != 0) {
        deg = (int) i;
      }
    }
    a_d->deg = deg;

    mpz_poly_resultant(norm_Z, a_Z, f_Z);
    mpz_abs(norm_Z, norm_Z);
    norm_d = fabs(double_poly_resultant(a_d, f_d));

    if (mpz_cmp_ui(norm_Z, 0) == 0) {
      size_Z = 1;
    } else {
      size_Z = (int) mpz_sizeinbase(norm_Z, 2);
    }
    if (norm_d == 0.0) {
      size_d = 1;
    } else {
      size_d = (int) ceil(log2(norm_d));
    }

    if (ABS(size_Z - size_d) > 1) {
      mpz_poly_fprintf(stdout, f_Z);
      mpz_poly_fprintf(stdout, a_Z);
      printf("Gap: %d; mpz: %u -- double: %u\n", ABS(size_Z - size_d), size_Z,
          size_d);
    }
  }

  mat_int64_clear(M);
  mpz_clear(norm_Z);
  sieving_bound_clear(H);
  int64_vector_clear(v);
  mpz_poly_clear(a_Z);
}

int main()
{
  resultant();

  exit(EXIT_SUCCESS);
}
