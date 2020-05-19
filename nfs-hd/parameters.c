#include "cado.h"
// TODO: did not work.

#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "sieving_bound.h"
#include "parameters.h"
#include "polyselect/rho.h"
#include "polynomials.h"
#include <stdint.h>

double size_sieving_region(sieving_bound_srcptr H)
{
  double res = 1;
  for (unsigned int i = 0; i < H->t - 1; i++) {
    res = res * 2 * (double) H->h[i];
  }
  res = res * (double) H->h[H->t - 1];

  ASSERT(res > 0.0);

  return res;
}

void find_sieving_region_square(sieving_bound_ptr H, double number,
    unsigned int t)
{
  H->t = t;
  unsigned int H0 = (unsigned int)
    ceil(pow((double)number / pow(2.0, (double)(t - 1)), 1.0 / (double) t));
  for (unsigned int i = 0; i < t; i++) {
    H->h[i] = H0;
  }
}

void rand_mpz_poly(mpz_poly_ptr a, sieving_bound_srcptr H)
{
  for (int i = 0; i < (int)H->t - 1; i++) {
    mpz_poly_setcoeff_si(a, i, rand() % (2 * (int)H->h[i]) - (int)H->h[i]);
  }
  mpz_poly_setcoeff_si(a, (int)H->t - 1, rand() % ((int)H->h[H->t - 1]));

#ifndef NDEBUG
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i < a->deg; i++) {
    mpz_poly_getcoeff(tmp, i, a);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  if (a->deg == (int)H->t - 1) {
    mpz_poly_getcoeff(tmp, H->t - 1, a);
    ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) <= 0);
    ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
    mpz_clear(tmp);
  }
#endif
}

void stat_approx_number(mpz_t * mean, mpz_t * max, mpz_poly * f,
    unsigned int nb_fields, uint64_t number_a, sieving_bound_srcptr H)
{
  for (unsigned int i = 0; i < nb_fields; i++) {
    if (mean != NULL) {
      mpz_set_ui(mean[i], 0);
    }
    if (max != NULL) {
      mpz_set_ui(max[i], 0);
    }
  }

  srand(time(NULL));
  mpz_poly a;
  mpz_poly_init(a, (int)H->t);

  mpz_t res;
  mpz_init(res);
  for (uint64_t i = 0; i < number_a; i++) {
    rand_mpz_poly(a, H);
    for(unsigned int j = 0; j < nb_fields; j++) {
      mpz_poly_resultant(res, f[j], a);
      mpz_abs(res, res);
      if (mean != NULL) {
        mpz_add(mean[j], mean[j], res);
      }
      if (max != NULL) {
        if (mpz_cmp(res, max[j]) > 0) {
          mpz_set(max[j], res);
        }
      }
    }
  }

  if (mean != NULL) {
    for(unsigned int i = 0; i < nb_fields; i++) {
      mpz_cdiv_q_ui(mean[i], mean[i], number_a);
    }
  }

  mpz_clear(res);
  mpz_poly_clear(a);
}

double compute_smoothness(mpz_t * max, unsigned int V, double lpb)
{
  double tmp = (double)mpz_sizeinbase(max[0], 2) / lpb;
  double res = 1.0;
  for (unsigned int i = 0; i < V; i++) {
    res = res * pow(tmp, -tmp);
  }
  return res;
}

void find_lpb(mpz_t * max, MAYBE_UNUSED unsigned int lpb, unsigned int V,
    MAYBE_UNUSED double smoothness, double error)
{
  size_t min = mpz_sizeinbase(max[0], 2);
  for (unsigned int i = 1; i < V; i++) {
    min = MIN(min, mpz_sizeinbase(max[i], 2));
  }

  double lpb_min = 0;
  double lpb_max = (double) min;

  while (lpb_max - lpb_min >= error) {
    double m = (lpb_min + lpb_max) / 2;
    if (compute_smoothness(max, V, lpb_min) * compute_smoothness(max, V, m)) {
      lpb_min = m;
    } else {
      lpb_max = m;
    }
  }
}

void find_parameters(mpz_poly * f, unsigned int V, MAYBE_UNUSED double size,
    unsigned int t, MAYBE_UNUSED double smoothness, unsigned int number_poly,
    MAYBE_UNUSED unsigned int found)
{
  sieving_bound_t H;
  sieving_bound_init(H, t);
  unsigned int * max_bit = (unsigned int * ) malloc(sizeof(unsigned int) * V);
  mpz_t * max = (mpz_t * ) malloc(sizeof(mpz_t) * V);
  for (unsigned int i = 0; i < V; i++) {
    mpz_init(max[i]);
  }
  stat_approx_number(NULL, max, f, V, number_poly, H);

  free(max_bit);
  free(max);
  sieving_bound_clear(H);
}
