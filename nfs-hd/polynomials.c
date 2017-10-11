#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "utils.h"
#include "polynomials.h"
#include "utils_cofactorisation.h"
#include "polyselect/auxiliary.h"

MAYBE_UNUSED static void extended_euclidean_algorithm(mpz_ptr r0, mpz_ptr t0,
    mpz_srcptr f, mpz_srcptr g)
{
  mpz_t s0;
  mpz_init(s0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t r2;
  mpz_init(r2);
  mpz_t s2;
  mpz_init(s2);
  mpz_t t2;
  mpz_init(t2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  unsigned int i = 1;
  while (mpz_cmp_ui(r1, 0) != 0) {
    mpz_fdiv_qr(q, r2, r0, r1);
    mpz_set(r2, r0);
    mpz_submul(r2, q, r1);
    mpz_set(s2, s0);
    mpz_submul(s2, q, s1);
    mpz_set(t2, t0);
    mpz_submul(t2, q, t1);

#ifndef NDEBUG
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set(tmp, f);
    mpz_mul(tmp, tmp, s1);
    mpz_addmul(tmp, t1, g);
    ASSERT(mpz_cmp(tmp, r1) == 0);
    mpz_clear(tmp);
#endif

    mpz_set(r0, r1);
    mpz_set(s0, s1);
    mpz_set(t0, t1);
    mpz_set(r1, r2);
    mpz_set(s1, s2);
    mpz_set(t1, t2);
    i++;
  }

  mpz_clear(s0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(r2);
  mpz_clear(s2);
  mpz_clear(t2);
  mpz_clear(q);
}

#ifdef EEA_BOUND
static void extended_euclidean_algorithm_bound(mpz_ptr r0,
    mpz_ptr t0, mpz_srcptr f, mpz_srcptr g, mpz_srcptr bound)
{
  mpz_t s0;
  mpz_init(s0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t r2;
  mpz_init(r2);
  mpz_t s2;
  mpz_init(s2);
  mpz_t t2;
  mpz_init(t2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  unsigned int i = 1;
  while (mpz_cmp_ui(r1, 0) != 0 && mpz_cmp(r1, bound) <= 0) {
    mpz_fdiv_qr(q, r2, r0, r1);
    mpz_set(r2, r0);
    mpz_submul(r2, q, r1);
    mpz_set(s2, s0);
    mpz_submul(s2, q, s1);
    mpz_set(t2, t0);
    mpz_submul(t2, q, t1);

#ifndef NDEBUG
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set(tmp, f);
    mpz_mul(tmp, tmp, s1);
    mpz_addmul(tmp, t1, g);
    ASSERT(mpz_cmp(tmp, r1) == 0);
    mpz_clear(tmp);
#endif

    mpz_set(r0, r1);
    mpz_set(s0, s1);
    mpz_set(t0, t1);
    mpz_set(r1, r2);
    mpz_set(s1, s2);
    mpz_set(t1, t2);
    i++;
  }

  mpz_clear(s0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(r2);
  mpz_clear(s2);
  mpz_clear(t2);
  mpz_clear(q);
}
#else // EEA_BOUND
static void extended_euclidean_algorithm_stop_2(mpz_ptr r2, mpz_ptr t2,
    mpz_srcptr f, mpz_srcptr g)
{
  mpz_t r0;
  mpz_init(r0);
  mpz_t s0;
  mpz_init(s0);
  mpz_t t0;
  mpz_init(t0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t s2;
  mpz_init(s2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  mpz_fdiv_qr(q, r2, r0, r1);
  mpz_set(r2, r0);
  mpz_submul(r2, q, r1);
  mpz_set(s2, s0);
  mpz_submul(s2, q, s1);
  mpz_set(t2, t0);
  mpz_submul(t2, q, t1);

  mpz_clear(r0);
  mpz_clear(s0);
  mpz_clear(t0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(s2);
  mpz_clear(q);
}
#endif // EEA_BOUND

unsigned int mpz_poly_is_reciprocal(mpz_poly_srcptr p)
{
  int nb = p->deg / 2;
  if (p->deg % 2 != 0) {
    nb++;
  }

  for (int i = 0; i < nb; i++) {
    if (mpz_cmp(p->coeff[i], p->coeff[p->deg - i]) != 0) {
      return 0;
    }
  }
  return 1;
}

static double find_epsilon(unsigned int q, mpz_srcptr p, unsigned int t)
{
  double q_d = pow(2.0, (double)q);
  double p_d = mpz_get_d(p);

  return log(q_d) / (2 * (double)(t - 1) * log(p_d));
}

//Generate a random number between [offset, offset + length[
static inline void sub_rand_mpz(mpz_ptr rand_Z, mpz_srcptr length,
    mpz_srcptr offset, gmp_randstate_t state)
{
  mpz_urandomm(rand_Z, state, length);
  mpz_add(rand_Z, rand_Z, offset);
}

static void rand_mpz(mpz_ptr rand_Z, mpz_srcptr min, mpz_srcptr max,
    gmp_randstate_t state)
{
  mpz_t length;
  mpz_init(length);
  mpz_sub(length, max, min);
  sub_rand_mpz(rand_Z, length, min, state);
  mpz_clear(length);
}

static void random_mpz_poly(mpz_poly_ptr g, mpz_srcptr min,
    mpz_srcptr max, int degree, int lc, gmp_randstate_t state)
{
  mpz_t rand_Z;
  mpz_init(rand_Z);
  for (int i = 0; i < degree; i++) {
    rand_mpz(rand_Z, min, max, state);
    mpz_poly_setcoeff(g, i, rand_Z);
  }
  if (lc == 0) {
    while(mpz_cmp_ui(rand_Z, 0) == 0) {
      rand_mpz(rand_Z, min, max, state);
    }
    mpz_poly_setcoeff(g, degree, rand_Z);
  } else {
    ASSERT(lc != 0);

    mpz_poly_setcoeff_si(g, degree, lc);
  }
  mpz_clear(rand_Z);
}

static void random_mpz_poly_reciprocal(mpz_poly_ptr g, mpz_srcptr min,
    mpz_srcptr max, int degree, int lc, gmp_randstate_t state)
{
  int i = 0;
  int nb = degree / 2;
  if (degree % 2 != 0) {
    nb++;
  }

  mpz_t rand_Z;
  mpz_init(rand_Z);
  for (i = 1; i < nb; i++) {
    rand_mpz(rand_Z, min, max, state);
    mpz_poly_setcoeff(g, i, rand_Z);
    mpz_poly_setcoeff(g, degree - i, rand_Z);
  }
  if (degree % 2 == 0) {
    rand_mpz(rand_Z, min, max, state);
    mpz_poly_setcoeff(g, i, rand_Z);
  }
  if (lc == 0) {
    while(mpz_cmp_ui(rand_Z, 0) == 0) {
      rand_mpz(rand_Z, min, max, state);
    }
    mpz_poly_setcoeff(g, degree, rand_Z);
    mpz_poly_setcoeff(g, 0, rand_Z);
  } else {
    ASSERT(lc != 0);

    mpz_poly_setcoeff_si(g, degree, lc);
    mpz_poly_setcoeff_si(g, 0, lc);
  }
  mpz_clear(rand_Z);

  ASSERT(mpz_poly_is_reciprocal(g));
}

static void random_mpz_poly_constraint(mpz_poly_ptr g, mpz_srcptr min,
    mpz_srcptr max, int degree, gmp_randstate_t state,
    mpz_poly_srcptr h)
{
  ASSERT(degree >= 1);
  ASSERT(degree <= h->deg);

  mpz_t rand_Z;
  mpz_init(rand_Z);
  mpz_t zero;
  mpz_init(zero);

  ASSERT(mpz_cmp_ui(zero, 0) == 0);

  for (int i = 0; i < degree; i++) {
    if (mpz_cmp_ui(h->coeff[i], 0) < 0) {
    /*if (mpz_cmp_ui(h->coeff[i], 0) > 0) {*/
      rand_mpz(rand_Z, zero, max, state);

      ASSERT(mpz_cmp_ui(rand_Z, 0) >= 0);
    } else {
      rand_mpz(rand_Z, min, zero, state);

      ASSERT(mpz_cmp_ui(rand_Z, 0) < 0);
    }
    mpz_poly_setcoeff(g, i, rand_Z);
  }

  if (mpz_cmp_ui(h->coeff[degree], 0) < 0) {
  /*if (mpz_cmp_ui(h->coeff[degree], 0) > 0) {*/
    rand_mpz(rand_Z, zero, max, state);

    ASSERT(mpz_cmp_ui(rand_Z, 0) >= 0);
  } else {
    rand_mpz(rand_Z, min, zero, state);

    ASSERT(mpz_cmp_ui(rand_Z, 0) < 0);
  }
  while(mpz_cmp_ui(rand_Z, 0) == 0) {
    if (mpz_cmp_ui(h->coeff[degree], 0) < 0) {
    /*if (mpz_cmp_ui(h->coeff[degree], 0) > 0) {*/
      rand_mpz(rand_Z, min, zero, state);
    } else {
      rand_mpz(rand_Z, zero, max, state);
    }
  }
  mpz_poly_setcoeff(g, degree, rand_Z);
  mpz_clear(rand_Z);
}

#ifndef ALPHA
static double double_mpz_poly_infinty_norm(mpz_poly_srcptr a)
{
  ASSERT(a->deg >= 1);

  double max = abs(mpz_get_d(a->coeff[0]));
  for (int i = 1; i <= a->deg; i++) {
    max = MAX(max, abs(mpz_get_d(a->coeff[i])));
  }

  return max;
}
#endif // ALPHA

static void find_c(mpz_ptr c, mpz_ptr a, mpz_ptr b, int c_tol, mpz_srcptr p) {
  mpz_t a_tmp;
  mpz_t b_tmp;
  mpz_t c_tmp;
  mpz_init(a_tmp);
  mpz_init(b_tmp);
  mpz_init(c_tmp);
  mpz_t c_tol_max;
  mpz_init(c_tol_max);
  mpz_set_si(c_tmp, -c_tol);
  mpz_add(c_tmp, c_tmp, c);
  mpz_set_si(c_tol_max, c_tol + 1);
  mpz_add(c_tol_max, c_tol_max, c);

  mpz_t min;
  mpz_init(min);
  mpz_set(min, c);
  mpz_t min_tmp;
  mpz_init(min_tmp);

  while (mpz_cmp(c_tmp, c_tol_max) <= 0) {
#ifdef EEA_BOUND
    extended_euclidean_algorithm_bound(a_tmp, b_tmp, p, c_tmp, c_tmp);
#else // EEA_BOUND
    extended_euclidean_algorithm_stop_2(a_tmp, b_tmp, p, c_tmp);
#endif // EEA_BOUND
    mpz_sub(min_tmp, c_tmp, a_tmp);
    mpz_abs(min_tmp, min_tmp);
    if (mpz_cmp(min, min_tmp) > 0) {
      mpz_set(c, c_tmp);
      mpz_set(min, min_tmp);
      mpz_set(a, a_tmp);
      mpz_set(b, b_tmp);
    }
    mpz_add_ui(c_tmp, c_tmp, 1);
  }

  mpz_clear(min_tmp);
  mpz_clear(min);
  mpz_clear(c_tol_max);
  mpz_clear(a_tmp);
  mpz_clear(b_tmp);
  mpz_clear(c_tmp);
}

static void rewrite_f(mpz_poly_ptr f, mpz_srcptr p) {
  mpz_t p2;
  mpz_init(p2);
  mpz_fdiv_q_ui(p2, p, 2);
  mpz_t p2n;
  mpz_init(p2n);
  mpz_set(p2n, p2);
  mpz_mul_si(p2n, p2n, -1);

  mpz_t Z_tmp;
  mpz_init(Z_tmp);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(Z_tmp, i, f);
    if (mpz_sgn(Z_tmp) == 1) {
      if (mpz_cmp(Z_tmp, p2) > 0) {
        mpz_sub(Z_tmp, Z_tmp, p);
        mpz_poly_setcoeff(f, i, Z_tmp);
      }
    } else if (mpz_cmp(Z_tmp, p2) < 0) {
      if (mpz_cmp(Z_tmp, p2n) < 0) {
        mpz_add(Z_tmp, Z_tmp, p);
        mpz_poly_setcoeff(f, i, Z_tmp);
      }
    }
  }
  mpz_clear(Z_tmp);
  mpz_clear(p2);
  mpz_clear(p2n);
}

static void gen_g_gal(mpz_poly_ptr g, mpz_poly_srcptr h,
    mpz_srcptr min, mpz_srcptr max, gmp_randstate_t state) {
  int degree = (rand() % h->deg) + 1;
  int degree_min = 0;
  if (h->deg % 2 == 0) {
    degree_min = h->deg / 2;
  } else {
    degree_min = h->deg / 2 + 1;
  }
  while (degree < degree_min) {
    degree = (rand() % h->deg) + 1;
  }

  ASSERT(degree >= 1);
  ASSERT(degree <= h->deg);
#ifndef NDEBUG
  if (h->deg % 2 == 0) {
    ASSERT(degree >= h->deg / 2);
  } else {
    ASSERT(degree >= h->deg / 2 + 1);
  }
#endif // NDEBUG

  mpz_t rand_Z;
  mpz_init(rand_Z);
  mpz_t zero;
  mpz_init(zero);
  ASSERT(mpz_cmp_ui(zero, 0) == 0);

  int nb = degree / 2;
  int i = 0;
  if (degree % 2 != 0) {
    nb++;
    i = nb;
  } else {
    i = nb;
    if (mpz_cmp_ui(h->coeff[i], 0) < 0) {
      /*if (mpz_cmp_ui(h->coeff[degree], 0) > 0) {*/
      rand_mpz(rand_Z, min, zero, state);
    } else {
      rand_mpz(rand_Z, zero, max, state);
    }
    mpz_poly_setcoeff(g, i, rand_Z);
    i++;
  }

  for (; i < degree; i++) {
    if (mpz_cmp_ui(h->coeff[i], 0) < 0) {
      /*if (mpz_cmp_ui(h->coeff[degree], 0) > 0) {*/
      rand_mpz(rand_Z, min, zero, state);
    } else {
      rand_mpz(rand_Z, zero, max, state);
    }
    mpz_poly_setcoeff(g, i, rand_Z);
    mpz_poly_setcoeff(g, h->deg - i, rand_Z);
  }

  while(mpz_cmp_ui(rand_Z, 0) == 0) {
    rand_mpz(rand_Z, min, max, state);
  }
  mpz_poly_setcoeff(g, degree, rand_Z);
  mpz_poly_setcoeff(g, h->deg - degree, rand_Z);

  mpz_clear(rand_Z);
  mpz_clear(zero);
}

static void mpz_poly_erase_coefficients(mpz_poly_ptr g) {
  for (int i = 0; i <= g->deg; i++) {
    mpz_set_ui(g->coeff[i], 0);
  }
  g->deg = -1;
}

double function_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1,
    mpz_poly_ptr g, mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p,
    mpz_poly_ptr h, int * coeff, unsigned int q, unsigned int t,
    unsigned int nb_times, double * weight, int c_tol, gmp_randstate_t state,
    unsigned int h_set, unsigned int gal)
{
  ASSERT(g->alloc > 0);

  mpz_t * coeff_Z = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (unsigned int i = 0; i < 2; i++) {
    mpz_init(coeff_Z[i]);
    mpz_set_si(coeff_Z[i], coeff[i]);
  }

  double epsilon = find_epsilon(q, p, t);
  mpz_set_d(c, pow(mpz_get_d(p), 0.5 - epsilon) + 1.0);

  //Select the c=a/b[p] such that |c - a| is minimal
  find_c(c, a, b, c_tol, p);

  double alpha0 = 0.0;
  double alpha1 = 0.0;
  double sum_alpha = DBL_MAX;
  mpz_poly f0_tmp;
  mpz_poly_init(f0_tmp, -1);
  mpz_poly f1_tmp;
  mpz_poly_init(f1_tmp, -1);
  mpz_poly g_tmp;
  mpz_poly_init(g_tmp, -1);
  mpz_poly h_tmp;
  mpz_poly_init(h_tmp, -1);
  if (h_set) {
    mpz_poly_set(h_tmp, h);
  }

  unsigned int k = 0;
  while (k < nb_times) {
    if (h_set) {
      if (gal == 0) {
        random_mpz_poly_constraint(g_tmp, coeff_Z[0], coeff_Z[1],
            (rand() % h->deg) + 1, state, h);

      } else if (gal == 2) {
        gen_g_gal(g_tmp, h, coeff_Z[0], coeff_Z[1], state);
      }
    } else {
      if (gal == 0) {
        random_mpz_poly(h_tmp, coeff_Z[0], coeff_Z[1], h->deg, 0, state);
        random_mpz_poly_constraint(g_tmp, coeff_Z[0], coeff_Z[1],
            (rand() % h->deg) + 1, state, h);
      } else if (gal == 2) {
        random_mpz_poly_reciprocal(h_tmp, coeff_Z[0], coeff_Z[1], h->deg, 0,
            state);
        gen_g_gal(g_tmp, h, coeff_Z[0], coeff_Z[1], state);
      }
    }
    mpz_poly_set(f0_tmp, h_tmp);
    mpz_poly_mul_mpz(f0_tmp, f0_tmp, c);
    mpz_poly_add(f0_tmp, f0_tmp, g_tmp);
#ifdef ALPHA
    alpha0 = get_alpha(f0_tmp, ALPHA_BOUND);
#else // ALPHA
    alpha0 = log2(double_mpz_poly_infinty_norm(f0_tmp));
#endif // ALPHA
    if (mpz_poly_is_irreducible(f0_tmp, p)) {
      k++;
      mpz_poly_mul_mpz(f1_tmp, g_tmp, b);
      mpz_poly f_tmp;
      mpz_poly_init(f_tmp, -1);
      mpz_poly_mul_mpz(f_tmp, h_tmp, a);
      mpz_poly_add(f1_tmp, f1_tmp, f_tmp);
      mpz_poly_clear(f_tmp);

      ASSERT(mpz_poly_is_irreducible(f1_tmp, p));

      rewrite_f(f0_tmp, p);
      rewrite_f(f1_tmp, p);

      if (gal == 2) {
        ASSERT(mpz_poly_is_reciprocal(f0_tmp));
        ASSERT(mpz_poly_is_reciprocal(f1_tmp));
      }

#ifdef ALPHA
      alpha1 = get_alpha(f1_tmp, ALPHA_BOUND);
#else // ALPHA
      alpha1 = log2(double_mpz_poly_infinty_norm(f1_tmp));
#endif // ALPHA
      if (weight[0] * alpha0 + weight[1] * alpha1 < sum_alpha) {
        mpz_poly_set(f0, f0_tmp);
        mpz_poly_set(f1, f1_tmp);
        mpz_poly_set(g, g_tmp);
        mpz_poly_set(h, h_tmp);
        sum_alpha = weight[0] * alpha0 + weight[1] * alpha1;
      }
    }
    mpz_poly_erase_coefficients(g_tmp);
  }

  mpz_poly_clear(f0_tmp);
  mpz_poly_clear(f1_tmp);
  mpz_poly_clear(g_tmp);
  mpz_poly_clear(h_tmp);

  for (unsigned int i = 0; i < 2; i++) {
    mpz_clear(coeff_Z[i]);
  }
  free(coeff_Z);

  return epsilon;
}

void function_classical(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_srcptr p,
    unsigned int n, int * coeff, unsigned int nb_times,
    double * weight, gmp_randstate_t state)
{
  mpz_t * coeff_Z = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (unsigned int i = 0; i < 2; i++) {
    mpz_init(coeff_Z[i]);
    mpz_set_si(coeff_Z[i], coeff[i]);
  }

  mpz_poly f0_tmp;
  mpz_poly_init(f0_tmp, -1);
  mpz_poly f1_tmp;
  mpz_poly_init(f1_tmp, -1);

  double sum_alpha = DBL_MAX;
  double alpha0 = 0.0;
  double alpha1 = 0.0;
  unsigned int k = 0;
  while (k < nb_times) {
    random_mpz_poly(f0_tmp, coeff_Z[0], coeff_Z[1], n, 1, state);
    if (mpz_poly_is_irreducible(f0_tmp, p)) {
      k++;
      mpz_poly_set(f1_tmp, f0_tmp);
      if (mpz_cmp_ui(f1_tmp->coeff[0], 0) < 0) {
        mpz_add(f1_tmp->coeff[0], f1_tmp->coeff[0], p);
      } else {
        mpz_sub(f1_tmp->coeff[0], f1_tmp->coeff[0], p);
      }
      alpha0 = get_alpha(f0_tmp, ALPHA_BOUND);
      alpha1 = get_alpha(f1_tmp, ALPHA_BOUND);
      if (weight[0] * alpha0 + weight[1] * alpha1 < sum_alpha) {
        mpz_poly_set(f0, f0_tmp);
        mpz_poly_set(f1, f1_tmp);
        sum_alpha = weight[0] * alpha0 + weight[1] * alpha1;
      }
    }
  }

  mpz_poly_clear(f0_tmp);
  mpz_poly_clear(f1_tmp);

  for (unsigned int i = 0; i < 2; i++) {
    mpz_clear(coeff_Z[i]);
  }
  free(coeff_Z);
}
