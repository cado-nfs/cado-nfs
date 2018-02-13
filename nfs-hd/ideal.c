#include "cado.h"
#include "ideal.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include "macros.h"
#include <inttypes.h>

/* Ideal */

void ideal_init(ideal_ptr ideal)
{
  mpz_poly_init(ideal->h, -1);
  ideal->r = 0;
}

void ideal_set_part(ideal_ptr ideal, uint64_t r, mpz_poly_srcptr h)
{
  mpz_poly_set(ideal->h, h);
  ideal->r = r;
}

void ideal_clear(ideal_ptr ideal)
{
  mpz_poly_clear(ideal->h);
  ideal->r = 0;
}

void ideal_fprintf(FILE * file, ideal_srcptr ideal)
{
  fprintf(file, "r: %" PRIu64 ", h: ", ideal->r);
  mpz_poly_fprintf(file, ideal->h);
}

/* Ideal with h of degree 1 */

void ideal_1_init(ideal_1_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0.0;
}

void ideal_1_set_part(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    unsigned int t)
{
  ASSERT(h->deg == 1);
  ASSERT(mpz_cmp_ui(mpz_poly_lc(h), 1) == 0);

  mpz_poly_set(ideal->ideal->h, h);

  //if r == 0, Tr is not created.
  if (ideal->ideal->r == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  mpz_set(ideal->Tr[0], h->coeff[0]);

  for (unsigned int i = 1; i < t - 1; i++) {
    mpz_set_si(ideal->Tr[i], -1);
    mpz_mul(ideal->Tr[i], ideal->Tr[i], ideal->Tr[i - 1]);
    mpz_mul(ideal->Tr[i], ideal->Tr[i], ideal->Tr[0]);
  }

  ideal->ideal->r = r;
  ideal->log = log2((double) r);
}

void ideal_1_set_element(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    mpz_t * Tr, double log, unsigned int t)
{
  mpz_poly_set(ideal->ideal->h, h);

  if (ideal->ideal->r == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set(ideal->Tr[i], Tr[i]);
  }

  ideal->log = log;

  ideal->ideal->r = r;
}

void ideal_1_set(ideal_1_ptr ideal_new, ideal_1_srcptr ideal_old,
    unsigned int t)
{
  mpz_poly_set(ideal_new->ideal->h, ideal_old->ideal->h);
  if (ideal_new->ideal->r == 0) {
    ideal_new->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal_new->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set(ideal_new->Tr[i], ideal_old->Tr[i]);
  }
  ideal_new->log = ideal_old->log;
  ideal_new->ideal->r = ideal_old->ideal->r;
}

void ideal_1_clear(ideal_1_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg == 1 && ideal->ideal->r != 0);

  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(ideal->Tr[i]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0.0;
}

void ideal_1_fprintf(FILE * file, ideal_1_srcptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg == 1);

  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %f\n", ideal->log);
  fprintf(file, "Tr: [");
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd]\n", ideal->Tr[t - 2]);
}

/* Ideal with h of degree u > 1 */

void ideal_u_init(ideal_u_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0.0;
}

/*
   If you want to line sieve, you can need to have Tr only with the coefficient
   modulo r.
   */
void ideal_u_set_part(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    unsigned int t)
{
  ASSERT(h->deg > 1);
  ASSERT(mpz_cmp_ui(mpz_poly_lc(h), 1) == 0);

  mpz_poly_set(ideal->ideal->h, h);
  //if r == 0, Tr is not created.
  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t ** ) malloc(sizeof(mpz_t * ) * (unsigned int)(h->deg));
    for (int row = 0; row < h->deg; row++) {
      ideal->Tr[row] = (mpz_t * ) malloc(sizeof(mpz_t) * (t - (unsigned
              int)h->deg));
    }
    for (int row = 0; row < h->deg; row++) {
      for (int col = 0; col < ((int)t - h->deg); col++) {
        mpz_init(ideal->Tr[row][col]);
      }
    }
  }

#ifndef NDEBUG
  else {
    ASSERT(h->deg == ideal->ideal->h->deg);
  }
#endif // NDEBUG

  //Initialize the elements of Tr.
  for (int col = 0; col < (int) t - h->deg; col++) {
    for (int row = col; row < h->deg; row++) {
      mpz_set(ideal->Tr[row][col], h->coeff[row-col]);
    }
  }

  //Perform the linear combination.
  for (int col = 1; col < (int) t - h->deg; col++) {
    for (int row = 0; row < h->deg; row++) {
      for (int k = 0; k < h->deg; k++) {
        if (col - h->deg + k >= 0) {
          mpz_submul(ideal->Tr[row][col], h->coeff[k],
              ideal->Tr[row][col - h->deg + k]);
        }
      }
    }
  }

  ideal->ideal->r = r;
  ideal->log = log2((double) r) * (double)h->deg;
}

void ideal_u_set_element(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    mpz_t * Tr, double log, unsigned int t)
{
  mpz_poly_set(ideal->ideal->h, h);
  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t ** ) malloc(sizeof(mpz_t * ) * h->deg);
    for (int row = 0; row < h->deg; row++) {
      ideal->Tr[row] = (mpz_t * ) malloc(sizeof(mpz_t) * ((int)t - h->deg));
      for (int col = 0; col < ((int)t - h->deg); col++) {
        mpz_init(ideal->Tr[row][col]);
      }
    }
  }
  for (int row = 0; row < h->deg; row++) {
    for (int col = 0; col < ((int)t - h->deg); col++) {
      mpz_set(ideal->Tr[row][col], Tr[((int)t - h->deg) * row + col]);
    }
  }
  ideal->log = log;
  ideal->ideal->r = r;
}

void ideal_u_set(ideal_u_ptr ideal_new, ideal_u_srcptr ideal_old,
    unsigned int t)
{
  if (ideal_new->ideal->r == 0) {
    ideal_new->Tr =
      (mpz_t ** ) malloc(sizeof(mpz_t * ) * ideal_old->ideal->h->deg);
    for (int row = 0; row < ideal_old->ideal->h->deg; row++) {
      ideal_new->Tr[row] =
        (mpz_t * ) malloc(sizeof(mpz_t) * ((int)t - ideal_old->ideal->h->deg));
      for (int col = 0; col < ((int)t - ideal_old->ideal->h->deg); col++) {
        mpz_init(ideal_new->Tr[row][col]);
      }
    }
  }
  for (int row = 0; row < ideal_old->ideal->h->deg; row++) {
    for (int col = 0; col < ((int)t - ideal_old->ideal->h->deg); col++) {
      mpz_set(ideal_new->Tr[row][col], ideal_old->Tr[row][col]);
    }
  }
  ideal_new->ideal->r = ideal_old->ideal->r;
  ideal_new->log = ideal_old->log;
  mpz_poly_set(ideal_new->ideal->h, ideal_old->ideal->h);
}

void ideal_u_clear(ideal_u_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg > 1 && ideal->ideal->r != 0);

  for (int row = 0; row < ideal->ideal->h->deg; row++) {
    for (int col = 0; col < ((int)t - ideal->ideal->h->deg); col++) {
      mpz_clear(ideal->Tr[row][col]);
    }
  }
  for (int row = 0; row < ideal->ideal->h->deg; row++) {
    free(ideal->Tr[row]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0.0;
}

void ideal_u_fprintf(FILE * file, ideal_u_srcptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg > 1);

  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %f\n", ideal->log);

  fprintf(file, "Tr: [");
  for (int row = 0; row < ideal->ideal->h->deg - 1; row++) {
    fprintf(file, "[");
    for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
      gmp_fprintf(file, "%Zd, ", ideal->Tr[row][col]);
    }
    gmp_fprintf(file, "%Zd],\n", ideal->Tr[row]
        [t - (unsigned int)ideal->ideal->h->deg - 1]);
  }
  fprintf(file, "[");
  for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[ideal->ideal->h->deg - 1][col]);
  }
  gmp_fprintf(file, "%Zd]]\n", ideal->Tr[ideal->ideal->h->deg - 1]
      [t - (unsigned int)ideal->ideal->h->deg - 1]);
}

/* Ideal projective root and Tr */

void ideal_pr_init(ideal_pr_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0.0;
}

void ideal_pr_set_part(ideal_pr_ptr ideal, uint64_t r, unsigned int t)
{
  mpz_poly h;
  mpz_poly_init(h, -1);
  mpz_poly_setcoeff_si(h, 1, 1);
  mpz_poly_set(ideal->ideal->h, h);
  mpz_poly_clear(h);

  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set_ui(ideal->Tr[i], 0);

  }
  ideal->ideal->r = r;
  ideal->log = log2((double) r);
}

void ideal_pr_set_element(ideal_pr_ptr ideal, uint64_t r, double log,
    unsigned int t)
{
  mpz_poly h;
  mpz_poly_init(h, -1);
  mpz_poly_setcoeff_si(h, 1, 1);
  mpz_poly_set(ideal->ideal->h, h);
  mpz_poly_clear(h);

  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set_ui(ideal->Tr[i], 0);
  }
  ideal->ideal->r = r;
  ideal->log = log;
}

void ideal_pr_set(ideal_pr_ptr ideal_new, ideal_pr_srcptr ideal_old,
    unsigned int t)
{
  mpz_poly_set(ideal_new->ideal->h, ideal_old->ideal->h);
  if (ideal_new->ideal->r == 0) {
    ideal_new->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal_new->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set(ideal_new->Tr[i], ideal_old->Tr[i]);
  }
  ideal_new->log = ideal_old->log;
  ideal_new->ideal->r = ideal_old->ideal->r;
}

void ideal_pr_clear(ideal_pr_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->r != 0);

  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(ideal->Tr[i]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0.0;
}

void ideal_pr_fprintf(FILE * file, ideal_pr_srcptr ideal, unsigned int t)
{
  fprintf(file, "Projective root, ");
  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %f\n", ideal->log);
  fprintf(file, "Tr: [");
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd]\n", ideal->Tr[t - 2]);
}

/* Ideal special-q. */

void ideal_spq_init(ideal_spq_ptr ideal)
{
  ideal->type = -1;
}

void ideal_spq_set_part(ideal_spq_ptr ideal, uint64_t q, mpz_poly_srcptr g,
    unsigned int t, int type)
{
  ASSERT(type >= 0 && type < 3);
#ifndef NDEBUG
  if (g->deg == 1) {
    ASSERT(type == 0 || type == 2);
  } else {
    ASSERT(type == 1);
  }
#endif // NDEBUG

  //TODO: clear only if type != ideal->type and then, init only if necessarily.
  ideal_spq_clear(ideal, t);

  ideal->type = type;
  if (type == 0) {
    ideal_1_init(ideal->ideal_1);
    ideal_1_set_part(ideal->ideal_1, q, g, t);
  } else if (type == 1) {
    ideal_u_init(ideal->ideal_u);
    ideal_u_set_part(ideal->ideal_u, q, g, t);
  } else {
    ASSERT(type == 2);
    ASSERT(g->deg == 1);
#ifndef NDEBUG
    mpz_t tmp;
    mpz_init(tmp);
    mpz_poly_getcoeff(tmp, 0, g);
    ASSERT(mpz_cmp_ui(tmp, 0) == 1);
    mpz_poly_getcoeff(tmp, 1, g);
    ASSERT(mpz_cmp_ui(tmp, 1) == 1);
    mpz_clear(tmp);
#endif // NDEBUG

    ideal_pr_init(ideal->ideal_pr);
    ideal_pr_set_part(ideal->ideal_pr, q, t);
  }
}

void ideal_spq_set(ideal_spq_ptr ideal, ideal_spq_srcptr ideal_old,
    unsigned int t)
{
  ASSERT(ideal_old->type >= 0 && ideal_old->type < 3);
  ASSERT(ideal->type == -1);

  //TODO: clear only if type != ideal->type and then, init only if necessarily.
  ideal_spq_clear(ideal, t);

  ideal->type = ideal_old->type;
  if (ideal_old->type == 0) {
    ideal_1_init(ideal->ideal_1);
    ideal_1_set(ideal->ideal_1, ideal_old->ideal_1, t);
  } else if (ideal_old->type == 1) {
    ideal_u_init(ideal->ideal_u);
    ideal_u_set(ideal->ideal_u, ideal_old->ideal_u, t);
  } else {
    ASSERT(ideal_old->type == 2);

    ideal_pr_init(ideal->ideal_pr);
    ideal_pr_set(ideal->ideal_pr, ideal_old->ideal_pr, t);
  }
}

void ideal_spq_clear(ideal_spq_ptr ideal, unsigned int t)
{
  if (ideal->type == 0) {
    ideal_1_clear(ideal->ideal_1, t);
  } else if (ideal->type == 1) {
    ideal_u_clear(ideal->ideal_u, t);
  } else if (ideal->type == 2) {
    ideal_pr_clear(ideal->ideal_pr, t);
  } else {
    ASSERT(ideal->type == -1);
  }

  ideal->type = -1;
}

void ideal_spq_fprintf(FILE * file, ideal_spq_srcptr ideal, unsigned int t)
{
  if (ideal->type == -1) {
    fprintf(file, "No special-q ideal.\n");
  }
  if (ideal->type == 0) {
    ideal_1_fprintf(file, ideal->ideal_1, t);
  } else if (ideal->type == 1) {
    ideal_u_fprintf(file, ideal->ideal_u, t);
  } else {
    ASSERT(ideal->type == 2);

    ideal_pr_fprintf(file, ideal->ideal_pr, t);
  }
}

void ideal_spq_fprintf_q_g(FILE * file, ideal_spq_srcptr ideal)
{
  if (ideal->type == -1) {
    fprintf(file, "No special-q ideal.\n");
  }
  if (ideal->type == 0) {
    ideal_fprintf(file, ideal->ideal_1->ideal);
  } else if (ideal->type == 1) {
    ideal_fprintf(file, ideal->ideal_u->ideal);
  } else {
    ASSERT(ideal->type == 2);

    ideal_fprintf(file, ideal->ideal_pr->ideal);
  }
}

double ideal_spq_get_log(ideal_spq_srcptr ideal)
{
  if (ideal->type == 0) {
    return ideal->ideal_1->log;
  } else if (ideal->type == 1) {
    return ideal->ideal_u->log;
  } else {
    ASSERT(ideal->type == 2);

    return ideal->ideal_pr->log;
  }
}

uint64_t ideal_spq_get_q(ideal_spq_srcptr ideal)
{
  uint64_t q = 0;
  if (ideal->type == 0) {
    q = ideal->ideal_1->ideal->r;
  } else if (ideal->type == 1) {
    q = ideal->ideal_u->ideal->r;
  } else {
    ASSERT(ideal->type == 2);

    q = ideal->ideal_pr->ideal->r;
  }
  return q;
}

int ideal_spq_get_deg_g(ideal_spq_srcptr ideal)
{
  int deg_g = 0;
  if (ideal->type == 0) {
    deg_g = ideal->ideal_1->ideal->h->deg;
  } else if (ideal->type == 1) {
    deg_g = ideal->ideal_u->ideal->h->deg;
  } else {
    ASSERT(ideal->type == 2);

    deg_g = ideal->ideal_pr->ideal->h->deg;
  }
  return deg_g;
}

void ideal_spq_get_g(mpz_poly_ptr g, ideal_spq_srcptr ideal)
{
  if (ideal->type == 0) {
    mpz_poly_set(g, ideal->ideal_1->ideal->h);
  } else if (ideal->type == 1) {
    mpz_poly_set(g, ideal->ideal_u->ideal->h);
  } else {
    ASSERT(ideal->type == 2);

    mpz_poly_set(g, ideal->ideal_pr->ideal->h);
  }
}
