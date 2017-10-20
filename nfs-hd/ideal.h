#ifndef IDEAL_H
#define IDEAL_H

#include <stdint.h>
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Represent an ideal.
 */
typedef struct {
  uint64_t r;
  mpz_poly h;
} s_ideal_t;

typedef s_ideal_t ideal_t[1];
typedef s_ideal_t * ideal_ptr;
typedef const s_ideal_t * ideal_srcptr;

typedef struct {
  ideal_t ideal;
  mpz_t * Tr; // Keep mpz_t because of the explosion of g0^i.
  double log;
} s_ideal_1_t;

//1 for ideal (r, h) with deg(h) = 1.
typedef s_ideal_1_t ideal_1_t[1];
typedef s_ideal_1_t * ideal_1_ptr;
typedef const s_ideal_1_t * ideal_1_srcptr;

//u for upper deg (h) > 1.
typedef struct {
  ideal_t ideal;
  mpz_t ** Tr; // Keep mpz_t because of big coefficients.
  double log;
} s_ideal_u_t;

typedef s_ideal_u_t ideal_u_t[1];
typedef s_ideal_u_t * ideal_u_ptr;
typedef const s_ideal_u_t * ideal_u_srcptr;

//pr for projective root.
typedef struct {
  ideal_t ideal;
  mpz_t * Tr; // Keep mpz_t because of big coefficients.
  double log;
} s_ideal_pr_t;

typedef s_ideal_pr_t ideal_pr_t[1];
typedef s_ideal_pr_t * ideal_pr_ptr;
typedef const s_ideal_pr_t * ideal_pr_srcptr;

// spq for special-q.
typedef struct {
  ideal_1_t ideal_1;
  ideal_u_t ideal_u;
  ideal_pr_t ideal_pr;
  // 0 if ideal_1, 1 if ideal_u, 2 if ideal_pr.
  int type;
} s_ideal_spq_t;

typedef s_ideal_spq_t ideal_spq_t[1];
typedef s_ideal_spq_t * ideal_spq_ptr;
typedef const s_ideal_spq_t * ideal_spq_srcptr;

/* Ideal: (r, h) */

/*
 * Initialize an ideal.
 *
 * ideal: an ideal.
 */
void ideal_init(ideal_ptr ideal);

/*
 * Set an ideal. Compute Tr and log.
 *
 * ideal: the ideal we want to set.
 * r: r is equal to the r of the ideal.
 * h: h is equal to the h of the ideal.
 */
void ideal_set_part(ideal_ptr ideal, uint64_t r, mpz_poly_srcptr h);

/*
 * Delete an ideal.
 *
 * ideal: the ideal we want to delete.
 */
void ideal_clear(ideal_ptr ideal);

/*
 * Write an ideal in a file.
 *
 * file: the file.
 * ideal: the ideal we want to write.
 */
void ideal_fprintf(FILE * file, ideal_srcptr ideal);

/* Ideal with h of degree 1 and Tr */

/*
 * Initialize an ideal_1.
 *
 * ideal: an ideal.
 */
void ideal_1_init(ideal_1_ptr ideal);

/*
 * Set an ideal_1, compute Tr and the logarithm.
 *
 * ideal: the ideal we want to set.
 * r: r is equal to the r of the ideal.
 * h: h is equal to the h of the ideal.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_1_set_part(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    unsigned int t);

/*
 * Set an ideal.
 *
 * ideal: the ideal we want to set.
 * r: r of the ideal.
 * h: h of the ideal.
 * Tr: Tr matrix.
 * log: frequently equal to log2(r).
 */
void ideal_1_set_element(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    mpz_t * Tr, double log, unsigned int t);

/*
 * Copy an ideal_1 in an other.
 *
 * ideal_new: the new ideal.
 * ideal_old: the old ideal.
 */
void ideal_1_set(ideal_1_ptr ideal_new, ideal_1_srcptr ideal_old,
    unsigned int t);

/*
 * Delete an ideal_1.
 *
 * ideal: the ideal we want to delete.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_1_clear(ideal_1_ptr ideal, unsigned int t);

/*
 * Write an ideal_1.
 *
 * ideal: the ideal we want to delete.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_1_fprintf(FILE * file, ideal_1_srcptr ideal, unsigned int t);

/* Ideal with h of degree u > 1 and Tr */

/*
 * Initialize an ideal_u.
 *
 * ideal: an ideal.
 */
void ideal_u_init(ideal_u_ptr ideal);

/*
 * Set an ideal_u, compute Tr and the logarithm.
 *
 * ideal: the ideal we want to set.
 * r: r is equal to the r of the ideal.
 * h: h is equal to the h of the ideal.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_u_set_part(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
                      unsigned int t);

/*
 * Set an ideal.
 *
 * ideal: the ideal we want to set.
 * r: r of the ideal.
 * h: h of the ideal.
 * Tr: Tr matrix.
 * log: frequently equal to log2(r^deg(h)).
 */
void ideal_u_set_element(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    mpz_t * Tr, double log, unsigned int t);

/*
 * Copy ideal_old in ideal_new.
 *
 * ideal_new: the new ideal.
 * ideal_old: the old ideal.
 * t: dimension of the lattice.
 */
void ideal_u_set(ideal_u_ptr ideal_new, ideal_u_srcptr ideal_old,
    unsigned int t);

/*
 * Delete an ideal_u.
 *
 * ideal: ideal_u_srcptr, the ideal we want to delete.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_u_clear(ideal_u_ptr ideal, unsigned int t);

/*
 * Write an ideal_u.
 *
 * file: the file in which we write.
 * ideal: ideal_u_srcptr, the ideal we want to write.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_u_fprintf(FILE * file, ideal_u_srcptr ideal, unsigned int t);

/* Ideal projective root and Tr */

/*
 * Init an ideal_pr.
 *
 * ideal: the ideal we initialize.
 */
void ideal_pr_init(ideal_pr_ptr ideal);

/*
 * Set an ideal_pr, compute Tr and the logarithm. The polynomial is x.
 *
 * ideal: the ideal we want to set.
 * r: r is equal to the r of the ideal.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_pr_set_part(ideal_pr_ptr ideal, uint64_t r, unsigned int t);

/*
 * Set an ideal_pr.
 *
 * ideal: the ideal we want to set.
 * r: r of the ideal.
 * Tr: Tr matrix.
 * log: frequently equal to log2(r).
 */
void ideal_pr_set_element(ideal_pr_ptr ideal, uint64_t r, double log,
    unsigned int t);

void ideal_pr_set(ideal_pr_ptr ideal_new, ideal_pr_srcptr ideal_old,
    unsigned int t);
/*
 * Delete an ideal_pr.
 *
 * ideal: ideal_pr_srcptr, the ideal we want to delete.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_pr_clear(ideal_pr_ptr ideal, unsigned int t);

/*
 * Write an ideal_pr.
 *
 * file: the file in which we write.
 * ideal: ideal_pr_srcptr, the ideal we want to write.
 * t: t is the dimension of the lattice we explore.
 */
void ideal_pr_fprintf(FILE * file, ideal_pr_srcptr ideal, unsigned int t);

/* Special-q ideal: (q, g) */
/* type: 0 if ideal_1, 1 if ideal_u, 2 if ideal_pr. */

/*
 * Init an ideal_spq.
 *
 * ideal: the ideal.
 */
void ideal_spq_init(ideal_spq_ptr ideal);

/*
 * Set an ideal spq.
 *
 * r: r of the ideal.
 * h: h of the ideal.
 * t: dimension of the lattice.
 * type: type of the ideal_spq.
 */
void ideal_spq_set_part(ideal_spq_ptr ideal, uint64_t q, mpz_poly_srcptr g,
    unsigned int t, int type);

void ideal_spq_set(ideal_spq_ptr ideal, ideal_spq_srcptr ideal_old, unsigned int t);

/*
 * Clear an ideal_spq.
 *
 * ideal: the ideal.
 * t: dimension of the lattice.
 */
void ideal_spq_clear(ideal_spq_ptr ideal, unsigned int t);

/*
 * Write an ideal_spq.
 *
 * file: file in which we print.
 * ideal: the ideal.
 * t: dimension of the lattice.
 */
void ideal_spq_fprintf(FILE * file, ideal_spq_srcptr ideal, unsigned int t);

void ideal_spq_fprintf_q_g(FILE * file, ideal_spq_srcptr ideal);

/*
 * Return the log of the ideal_spq.
 *
 * ideal: the ideal.
 */
double ideal_spq_get_log(ideal_spq_srcptr ideal);

/*
 * Return q of an ideal_spq.
 *
 * ideal: the ideal.
 */
uint64_t ideal_spq_get_q(ideal_spq_srcptr ideal);

/*
 * Return g of an ideal_spq.
 *
 * ideal: the ideal.
 */
int ideal_spq_get_deg_g(ideal_spq_srcptr ideal);

/*
 * Get g of an ideal_spq.
 */
void ideal_spq_get_g(mpz_poly_ptr g, ideal_spq_srcptr ideal);

#ifdef __cplusplus
}
#endif

#endif /* IDEAL_H */
