#ifndef ROPT_ARITH_H
#define ROPT_ARITH_H

#include <gmp.h>

#define ROPT_NPRIMES 46

/* declarations */
#ifdef __cplusplus
extern "C" {
#endif

void compute_fuv_mp ( mpz_t *fuv,
                      mpz_t *f,
                      mpz_t *g,
                      int d,
                      mpz_t u,
                      mpz_t v );

void compute_fuv_ui ( unsigned int *fuv_ui,
                      unsigned int *f_ui,
                      unsigned int *g_ui,
                      int d,
                      unsigned int u,
                      unsigned int v,
                      unsigned int p );

unsigned int eval_poly_ui_mod ( unsigned int *f,
                                int d,
                                unsigned int r,
                                unsigned int pe );

void Lemma21 ( mpz_t *a,
               mpz_t N,
               int d,
               mpz_t p,
               mpz_t m,
               mpz_t res );

void eval_polys ( mpz_t *f,
                  mpz_t *g,
                  mpz_t *fr,
                  mpz_t *gr,
                  mpz_t *numerator,
                  const unsigned int *primes,
                  int d );

unsigned long
solve_lineq ( unsigned long a,
              unsigned long b,
              unsigned long c,
              unsigned long p );

unsigned int compute_v_ui ( unsigned int fx,
                            unsigned int gx,
                            unsigned int r,
                            unsigned int u,
                            unsigned int p);


void reduce_poly_ul ( unsigned int *f_ui,
                      mpz_t *f,
                      int d,
                      unsigned int pe );


void crt_pair_mp ( mpz_t a,
                   mpz_t p1,
                   mpz_t b,
                   mpz_t p2,
                   mpz_t re );


void ab2uv ( mpz_t A,
             mpz_t MOD,
             long a,
             mpz_t u );


long ab2ij ( long Amin,
             long a );


void ij2uv ( mpz_t A,
             mpz_t MOD,
             long Amin,
             long i,
             mpz_t u );


long uv2ij_mod ( mpz_t A,
                 long Amin,
                 mpz_t MOD,
                 unsigned int U,
                 unsigned int p );


#ifdef __cplusplus
}
#endif

#endif
