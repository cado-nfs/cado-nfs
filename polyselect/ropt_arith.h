#ifndef CADO_ROPT_ARITH_H
#define CADO_ROPT_ARITH_H

#include <gmp.h>
#include "mpz_poly.h"
#include "ropt_str.h"

#define ROPT_NPRIMES 46

/* declarations */
#ifdef __cplusplus
extern "C" {
#endif

void compute_fuv_mp ( mpz_poly_ptr fuv,
                      mpz_poly_srcptr f,
                      mpz_poly_srcptr g,
                      mpz_srcptr u,
                      mpz_srcptr v );

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

void Lemma21 ( ropt_poly_ptr poly,
               mpz_t N,
               int d,
               mpz_srcptr ad,
               mpz_srcptr p,
               mpz_srcptr m,
               mpz_ptr res);

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


void reduce_poly_uint ( unsigned int *f_ui,
                      mpz_poly_srcptr f,
                      unsigned int pe );


void crt_pair_mp ( mpz_t a,
                   mpz_t p1,
                   mpz_t b,
                   mpz_t p2,
                   mpz_t re );


void ab2uv ( mpz_srcptr A,
             mpz_srcptr MOD,
             long a,
             mpz_ptr u );


long ab2ij ( long Amin, long a );


void ij2uv ( mpz_srcptr A,
             mpz_srcptr MOD,
             long Amin,
             long i,
             mpz_ptr u );


long uv2ij_mod ( mpz_srcptr A,
                 long Amin,
                 mpz_srcptr MOD,
                 unsigned int U,
                 unsigned int p );


#ifdef __cplusplus
}
#endif

#endif
