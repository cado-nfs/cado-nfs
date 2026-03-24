#ifndef CADO_ROPT_STR_H
#define CADO_ROPT_STR_H
#include <gmp.h>
#include <stdbool.h>    // for bool (in C)
#include "cado_poly.h"
#include "mpz_poly.h"

/* --- structs for ropt --- */

/**
 * Struct for the polynomial currently being ropt-ed.
 *
 * TODO: this should at least contain a cado_poly. The only reason not to
 * do it would be to save a pointer indirection for the pols member.
 */
struct ropt_poly_s {
  /* polynomial */
  double alpha_proj;
  mpz_t m;
  cado_poly cpoly;

  /* polynomial values */
  mpz_t *fx;
  mpz_t *gx;
  mpz_t *numerator;
}; 
typedef struct ropt_poly_s ropt_poly[1];
typedef struct ropt_poly_s * ropt_poly_ptr;
typedef const struct ropt_poly_s * ropt_poly_srcptr;


/**
 * Struct for size bound.
 */
struct ropt_bound_s {
  /* global bounds */
  long global_w_boundl;
  long global_w_boundr;
  long global_u_boundl;
  long global_u_boundr;
  mpz_t global_v_boundl;
  mpz_t global_v_boundr;

  /* used to decide global bounds */
  double init_lognorm;
  double bound_lognorm;
  double bound_E;
  double exp_min_alpha;
};
typedef struct ropt_bound_s ropt_bound[1];
typedef struct ropt_bound_s * ropt_bound_ptr;
typedef const struct ropt_bound_s * ropt_bound_srcptr;


/**
 * Struct for stage 1 parameters.
 */
struct ropt_s1param_s {
  unsigned int len_e_sl;   /* length */
  unsigned int tlen_e_sl;  /* true length (excluding zeros), <= len_e_sl */

  /* num of sublattices */
  unsigned int nbest_sl;

  /* num of sublattices for tune */
  unsigned int nbest_sl_tune;
  
  /* num of all sublattices */
  unsigned int nbest_sieve;

  /* tune mode 1 */
  unsigned int nbest_sl_tunemode;

  /* num of individual sublattices: for each p_i^{e_i} */
  unsigned int *individual_nbest_sl;

  /* upper bound for e_i */
  unsigned int *e_sl;

  // mpz_t modulus;

  unsigned long modbound;
};
typedef struct ropt_s1param_s ropt_s1param[1];
typedef struct ropt_s1param_s * ropt_s1param_ptr;
typedef const struct ropt_s1param_s * ropt_s1param_srcptr;

/**
 * Struct for stage 2 parameters:
 * sieve bound for (A + MOD*i)*x + (B + MOD*j).
 */
struct ropt_s2param_s {
  /* maybe different to global_*_bound */
  mpz_t Umax;
  mpz_t Umin;
  mpz_t Vmax;
  mpz_t Vmin;

  /* sieve length */
  long Amax;
  long Amin;
  long Bmax;
  long Bmin;
  long Blocksize;

  /* sublattice */
  mpz_t A;
  mpz_t B;
  mpz_t MOD;
  unsigned int len_p_rs;

  /* polynomial */
  mpz_poly f;
  mpz_poly g;
};
typedef struct ropt_s2param_s ropt_s2param[1];
typedef struct ropt_s2param_s * ropt_s2param_ptr;
typedef const struct ropt_s2param_s * ropt_s2param_srcptr;


/**
 * Struct for manually-input parameters:
 * it has miscellaneous parameters filled from stdin.
 */
struct ropt_param_s {
  /* for msieve format, we need n, d */
  mpz_t n;
  int d;

  /* quadratic rotation bound */
  int w_left_bound;
  int w_length;

  /* stage 1 parameters */
  unsigned short s1_num_e_sl;
  unsigned short *s1_e_sl;

  /* sublattice parameters */
  int s2_w;
  mpz_t s2_u;
  mpz_t s2_v;
  mpz_t s2_mod;

  /* sieve bound parameters */
  long s2_Amax;
  long s2_Bmax;

  /* bound */
  double bound_lognorm;

  /* flag == 0, only use -w and -l
     flag == 1, manually input stage 1 params
     flag == 2, manually input stage 2 params */
  int stage_flag;

  /* rootsieve effort (default is DEFAULT_ROPTEFFORT) */
  double effort;

  /* skip ropt */
  int skip_ropt;

  /* skip ropt and output polynomial with short c5 */
  int gen_raw;

  /* do sopt if true; default not */
  int sopt;

  /* verbose */
  int verbose;
};
typedef struct ropt_param_s ropt_param[1];
typedef struct ropt_param_s * ropt_param_ptr;
typedef const struct ropt_param_s * ropt_param_srcptr;


/**
 * Struct for top polynomials.
 */
struct ropt_bestpoly_s {
  mpz_poly f;
  mpz_poly g;
};
typedef struct ropt_bestpoly_s ropt_bestpoly[1];
typedef struct ropt_bestpoly_s * ropt_bestpoly_ptr;
typedef const struct ropt_bestpoly_s * ropt_bestpoly_srcptr;


/**
 * Struct for info (e.g. miscellaneous and for extension purpose).
 */
#define ROPT_MODE_INIT 0
#define ROPT_MODE_TUNE 1
struct ropt_info_s {
  double ave_MurphyE;
  double best_MurphyE;
  int mode; /* ROPT_MODE_INIT or ROPT_MODE_TUNE */
  /* record quadratic rotation information */
  int w; 
  double ropt_time_stage1;
  double ropt_time_tuning;
  double ropt_time_stage2;
};
typedef struct ropt_info_s ropt_info[1];
typedef struct ropt_info_s * ropt_info_ptr;
typedef const struct ropt_info_s * ropt_info_srcptr;


/* --- declarations --- */


#ifdef __cplusplus
extern "C" {
#endif

/* ropt_poly_t */
void ropt_poly_init ( ropt_poly_ptr );

void ropt_poly_refresh ( ropt_poly_ptr );

void ropt_poly_setup ( ropt_poly_ptr );

bool ropt_poly_setup_check ( ropt_poly_ptr );

void ropt_poly_clear ( ropt_poly_ptr );


/* ropt_bound_t */
void ropt_bound_init ( ropt_bound_ptr );

void ropt_bound_setup ( ropt_poly_srcptr poly,
                        ropt_bound_ptr bound,
                        ropt_param_ptr param,
                        double incr );

void ropt_bound_setup_incr ( ropt_poly_srcptr poly,
                             ropt_bound_ptr bound,
                             ropt_param_ptr param,
                             double incr );

#if 0
double ropt_bound_expected_E (mpz_poly F, mpz_poly G);
#endif

void ropt_bound_reset ( ropt_poly_srcptr poly,
                        ropt_bound_ptr bound,
                        ropt_param_srcptr param );

void ropt_bound_clear ( ropt_bound_ptr bound );


/* ropt_s1param */
void ropt_s1param_init ( ropt_s1param_ptr s1param);

void ropt_s1param_setup_individual_nbest_sl (ropt_s1param_ptr s1param);

void ropt_s1param_setup_individual_nbest_sl_tune (ropt_s1param_ptr s1param);

void ropt_s1param_setup ( ropt_poly_srcptr poly,
                          ropt_s1param_ptr s1param,
                          ropt_bound_srcptr bound,
                          ropt_param_srcptr param);

void ropt_s1param_resetup ( ropt_poly_srcptr poly,
                            ropt_s1param_ptr s1param,
                            ropt_bound_srcptr bound,
                            ropt_param_srcptr param,
                            unsigned int nbest );

void
ropt_s1param_resetup_modbound ( ropt_poly_srcptr poly,
                                ropt_s1param_ptr s1param,
                                ropt_bound_srcptr bound,
                                ropt_param_srcptr param,
                                unsigned int nbest,
                                unsigned long modbound);

void ropt_s1param_clear ( ropt_s1param_ptr s1param );


/* ropt_s2param_t */
void ropt_s2param_init ( ropt_s2param_ptr s2param );

void ropt_s2param_clear (ropt_s2param_ptr s2param );

void ropt_s2param_setup ( ropt_bound_srcptr bound,
                          ropt_s1param_srcptr s1param,
                          ropt_s2param_ptr s2param,
                          ropt_param_srcptr param,
                          mpz_srcptr A,
                          mpz_srcptr B,
                          mpz_srcptr MOD );

void ropt_s2param_setup_stage2_only ( ropt_bound_srcptr bound,
                                      ropt_s2param_ptr s2param,
                                      ropt_param_srcptr param,
                                      mpz_srcptr A,
                                      mpz_srcptr B,
                                      mpz_srcptr MOD );

void ropt_s2param_setup_tune ( ropt_s1param_srcptr s1param,
                               ropt_s2param_ptr s2param,
                               mpz_srcptr A,
                               mpz_srcptr B,
                               mpz_srcptr MOD,
                               unsigned long Amax,
                               unsigned long Bmax,
                               unsigned int len_p_rs );

void ropt_s2param_print ( ropt_s2param_srcptr s2param );


/* bestpoly */
void ropt_bestpoly_init ( ropt_bestpoly_ptr bestpoly);
                      
void ropt_bestpoly_set ( ropt_bestpoly_ptr bestpoly,
                           mpz_poly_srcptr f,
                           mpz_poly_srcptr g);

void ropt_bestpoly_clear ( ropt_bestpoly_ptr bestpoly);


/* param */
void ropt_param_init ( ropt_param_ptr param );

void ropt_param_clear ( ropt_param_ptr param );


/* info */
void ropt_info_init ( ropt_info_ptr info );

void ropt_info_clear ( ropt_info_ptr info );


#ifdef __cplusplus
}
#endif

#endif
