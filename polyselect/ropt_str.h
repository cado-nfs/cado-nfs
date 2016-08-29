#ifndef ROPT_STR_H
#define ROPT_STR_H


#include "ropt_arith.h"
#include "ropt_param.h"


/* --- structs for ropt --- */


/**
 * Struct for the polynomial currently being ropt-ed.
 */
typedef struct {
  /* polynomial */
  int d;
  double skew;
  double alpha_proj;
  mpz_t n;
  mpz_t m;
  mpz_t *f;
  mpz_t *g;

  /* polynomial values */
  mpz_t *fx;
  mpz_t *gx;
  mpz_t *numerator;
} _ropt_poly_t; 
typedef _ropt_poly_t ropt_poly_t[1];


/**
 * Struct for size bound.
 */
typedef struct {
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
} _ropt_bound_t;
typedef _ropt_bound_t ropt_bound_t[1];


/**
 * Struct for stage 1 parameters.
 */
typedef struct {
  unsigned int len_e_sl;
  unsigned int tlen_e_sl;

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

  mpz_t modulus;

  unsigned long modbound;
} _ropt_s1param_t;
typedef _ropt_s1param_t ropt_s1param_t[1];


/**
 * Struct for stage 2 parameters:
 * sieve bound for (A + MOD*i)*x + (B + MOD*j).
 */
typedef struct {
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
  mpz_t *f;
  mpz_t *g;
} _ropt_s2param_t;
typedef _ropt_s2param_t ropt_s2param_t[1];


/**
 * Struct for manually-input parameters:
 * it has miscellaneous parameters filled from stdin.
 */
typedef struct {
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
} _ropt_param_t;
typedef _ropt_param_t ropt_param_t[1];


/**
 * Struct for top polynomials.
 */
typedef struct {
  mpz_t *f;
  mpz_t *g;
} _ropt_bestpoly_t;
typedef _ropt_bestpoly_t ropt_bestpoly_t[1];


/**
 * Struct for info (e.g. miscellaneous and for extension purpose).
 */
#define ROPT_MODE_INIT 0
#define ROPT_MODE_TUNE 1
typedef struct {
  double ave_MurphyE;
  double best_MurphyE;
  int mode; /* ROPT_MODE_INIT or ROPT_MODE_TUNE */
  /* record quadratic rotation information */
  int w; 
  double ropt_time_stage1;
  double ropt_time_tuning;
  double ropt_time_stage2;
} _ropt_info_t;
typedef _ropt_info_t ropt_info_t[1];


/* --- declarations --- */


/* ropt_poly_t */
void ropt_poly_init ( ropt_poly_t );

void ropt_poly_set ( ropt_poly_t curr_poly, ropt_poly_t poly );

void ropt_poly_refresh ( ropt_poly_t );

void ropt_poly_setup ( ropt_poly_t );

bool ropt_poly_setup_check ( ropt_poly_t );

void ropt_poly_free ( ropt_poly_t );


/* ropt_bound_t */
void ropt_bound_init ( ropt_bound_t );

void ropt_bound_setup ( ropt_poly_t poly,
                        ropt_bound_t bound,
                        ropt_param_t param,
                        double incr );

void ropt_bound_setup_incr ( ropt_poly_t poly,
                             ropt_bound_t bound,
                             ropt_param_t param,
                             double incr );

#if 0
double ropt_bound_expected_E (mpz_poly F, mpz_poly G);
#endif

void ropt_bound_reset ( ropt_poly_t poly,
                        ropt_bound_t bound,
                        ropt_param_t param );

void ropt_bound_free ( ropt_bound_t bound );


/* ropt_s1param_t */
void ropt_s1param_init ( ropt_s1param_t s1param);

void ropt_s1param_setup_individual_nbest_sl (ropt_s1param_t s1param);

void ropt_s1param_setup_individual_nbest_sl_tune (ropt_s1param_t s1param);

void ropt_s1param_setup ( ropt_poly_t poly,
                          ropt_s1param_t s1param,
                          ropt_bound_t bound,
                          ropt_param_t param);

void ropt_s1param_resetup ( ropt_poly_t poly,
                            ropt_s1param_t s1param,
                            ropt_bound_t bound,
                            ropt_param_t param,
                            unsigned int nbest );

void
ropt_s1param_resetup_modbound ( ropt_poly_t poly,
                                ropt_s1param_t s1param,
                                ropt_bound_t bound,
                                ropt_param_t param,
                                unsigned int nbest,
                                unsigned long modbound);

void ropt_s1param_free ( ropt_s1param_t s1param );


/* ropt_s2param_t */
void ropt_s2param_init ( ropt_poly_t poly,
                         ropt_s2param_t s2param );

void ropt_s2param_free ( ropt_poly_t poly,
                         ropt_s2param_t s2param );

void ropt_s2param_setup ( ropt_bound_t bound,
                          ropt_s1param_t s1param,
                          ropt_s2param_t s2param,
                          ropt_param_t param,
                          mpz_t A,
                          mpz_t B,
                          mpz_t MOD );

void ropt_s2param_setup_stage2_only ( ropt_bound_t bound,
                                      ropt_s2param_t s2param,
                                      ropt_param_t param,
                                      mpz_t A,
                                      mpz_t B,
                                      mpz_t MOD );

void ropt_s2param_setup_tune ( ropt_s1param_t s1param,
                               ropt_s2param_t s2param,
                               mpz_t A,
                               mpz_t B,
                               mpz_t MOD,
                               unsigned long Amax,
                               unsigned long Bmax,
                               unsigned int len_p_rs );

void ropt_s2param_print ( ropt_s2param_t s2param );


/* bestpoly */
void ropt_bestpoly_init ( ropt_bestpoly_t bestpoly,
                          int d );
                      
void ropt_bestpoly_setup ( ropt_bestpoly_t bestpoly,
                           mpz_t *f,
                           mpz_t *g,
                           int d );

void ropt_bestpoly_free ( ropt_bestpoly_t bestpoly,
                          int d );


/* param */
void ropt_param_init ( ropt_param_t param );

void ropt_param_free ( ropt_param_t param );


/* info */
void ropt_info_init ( ropt_info_t info );

void ropt_info_set ( ropt_info_t curr_info,  ropt_info_t info );
  
void ropt_info_free ( ropt_info_t info );

#endif
