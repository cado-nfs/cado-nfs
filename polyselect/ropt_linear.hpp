#ifndef CADO_ROPT_LINEAR_H
#define CADO_ROPT_LINEAR_H

#include "ropt_str.hpp"
#include "ropt_tree.h"
#include "ropt_param.h"
#include "gmp_aux.h"

/* -- declarations -- */

void ropt_linear ( ropt_poly & poly,
                   ropt_bestpoly_ptr bestpoly,
                   ropt_param_ptr param,
                   ropt_info_ptr info);

double
ropt_tune_stage2_fast ( ropt_poly & poly,
                        ropt_s1param_srcptr s1param,
                        ropt_param_srcptr param,
                        ropt_info_ptr info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune );

void
ropt_tune_stage2_slow ( ropt_poly & poly,
                        ropt_bound_srcptr bound,
                        ropt_s1param_ptr s1param,
                        ropt_param_srcptr param,
                        ropt_info_ptr info,
                        alpha_pq *alpha_pqueue,
                        MurphyE_pq *global_E_pqueue,
                        unsigned int curr_size_tune,
                        int tune_round,
                        unsigned int curr_nbest );

void
ropt_tune_stage2 ( ropt_poly & poly,
                   ropt_bound_srcptr bound,
                   ropt_s1param_ptr s1param,
                   ropt_param_srcptr param,
                   ropt_info_ptr info,
                   alpha_pq *alpha_pqueue,
#if TUNE_LOGNORM_INCR
                   alpha_pq *tune_E_pqueue,
#endif
                   MurphyE_pq *global_E_pqueue );


#if TUNE_LOGNORM_INCR
double
ropt_linear_tune_stage1 ( ropt_poly & poly,
                          ropt_s1param_srcptr s1param,
                          ropt_param_ptr param,
                          alpha_pq *tune_E_pqueue,
                          alpha_pq *alpha_pqueue,
                          ropt_info_ptr info,
                          MurphyE_pq *global_E_pqueue,
                          int quad);
#endif

void
ropt_call_sieve ( ropt_poly & poly,
                  ropt_bound_srcptr bound,
                  ropt_s1param_srcptr s1param,
                  ropt_param_srcptr param,
                  ropt_info_ptr info,
                  alpha_pq *alpha_pqueue,
                  MurphyE_pq *global_E_pqueue );

void
ropt_MurphyE_to_alpha ( MurphyE_pq *E_pqueue,
                        alpha_pq *alpha_pqueue );


#endif /* CADO_ROPT_LINEAR_H */
