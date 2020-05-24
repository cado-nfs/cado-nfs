#ifndef ROPT_STAGE1_H
#define ROPT_STAGE1_H

#include "ropt_str.h"    // ropt_param_t ropt_bound_t ...
#include "ropt_tree.h" // alpha_pq


#ifdef __cplusplus
extern "C" {
#endif

int ropt_stage1 ( ropt_poly_t poly,
                  ropt_bound_t bound,
                  ropt_s1param_t s1param,
                  ropt_param_t param,
                  alpha_pq *alpha_pqueue,
                  int current_w );


#ifdef __cplusplus
}
#endif


#endif /* ROPT_STAGE1_H */
