#ifndef CADO_ROPT_STAGE1_H
#define CADO_ROPT_STAGE1_H

#include "ropt_str.h"    // ropt_param_t ropt_bound_t ...
#include "ropt_tree.h" // alpha_pq


#ifdef __cplusplus
extern "C" {
#endif

int ropt_stage1 ( ropt_poly_srcptr poly,
                  ropt_bound_ptr bound,
                  ropt_s1param_ptr s1param,
                  ropt_param_srcptr param,
                  alpha_pq *alpha_pqueue,
                  int current_w );


#ifdef __cplusplus
}
#endif


#endif /* CADO_ROPT_STAGE1_H */
