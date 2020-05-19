#ifndef ROPT_STAGE1_H
#define ROPT_STAGE1_H

#include "ropt_param.h"
#include "ropt_tree.h"
#include "ropt_arith.h"


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
