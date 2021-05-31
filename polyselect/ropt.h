#ifndef ROPT_H
#define ROPT_H

#include "cado_poly.h"
#include "ropt_str.h"
#include "ropt_tree.h"  // MurphyE_pq
#include "gmp_aux.h"           // for gmp_randstate_ptr


/* timing structure for ropt */
struct ropt_sime_struct {
  double ropt_time;
  double ropt_time_stage1;
  double ropt_time_tuning;
  double ropt_time_stage2;
};
typedef struct ropt_sime_struct ropt_time_t[1];


/* -- declarations -- */
#ifdef __cplusplus
extern "C" {
#endif

void ropt ( ropt_poly_t poly,
            ropt_bestpoly_t bestpoly,
            ropt_param_t param,
            ropt_info_t info);

void ropt_get_bestpoly ( ropt_poly_t poly,
                         MurphyE_pq *global_E_pqueue,
                         ropt_bestpoly_t bestpoly);

void ropt_polyselect (cado_poly_ptr output_poly, cado_poly_ptr input_poly,
                      ropt_param_t param, ropt_time_t thr);

#ifdef __cplusplus
}
#endif


#endif /* ROPT_H */
