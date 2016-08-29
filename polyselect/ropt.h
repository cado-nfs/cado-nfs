#ifndef ROPT_H
#define ROPT_H

#include "cado.h"
#include "ropt_linear.h"
#include "ropt_quadratic.h"
#include "ropt_param.h"
#include "ropt_arith.h"
#include "ropt_io.h"
#include "ropt_str.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "portability.h"


/* timing structure for ropt */
typedef struct
{
  double ropt_time;
  double ropt_time_stage1;
  double ropt_time_tuning;
  double ropt_time_stage2;
} __ropt_time_struct;
typedef __ropt_time_struct ropt_time_t[1];


/* -- declarations -- */
void ropt ( ropt_poly_t poly,
            ropt_bestpoly_t bestpoly,
            ropt_param_t param,
            ropt_info_t info );

void ropt_get_bestpoly ( ropt_poly_t poly,
                         MurphyE_pq *global_E_pqueue,
                         ropt_bestpoly_t bestpoly,
                         ropt_param_t param );

void ropt_polyselect (cado_poly_ptr output_poly, cado_poly_ptr input_poly,
                      ropt_param_t param, ropt_time_t thr);

#endif /* ROPT_H */
