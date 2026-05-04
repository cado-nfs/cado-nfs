#ifndef CADO_ROPT_HPP
#define CADO_ROPT_HPP

#include "cado_poly.hpp"
#include "ropt_str.hpp"
#include "ropt_tree.h"  // MurphyE_pq
#include "gmp_aux.h"           // for gmp_randstate_ptr


/* timing structure for ropt */
struct ropt_sime_struct {
  double ropt_time;
  double ropt_time_stage1;
  double ropt_time_tuning;
  double ropt_time_stage2;
};
typedef struct ropt_sime_struct ropt_time[1];
typedef struct ropt_sime_struct * ropt_time_ptr;
typedef const struct ropt_sime_struct * ropt_time_srcptr;


/* -- declarations -- */
void ropt ( ropt_poly & poly,
            ropt_bestpoly & bestpoly,
            ropt_param_ptr param,
            ropt_info_ptr info);

void ropt_get_bestpoly ( ropt_poly const & poly,
                         MurphyE_pq *global_E_pqueue,
                         ropt_bestpoly & bestpoly);

void ropt_polyselect (cxx_cado_poly & output_poly, cxx_cado_poly const & input_poly,
                      ropt_param_ptr param, ropt_time_ptr thr);


#endif /* CADO_ROPT_HPP */
