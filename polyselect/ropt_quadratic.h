#ifndef CADO_ROPT_QUADRATIC_H
#define CADO_ROPT_QUADRATIC_H

#include "ropt_str.h"   // ropt_bestpoly_t ropt_poly_t ropt_param_t ropt_info_t
#include "gmp_aux.h"           // for gmp_randstate_ptr

/* -- declarations -- */

#ifdef __cplusplus
extern "C" {
#endif

void ropt_quadratic ( ropt_poly_ptr rs,
                      ropt_bestpoly_ptr bestpoly,
                      ropt_param_ptr param,
                      ropt_info_ptr info);

#ifdef __cplusplus
}
#endif



#endif /* CADO_ROPT_QUADRATIC_H */
