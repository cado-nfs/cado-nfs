#ifndef ROPT_QUADRATIC_H
#define ROPT_QUADRATIC_H

#include "ropt_str.h"   // ropt_bestpoly_t ropt_poly_t ropt_param_t ropt_info_t
#include "gmp_aux.h"           // for gmp_randstate_ptr

/* -- declarations -- */

#ifdef __cplusplus
extern "C" {
#endif

void ropt_quadratic ( ropt_poly_t rs,
                      ropt_bestpoly_t bestpoly,
                      ropt_param_t param,
                      ropt_info_t info);

#ifdef __cplusplus
}
#endif



#endif /* ROPT_QUADRATIC_H */
