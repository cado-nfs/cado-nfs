#ifndef CADO_ROPT_QUADRATIC_H
#define CADO_ROPT_QUADRATIC_H

#include "ropt_str.hpp"   // ropt_bestpoly_t ropt_poly_t ropt_param_t ropt_info_t
#include "gmp_aux.h"           // for gmp_randstate_ptr

/* -- declarations -- */

void ropt_quadratic ( ropt_poly & rs,
                      ropt_bestpoly & bestpoly,
                      ropt_param_ptr param,
                      ropt_info_ptr info);



#endif /* CADO_ROPT_QUADRATIC_H */
