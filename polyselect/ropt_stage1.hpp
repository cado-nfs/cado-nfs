#ifndef CADO_ROPT_STAGE1_HPP
#define CADO_ROPT_STAGE1_HPP

#include "ropt_str.hpp"    // ropt_param_t ropt_bound_t ...
#include "ropt_tree.h" // alpha_pq


int ropt_stage1 ( ropt_poly const & poly,
                  ropt_bound_ptr bound,
                  ropt_s1param_ptr s1param,
                  ropt_param_srcptr param,
                  alpha_pq *alpha_pqueue,
                  int current_w );

#endif /* CADO_ROPT_STAGE1_HPP */
