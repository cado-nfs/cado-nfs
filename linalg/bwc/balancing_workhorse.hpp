#ifndef BALANCING_WORKHORSE_HPP_
#define BALANCING_WORKHORSE_HPP_

#include "params.h"              // for param_list_ptr, param_list
#include "parallelizing_info.hpp"
#include "raw_matrix_u32.h"

void balancing_decl_usage(param_list_ptr pl);
void balancing_lookup_parameters(param_list_ptr pl);

void balancing_get_matrix_u32(parallelizing_info_ptr pi, param_list_ptr pl, matrix_u32_ptr arg);

#endif	/* BALANCING_WORKHORSE_HPP_ */
