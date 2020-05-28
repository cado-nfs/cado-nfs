#ifndef BALANCING_WORKHORSE_H_
#define BALANCING_WORKHORSE_H_

#include "params.h"              // for param_list_ptr, param_list
#include "parallelizing_info.h"
#include "raw_matrix_u32.h"

#ifdef __cplusplus
extern "C" {
#endif

void balancing_decl_usage(param_list_ptr pl);
void balancing_lookup_parameters(param_list_ptr pl);

void balancing_get_matrix_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg);


#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_WORKHORSE_H_ */
