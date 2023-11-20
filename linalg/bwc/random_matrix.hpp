#ifndef RANDOM_MATRIX_HPP_
#define RANDOM_MATRIX_HPP_

#include "parallelizing_info.hpp"
#include "raw_matrix_u32.h"
#include "balancing.hpp"
#include "params.h"

/* fills arg with simulated data, which should correspond to the inner
 * blocks of a matrix split according to pi. The parameter list pl is
 * used to read the random_matrix= parameter, which corresponds to the
 * full matrix characteristics: the format is a comma separated list of
 * settings, which then get parsed as a parameter list in the same manner
 * as for the standalone "random_matrix" program. Therefore, the easiest
 * way to state this argument is for example
 * random_matrix=2000,density=4,seed=1
 */
void random_matrix_get_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg, unsigned long data_nrows, unsigned long data_ncols, unsigned long padded_nrows, unsigned long padded_ncols);

void random_matrix_fill_fake_balancing_header(balancing & bal, parallelizing_info_ptr pi, const char * rtmp);

#endif	/* RANDOM_MATRIX_HPP_ */
