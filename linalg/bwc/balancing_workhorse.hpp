#ifndef BALANCING_WORKHORSE_HPP_
#define BALANCING_WORKHORSE_HPP_

#include "params.h"
#include "parallelizing_info.hpp"
#include "matrix_u32.hpp"

void balancing_decl_usage(cxx_param_list & pl);
void balancing_lookup_parameters(cxx_param_list & pl);

matrix_u32 balancing_get_matrix_u32(
        parallelizing_info_ptr pi,
        cxx_param_list & pl,
        std::string const & mfile,
        std::string const & bfile,
        bool withcoeffs,
        bool transpose_while_dispatching);

#endif	/* BALANCING_WORKHORSE_HPP_ */
