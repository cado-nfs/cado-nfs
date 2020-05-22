#ifndef LINGEN_TUNING_HPP_
#define LINGEN_TUNING_HPP_

#include <map>
#include "params.h"
#include "lingen.hpp"
#include "select_mpi.h"
#include "timing.h"     /* for weighted_double */
#include "tree_stats.hpp"
#include "lingen_hints.hpp"

void lingen_tuning_decl_usage(cxx_param_list & pl);
void lingen_tuning_lookup_parameters(cxx_param_list & pl);
lingen_hints lingen_tuning(bw_dimensions & d, size_t, MPI_Comm comm, cxx_param_list & pl);

#endif	/* LINGEN_TUNING_HPP_ */
