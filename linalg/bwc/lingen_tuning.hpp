#ifndef LINGEN_TUNING_HPP_
#define LINGEN_TUNING_HPP_

#include <stddef.h>          // for size_t
#include "lingen_hints.hpp"  // for lingen_hints
#include "select_mpi.h"      // for MPI_Comm
struct bw_dimensions;
struct cxx_param_list;


void lingen_tuning_decl_usage(cxx_param_list & pl);
void lingen_tuning_lookup_parameters(cxx_param_list & pl);
lingen_hints lingen_tuning(bw_dimensions & d, size_t, MPI_Comm comm, cxx_param_list & pl);

#endif	/* LINGEN_TUNING_HPP_ */
