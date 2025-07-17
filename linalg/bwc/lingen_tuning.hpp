#ifndef CADO_LINGEN_TUNING_HPP
#define CADO_LINGEN_TUNING_HPP

#include <cstddef>

#include "lingen_hints.hpp"
#include "select_mpi.h"

template<bool is_binary>
struct bw_dimensions;

struct cxx_param_list;

void lingen_tuning_decl_usage(cxx_param_list & pl);
void lingen_tuning_lookup_parameters(cxx_param_list & pl);

template<bool is_binary>
lingen_hints lingen_tuning(bw_dimensions<is_binary> & d, size_t, MPI_Comm comm, cxx_param_list & pl);

#endif	/* LINGEN_TUNING_HPP_ */
