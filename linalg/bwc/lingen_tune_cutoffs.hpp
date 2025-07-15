#ifndef CADO_LINGEN_TUNE_CUTOFFS_HPP
#define CADO_LINGEN_TUNE_CUTOFFS_HPP

#include "select_mpi.h"

template<bool is_binary>
struct bw_dimensions;

struct cxx_param_list;

extern void lingen_tune_cutoffs_decl_usage(cxx_param_list & pl);
extern void lingen_tune_cutoffs_lookup_parameters(cxx_param_list & pl);

template<bool is_binary>
extern void lingen_tune_cutoffs(bw_dimensions<is_binary> & d, MPI_Comm comm, cxx_param_list & pl);


#endif	/* LINGEN_TUNE_CUTOFFS_HPP_ */
