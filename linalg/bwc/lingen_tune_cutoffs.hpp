#ifndef LINGEN_TUNE_CUTOFFS_HPP_
#define LINGEN_TUNE_CUTOFFS_HPP_

#include "select_mpi.h"
struct bw_dimensions;
struct cxx_param_list;

extern void lingen_tune_cutoffs_decl_usage(cxx_param_list & pl);
extern void lingen_tune_cutoffs_lookup_parameters(cxx_param_list & pl);
extern void lingen_tune_cutoffs(bw_dimensions & d, MPI_Comm comm, cxx_param_list & pl);

#endif	/* LINGEN_TUNE_CUTOFFS_HPP_ */
