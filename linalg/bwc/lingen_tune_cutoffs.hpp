#ifndef LINGEN_TUNE_CUTOFFS_HPP_
#define LINGEN_TUNE_CUTOFFS_HPP_

#include "params.h"
#include "lingen.hpp"
#include "select_mpi.h"

extern void lingen_tune_cutoffs_decl_usage(cxx_param_list & pl);
extern void lingen_tune_cutoffs_lookup_parameters(cxx_param_list & pl);
extern void lingen_tune_cutoffs(bw_dimensions & d, MPI_Comm comm, cxx_param_list & pl);

#endif	/* LINGEN_TUNE_CUTOFFS_HPP_ */
