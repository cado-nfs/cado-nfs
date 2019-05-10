#ifndef PLINGEN_TUNE_CUTOFFS_HPP_
#define PLINGEN_TUNE_CUTOFFS_HPP_

#include "params.h"
#include "plingen.hpp"
#include "select_mpi.h"

extern void plingen_tune_cutoffs_decl_usage(cxx_param_list & pl);
extern void plingen_tune_cutoffs_lookup_parameters(cxx_param_list & pl);
extern void plingen_tune_cutoffs(bw_dimensions & d, MPI_Comm comm, cxx_param_list & pl);

#endif	/* PLINGEN_TUNE_CUTOFFS_HPP_ */
