#ifndef PLINGEN_TUNING_HPP_
#define PLINGEN_TUNING_HPP_

#include "params.h"
#include "plingen.h"
#include "select_mpi.h"

void plingen_tuning_decl_usage(cxx_param_list & pl);
void plingen_tuning_lookup_parameters(cxx_param_list & pl);
void plingen_tuning(dims * d, MPI_Comm comm, cxx_param_list & pl);

#endif	/* PLINGEN_TUNING_HPP_ */
