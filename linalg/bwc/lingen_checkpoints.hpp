#ifndef LINGEN_CHECKPOINTS_HPP_
#define LINGEN_CHECKPOINTS_HPP_

#include "lingen.hpp"
#include "lingen_bmstatus.hpp"
#include "params.h"

#include "lingen_matpoly_select.hpp"
#include "lingen_bigmatpoly.hpp"

enum cp_which {
    LINGEN_CHECKPOINT_E,
    LINGEN_CHECKPOINT_PI,
};

template<typename matpoly_type>
int load_checkpoint_file(bmstatus & bm, cp_which which, matpoly_type &, unsigned int t0, unsigned int t1);
template<typename matpoly_type>
int save_checkpoint_file(bmstatus & bm, cp_which which, matpoly_type const &, unsigned int t0, unsigned int t1);

template<>
int load_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly &, unsigned int t0, unsigned int t1);
template<>
int save_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly const &, unsigned int t0, unsigned int t1);

template<>
int load_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly &, unsigned int t0, unsigned int t1);
template<>
int save_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly const &, unsigned int t0, unsigned int t1);

void lingen_checkpoints_decl_usage(cxx_param_list & pl);
void lingen_checkpoints_lookup_parameters(cxx_param_list & pl);
void lingen_checkpoints_interpret_parameters(cxx_param_list & pl);

#endif	/* LINGEN_CHECKPOINTS_HPP_ */
