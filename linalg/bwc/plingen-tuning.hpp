#ifndef PLINGEN_TUNING_HPP_
#define PLINGEN_TUNING_HPP_

#include "params.h"
#include "plingen.hpp"
#include "select_mpi.h"
#include "lingen-substep-schedule.h"
#include <map>

size_t plingen_round_operand_size(size_t x, int bits=2);

/* This object is passed as a companion info to a call
 * of bw_biglingen_recursive ; it is computed by the code in
 * plingen-tuning.cpp
 * but once tuning is over, it is essentially fixed.
 */

struct lingen_call_companion {
    bool recurse;
    bool go_mpi;
    unsigned int weight;
    double ttb;
    struct mul_or_mp_times {
        lingen_substep_schedule S;
        double tt;      /* time per call to the mul operation */
        double t_dft_A; /* time per dft of the first operand, and so on */
        double t_dft_A_comm;
        double t_dft_B;
        double t_dft_B_comm;
        double t_addmul;
        double t_ift_C;
        size_t reserved_ram;
        size_t ram;
    };
    mul_or_mp_times mp, mul;
    struct key {
        int depth;
        size_t L;
        bool operator<(key const& a) const {
            if (depth < a.depth) return true;
            if (depth > a.depth) return false;
            return plingen_round_operand_size(L) < plingen_round_operand_size(a.L);
        }
    };
};

struct lingen_hints_t : public std::map<lingen_call_companion::key, lingen_call_companion> {
    typedef std::map<lingen_call_companion::key, lingen_call_companion> super;
    double tt_scatter_per_unit;
    double tt_gather_per_unit;
    int ipeak;
    size_t peak;
    void share(int root, MPI_Comm comm);
};

void plingen_tuning_decl_usage(cxx_param_list & pl);
void plingen_tuning_lookup_parameters(cxx_param_list & pl);
lingen_hints_t plingen_tuning(bw_dimensions & d, size_t, MPI_Comm comm, cxx_param_list & pl);

#endif	/* PLINGEN_TUNING_HPP_ */
