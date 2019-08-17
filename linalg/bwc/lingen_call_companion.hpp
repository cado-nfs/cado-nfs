#ifndef LINGEN_CALL_COMPANION_HPP_
#define LINGEN_CALL_COMPANION_HPP_

#include "lingen_substep_schedule.hpp"
#include "timing.h"     /* weighted_double */
#include "lingen_round_operand_size.hpp"

/* This object is passed as a companion info to a call
 * of bw_biglingen_recursive ; it is computed by the code in
 * plingen-tuning.cpp
 * but once tuning is over, it is essentially fixed.
 */

struct lingen_call_companion {
    bool recurse;
    bool go_mpi;
    double ttb;
    /* total_ncalls is a priori a power of two, but not always.
     * It is the number of calls that correspond to identical
     * lingen_call_companion::key keys.  In particular, since comparison
     * of ::key types is coarse, this means that total_ncalls is the
     * addition of the number of calls for two possibly different input
     * lengths.
     */
    size_t total_ncalls;
    struct mul_or_mp_times {
        /* XXX This must be trivially copyable because we share it via
         * MPI... ! */
        lingen_substep_schedule S;
        weighted_double
            tt,         /* 1, time per call to the mul operation */
            /* For the following, we have both the number of times the
             * operation is done within 1 call of the mul (or mp)
             * operation, plus the time of _each individual call_.
             */
            t_dft_A,    /* time per dft of the first operand, and so on */
            t_dft_A_comm,
            t_dft_B,
            t_dft_B_comm,
            t_conv,
            t_ift_C;
        size_t reserved_ram;
        size_t ram;

        /* we store the per-transform ram here, so that we can act
         * appropriately if we ever detect that it changes for one
         * specific call */
        size_t per_transform_ram;

        size_t asize, bsize, csize;
    };
    mul_or_mp_times mp, mul;
    struct key {
        int depth;
        size_t L;
        bool operator<(key const& a) const {
            if (depth < a.depth) return true;
            if (depth > a.depth) return false;
            return lingen_round_operand_size(L) < lingen_round_operand_size(a.L);
        }
    };
};


#endif	/* LINGEN_CALL_COMPANION_HPP_ */
