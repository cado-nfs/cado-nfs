#ifndef LINGEN_CALL_COMPANION_HPP_
#define LINGEN_CALL_COMPANION_HPP_

#include <istream>
#include <ostream>
#include "lingen_substep_schedule.hpp"
#include "timing.h"     /* weighted_double */
#include "lingen_round_operand_size.hpp"

/* This object is passed as a companion info to a call of
 * bw_biglingen_recursive ; it is computed by the code in
 * plingen-tuning.cpp but once tuning is over, it is essentially fixed.
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
    struct mul_or_mp_times {/*{{{*/
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
        
        /* This unserializes only part of the data: the schedule S/
         * The rest is always recomputed. Even for operator== (which we
         * chiefly use for compatibility checking), only the schedule
         * matters.
         */
        std::istream& unserialize(std::istream& is);
        std::ostream& serialize(std::ostream& os) const;
        bool operator==(mul_or_mp_times const & o) const {
            return S == o.S;
        }
        inline bool operator!=(mul_or_mp_times const & o) const { return !(*this == o); }
    };/*}}}*/
    mul_or_mp_times mp, mul;

    /* This unserializes only part of the data -- recurse, go_mpi,
     * and the schedules. The rest is always recomputed.
     */
    private:
    static constexpr const char * io_token_recursive = "recursive";
    static constexpr const char * io_token_quadratic = "quadratic";
    static constexpr const char * io_token_collective = "collective";
    static constexpr const char * io_token_single = "single";
    static constexpr const char * io_token_MP = "MP";
    static constexpr const char * io_token_MUL = "MUL";
    static constexpr const char * io_token_ignored = "-";
    public:
    std::istream& unserialize(std::istream& is);
    std::ostream& serialize(std::ostream& os) const;
    bool operator==(lingen_call_companion const & o) const;
    inline bool operator!=(lingen_call_companion const & o) const { return !(*this == o); }
    bool check() const {
        return mul.S.check() && mp.S.check();
    }
    struct key {
        int depth;
        size_t L;
        std::istream& unserialize(std::istream& is) {
            return is >> depth >> L;
        }
        std::ostream& serialize(std::ostream& os) const {
            return os << " " << depth << " " << L;
        }
        bool operator==(key const & o) const {
            return depth == o.depth && L == o.L;
        }
        bool operator<(key const& a) const {
            if (depth < a.depth) return true;
            if (depth > a.depth) return false;
            return lingen_round_operand_size(L) < lingen_round_operand_size(a.L);
        }
    };
};

inline std::ostream& operator<<(std::ostream& os, lingen_call_companion const & c) {
    return c.serialize(os);
}

inline std::istream& operator>>(std::istream& is, lingen_call_companion & c) {
    return c.unserialize(is);
}

inline std::ostream& operator<<(std::ostream& os, lingen_call_companion::key const & c) {
    return c.serialize(os);
}

inline std::istream& operator>>(std::istream& is, lingen_call_companion::key & c) {
    return c.unserialize(is);
}

#endif	/* LINGEN_CALL_COMPANION_HPP_ */
