#ifndef LINGEN_CALL_COMPANION_HPP_
#define LINGEN_CALL_COMPANION_HPP_

#include <istream>
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
        std::istream& unserialize(std::istream& is) {
            S.unserialize(is);
            is  >> tt.n >> tt.t
                >> t_dft_A.n >> t_dft_A.t
                >> t_dft_A_comm.n >> t_dft_A_comm.t
                >> t_dft_B.n >> t_dft_B.t
                >> t_dft_B_comm.n >> t_dft_B_comm.t
                >> t_conv.n >> t_conv.t
                >> t_ift_C.n >> t_ift_C.t
                >> reserved_ram
                >> ram
                >> per_transform_ram
                >> asize
                >> bsize
                >> csize;
            return is;
        }
        std::ostream& serialize(std::ostream& os) const {
            S.serialize(os);
            os  << " " << tt.n << " " << tt.t
                << " " << t_dft_A.n << " " << t_dft_A.t
                << " " << t_dft_A_comm.n << " " << t_dft_A_comm.t
                << " " << t_dft_B.n << " " << t_dft_B.t
                << " " << t_dft_B_comm.n << " " << t_dft_B_comm.t
                << " " << t_conv.n << " " << t_conv.t
                << " " << t_ift_C.n << " " << t_ift_C.t
                << " " << reserved_ram
                << " " << ram
                << " " << per_transform_ram
                << " " << asize
                << " " << bsize
                << " " << csize;
            return os;
        }
    };/*}}}*/
    mul_or_mp_times mp, mul;
    std::istream& unserialize(std::istream& is) {
        is >> recurse
            >> go_mpi
            >> ttb
            >> total_ncalls;
        for(int i = 2 ; is && i-- ; ) {
            std::string s;
            is >> s;
            if (s == "MP") {
                mp.unserialize(is);
            } else if (s == "MUL") {
                mul.unserialize(is);
            } else {
                is.setstate(std::ios::failbit);
            }
        }
        return is;
    }
    std::ostream& serialize(std::ostream& os) const {
        os << " " << recurse
            << " " << go_mpi
            << " " << ttb
            << " " << total_ncalls;
        os << "\n";
        os << "\t" "MP"; mp.serialize(os); os << "\n";
        os << "\t" "MUL"; mul.serialize(os); os << "\n";
        return os;
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
