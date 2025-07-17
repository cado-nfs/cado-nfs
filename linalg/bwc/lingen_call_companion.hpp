#ifndef CADO_LINGEN_CALL_COMPANION_HPP
#define CADO_LINGEN_CALL_COMPANION_HPP

#include <cstddef>

#include <stdexcept>
#include <string>
#include <iosfwd>
#include <istream>
#include <array>

#include "lingen_substep_schedule.hpp"
#include "timing.h"
#include "lingen_round_operand_size.hpp"
#include "lingen_mul_substeps_base.hpp"
#include "fmt/ostream.h"

/* This object is passed as a companion info to a call of
 * bw_biglingen_recursive ; it is computed by the code in
 * plingen-tuning.cpp but once tuning is over, it is essentially fixed.
 */

struct lingen_call_companion {
    /* This is a safeguard. Only the tuning code can create complete call
     * companions.
     */
    bool complete = false;

    unsigned int mesh = 0;
    bool recurse() const { return mesh > 0; }
    bool go_mpi() const { return mesh > 1; }

    double ttb = 0;
    /* total_ncalls is a priori a power of two, but not always.
     * It is the number of calls that correspond to identical
     * lingen_call_companion::key keys.  In particular, since comparison
     * of ::key types is coarse, this means that total_ncalls is the
     * addition of the number of calls for two possibly different input
     * lengths.
     */
    size_t total_ncalls = 0;
    struct mul_or_mp_times {/*{{{*/
        op_mul_or_mp_base::op_type_t op_type;
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

        mul_or_mp_times(op_mul_or_mp_base::op_type_t op_type) : op_type(op_type) {}
        const char * fft_name() const { return S.fft_name(); }
        std::string step_name() const {
            std::string s = op_mul_or_mp_base::op_name(op_type);
            s += ';';
            s += fft_name();
            return s;
        }
        /* we store the per-transform ram here, so that we can act
         * appropriately if we ever detect that it changes for one
         * specific call. This is supposed to be the "peak" of the
         * recorded sizes.
         */
        std::array<size_t, 3> fft_alloc_sizes;
        std::array<std::array<unsigned int, 3>, 2> peak_ram_multipliers;
        size_t ram(std::array<size_t, 3> fft_alloc_sizes) const {
            size_t rpeak = 0;
            for(auto const & M : peak_ram_multipliers) {
                size_t r = 0;
                for(unsigned int i = 0 ; i < 3 ; i++)
                    r += M[i] * fft_alloc_sizes[i];
                if (r > rpeak) rpeak = r;
            }
            return rpeak;
        }
        size_t ram() const {
            return ram(fft_alloc_sizes);
        }
        size_t ram_total() const {
            return ram(fft_alloc_sizes) + reserved_ram;
        }

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
    mul_or_mp_times mp  { op_mul_or_mp_base::OP_MP };
    mul_or_mp_times mul { op_mul_or_mp_base::OP_MUL };

    mul_or_mp_times operator[](op_mul_or_mp_base::op_type_t op_type) const {
        switch(op_type) {
            case op_mul_or_mp_base::OP_MP: return mp;
            case op_mul_or_mp_base::OP_MUL: return mul;
            default: throw std::runtime_error("bad op");
        }
    }


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
            return lingen_round_operand_size(L) > lingen_round_operand_size(a.L);
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

namespace fmt {
    template <> struct formatter<lingen_call_companion::key>: ostream_formatter {};
}


#endif	/* LINGEN_CALL_COMPANION_HPP_ */
