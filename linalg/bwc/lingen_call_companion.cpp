#include "cado.h"
#include "lingen_call_companion.hpp"

std::istream& lingen_call_companion::unserialize(std::istream& is) {
    std::string s;
    is >> s;
    if (s == io_token_recursive) {
        recurse = true;
    } else if (s == io_token_quadratic) {
        recurse = false;
    } else {
        is.setstate(std::ios::failbit);
        return is;
    }
    if (recurse) {
        is >> s;
        if (s == io_token_collective) {
            go_mpi = true;
        } else if (s == io_token_single) {
            go_mpi = false;
        } else {
            is.setstate(std::ios::failbit);
            return is;
        }
        for(int i = 2 ; is && i-- ; ) {
            std::string s;
            is >> s;
            if (s == io_token_MP) {
                mp.unserialize(is);
            } else if (s == io_token_MUL) {
                mul.unserialize(is);
            } else {
                is.setstate(std::ios::failbit);
            }
        }
    }
    /*
       >> ttb
       >> total_ncalls
       */
    ;
    return is;
}
std::ostream& lingen_call_companion::serialize(std::ostream& os) const {
    os << " " << (recurse ? io_token_recursive : io_token_quadratic);
    if (recurse) {
        os << " " << (go_mpi ? io_token_collective : io_token_single);
        os << "\n";
        os << "\t" << io_token_MP;  mp.serialize(os); os << "\n";
        os << "\t" << io_token_MUL; mul.serialize(os); os << "\n";
    }
    /*
       << " " << ttb
       << " " << total_ncalls;
       */
    return os;
}

/* This unserializes only part of the data: the schedule S/
 * The rest is always recomputed.
 */
std::istream& lingen_call_companion::mul_or_mp_times::unserialize(std::istream& is) {
    S.unserialize(is);
    /*
       is
       >> tt.n >> tt.t
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
       */
    return is;
}
std::ostream& lingen_call_companion::mul_or_mp_times::serialize(std::ostream& os) const {
    S.serialize(os);
    /*
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
       */
    return os;
}
