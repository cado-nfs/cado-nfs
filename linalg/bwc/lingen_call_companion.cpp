#include "cado.h" // IWYU pragma: keep

#include <ios>
#include <ostream>
#include <istream>
#include <string>

#include "lingen_substep_schedule.hpp"
#include "lingen_call_companion.hpp"

std::istream& lingen_call_companion::unserialize(std::istream& is) {
    std::string s;
    is >> s;
    if (s == io_token_recursive) {
        mesh = 1;
    } else if (s == io_token_quadratic) {
        mesh = 0;
    } else {
        is.setstate(std::ios::failbit);
        return is;
    }
    if (recurse()) {
        is >> s;
        if (s == io_token_collective) {
            is >> mesh;
        } else if (s == io_token_single) {
            mesh = 1;
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
    os << (recurse() ? io_token_recursive : io_token_quadratic);
    if (recurse()) {
        if (go_mpi()) {
            os << " " << io_token_collective << " " << mesh;
        } else {
            os << " " << io_token_single;
        }
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
bool lingen_call_companion::operator==(lingen_call_companion const & o) const
{
    if (mesh != o.mesh) return false;
    if (mp != o.mp) return false;
    if (mul != o.mul) return false;
    return true;
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
