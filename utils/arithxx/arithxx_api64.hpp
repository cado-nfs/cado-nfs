#ifndef UTILS_ARITHXX_API64_HPP_
#define UTILS_ARITHXX_API64_HPP_

#include "arithxx_api.hpp"

namespace arithxx_details {
template <typename layer>
struct api64
: public api<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    bool div3(Residue &, Residue const &) const;
    bool div5(Residue &, Residue const &) const;
    bool div7(Residue &, Residue const &) const;
    bool div11(Residue &, Residue const &) const;
    bool div13(Residue &, Residue const &) const;

    void gcd (Integer &g, const Residue &r) const;
    int jacobi(Residue const &) const;

    bool sprp2_is_enough() const { return downcast().m < 2047; }
};
}

#endif	/* UTILS_ARITHXX_API64_HPP_ */
