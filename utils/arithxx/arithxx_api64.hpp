#ifndef UTILS_ARITHXX_API64_HPP_
#define UTILS_ARITHXX_API64_HPP_

#include <cstdint>
#include <cstddef>

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

    /* We have overrides for all of these that happen to be
     * common to both mod64 and modredc64
     */
    bool sprp(Residue const &) const;
    bool sprp2() const;
    bool isprime() const;

    bool div3(Residue &, Residue const &) const;
    bool div5(Residue &, Residue const &) const;
    bool div7(Residue &, Residue const &) const;
    bool div11(Residue &, Residue const &) const;
    bool div13(Residue &, Residue const &) const;

    void gcd (Integer &g, const Residue &r) const;
    int jacobi(Residue const &) const;

    template <typename value_type>
        void
        pow2_oneWord(value_type mask, value_type word, Residue &t) const;

    void pow2 (Residue &r, uint64_t e) const;
    void pow2 (Residue &r, const uint64_t *e, size_t e_nrwords) const;
    void pow2 (Residue &r, const Integer &e) const;
    void pow3 (Residue &r, uint64_t e) const;
};
}

#endif	/* UTILS_ARITHXX_API64_HPP_ */
