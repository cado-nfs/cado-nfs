#ifndef UTILS_ARITHXX_ARITHXX_API128_HPP_
#define UTILS_ARITHXX_ARITHXX_API128_HPP_

#include <type_traits>

#include "arithxx_api.hpp"
#include "modint.hpp"

namespace arithxx_details {
template <typename layer>
struct api128
: public api<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    static_assert(std::is_same<Integer, Integer128>::value,
            "api128 assumes that our underlying integer type is Integer128");

    void gcd (Integer &g, const Residue &r) const;
    int jacobi(Residue const &) const;

    bool inv(Residue &, Residue const &) const;

};


template<typename parent>
struct redc_interface
: public parent
{
    using typename parent::Modulus;
    using typename parent::Residue;
    using typename parent::Integer;
    using parent::downcast;

    protected:

    /* r = s * 2^64 mod m */
    void tomontgomery1(Residue & r, Residue const & s) const
    {
        auto const & me = downcast();
        int const shift = me.m.clz();
        Integer ml = me.m << shift;
        uint64_t dummy;
        u64arith_divqr_3_2_1_recip_precomp(&dummy, r.r.data(), r.r.data()+1,
                0, s.r[0], s.r[1], ml[0], ml[1], me.mrecip, shift);
    }

    void redc1_wide_inplace(uint64_t *t) const
    {
        auto const & me = downcast();
        uint64_t k;
        uint64_t pl, ph;

        k = t[0] * me.invm;
        u64arith_mul_1_1_2 (&pl, &ph, k, me.m[0]);
        ph += (t[0] != 0UL);        // t[0] == 0
        u64arith_add_1_2 (&(t[1]), &(t[2]), ph);
        u64arith_mul_1_1_2 (&pl, &ph, k, me.m[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[1]), &(t[2]), pl, ph);
    }

    void redc1(uint64_t * r, uint64_t const * s) const
    {
        uint64_t t[3] = { s[0], s[1], 0 };
        redc1_wide_inplace(t);
        r[0] = t[1];
        r[1] = t[2];
    }

    /* Do a one-word REDC, i.e., r == s / w (mod m), w = 2^64.
       If m > w, r < 2m. If s < m, then r < m */
    void redc1(Residue & r, Residue const & s) const { redc1(r.r.data(), s.r.data()); }
    void redc1(Integer & r, Integer const & s) const { redc1(r.data(), s.data()); }



    /* r = s * 2^128 mod m */
    void tomontgomery(Residue & r, Residue const & s) const
    {
        tomontgomery1(r, s);
        tomontgomery1(r, r);
    }

    /* Converts s out of Montgomery form by dividing by 2^(2*ULONG_BITS).
       Requires s < m. */
    void frommontgomery(Residue & r, Residue const & s) const
    {
        /* Do two REDC steps */
        redc1(r, s);
        redc1(r, r);
    }
};

}



#endif	/* UTILS_ARITHXX_ARITHXX_API128_HPP_ */
