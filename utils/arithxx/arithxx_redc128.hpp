#ifndef CADO_UTILS_ARITHXX_REDC128_HPP
#define CADO_UTILS_ARITHXX_REDC128_HPP

#include <cstdint>

#include "u64arith.h"

#include "arithxx_common.hpp"
#include "arithxx_api.hpp"
#include "arithxx_redc.hpp"

namespace arithxx_details {
    template<typename layer> struct redc128;
}

template<typename layer>
struct arithxx_details::redc128
: public redc<layer>
{
    /* this gets instantiated as a base class of layer::Modulus, and
     * relies on layer::Modulus having the following data members.
     */
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    using redc<layer>::redc;

    Residue one;
    uint64_t invm;
    uint64_t mrecip;

    explicit redc128(Integer const & m)
        : redc<layer>(m)
        , one(downcast())
        , invm(-u64arith_invmod(m[0]))
    {
        int const shift = u64arith_clz(m[1]);
        uint64_t dummy, ml[2] = {m[0], m[1]};
        u64arith_shl_2(&ml[0], &ml[1], shift);
        mrecip = u64arith_reciprocal_for_div_3by2(ml[0], ml[1]);
        u64arith_divqr_3_2_1_recip_precomp(&dummy, one.r.data(), one.r.data()+1, 0, 0,
                                           1, ml[0], ml[1], mrecip, shift);
    }

    protected:
    /* r = s * 2^64 mod m */
    void unredc1(Integer & r, Integer const & s) const
    {
        auto const & me = arithxx_details::api<layer>::downcast();
        int const shift = int(me.m.clz());
        typename layer::Integer ml = me.m << shift;
        uint64_t dummy;
        u64arith_divqr_3_2_1_recip_precomp(&dummy, r.data(), r.data()+1,
                0, s[0], s[1], ml[0], ml[1], me.mrecip, shift);
    }

    void redc1_wide_inplace(uint64_t *t) const
    {
        auto const & me = arithxx_details::api<layer>::downcast();
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

    /* redc1 has overloads in redc<layer>, we need to take them into
     * account */
    using redc<layer>::redc1;

  private:
    /* Computes r = (a * b * 2^-64) mod m, where a is in REDC
     * representation */
    void mul_ul(Residue & r, Residue const & a, uint64_t b) const;
  protected:
    std::vector<Integer> batchinv_redc(std::vector<uint64_t> const & a, Integer const & c) const;
    friend struct arithxx_details::batch_Q_to_Fp_context<layer>;
};
#endif	/* UTILS_ARITHXX_REDC128_HPP_ */
