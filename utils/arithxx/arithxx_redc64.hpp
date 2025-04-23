#ifndef CADO_UTILS_ARITHXX_REDC64_HPP
#define CADO_UTILS_ARITHXX_REDC64_HPP

#include <cstdint>

#include <vector>

#include "u64arith.h"

#include "arithxx_common.hpp"
#include "arithxx_api.hpp"
#include "arithxx_redc.hpp"

namespace arithxx_details {
    template<typename layer> struct redc64;
}

template<typename layer>
struct arithxx_details::redc64
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

    explicit redc64(Integer const & m)
        : redc<layer>(m)
        , one(downcast())
        , invm(-u64arith_invmod(m[0]))
    {
        int const shift = u64arith_clz(m[0]);
        uint64_t dummy;
        uint64_t const ml = m[0] << shift;
        mrecip = u64arith_reciprocal_for_div(ml);
        u64arith_divqr_2_1_1_recip_precomp(&dummy, one.r.data(), 0, 1, ml, mrecip,
                                           shift);
    }

    /* r = s * 2^64 mod m */
    void unredc1(Integer & r, Integer const & a) const
    {
        auto const & me = downcast();
        const int shift = u64arith_clz(me.m[0]);
        const uint64_t ml = me.m[0] << shift;
        uint64_t dummy;
        u64arith_divqr_2_1_1_recip_precomp(&dummy, r.data(), 0, a[0], ml, mrecip,
                                           shift);
    }

#if 0
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
#endif

    void redc1(uint64_t * r, uint64_t const * a) const
    {
        auto const & me = downcast();
        uint64_t tlow, thigh;
        tlow = a[0] * invm;
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m[0]);
        r[0] = thigh + (a[0] != 0);
    }

    /* redc1 has overloads in redc<layer>, we need to take them into
     * account */
    using redc<layer>::redc1;

  protected:
    // return c/a[i] mod N, or an empty vector if one of the a[i] is not
    // invertible
    std::vector<Integer> batchinv_redc(std::vector<uint64_t> const & a, Integer const & c) const;
    friend struct arithxx_details::batch_Q_to_Fp_context<layer>;
};


#endif	/* UTILS_ARITHXX_REDC64_HPP_ */
