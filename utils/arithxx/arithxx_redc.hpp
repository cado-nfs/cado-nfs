#ifndef CADO_UTILS_ARITHXX_REDC_INTERFACE_HPP
#define CADO_UTILS_ARITHXX_REDC_INTERFACE_HPP

#include <cstdint>
#include <cstddef>

#include "arithxx_common.hpp"
#include "macros.h"

template<typename layer>
struct arithxx_details::redc
: public api_bysize<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;
    using api_bysize<layer>::api_bysize;

    static_assert(!layer::even_moduli_allowed::value, "redc wants odd moduli");

    /* Do a one-word REDC, i.e., r == s / w (mod m), w = 2^64.
       If m > w, r < 2m. If s < m, then r < m */
    void redc1(Residue & r, Residue const & s) const { 
        downcast().redc1(r.r.data(), s.r.data());
    }
    void redc1(Integer & r, Integer const & s) const {
        downcast().redc1(r.data(), s.data());
    }

    /* r = s * 2^128 mod m */
    void tomontgomery(Integer & r, Integer const & s) const
    {
        downcast().unredc1(r, s);
        /* of course this is optimized out if max_size_in_words==1 */
        for(size_t i = 1 ; i < Integer::max_size_in_words ; i++)
            downcast().unredc1(r, r);
    }

    /* Converts s out of Montgomery form by dividing by 2^(2*ULONG_BITS).
       Requires s < m. */
    void frommontgomery(Integer & r, Integer const & s) const
    {
        downcast().redc1(r, s);
        /* of course this is optimized out if max_size_in_words==1 */
        for(size_t i = 1 ; i < Integer::max_size_in_words ; i++)
            downcast().redc1(r, r);
    }

    /* {{{ convenient overloads. */
    void tomontgomery(Residue & r, Residue const & a) const
    {
        downcast().assertValid(a);
        downcast().tomontgomery(r.r, a.r);
    }
    Integer tomontgomery(Integer const & a) const
    {
        Integer r;
        downcast().tomontgomery(r, a);
        return r;
    }

    void frommontgomery(Residue & r, Residue const & a) const
    {
        downcast().assertValid(a);
        downcast().frommontgomery(r.r, a.r);
    }
    Integer frommontgomery(Integer const & a) const
    {
        Integer r;
        downcast().frommontgomery(r, a);
        return r;
    }
    /* }}} */

    /* the overloads here are *in addition* to the ones that exist by
     * default */
    using api<layer>::set;
    using api<layer>::set_reduced;

    /* {{{ set(*2), set_reduced(*1), set1, get, is1, add1, sub1: dependent on redc */
    void set(Residue & r, uint64_t const s) const
    {
        /* s is only one word, so it's of course reduced.
         */
        r.r = tomontgomery(Integer(s));
    }
    void set(Residue & r, Integer const & s) const
    {
        auto const & me = downcast();
        if (s < me.m) {
            r.r = s;
        } else {
            r.r = s % me.m;
        }
        tomontgomery(r, r);
    }
    void set_reduced(Residue & r, Integer const & s) const
    {
        ASSERT(s < downcast().m);
        r.r = s;
        tomontgomery(r, r);
    }
    void set1(Residue & r) const {
        auto const & me = downcast();
        me.set(r, me.one);
    }

    /* Returns the residue as an integer */
    Integer get(Residue const & s) const
    {
        return frommontgomery(s.r);
    }

    bool is1(Residue const & a) const {
        auto const & me = downcast();
        return me.equal(a, me.one);
    }
    void add1(Residue & r, Residue const & a) const {
        auto const & me = downcast();
        me.add(r, a, me.one);
    }
    void sub1(Residue & r, Residue const & a) const {
        auto const & me = downcast();
        me.sub(r, a, me.one);
    }
    /* }}} */

    public:

    /* compute num/den/2^k mod p[i] for all i. den must be odd. p[i] must
     * be odd (unless k=0, in which case even moduli are ok).
     *
     * This uses the modulus layer only under the hood, and relies on the
     * fact that den is an acceptable modulus.
     */
    static std::vector<uint64_t> batch_Q_to_Fp(Integer const & num,
            Integer const & den, int k,
            std::vector<uint64_t> const & p);
};




#endif	/* UTILS_ARITHXX_REDC_INTERFACE_HPP_ */
