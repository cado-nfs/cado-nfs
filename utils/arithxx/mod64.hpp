#ifndef CADO_UTILS_ARITHXX_MOD64_HPP
#define CADO_UTILS_ARITHXX_MOD64_HPP

/* A class for modular arithmetic with residues and modulus of up to 64
 * bits. */

#include <cstdint>
#include <cstdlib>

#include <gmp.h>

#include "macros.h"
#include "modint.hpp"
#include "misc.h"
#include "u64arith.h"
#include "arithxx_common.hpp"
#include "cxx_mpz.hpp"

struct arithxx_mod64 {
    class Modulus;
    class Residue;
    typedef Integer64 Integer;
};

class arithxx_mod64::Residue : public arithxx_details::Residue_base<arithxx_mod64>
{
    typedef arithxx_mod64 layer;
    friend class layer::Modulus;
    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api64<layer>;

    protected:
    uint64_t r;

    public:
    explicit Residue(Modulus const & m MAYBE_UNUSED)
        : r(0)
    { }
    Residue(Modulus const & m MAYBE_UNUSED, Residue const & s)
        : r(s.r)
    { }

    Residue() = delete;

    protected:
    /* These two prototypes are absent in arithxx_mod_mpz_new::Residue,
     * so what gives? We'll mark them as deprecated for the time being.
     */
    Residue & operator=(Integer const & s) ATTRIBUTE_DEPRECATED
    {
        r = s[0];
        return *this;
    }
    Residue & operator=(uint64_t const s) ATTRIBUTE_DEPRECATED
    {
        r = s;
        return *this;
    }
};

class arithxx_mod64::Modulus
    : public arithxx_details::api64<arithxx_mod64>
{
    typedef arithxx_mod64 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api64<layer>;

  protected:
    /* Data members */
    uint64_t m;

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m) { return m % 2 == 1; }
    static bool valid(cxx_mpz const & m) { return m % 2 == 1 && m > 0 && mpz_sizeinbase(m, 2) <= Integer::max_bits; }

    explicit Modulus(uint64_t const s)
        : m(s)
    { }
    explicit Modulus(Integer const & s)
        : m(s)
    { }
    Integer getmod() const { return Integer(m); }

  protected:
    /* Methods used internally */
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(Residue const & a MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(a.r < m);
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(uint64_t const a MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(a < m);
    }
    /* }}} */

  protected:
#if 0
    /* Computes (a * 2^64) % m ; do we need it?*/
    void tomontgomery(Residue & r, Residue const & a) const ATTRIBUTE_DEPRECATED
    {
        assertValid(a);
        u64arith_divr_2_1_1(&r.r, 0, a.r, m);
    }

    /* Computes (a / 2^64) % m */
    void frommontgomery(Residue & r, Residue const & a,
                        uint64_t const invm) const ATTRIBUTE_DEPRECATED
    {
        uint64_t tlow, thigh;
        assertValid(a);
        tlow = a.r * invm;
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, m);
        r.r = thigh + (a.r != 0 ? 1 : 0);
    }
#endif

    uint64_t get_u64(Residue const & s) const
    {
        assertValid(s);
        return s.r;
    }

    uint64_t getmod_u64() const { return m; }

  public:

    /* {{{ set(*4), set_reduced(*2), set0, set1 */
    void set(Residue & r, Residue const & s) const
    {
        assertValid(s);
        r = s;
    }

    void set(Residue & r, uint64_t const s) const { r.r = s % m; }

    void set(Residue & r, Integer const & s) const { set(r, s.getWord(0)); }


    /* Sets the residue to the class represented by the integer s.
     * Assumes that s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced(Residue & r, uint64_t const s) const
    {
        assertValid(s);
        r.r = s;
    }
    void set_reduced(Residue & r, Integer const & s) const
    {
        set_reduced(r, s.getWord(0));
    }
    void set(Residue & r, int64_t const s) const
    {
        set(r, safe_abs64(s));
        if (s < 0)
            neg(r, r);
    }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void set0(Residue & r) const { r.r = 0; }

    void set1(Residue & r) const { r.r = (m != 1); }
    /* }}} */

    /* {{{ get equal is0 is1 */
    Integer get(Residue const & s) const
    {
        assertValid(s);
        return Integer(s.r);
    }
    bool equal(Residue const & a, Residue const & b) const
    {
        assertValid(a);
        assertValid(b);
        return (a.r == b.r);
    }
    bool is0(Residue const & a) const
    {
        assertValid(a);
        return (a.r == 0);
    }
    bool is1(Residue const & a) const
    {
        assertValid(a);
        return (a.r == 1);
    }
    /* }}} */

    /* {{{ neg add(*2) add1 sub(*2) sub1 div2 */
    void neg(Residue & r, Residue const & a) const
    {
        assertValid(a);
        if (a.r == 0)
            r.r = a.r;
        else
            r.r = m - a.r;
    }
    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_addmod_1_1(&r.r, a.r, b.r, m);
    }


    void add1(Residue & r, Residue const & a) const
    {
        assertValid(a);
        r.r = a.r + 1;
        if (r.r == m)
            r.r = 0;
    }

    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        u64arith_addmod_1_1(&r.r, a.r, b % m, m);
    }
    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_submod_1_1(&r.r, a.r, b.r, m);
    }
    void sub1(Residue & r, Residue const & a) const
    {
        u64arith_submod_1_1(&r.r, a.r, 1, m);
    }
    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        u64arith_submod_1_1(&r.r, a.r, b % m, m);
    }
    bool div2(Residue & r, Residue const & a) const
    {
        if (m % 2 == 0)
            return false;
        else {
            r.r = u64arith_div2mod(a.r, m);
            return true;
        }
    }

    /* }}} */

    /* {{{ mul sqr */
    void mul(Residue & r, Residue const & a, Residue const & b) const
    {
        uint64_t t1, t2;
        assertValid(a);
        assertValid(b);
        u64arith_mul_1_1_2(&t1, &t2, a.r, b.r);
        u64arith_divr_2_1_1(&r.r, t1, t2, m);
    }
    void sqr(Residue & r, Residue const & a) const
    {
        uint64_t t1, t2;
        assertValid(a);
        u64arith_mul_1_1_2(&t1, &t2, a.r, a.r);
        u64arith_divr_2_1_1(&r.r, t1, t2, m);
    }
    /* }}} */

    /* Computes (a / 2^wordsize) % m, but result can be r = m.
       Input a must not be equal 0 */
    void redcsemi_u64_not0(Residue & r, uint64_t const a,
                           uint64_t const invm) const
        ATTRIBUTE_DEPRECATED    /* never used, never tested! */
    {
        uint64_t tlow, thigh;
        ASSERT(a != 0);
        tlow = a * invm; /* tlow <= 2^w-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, m);
        /* thigh:tlow <= (2^w-1) * m */
        r.r = thigh + 1;
        /* (thigh+1):tlow <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m
                          <= 2^w * (m + 1) - m */
        /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
             <= m */
    }
};
#endif /* CADO_UTILS_ARITHXX_MOD64_HPP */
