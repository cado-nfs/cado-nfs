#ifndef CADO_UTILS_ARITHXX_MOD64_HPP
#define CADO_UTILS_ARITHXX_MOD64_HPP

/* A class for modular arithmetic with residues and modulus of up to 64
 * bits. */

#include <cstdint>
#include <cstdlib>

#include <type_traits>

#include <gmp.h>

#include "macros.h"
#include "modint.hpp"
#include "u64arith.h"
#include "arithxx_common.hpp"
#include "cxx_mpz.hpp"

struct arithxx_mod64 {
    class Modulus;
    class Residue;
    typedef Integer64 Integer;

    /* When we multiply by a small constant, we use a left-to-right
     * binary method. So we typically have log(n) shifts and log(n)/2
     * additions, which should be compared to the cost of a runtime
     * multiplication.
     */
    typedef std::integral_constant<int, 8> mul_c_cutoff;

    /* this gives k such that 2^k*modulus-1 <= Integer::max_value
     */
    typedef std::integral_constant<int, 0> overflow_bits;

    typedef std::false_type uses_montgomery_representation;

    /* IDK what this is good for, but as a matter fact, this code does
     * seem to grok even moduli.
     */
    typedef std::true_type even_moduli_allowed;
};

/* okay, at this point it's a typedef, really... */
class arithxx_mod64::Residue
    : public arithxx_details::Residue_base<arithxx_mod64>
{
    using Residue_base::Residue_base;
};

class arithxx_mod64::Modulus
    : public arithxx_details::api_bysize<arithxx_mod64>
{
    typedef arithxx_mod64 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api_bysize<layer>;

  protected:

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const &) { return true; }
    static bool valid(cxx_mpz const & m) { return m > 0 && mpz_sizeinbase(m, 2) <= Integer::max_bits; }

    explicit Modulus(uint64_t const s)
        : Modulus(Integer(s))
    {
    }
    explicit Modulus(Integer const & s)
        : api_bysize(s)
    {
        ASSERT(valid(m));
    }

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
    uint64_t getmod_u64() const { return m[0]; }

  public:

    /* {{{ add(*2) add1 sub(*2) sub1 div2 */
    void add1(Residue & r, Residue const & a) const
    {
        assertValid(a);
        r.r = a.r + 1;
        if (r.r == m)
            r.r = 0;
    }

    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_addmod_1_1(r.r.data(), a.r[0], b.r[0], m[0]);
    }
    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        u64arith_addmod_1_1(r.r.data(), a.r[0], b % m[0], m[0]);
    }
    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_submod_1_1(r.r.data(), a.r[0], b.r[0], m[0]);
    }
    void sub1(Residue & r, Residue const & a) const
    {
        u64arith_submod_1_1(r.r.data(), a.r[0], 1, m[0]);
    }
    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        u64arith_submod_1_1(r.r.data(), a.r[0], b % m[0], m[0]);
    }
    bool div2(Residue & r, Residue const & a) const
    {
        if (m % 2 == 0)
            return false;
        else {
            r.r = u64arith_div2mod(a.r[0], m[0]);
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
        u64arith_mul_1_1_2(&t1, &t2, a.r[0], b.r[0]);
        u64arith_divr_2_1_1(r.r.data(), t1, t2, m[0]);
    }
    void sqr(Residue & r, Residue const & a) const
    {
        uint64_t t1, t2;
        assertValid(a);
        u64arith_mul_1_1_2(&t1, &t2, a.r[0], a.r[0]);
        u64arith_divr_2_1_1(r.r.data(), t1, t2, m[0]);
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
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, m[0]);
        /* thigh:tlow <= (2^w-1) * m */
        r.r = thigh + 1;
        /* (thigh+1):tlow <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m
                          <= 2^w * (m + 1) - m */
        /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
             <= m */
    }
};
#endif /* CADO_UTILS_ARITHXX_MOD64_HPP */
