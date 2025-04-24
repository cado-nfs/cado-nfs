#ifndef CADO_UTILS_ARITHXX_MODREDC64_HPP
#define CADO_UTILS_ARITHXX_MODREDC64_HPP

/* A class for modular arithmetic with residues and modulus of
 * up to 64 bits. Residues are stored in Montgomery form,
 * reduction after multiplication is done with REDC.
 */

#include <cstdint>
#include <cstdlib>

#include <type_traits>
#include <vector>

#include <gmp.h>

#include "macros.h"
#include "misc.h"
#include "modint.hpp"
#include "u64arith.h"
#include "arithxx_common.hpp"
#include "cxx_mpz.hpp"
#include "arithxx_redc64.hpp"

struct arithxx_modredc64 {
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

    typedef std::true_type uses_montgomery_representation;

    typedef std::false_type even_moduli_allowed;
};

/* okay, at this point it's a typedef, really... */
class arithxx_modredc64::Residue
    : public arithxx_details::Residue_base<arithxx_modredc64>
{
    using Residue_base::Residue_base;
};

class arithxx_modredc64::Modulus
    : public arithxx_details::redc64<arithxx_modredc64>
{
    typedef arithxx_modredc64 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api_bysize<layer>;

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m) { return m % 2 == 1; }
    static bool valid(cxx_mpz const & m) { return m % 2 == 1 && m > 0 && mpz_sizeinbase(m, 2) <= Integer::max_bits; }

    explicit Modulus(uint64_t s)
        : Modulus(Integer(s))
    { }

    explicit Modulus(Integer const & s)
        : redc64(s)
    { }

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

  public:

    /* {{{ add(*2) add1 sub(*2) sub1 div2 */
    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_addmod_1_1(r.r.data(), a.r[0], b.r[0], m[0]);
    }

    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_submod_1_1(r.r.data(), a.r[0], b.r[0], m[0]);
    }

    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        add(r, a, downcast()(b));
    }

    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        sub(r, a, downcast()(b));
    }

    bool div2(Residue & r, Residue const & a) const
    {
        r.r = u64arith_div2mod(a.r[0], m[0]);
        return true;
    }
    /* }}} */

  private:
    /* This is part of the mul interface, but not everything is taken as
     * Residue types. Here we compute a*b mod m, with a and b in
     * montgomery form.
     */
    void mul_u64_u64(uint64_t & r, Residue const & a, uint64_t const b) const
    {
        uint64_t plow, phigh;

        ASSERT_EXPENSIVE(m[0] % 2 != 0);
        assertValid(a);

        u64arith_mul_1_1_2(&plow, &phigh, a.r[0], b);
        /* We have a <= m-1, b <= 2^64 - 1. Thus the product
           phigh:plow <= (m-1)*(2^64 - 1) = m*2^64 - 2^64 - m + 1,
           and with m >= 1,
           phigh:plow <= m*2^64 - 2^64, so phigh < m. */
        u64arith_redc(&r, plow, phigh, m[0], invm);
    }
    friend struct arithxx_details::redc64<layer>;


  public:
    /* {{{ mul sqr */
    /** XXX This is specific to Montgomery form */
    void mul(Residue & r, Residue const & a, Residue const & b) const
    {
        assertValid(b);
        mul_u64_u64(r.r[0], a, b.r[0]);
    }

    /** XXX This is specific to Montgomery form */
    void sqr(Residue & r, Residue const & a) const
    {
        uint64_t plow, phigh;

        assertValid(a);

        u64arith_sqr_1_2(&plow, &phigh, a.r[0]);
        u64arith_redc(r.r.data(), plow, phigh, m[0], invm);
    }
    /* }}} */
};

#endif /* CADO_UTILS_ARITHXX_MODREDC64_HPP */
