#ifndef CADO_UTILS_ARITHXX_MODREDC64_HPP
#define CADO_UTILS_ARITHXX_MODREDC64_HPP

/* A class for modular arithmetic with residues and modulus of
 * up to 64 bits. Residues are stored in Montgomery form,
 * reduction after multiplication is done with REDC.
 */

#include <cstdint>
#include <cstdlib>

#include <vector>

#include "macros.h"
#include "misc.h"
#include "modint.hpp"
#include "u64arith.h"
#include "arithxx_common.hpp"

struct arithxx_modredc64 {
    class Modulus;
    class Residue;
    typedef Integer64 Integer;
};

class arithxx_modredc64::Residue : public arithxx_details::Residue_base<arithxx_modredc64>
{
    typedef arithxx_modredc64 layer;
    friend class layer::Modulus;
    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api64<layer>;

    protected:
    Integer r;

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

class arithxx_modredc64::Modulus
    : public arithxx_details::api64<arithxx_modredc64>
{
    typedef arithxx_modredc64 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api64<layer>;

  protected:
    /* Data members */
    uint64_t m;
    uint64_t invm;
    uint64_t mrecip;
    Residue one;


    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m) { return m % 2 == 1; }

    explicit Modulus(uint64_t const s)
        : m(s)
        , invm(-u64arith_invmod(s))
        , one(*this)
    {
        int const shift = u64arith_clz(m);
        uint64_t const ml = m << shift;
        uint64_t dummy;
        mrecip = u64arith_reciprocal_for_div(ml);
        u64arith_divqr_2_1_1_recip_precomp(&dummy, one.r.data(), 0, 1, ml, mrecip,
                                           shift);
    }
    explicit Modulus(Integer const & s)
        : Modulus(s.getWord(0))
    {
    }
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
    /* Computes (a * 2^64) % m */
    void tomontgomery(Residue & r, Residue const & a) const
    {
        assertValid(a);
        const int shift = u64arith_clz(m);
        const uint64_t ml = m << shift;
        uint64_t dummy;
        u64arith_divqr_2_1_1_recip_precomp(&dummy, r.r.data(), 0, a.r[0], ml, mrecip,
                                           shift);
    }

    /* Computes (a / 2^64) % m. Assumes a < m */
    void frommontgomery(uint64_t & r, uint64_t const a) const
    {
        uint64_t tlow, thigh;
        assertValid(a);
        tlow = a * invm;
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, m);
        r = thigh + ((a != 0) ? 1 : 0);
    }


    uint64_t get_u64(Residue const & s) const
    {
        uint64_t r;
        assertValid(s);
        frommontgomery(r, s.r[0]);
        return r;
    }

    uint64_t getmod_u64() const { return m; }

  public:

    /* {{{ set(*4), set_reduced(*2), set0, set1 */
    void set(Residue & r, Residue const & s) const
    {
        assertValid(s);
        r = s;
    }

    /** XXX This is specific to Montgomery form */
    /* Puts in r the value of s * beta mod m, where beta is the word base.
       Note: s can be any uint64_t, in particular can be larger than m.
       When 0 <= s < m, use set_reduced for better efficiency. */
    void set(Residue & r, uint64_t const s) const
    {
        uint64_t plow, phigh;

        u64arith_mul_1_1_2(&plow, &phigh, s, one.r[0]);
        u64arith_redc(r.r.data(), plow, phigh, m, invm);
        tomontgomery(r, r);
    }

    void set(Residue & r, Integer const & s) const { set(r, s.getWord(0)); }

    /** XXX This is specific to Montgomery form */
    /* Sets the residue_t to the class represented by the integer s.
     * Assumes that s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced(Residue & r, uint64_t const s) const
    {
        assertValid(s);
        r.r = s;
        tomontgomery(r, r);
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

    /** XXX This is specific to Montgomery form */
    void set1(Residue & r) const { r = one; }
    /* }}} */

    /* {{{ get equal is0 is1 */
    Integer get(Residue const & s) const
    {
        assertValid(s);
        uint64_t t;
        frommontgomery(t, s.r[0]);
        return Integer(t);
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
        return equal(a, one);
    }
    /* }}} */

    /* {{{ neg add(*2) add1 sub(*2) sub1 div2 */
    void neg(Residue & r, Residue const & a) const
    {
        assertValid(a);
        if (a.r == 0)
            r.r = a.r;
        else
            r.r[0] = m - a.r[0];
    }
    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_addmod_1_1(r.r.data(), a.r[0], b.r[0], m);
    }

    /** XXX This is specific to Montgomery form */
    void add1(Residue & r, Residue const & a) const { add(r, a, one); }

    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        Residue t(*this);

        assertValid(a);
        set(t, b);
        add(r, a, t);
    }
    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
        u64arith_submod_1_1(r.r.data(), a.r[0], b.r[0], m);
    }

    /** XXX This is specific to Montgomery form */
    void sub1(Residue & r, Residue const & a) const { sub(r, a, one); }

    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        Residue t(*this);

        assertValid(a);
        set(t, b);
        sub(r, a, t);
    }

    bool div2(Residue & r, Residue const & a) const
    {
        r.r = u64arith_div2mod(a.r[0], m);
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

        ASSERT_EXPENSIVE(m % 2 != 0);
        assertValid(a);

        u64arith_mul_1_1_2(&plow, &phigh, a.r[0], b);
        /* We have a <= m-1, b <= 2^64 - 1. Thus the product
           phigh:plow <= (m-1)*(2^64 - 1) = m*2^64 - 2^64 - m + 1,
           and with m >= 1,
           phigh:plow <= m*2^64 - 2^64, so phigh < m. */
        u64arith_redc(&r, plow, phigh, m, invm);
    }
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
        u64arith_redc(r.r.data(), plow, phigh, m, invm);
    }
    /* }}} */

    /* {{{ V_dadd and V_dbl for Lucas sequences */
    /* Given a = V_n (x), b = V_m (x) and d = V_{n-m} (x), compute V_{m+n} (x).
     * r can be the same variable as a or b but must not be the same variable as
     * d.
     */
    void V_dadd(Residue & r, Residue const & a, Residue const & b,
                Residue const & d) const
    {
        ASSERT(&r != &d);
        mul(r, a, b);
        sub(r, r, d);
    }

    /* Given a = V_n (x) and two = 2, compute V_{2n} (x).
     * r can be the same variable as a but must not be the same variable as two.
     */
    void V_dbl(Residue & r, Residue const & a, Residue const & two) const
    {
        ASSERT(&r != &two);
        sqr(r, a);
        sub(r, r, two);
    }
    /* }}} */

  private:
    // return c/a[i] mod N, or an empty vector if one of the a[i] is not
    // invertible
    std::vector<Integer> batchinv_redc(std::vector<uint64_t> const & a, Integer c) const;

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

    struct batch_Q_to_Fp_context;
};

struct arithxx_modredc64::Modulus::batch_Q_to_Fp_context {
    Integer remainder, quotient;
    Modulus D;
    batch_Q_to_Fp_context(Integer const & num, Integer const & den);
    std::vector<uint64_t> operator()(std::vector<uint64_t> const & p, int k=0) const;
};

#endif /* CADO_UTILS_ARITHXX_MODREDC64_HPP */
