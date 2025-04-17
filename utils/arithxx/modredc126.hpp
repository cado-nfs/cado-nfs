#ifndef CADO_UTILS_ARITHXX_MODREDC126_HPP
#define CADO_UTILS_ARITHXX_MODREDC126_HPP

/* A class for modular arithmetic with modulus in [2^64+1, 2^126-1].
 * Moduli must be odd. Residues are stored in Montgomery form, reduction
 * after multiplication is done with REDC. */


#include "cado_config.h" // for HAVE_GCC_STYLE_AMD64_INLINE_ASM

#include <cstdint>

#include <array>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "macros.h"
#include "misc.h"
#include "modint.hpp"
#include "u64arith.h"
#include "arithxx_common.hpp"

struct arithxx_modredc126 {
    class Modulus;
    class Residue;
    typedef Integer128 Integer;
};

class arithxx_modredc126::Residue : public arithxx_details::Residue_base<arithxx_modredc126>
{
    typedef arithxx_modredc126 layer;
    friend class layer::Modulus;
    friend struct arithxx_details::api<layer>;

    protected:
    Integer r;

    public:
    explicit Residue(Modulus const & m MAYBE_UNUSED)
    { }
    Residue(Modulus const & m MAYBE_UNUSED, Residue const & s)
        : r {s.r}
    { }

    Residue() = delete;

private:
    /* This is only used in batchinv_redc, where we play tricks
     * with the Montgomery representation.
     */
    Residue(Modulus const & m MAYBE_UNUSED, Integer const & s)
        : r(s)
    {}
};

class arithxx_modredc126::Modulus
    : public arithxx_details::api<arithxx_modredc126>
{
    typedef arithxx_modredc126 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;

  protected:
    /* Data members */
    std::array<uint64_t, 2> m;
    uint64_t invm;
    uint64_t mrecip;
    Residue one;

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m)
    {
        return Integer(1, 1) <= m && m % 2 == 1 && !(m[1] >> 62);
    }
    static bool valid(cxx_mpz const & m) {
        return m % 2 == 1 && m > 0 && mpz_sizeinbase(m, 2) > 64 && mpz_sizeinbase(m, 2) <= 126;
    }

    explicit Modulus(Integer const s)
        : m(s)
        , one(*this)
    {
        ASSERT(valid(s));
        invm = -u64arith_invmod(m[0]);

        int const shift = u64arith_clz(m[1]);
        uint64_t dummy, ml[2] = {m[0], m[1]};
        u64arith_shl_2(&ml[0], &ml[1], shift);
        mrecip = u64arith_reciprocal_for_div_3by2(ml[0], ml[1]);
        u64arith_divqr_3_2_1_recip_precomp(&dummy, one.r.data(), one.r.data()+1, 0, 0,
                                           1, ml[0], ml[1], mrecip, shift);
    }

    /* Returns the modulus as an Integer. */
    // void getmod(Integer & r) const { r.set(m); }
    Integer getmod() const { return Integer(m); }

  protected:
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(Residue const & r MAYBE_UNUSED) const
    {
        ASSERT(u64arith_lt_2_2(r.r[0], r.r[1], m[0], m[1]));
    }

    /* r = s * 2^64 mod m */
    void tomontgomery1(Residue & r, Residue const & s) const
    {
        const int shift = u64arith_clz(m[1]);
        uint64_t dummy, ml[2] = {m[0], m[1]};
        u64arith_shl_2(&ml[0], &ml[1], shift);
        u64arith_divqr_3_2_1_recip_precomp(&dummy, r.r.data(), r.r.data()+1, 0, s.r[0],
                                           s.r[1], ml[0], ml[1], mrecip, shift);
    }
    /* }}} */

  protected:
    /* r = s * 2^128 mod m */
    void tomontgomery(Residue & r, Residue const & s) const
    {
        tomontgomery1(r, s);
        tomontgomery1(r, r);
    }

    void redc1_wide_inplace(uint64_t *t) const
    {
        uint64_t k;
        uint64_t pl, ph;

        k = t[0] * invm;
        u64arith_mul_1_1_2 (&pl, &ph, k, m[0]);
        ph += (t[0] != 0UL);        // t[0] == 0
        u64arith_add_1_2 (&(t[1]), &(t[2]), ph);
        u64arith_mul_1_1_2 (&pl, &ph, k, m[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[1]), &(t[2]), pl, ph);
    }

    void redc1(uint64_t * r, uint64_t const * s) const
    {
#if 1
        uint64_t t[3] = { s[0], s[1], 0 };
        redc1_wide_inplace(t);
        r[0] = t[1];
        r[1] = t[2];
#else
        uint64_t t[4], k;

        k = s[0] * invm;
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), k, m[0]);
        if (s[0] != 0)
            t[1]++;
        t[2] = 0;
        u64arith_add_1_2(&(t[1]), &(t[2]), s[1]);       /* t[2] <= 1 */
        u64arith_mul_1_1_2(&(t[0]), &(t[3]), k, m[1]);  /* t[3] < 2^w-1 */
        u64arith_add_2_2(&(t[1]), &(t[2]), t[0], t[3]); /* t[2] < 2^w */

        /* r = (k*m + s) / w, k <= w-1. If s < m, then r < m */
        r[0] = t[1];
        r[1] = t[2];
#endif
    }

    /* Do a one-word REDC, i.e., r == s / w (mod m), w = 2^64.
       If m > w, r < 2m. If s < m, then r < m */
    void redc1(Residue & r, Residue const & s) const { redc1(r.r.data(), s.r.data()); }
    void redc1(Integer & r, Integer const & s) const
    {
        std::array<uint64_t, 2> t = s.get();
        redc1(t.data(), t.data());
        r = Integer(t);
    }

    /* Converts s out of Montgomery form by dividing by 2^(2*ULONG_BITS).
       Requires s < m. */
    void frommontgomery(Residue & r, Residue const & s) const
    {
        /* Do two REDC steps */
        redc1(r, s);
        redc1(r, r);
    }

  public:
    /* {{{ set(*4), set_reduced(*2), set0, set1 */

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void set(Residue & r, Residue const & s) const { r = s; }
    void set(Residue & r, uint64_t const s) const
    {
        r.r[0] = s;
        r.r[1] = 0;
        tomontgomery(r, r);
    }
    void set(Residue & r, int64_t const s) const
    {
        set(r, safe_abs64(s));
        if (s < 0)
            neg(r, r);
    }

    void set(Residue & r, Integer const & s) const
    {
        if (s < Integer(m)) {
            r.r = s;
        } else {
            r.r = s % Integer(m);
        }
        tomontgomery(r, r);
    }
    /* Sets the residueredc2ul2_t to the class represented by the integer s.
       Assumes that s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced(Residue & r, uint64_t const s) const { set(r, s); }
    void set_reduced(Residue & r, Integer s) const
    {
        ASSERT(s < Integer(m));
        r.r = s;
        tomontgomery(r, r);
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void set0(Residue & r) const
    {
        r.r[0] = 0;
        r.r[1] = 0;
    }
    void set1(Residue & r) const { set(r, one); }
    /* }}} */

    /* {{{ get equal is0 is1 */
    /* Returns the residue as a modintredc2ul2_t */
    Integer get(Residue const & s) const
    {
        Residue t(*this);
        frommontgomery(t, s);
        return t.r;
    }

    /* do we really want to keep these two, or should we use operator== ?
     * comparison to 1 in montgomery form is tricky, though.
     */
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    bool equal(Residue const & a, Residue const & b) const
    {
        return (a.r[0] == b.r[0] && a.r[1] == b.r[1]);
    }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    bool is0(Residue const & a) const { return (a.r[0] == 0 && a.r[1] == 0); }

    bool is1(Residue const & a) const { return (equal(a, one)); }
    /* }}} */

    /* {{{ neg add(*2) add1 sub(*2) sub1 div2 */
    void neg(Residue & r, Residue const & a) const
    {
        if (is0(a)) {
            set0(r);
        } else {
            uint64_t t0 = m[0], t1 = m[1];
            u64arith_sub_2_2(&t0, &t1, a.r[0], a.r[1]);
            r.r[0] = t0;
            r.r[1] = t1;
        }
    }

    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        uint64_t const t0 = b.r[0],
                       t1 = b.r[1]; /* r, a, and/or b may overlap */
        set(r, a);
        u64arith_add_2_2(r.r.data(), r.r.data() + 1, t0, t1);
        u64arith_sub_2_2_ge(r.r.data(), r.r.data() + 1, m[0], m[1]);
    }

    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        {
            uint64_t s1 = m[0]; // NOLINT(misc-const-correctness)
            uint64_t s2 = m[1]; // NOLINT(misc-const-correctness)
            uint64_t t1 = a.r[0];
            uint64_t t2 = a.r[1];

            __asm__ __VOLATILE("subq %4, %0\n\t"
                               "sbbq %5, %1\n\t"    /* t -= b */
                               "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
                               "cmovncq %6, %3\n"
                               : "+&r"(t1), "+&r"(t2), "+&r"(s1), "+r"(s2)
                               : "g"(b.r[0]), "g"(b.r[1]), "rm"(0UL)
                               : "cc");
            u64arith_add_2_2(&t1, &t2, s1, s2);
            r.r[0] = t1;
            r.r[1] = t2;
        }
#else
        {
            uint64_t t1 = a.r[0], t2 = a.r[1];
            u64arith_sub_2_2(&t1, &t2, b.r[0], b.r[1]);

            if (t2 > a.r[1] || (t2 == a.r[1] && t1 > a.r[0]))
                u64arith_add_2_2(&t1, &t2, m[0], m[1]);

            r.r[0] = t1;
            r.r[1] = t2;
        }
#endif
    }

    void add1(Residue & r, Residue const & a) const { add(r, a, one); }

    void sub1(Residue & r, Residue const & a) const { sub(r, a, one); }

    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        Residue t(*this);

        set(t, b);
        add(r, a, t);
    }

    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        Residue t(*this);

        set(t, b);
        sub(r, a, t);
    }

    bool div2(Residue & r, Residue const & a) const
    {
        set(r, a);
        if (r.r[0] % 2 == 1)
            u64arith_add_2_2(r.r.data(), r.r.data() + 1, m[0], m[1]);
        u64arith_shr_2(r.r.data(), r.r.data() + 1, 1);
        return true;
    }
    /* }}} */

    /* {{{ mul sqr */
#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY                                                            \
    "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall "               \
    "abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif

  private:
    void mul(Integer & r, const Integer & a, const Integer & b) const
    {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
        uint64_t dummy;
        __asm__ __VOLATILE(
            /* Product of low words */
            "movq %[a0], %%rax\n\t"
            "mulq %[b0]\n\t" /* rdx:rax = a0*b0 <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t"
            /* Compute u0*m, add to t0:rax */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t2]\n\t" /* t2 = u0 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"        /* rdx:rax = u0*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"        /* if low word != 0, carry to high word */
            "movq %[t2], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t" /* t1:t0 = (a0*b0 + u0*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"
            "addq %%rax, %[t0]\n\t"
            "movq %[a0], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0*b0+u0*m)/2^64 */
            ABORT_IF_CY             /* <= 2^126 - 2^62 */

            /* 2 products of low and high word */
            "xorl %k[t2], %k[t2]\n\t"
            "mulq %[b1]\n\t" /* rdx:rax = a0*b1 <= (2^64-1)*(2^63-2^32-1) */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 */
            ABORT_IF_CY             /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
                               = 2^127 + 2^126 - 2^96 ... */

            /* Free slot here */
            "mulq %[b0]\n\t" /* rdx:rax = a1*b0 <= (2^63-2^32-1)*(2^64-1) */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"
            "setc %b[t2]\n\t" /* t2:t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 */
            /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                       = 2^128 + 2^126 - 2*2^96 ... */
            /* Product of high words */
            "mulq %[b1]\n\t" /* rdx:rax = a1*b1 <= (2^63-2^32-1)^2 */
            "addq %%rax, %[t1]\n\t"
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1:t0 = (a*b+u0*m)/2^64 */
            ABORT_IF_CY /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64
                   = 2^190 - 2^160 ... */
            /* Free slot here */
            /* Compute u1*m, add to t2:t1:t0 */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t0]\n\t" /* t0 = u1 */
            /* Free slot here */
            "mulq %[m0]\n\t" /* rdx:rax = m0*u1 <= (2^64-1)^2 */
            "negq %%rax\n\t" /* if low word != 0, carry to high word */
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "adcq $0,%[t2]\n\t" /* t2:t1:0 = (a*b+u0*m)/2^64 + u1*m0 */
            ABORT_IF_CY         /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */

            "mulq %[m1]\n\t" /* rdx:rax = u1*m1 */
            "addq %%rax, %[t1]\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*b+u0*m)/2^64 + u1*m)/2^64 */
            ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */

            "movq %[t1], %%rax\n\t" /* See if result > m */
            "movq %[t2], %%rdx\n\t"
            "subq %[m0], %[t1]\n\t"
            "sbbq %[m1], %[t2]\n\t"
            "cmovc %%rax, %[t1]\n\t" /* Carry -> restore old result */
            "cmovc %%rdx, %[t2]\n\t"
            : [t0] "=&r"(dummy), [t1] "=&r"(r[0]), [t2] "=&r"(r[1])
            : [a0] "g"(a[0]), [a1] "g"(a[1]), [b0] "rm"(b[0]),
              [b1] "rm"(b[1]), [m0] "rm"(m[0]), [m1] "rm"(m[1]),
              [invm] "rm"(invm)
            : "%rax", "%rdx", "cc");
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */
        uint64_t pl, ph, t[3];

        /* m < 1/4 W^2,  a,b < m */

        /* Product of the two low words */
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), a[0], b[0]);

        /* One REDC step */
        redc1(t, t); /* t < 2m < 1/2 W^2 */

        /* Products of one low and one high word  */
        u64arith_mul_1_1_2(&pl, &ph, a[1], b[0]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < 3/4 W^2 */
        u64arith_mul_1_1_2(&pl, &ph, a[0], b[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < W^2 */

        /* Product of the two high words */
        u64arith_mul_1_1_2(&pl, &(t[2]), a[1], b[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2(&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        redc1_wide_inplace(t);
        r[0] = t[1];
        r[1] = t[2];

        /* Result may be larger than m, but is < 2*m */
        u64arith_sub_2_2_ge(&(r[0]), &(r[1]), m[0], m[1]);
#endif
    }
    void mul(Residue & r, const Residue & a, const Integer & b) const
    {
        mul(r.r, a.r, b);
    }
  public:
    /* of course a Residue is just an Integer with some added type info.
     * We have a few distinct overloads of the mul code for convenience,
     * mostly because of the batchinv_redc code which plays a few tricks.
     */
    void mul(Residue & r, const Residue & a, const Residue & b) const
    {
        mul(r.r, a.r, b.r);
    }

    void sqr(Residue & r, Residue const & a) const
    {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM

        /* m <= 2^126-1
           Since m1>0, m*u is maximal for m0=1 and u=2^64-1, so
           u*m is bounded by (2^126 - 2^64 + 1)*(2^64 - 1) =
           2^190 - 2^128 - 2^126 + 2*2^64 - 1.
           If a,b <= 2^127-2^96-1, then
           ((a*b+u0*m)/2^64 + u1*m)/2^64 <=  2^127-2^96-1
           If we allow non-canonical residues up to 2^127-2^96-1, we can skip
           the final conditional subtraction. These residues are still < 2^127,
           so an addition does not overflow */

        uint64_t dummy;
        __asm__ __VOLATILE(
            /* Product of low words */
            "movq %[a0], %%rax\n\t"
            "mulq %[a0]\n\t" /* rdx:rax = a0^2 <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t"
            /* Compute u0*m, add to t0:rax */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t2]\n\t" /* t2 = u0 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"        /* rdx:rax = u0*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"        /* if low word != 0, carry to high word */
            "movq %[t2], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t" /* t1:t0 = (a0^2 + u0*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"
            "addq %%rax, %[t0]\n\t"
            "movq %[a0], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0^2+u0*m)/2^64 */
            ABORT_IF_CY             /* <= 2^126 - 2^62 */

            /* 2 products of low and high word */
            "xorl %k[t2], %k[t2]\n\t"
            "mulq %[a1]\n\t" /* rdx:rax = a0*a1 <= (2^64-1)*(2^63-2^32-1) */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0^2+u0*m)/2^64 + a0*a1 */
            ABORT_IF_CY             /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
                               = 2^127 + 2^126 - 2^96 ... */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "movq %[a1], %%rax\n\t" /* independent, goes in pipe 0 */
            "setc %b[t2]\n\t"       /* t2:t1:t0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 */
            /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                       = 2^128 + 2^126 - 2*2^96 ... */

            /* Product of high words */
            "mulq %[a1]\n\t" /* rdx:rax = a1^2 <= (2^63-2^32-1)^2 */
            "addq %%rax, %[t1]\n\t"
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1:t0 = (a^2+u0*m)/2^64 */
            ABORT_IF_CY /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64
                   = 2^190 - 2^160 ... */
            /* Free slot here */
            /* Compute u1*m, add to t2:t1:t0 */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t0]\n\t" /* t0 = u1 */
            /* Free slot here */
            "mulq %[m0]\n\t" /* rdx:rax = m0*u1 <= (2^64-1)^2 */
            "negq %%rax\n\t" /* if low word != 0, carry to high word */
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "adcq $0,%[t2]\n\t" /* t2:t1:0 = (a*a+u0*m)/2^64 + u1*m0 */
            ABORT_IF_CY         /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */

            "mulq %[m1]\n\t" /* rdx:rax = u1*m1 */
            "addq %%rax, %[t1]\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*a+u0*m)/2^64 + u1*m)/2^64 */
            ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */

            "movq %[t1], %%rax\n\t" /* See if result > m */
            "movq %[t2], %%rdx\n\t"
            "subq %[m0], %[t1]\n\t"
            "sbbq %[m1], %[t2]\n\t"
            "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
            "cmovc %%rdx, %[t2]\n\t"
            : [t0] "=&r"(dummy), [t1] "=&r"(r.r[0]), [t2] "=&r"(r.r[1])
            : [a0] "g"(a.r[0]), [a1] "g"(a.r[1]), [m0] "rm"(m[0]),
              [m1] "rm"(m[1]), [invm] "rm"(invm)
            : "%rax", "%rdx", "cc");
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */

        uint64_t pl, ph, t[3];

        /* m < 1/4 W^2,  a < m */

        /* Square of the low word */
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), a.r[0], a.r[0]);

        /* One REDC step */
        redc1(t, t); /* t < 2m < 1/2 W^2 */

        /* Products of low and high word  */
        u64arith_mul_1_1_2(&pl, &ph, a.r[1], a.r[0]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < 3/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < W^2 */

        /* Square of high word */
        u64arith_mul_1_1_2(&pl, &(t[2]), a.r[1], a.r[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2(&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        redc1_wide_inplace(t);
        r.r[0] = t[1];
        r.r[1] = t[2];

        /* Result may be larger than m, but is < 2*m */
        u64arith_sub_2_2_ge(r.r.data(), r.r.data() + 1, m[0], m[1]);
#endif
    }
    /* }}} */

  protected:
    /* Computes r = (a * b * 2^-64) mod m, where a is in REDC
     * representation */
    void mul_ul(Residue & r, Residue const & a, uint64_t const b) const
    {
        assertValid(a);

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        uint64_t u = invm;    // NOLINT(misc-const-correctness)
        uint64_t a0 = a.r[0]; // NOLINT(misc-const-correctness)
        __asm__ __VOLATILE(
            /* Product of low words */
            "mulq %[b]\n\t"         /* rdx:rax = a0*b <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t" /* t0:rax = a0*b <= (2^64-1)^2 */
            /* Compute u such that a0*b + u*m == 0 (mod 2^64), compute u*m, add
               to t0:rax */
            "imulq %[u], %%rax\n\t"
            "movq %%rax, %[u]\n\t" /* u <= 2^64-1 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"       /* rdx:rax = u*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"       /* if low word != 0, carry to high word */
            "movq %[u], %%rax\n\t" /* rax = u, independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t" /* t1:t0 = (a0*b + u*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"  /* rdx:rax = u*m1 */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t" /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = (a0*b+u*m)/2^64 <= 2^126 - 2^62
                                     */
            /* Free slot in pipe 2 here */
            ABORT_IF_CY

            /* Product of low and high word */
            "mulq %[b]\n\t" /* rdx:rax = a1*b <= (2^63-2^32-1)*(2^64-1) */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t" /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
                                    <= ((2^126-1)*(2^64-1) + (2^64-1)*(2^126-1))
                                    / 2^64 < 2^127 - 2^63 - 1, thus no carry */
            ABORT_IF_CY
            /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
                    <= ((m-1)*(2^64-1) + (2^64-1)*m) / 2^64
                     = 2*m*(1-1/2^64) - 1*(1-1/2^64). May need to subtract m. */
            "movq %[t0], %%rax\n\t" /* See if result > m */
            "movq %[t1], %%rdx\n\t"
            "subq %[m0], %[t0]\n\t" /* Try subtracting m, see if there's a carry
                                     */
            "sbbq %[m1], %[t1]\n\t"
            "cmovc %%rax, %[t0]\n\t" /* Carry -> restore old result */
            "cmovc %%rdx, %[t1]\n\t"
            :
            [u] "+&r"(u), [t0] "=&r"(r.r[0]), [t1] "=&r"(r.r[1]), [a0] "+&a"(a0)
            : [a1] "g"(a.r[1]), [b] "rm"(b), [m0] "rm"(m[0]), [m1] "rm"(m[1])
            : "%rdx", "cc");
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */
        uint64_t pl, ph, t[2];

        /* m < 1/4 W^2,  a < m, b < W */

        /* Product of b and low word */
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), a.r[0],
                           b); /* t1:t0 = a0*b < W^2 */

        /* One REDC step */
        redc1(t, t); /* t1:t0 = (a0*b + k*m) / W < m + W < 1/4 W^2 + W */

        /* Product of b and high word  */
        u64arith_mul_1_1_2(&pl, &ph, a.r[1], b);    /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 1/2 W^2 + W */

        /* Result may be larger than m, but is < 2*m */
        u64arith_sub_2_2_ge(&(t[0]), &(t[1]), m[0], m[1]);

        r.r[0] = t[0];
        r.r[1] = t[1];
#endif
        assertValid(r);
    }

  protected:
    /* These functions are used as building blocks in the specific
     * instantiations of the api layer routines.
     */
    bool divn(Residue & r, Residue const & a, uint64_t n,
              uint64_t w_mod_n,
              uint64_t const * inv_n,
              uint64_t c) const;

  private:
    std::vector<Integer> batchinv_redc(std::vector<uint64_t> const & a, Integer c) const;

  public:
    static std::vector<uint64_t> batch_Q_to_Fp(Integer const & num,
            Integer const & den, int k,
            std::vector<uint64_t> const & p);
    struct batch_Q_to_Fp_context;
};

struct arithxx_modredc126::Modulus::batch_Q_to_Fp_context {
    Integer remainder, quotient;
    Modulus D;
    batch_Q_to_Fp_context(Integer const & num, Integer const & den);
    std::vector<uint64_t> operator()(std::vector<uint64_t> const & p, int k=0) const;
};




#endif /* CADO_UTILS_ARITHXX_MODREDC126_HPP */
