#ifndef CADO_UTILS_ARITHXX_MODREDC126_HPP
#define CADO_UTILS_ARITHXX_MODREDC126_HPP

/* A class for modular arithmetic with modulus in [2^64+1, 2^126-1].
 * Moduli must be odd. Residues are stored in Montgomery form, reduction
 * after multiplication is done with REDC. */


#include "cado_config.h" // for HAVE_GCC_STYLE_AMD64_INLINE_ASM

#include <cstdint>

#include <array>
#include <type_traits>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "macros.h"
#include "modint.hpp"
#include "u64arith.h"
#include "arithxx_common.hpp"
#include "arithxx_redc.hpp"
#include "arithxx_redc128.hpp"

struct arithxx_modredc126 {
    class Modulus;
    class Residue;
    typedef Integer128 Integer;

    /* When we multiply by a small constant, we use a left-to-right
     * binary method. So we typically have log(n) shifts and log(n)/2
     * additions, which should be compared to the cost of a runtime
     * multiplication.
     */
    typedef std::integral_constant<int, 8> mul_c_cutoff;

    /* this gives k such that 2^k*modulus-1 <= Integer::max_value
     */
    typedef std::integral_constant<int, 2> overflow_bits;

    typedef std::true_type uses_montgomery_representation;

    typedef std::false_type even_moduli_allowed;
};

/* okay, at this point it's a typedef, really... */
class arithxx_modredc126::Residue
    : public arithxx_details::Residue_base<arithxx_modredc126>
{
    using Residue_base::Residue_base;
};

class arithxx_modredc126::Modulus
    : public arithxx_details::redc128<arithxx_modredc126>
{
    typedef arithxx_modredc126 layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<layer>;
    friend struct arithxx_details::api_bysize<layer>;
    friend struct arithxx_details::redc<layer>;
    friend struct arithxx_details::redc128<layer>;

  protected:

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m)
    {
        return Integer(1, 1) <= m && m % 2 == 1 && !(m[1] >> (64 - layer::overflow_bits::value));
    }
    static bool valid(cxx_mpz const & m) {
        return m % 2 == 1 && m > 0 && mpz_sizeinbase(m, 2) > 64 && mpz_sizeinbase(m, 2) <= (128 - layer::overflow_bits::value);
    }

    explicit Modulus(Integer const & s)
        : redc128(s)
    {
        ASSERT(valid(s));
    }

  protected:
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(Residue const & r MAYBE_UNUSED) const
    {
        ASSERT(r.r < m);
    }
    /* }}} */

  public:

    /* {{{ add(*2) add1 sub(*2) sub1 div2 */
    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        r.r = a.r + b.r;
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
        auto const & me = downcast();
        uint64_t pl, ph, t[3];

        /* m < 1/4 W^2,  a,b < m */

        /* Product of the two low words */
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), a[0], b[0]);

        /* One REDC step */
        me.redc1(t, t); /* t < 2m < 1/2 W^2 */

        /* Products of one low and one high word  */
        u64arith_mul_1_1_2(&pl, &ph, a[1], b[0]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < 3/4 W^2 */
        u64arith_mul_1_1_2(&pl, &ph, a[0], b[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < W^2 */

        /* Product of the two high words */
        u64arith_mul_1_1_2(&pl, &(t[2]), a[1], b[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2(&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        me.redc1_wide_inplace(t);
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
        auto const & me = downcast();
        uint64_t pl, ph, t[3];

        /* m < 1/4 W^2,  a < m */

        /* Square of the low word */
        u64arith_mul_1_1_2(&(t[0]), &(t[1]), a.r[0], a.r[0]);

        /* One REDC step */
        me.redc1(t, t); /* t < 2m < 1/2 W^2 */

        /* Products of low and high word  */
        u64arith_mul_1_1_2(&pl, &ph, a.r[1], a.r[0]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < 3/4 W^2 */
        u64arith_add_2_2(&(t[0]), &(t[1]), pl, ph);   /* t1:t0 < W^2 */

        /* Square of high word */
        u64arith_mul_1_1_2(&pl, &(t[2]), a.r[1], a.r[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2(&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        me.redc1_wide_inplace(t);
        r.r[0] = t[1];
        r.r[1] = t[2];

        /* Result may be larger than m, but is < 2*m */
        u64arith_sub_2_2_ge(r.r.data(), r.r.data() + 1, m[0], m[1]);
#endif
    }
    /* }}} */

#undef ABORT_IF_CY
};

#endif /* CADO_UTILS_ARITHXX_MODREDC126_HPP */
