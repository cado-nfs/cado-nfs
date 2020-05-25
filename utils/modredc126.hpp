#ifndef MODREDC126_HPP
#define MODREDC126_HPP

/* Some functions for modular arithmetic with modulus in [2^64+1, 2^126-1].
   Moduli must be odd. Residues are stored in Montgomery form, reduction
   after multiplication is done with REDC. */

/**********************************************************************/
#include "cado_config.h"  // for HAVE_GCC_STYLE_AMD64_INLINE_ASM
#include <cstdint>
#include <cstddef>       // for size_t, NULL
#include <new>            // for operator new
#include "macros.h"
#include "u64arith.h"
#include "modint.hpp"
#include "mod_stdop.hpp"

class ModulusREDC126 {
    /* Type definitions */
public:
    typedef Integer128 Integer;
    class Residue {
        friend class ModulusREDC126;
    protected:
        uint64_t r[2];
    public:
        typedef ModulusREDC126 Modulus;
        typedef Modulus::Integer Integer;
        typedef bool IsResidueType;
        Residue() = delete;
        Residue(const Modulus &m MAYBE_UNUSED) : r{0,0} {}
        Residue(const Modulus &m MAYBE_UNUSED, const Residue &s) : r{s.r[0], s.r[1]} {}
        Residue(const Residue &&s) : r{s.r[0], s.r[1]} {}
    protected:
        Residue &operator=(const Residue &s) {r[0] = s.r[0]; r[1] = s.r[1]; return *this;}
        Residue &operator=(const Integer &s) {r[1] = 0; s.get(r, 2); return *this;}
        Residue &operator=(const uint64_t s) {r[0] = s; r[1] = 0; return *this;}
    };

    typedef ResidueStdOp<Residue> ResidueOp;

protected:
    /* Data members */
    uint64_t m[2];
    uint64_t invm;
    uint64_t mrecip;
    Residue one;

public:
    static Integer getminmod() {return Integer(1,1);} /* 2^64+1 */
    static Integer getmaxmod() {return Integer(UINT64_MAX, UINT64_MAX >> 2);} /* 2^126-1 */
    static bool valid(const Integer &m) {
        return getminmod() <= m && m <= getmaxmod() && m % 2 == 1;
    }
    ModulusREDC126(const ModulusREDC126 &s) : m{s.m[0], s.m[1]}, invm(s.invm), mrecip(s.mrecip), one(s) {one = s.one;}
    ModulusREDC126 (const Integer s) : one(*this) {
        ASSERT (valid(s));
        s.get (m, 2);
        invm = -u64arith_invmod (m[0]);

        int shift = u64arith_clz(m[1]);
        uint64_t dummy, ml[2] = {m[0], m[1]};
        u64arith_shl_2(&ml[0], &ml[1], shift);
        mrecip = u64arith_reciprocal_for_div_3by2(ml[0], ml[1]);
        u64arith_divqr_3_2_1_recip_precomp(&dummy, &one.r[0], &one.r[1],
            0, 0, 1, ml[0], ml[1], mrecip, shift);
    }

    /* Returns the modulus as an Integer. */
    void getmod(Integer &r) const {
        r.set(m, 2);
    }
    
protected:
    void assertValid(const Residue &r MAYBE_UNUSED) const {ASSERT(u64arith_lt_2_2(r.r[0], r.r[1], m[0], m[1]));}
    
    /* r = s * 2^64 mod m */
    void
    tomontgomery1 (Residue &r, const Residue &s) const
    {
        int shift = u64arith_clz(m[1]);
        uint64_t dummy, ml[2] = {m[0], m[1]};
        u64arith_shl_2(&ml[0], &ml[1], shift);
        u64arith_divqr_3_2_1_recip_precomp(&dummy, &r.r[0], &r.r[1],
            0, s.r[0], s.r[1], ml[0], ml[1], mrecip, shift);
    }

    /* r = s * 2^128 mod m */
    void
    tomontgomery (Residue &r, const Residue &s) const
    {
        tomontgomery1(r, s);
        tomontgomery1(r, r);
    }

    void redc1 (uint64_t *r, const uint64_t *s) const
    {
        uint64_t t[4], k;

        k = s[0] * invm;
        u64arith_mul_1_1_2 (&(t[0]), &(t[1]), k, m[0]);
        if (s[0] != 0)
            t[1]++;
        t[2] = 0;
        u64arith_add_1_2 (&(t[1]), &(t[2]), s[1]);     /* t[2] <= 1 */
        u64arith_mul_1_1_2 (&(t[0]), &(t[3]), k, m[1]);  /* t[3] < 2^w-1 */
        u64arith_add_2_2 (&(t[1]), &(t[2]), t[0], t[3]); /* t[2] < 2^w */

        /* r = (k*m + s) / w, k <= w-1. If s < m, then r < m */
        r[0] = t[1];
        r[1] = t[2];
    }

    /* Do a one-word REDC, i.e., r == s / w (mod m), w = 2^64.
       If m > w, r < 2m. If s < m, then r < m */
    void redc1 (Residue &r, const Residue &s) const
    {
        redc1(r.r, s.r);
    }
    void redc1 (Integer &r, const Integer &s) const
    {
        uint64_t t[2];
        s.get(t, 2);
        redc1(t, t);
        r = Integer(t[0], t[1]);
    }

    /* Converts s out of Montgomery form by dividing by 2^(2*LONG_BIT).
       Requires s < m. */
    void frommontgomery (Residue &r, const Residue &s) const
    {
        /* Do two REDC steps */
        redc1 (r, s);
        redc1 (r, r);
    }

public:
    /** Allocate an array of len residues.
     *
     * Must be freed with deleteArray(), not with delete[].
     */
    Residue *newArray(const size_t len) const {
        void *t = operator new[](len * sizeof(Residue));
        if (t == NULL)
            return NULL;
        Residue *ptr = static_cast<Residue *>(t);
        for(size_t i = 0; i < len; i++) {
            new(&ptr[i]) Residue(*this);
        }
        return ptr;
    }

    void deleteArray(Residue *ptr, const size_t len) const {
        for(size_t i = len; i > 0; i++) {
            ptr[i - 1].~Residue();
        }
        operator delete[](ptr);
    }

    void set(Residue &r, const Residue &s) const {r = s;}
    void set (Residue &r, const uint64_t s) const {
        r.r[0] = s;
        r.r[1] = 0;
        tomontgomery (r, r);
    }
    void set (Residue &r, const Integer &s) const {
        if (s < Integer(m[0], m[1])) {
            s.get(r.r, 2);
        } else {
            Integer t = s;
            t = t % Integer(m[0], m[1]);
            t.get(r.r, 2);
        }
        tomontgomery (r, r);
    }
    /* Sets the residueredc2ul2_t to the class represented by the integer s.
       Assumes that s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced (Residue &r, const uint64_t s) const {set (r, s);}
    void set_reduced (Residue &r, Integer s) const {
        ASSERT (s < Integer(m[0], m[1]));
        s.get(r.r, 2);
        tomontgomery (r, r);
    }
    void set0(Residue &r) const {r.r[0] = 0; r.r[1] = 0;}
    void set1 (Residue &r) const {set (r, one);}
    /* Exchanges the values of the two arguments */

    void swap (Residue &a, Residue &b) const {
        Residue t(*this);
        set (t, a);
        set (a, b);
        set (b, t);
    }

    /* Returns the residue as a modintredc2ul2_t */
    void get (Integer &r, const Residue &s) const {
        Residue t(*this);
        frommontgomery (t, s);
        r.set(t.r, 2);
    }

    bool equal (const Residue &a, const Residue &b) const {
        return (a.r[0] == b.r[0] && a.r[1] == b.r[1]);
    }

    bool is0 (const Residue &a) const {
        return (a.r[0] == 0 && a.r[1] == 0);
    }

    bool is1 (const Residue &a) const {
        return (equal(a, one));
    }

    void add (Residue &r, const Residue &a, const Residue &b) const {
        const uint64_t t0 = b.r[0], t1 = b.r[1]; /* r, a, and/or b may overlap */
        set (r, a);
        u64arith_add_2_2 (&(r.r[0]), &(r.r[1]), t0, t1);
        u64arith_sub_2_2_ge (&(r.r[0]), &(r.r[1]), m[0], m[1]);
    }
    
    void sub (Residue &r, const Residue &a, const Residue &b) const {
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        {
            uint64_t s1 = m[0], s2 = m[1], t1 = a.r[0], t2 = a.r[1];

            __asm__ __VOLATILE (
                "subq %4, %0\n\t"
                "sbbq %5, %1\n\t"    /* t -= b */
                "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
                "cmovncq %6, %3\n"
                : "+&r" (t1), "+&r" (t2), "+&r" (s1), "+r" (s2)
                : "g" (b.r[0]), "g" (b.r[1]), "rm" (0UL)
                : "cc"
            );
            u64arith_add_2_2 (&t1, &t2, s1, s2);
            r.r[0] = t1;
            r.r[1] = t2;
        }
#else
        {
            uint64_t t1 = a.r[0], t2 = a.r[1];
            u64arith_sub_2_2 (&t1, &t2, b.r[0], b.r[1]);

            if (t2 > a.r[1] || (t2 == a.r[1] && t1 > a.r[0]))
                u64arith_add_2_2 (&t1, &t2, m[0], m[1]);

            r.r[0] = t1;
            r.r[1] = t2;
        }
#endif
    }

    void add1 (Residue &r, const Residue &a) const {
        add(r, a, one);
    }

    void sub1 (Residue &r, const Residue &a) const {
        sub(r, a, one);
    }

    void add (Residue &r, const Residue &a, const uint64_t b) const
    {
        Residue t(*this);

        set (t, b);
        add (r, a, t);
    }

    void sub (Residue &r, const Residue &a, const uint64_t b) const {
        Residue t(*this);

        set (t, b);
        sub (r, a, t);
    }

    void neg (Residue &r, const Residue &a) const {
        if (is0 (a)) {
            set0 (r);
        } else {
            uint64_t t0 = m[0], t1 = m[1];
            u64arith_sub_2_2(&t0, &t1, a.r[0], a.r[1]);
            r.r[0] = t0;
            r.r[1] = t1;
        }
    }

    bool div2 (Residue &r, const Residue &a) const {
        set (r, a);
        if (r.r[0] % 2 == 1)
            u64arith_add_2_2 (&r.r[0], &r.r[1], m[0], m[1]);
        u64arith_shr_2(&r.r[0], &r.r[1], 1);
        return true;
    }
    
#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif

    void mul (Residue &r, const Residue &a, const Residue &b) const {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
        uint64_t dummy;
        __asm__ __VOLATILE (
            /* Product of low words */
            "movq %[a0], %%rax\n\t"
            "mulq %[b0]\n\t"         /* rdx:rax = a0*b0 <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t"
            /* Compute u0*m, add to t0:rax */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"         /* if low word != 0, carry to high word */
            "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t"        /* t1:t0 = (a0*b0 + u0*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"
            "addq %%rax, %[t0]\n\t"
            "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 */
            ABORT_IF_CY              /* <= 2^126 - 2^62 */


            /* 2 products of low and high word */
            "xorl %k[t2], %k[t2]\n\t"
            "mulq %[b1]\n\t"         /* rdx:rax = a0*b1 <= (2^64-1)*(2^63-2^32-1) */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 */
            ABORT_IF_CY              /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
    			        = 2^127 + 2^126 - 2^96 ... */

            /* Free slot here */
            "mulq %[b0]\n\t"         /* rdx:rax = a1*b0 <= (2^63-2^32-1)*(2^64-1) */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"
            "setc %b[t2]\n\t"        /* t2:t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 */
            /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                       = 2^128 + 2^126 - 2*2^96 ... */
            /* Product of high words */
            "mulq %[b1]\n\t"         /* rdx:rax = a1*b1 <= (2^63-2^32-1)^2 */
            "addq %%rax, %[t1]\n\t"
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t2]\n\t"  /* t2:t1:t0 = (a*b+u0*m)/2^64 */
            ABORT_IF_CY              /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64
                                = 2^190 - 2^160 ... */
            /* Free slot here */
            /* Compute u1*m, add to t2:t1:t0 */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t0]\n\t" /* t0 = u1 */
            /* Free slot here */
            "mulq %[m0]\n\t"        /* rdx:rax = m0*u1 <= (2^64-1)^2 */
            "negq %%rax\n\t"        /* if low word != 0, carry to high word */
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "adcq $0,%[t2]\n\t"     /* t2:t1:0 = (a*b+u0*m)/2^64 + u1*m0 */
            ABORT_IF_CY             /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */

            "mulq %[m1]\n\t"        /* rdx:rax = u1*m1 */
            "addq %%rax, %[t1]\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*b+u0*m)/2^64 + u1*m)/2^64 */
            ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */

            "movq %[t1], %%rax\n\t" /* See if result > m */
            "movq %[t2], %%rdx\n\t"
            "subq %[m0], %[t1]\n\t"
            "sbbq %[m1], %[t2]\n\t"
            "cmovc %%rax, %[t1]\n\t" /* Carry -> restore old result */
            "cmovc %%rdx, %[t2]\n\t"
            : [t0] "=&r" (dummy), [t1] "=&r" (r.r[0]), [t2] "=&r" (r.r[1])
            : [a0] "g" (a.r[0]), [a1] "g" (a.r[1]), [b0] "rm" (b.r[0]), [b1] "rm" (b.r[1]),
            [m0] "rm" (m[0]), [m1] "rm" (m[1]), [invm] "rm" (invm)
            : "%rax", "%rdx", "cc"
        );
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */
        uint64_t pl, ph, t[4], k;

        /* m < 1/4 W^2,  a,b < m */

        /* Product of the two low words */
        u64arith_mul_1_1_2 (&(t[0]), &(t[1]), a.r[0], b.r[0]);

        /* One REDC step */
        redc1 (t, t); /* t < 2m < 1/2 W^2 */

        /* Products of one low and one high word  */
        u64arith_mul_1_1_2 (&pl, &ph, a.r[1], b.r[0]);   /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 3/4 W^2 */
        u64arith_mul_1_1_2 (&pl, &ph, a.r[0], b.r[1]);   /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < W^2 */

        /* Product of the two high words */
        u64arith_mul_1_1_2 (&pl, &(t[2]), a.r[1], b.r[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2 (&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        k = t[0] * invm;
        u64arith_mul_1_1_2 (&pl, &ph, k, m[0]);
        if (t[0] != 0UL)
            ph++; /* t[0] = 0 */
        u64arith_add_1_2 (&(t[1]), &(t[2]), ph);
        u64arith_mul_1_1_2 (&pl, &ph, k, m[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[1]), &(t[2]), pl, ph);
        /* t2:t1:0 < 1/16 W^3 + W^2 + 1/4 W^3 < 5/16 W^3 + W^2 */

        /* Result may be larger than m, but is < 2*m */

        u64arith_sub_2_2_ge (&(t[1]), &(t[2]), m[0], m[1]);

        r.r[0] = t[1];
        r.r[1] = t[2];
#endif
    }

    void sqr (Residue &r, const Residue &a) const {
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
        __asm__ __VOLATILE (
            /* Product of low words */
            "movq %[a0], %%rax\n\t"
            "mulq %[a0]\n\t"         /* rdx:rax = a0^2 <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t"
            /* Compute u0*m, add to t0:rax */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"         /* if low word != 0, carry to high word */
            "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t"        /* t1:t0 = (a0^2 + u0*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"
            "addq %%rax, %[t0]\n\t"
            "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0^2+u0*m)/2^64 */
            ABORT_IF_CY              /* <= 2^126 - 2^62 */

            /* 2 products of low and high word */
            "xorl %k[t2], %k[t2]\n\t"
            "mulq %[a1]\n\t"         /* rdx:rax = a0*a1 <= (2^64-1)*(2^63-2^32-1) */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0^2+u0*m)/2^64 + a0*a1 */
            ABORT_IF_CY              /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
    			        = 2^127 + 2^126 - 2^96 ... */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
            "setc %b[t2]\n\t"        /* t2:t1:t0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 */
            /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                       = 2^128 + 2^126 - 2*2^96 ... */

            /* Product of high words */
            "mulq %[a1]\n\t"         /* rdx:rax = a1^2 <= (2^63-2^32-1)^2 */
            "addq %%rax, %[t1]\n\t"
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t2]\n\t"  /* t2:t1:t0 = (a^2+u0*m)/2^64 */
            ABORT_IF_CY              /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64
                                = 2^190 - 2^160 ... */
            /* Free slot here */
            /* Compute u1*m, add to t2:t1:t0 */
            "imulq %[invm], %%rax\n\t"
            "movq %%rax, %[t0]\n\t" /* t0 = u1 */
            /* Free slot here */
            "mulq %[m0]\n\t"        /* rdx:rax = m0*u1 <= (2^64-1)^2 */
            "negq %%rax\n\t"        /* if low word != 0, carry to high word */
            "movq %[t0], %%rax\n\t"
            "adcq %%rdx, %[t1]\n\t"
            "adcq $0,%[t2]\n\t"     /* t2:t1:0 = (a*a+u0*m)/2^64 + u1*m0 */
            ABORT_IF_CY             /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */

            "mulq %[m1]\n\t"        /* rdx:rax = u1*m1 */
            "addq %%rax, %[t1]\n\t"
            "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*a+u0*m)/2^64 + u1*m)/2^64 */
            ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */

            "movq %[t1], %%rax\n\t" /* See if result > m */
            "movq %[t2], %%rdx\n\t"
            "subq %[m0], %[t1]\n\t"
            "sbbq %[m1], %[t2]\n\t"
            "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
            "cmovc %%rdx, %[t2]\n\t"
            : [t0] "=&r" (dummy), [t1] "=&r" (r.r[0]), [t2] "=&r" (r.r[1])
            : [a0] "g" (a.r[0]), [a1] "g" (a.r[1]),
            [m0] "rm" (m[0]), [m1] "rm" (m[1]), [invm] "rm" (invm)
            : "%rax", "%rdx", "cc"
        );
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */

        uint64_t pl, ph, t[4], k;

        /* m < 1/4 W^2,  a < m */

        /* Square of the low word */
        u64arith_mul_1_1_2 (&(t[0]), &(t[1]), a.r[0], a.r[0]);

        /* One REDC step */
        redc1 (t, t); /* t < 2m < 1/2 W^2 */

        /* Products of low and high word  */
        u64arith_mul_1_1_2 (&pl, &ph, a.r[1], a.r[0]);   /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 3/4 W^2 */
        u64arith_add_2_2 (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < W^2 */

        /* Square of high word */
        u64arith_mul_1_1_2 (&pl, &(t[2]), a.r[1], a.r[1]); /* t2:pl < 1/16 W^2 */
        u64arith_add_1_2 (&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

        /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
        k = t[0] * invm;
        u64arith_mul_1_1_2 (&pl, &ph, k, m[0]);
        if (t[0] != 0UL)
            ph++; /* t[0] = 0 */
        u64arith_add_1_2 (&(t[1]), &(t[2]), ph);
        u64arith_mul_1_1_2 (&pl, &ph, k, m[1]); /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[1]), &(t[2]), pl, ph);
        /* t2:t1:0 < 1/16 W^3 + W^2 + 1/4 W^3 < 5/16 W^3 + W^2 */

        /* Result may be larger than m, but is < 2*m */

        u64arith_sub_2_2_ge (&(t[1]), &(t[2]), m[0], m[1]);

        r.r[0] = t[1];
        r.r[1] = t[2];
#endif
    }

    bool next (Residue &r) const {
        u64arith_add_1_2 (&(r.r[0]), &(r.r[1]), 1);
        return (finished (r));
    }

    bool finished (const Residue &r) const {
        return (r.r[0] == m[0] && r.r[1] == m[1]);
    }

    /* Division by small integer n, where (n-1)*m may overflow the most
       significant word. Returns 1 if n is invertible modulo m, 0 if not.

       w_mod_n is word base (e.g., 2^32 or  2^64) mod n
       inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime
       to n, for 0 <= i < n
       c = n^(-1) (mod word base)
    */

    bool divn (Residue &r, const Residue &a,
                      const uint64_t n, const uint64_t w_mod_n,
                      const uint64_t *inv_n, const uint64_t c) const {
        const uint64_t an = ((a.r[1] % n) * w_mod_n + a.r[0] % n) % n;
        const uint64_t mn = ((m[1] % n) * w_mod_n + m[0] % n) % n;

        Residue t(*this), t2(*this);
        uint64_t k;
#ifdef WANT_ASSERT_EXPENSIVE
        Residue a_backup(*this);
#endif

        if (inv_n[mn] == 0)
            return false;

#ifdef WANT_ASSERT_EXPENSIVE
        set (a_backup, a);
#endif
        set (t, a);

        /* Make t[1]:t[0] == a+km (mod w^2) with a+km divisible by n */
        /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
        k = (inv_n[mn] * an) % n;
        ASSERT_EXPENSIVE ((an + k*mn) % n == 0);
        u64arith_mul_1_1_2 (&(t.r[0]), &(t.r[1]), m[0], k);
        t.r[1] += m[1] * k;
        u64arith_add_2_2 (&(t.r[0]), &(t.r[1]), a.r[0], a.r[1]);

        /* We want r = (a+km)/n. */

        /* May overwrite a */
        r.r[0] = t.r[0] * c;

        /* r0 == (a+km)/n (mod w)
           (r1*w + r0) * n = (a+km)
           (r1*w + r0) * n == t (mod w^2)
           r1*w*n == t - n*r0 (mod w^2)
                                  t - n*r0 == 0 (mod w), thus
           r1*n == (t - n*r0)/w (mod w) */

        u64arith_mul_1_1_2 (&(t2.r[0]), &(t2.r[1]), r.r[0], n);
        u64arith_sub_2_2 (&(t.r[0]), &(t.r[1]), t2.r[0], t2.r[1]);
        ASSERT_EXPENSIVE (t.r[0] == 0);
        r.r[1] = t.r[1] * c;

#ifdef WANT_ASSERT_EXPENSIVE
        {
            uint64_t i;
            set (t, r);
            for (i = 1; i < n; i++)
                add (t, t, r);
            ASSERT_EXPENSIVE (equal (t, a_backup));
        }
#endif

        return true;
    }

    /* Given a = V_n (x), b = V_m (x) and d = V_{n-m} (x), compute V_{m+n} (x).
     * r can be the same variable as a or b but must not be the same variable as d.
     */
    void V_dadd (Residue &r, const Residue &a, const Residue &b,
                     const Residue &d) const {
        ASSERT (&r != &d);
        mul (r, a, b);
        sub (r, r, d);
    }

    /* Given a = V_n (x) and two = 2, compute V_{2n} (x).
     * r can be the same variable as a but must not be the same variable as two.
     */
    void V_dbl (Residue &r, const Residue &a, const Residue &two) const {
        ASSERT (&r != &two);
        sqr (r, a);
        sub (r, r, two);
    }

    
    /* prototypes of non-inline functions */
    bool div3 (Residue &, const Residue &) const;
    bool div5 (Residue &, const Residue &) const;
    bool div7 (Residue &, const Residue &) const;
    bool div11 (Residue &, const Residue &) const;
    bool div13 (Residue &, const Residue &) const;
    void gcd (Integer &, const Residue &) const;
    void pow (Residue &, const Residue &, uint64_t) const;
    void pow (Residue &, const Residue &, const uint64_t *, size_t) const;
    void pow (Residue &r, const Residue &b, const Integer &e) const;
    void pow2 (Residue &, uint64_t) const;
    void pow2 (Residue &, const uint64_t *, size_t) const;
    void pow2 (Residue &r, const Integer &e) const;
    void V (Residue &, const Residue &, const uint64_t) const;
    void V (Residue &, const Residue &, const uint64_t *, size_t) const;
    void V (Residue &r, const Residue &b, const Integer &e) const;
    void V (Residue &r, Residue *rp1, const Residue &b,
            const uint64_t k) const;
    bool sprp (const Residue &) const;
    bool sprp2 () const;
    bool isprime () const;
    bool inv (Residue &, const Residue &) const;
    bool batchinv (Residue *, const Residue *, size_t, const Residue *) const;
    bool batchinv (Residue *, const uint64_t *, size_t, const Integer *) const;
    struct batch_Q_to_Fp_context_s {
        Integer c;
        uint64_t rem_ul, ratio_ul, den_inv;
        ModulusREDC126 &m;
    };
    typedef struct batch_Q_to_Fp_context_s batch_Q_to_Fp_context_t;

    batch_Q_to_Fp_context_t *
    batch_Q_to_Fp_init (const Integer &, const Integer &) const;
    void batch_Q_to_Fp_clear (batch_Q_to_Fp_context_t *) const;

    bool batch_Q_to_Fp (uint64_t *, const batch_Q_to_Fp_context_t *,
                        uint64_t, int, const uint64_t *, size_t) const;
    int jacobi (const Residue &) const;
protected:
    bool find_minus1 (Residue &, const Residue &, const int) const;
    /* Computes r = (a * b * 2^-64) mod m, where a is in REDC
     * representation */
    void
    mul_ul (Residue &r, const Residue &a, const uint64_t b) const
    {
        assertValid(a);

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        uint64_t u = invm, a0 = a.r[0];
        __asm__ __VOLATILE (
            /* Product of low words */
            "mulq %[b]\n\t"          /* rdx:rax = a0*b <= (2^64-1)^2 */
            "movq %%rdx, %[t0]\n\t"  /* t0:rax = a0*b <= (2^64-1)^2 */
            /* Compute u such that a0*b + u*m == 0 (mod 2^64), compute u*m, add to t0:rax */
            "imulq %[u], %%rax\n\t"
            "movq %%rax, %[u]\n\t"  /* u <= 2^64-1 */
            "xorl %k[t1], %k[t1]\n\t"
            "mulq %[m0]\n\t"         /* rdx:rax = u*m0 <= (2^64-1)^2 */
            "negq %%rax\n\t"         /* if low word != 0, carry to high word */
            "movq %[u], %%rax\n\t"  /* rax = u, independent, goes in pipe 0 */
            "adcq %%rdx, %[t0]\n\t"
            "setc %b[t1]\n\t"        /* t1:t0 = (a0*b + u*m0) / 2^64 <= 2*2^64 - 4 */
            "mulq %[m1]\n\t"         /* rdx:rax = u*m1 */
            "addq %%rax, %[t0]\n\t"
            "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b+u*m)/2^64 <= 2^126 - 2^62 */
            /* Free slot in pipe 2 here */
            ABORT_IF_CY

            /* Product of low and high word */
            "mulq %[b]\n\t"          /* rdx:rax = a1*b <= (2^63-2^32-1)*(2^64-1) */
            "addq %%rax, %[t0]\n\t"
            "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
                                     <= ((2^126-1)*(2^64-1) + (2^64-1)*(2^126-1)) / 2^64
                                      < 2^127 - 2^63 - 1, thus no carry */
            ABORT_IF_CY
            /* t1:t0 = ((a1*2^64 + a0)*b + u*m) / 2^64
                    <= ((m-1)*(2^64-1) + (2^64-1)*m) / 2^64
                     = 2*m*(1-1/2^64) - 1*(1-1/2^64). May need to subtract m. */
            "movq %[t0], %%rax\n\t" /* See if result > m */
            "movq %[t1], %%rdx\n\t"
            "subq %[m0], %[t0]\n\t" /* Try subtracting m, see if there's a carry */
            "sbbq %[m1], %[t1]\n\t"
            "cmovc %%rax, %[t0]\n\t" /* Carry -> restore old result */
            "cmovc %%rdx, %[t1]\n\t"
            : [u] "+&r" (u), [t0] "=&r" (r.r[0]), [t1] "=&r" (r.r[1]), [a0] "+&a" (a0)
            : [a1] "g" (a.r[1]), [b] "rm" (b),
            [m0] "rm" (m[0]), [m1] "rm" (m[1])
            : "%rdx", "cc"
        );
#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */
        uint64_t pl, ph, t[2];

        /* m < 1/4 W^2,  a < m, b < W */

        /* Product of b and low word */
        u64arith_mul_1_1_2 (&(t[0]), &(t[1]), a.r[0], b); /* t1:t0 = a0*b < W^2 */

        /* One REDC step */
        redc1 (t, t); /* t1:t0 = (a0*b + k*m) / W < m + W < 1/4 W^2 + W */

        /* Product of b and high word  */
        u64arith_mul_1_1_2 (&pl, &ph, a.r[1], b);   /* ph:pl < 1/4 W^2 */
        u64arith_add_2_2 (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 1/2 W^2 + W */

        u64arith_sub_2_2 (&(t[0]), &(t[1]), m[0], m[1]);

        r.r[0] = t[0];
        r.r[1] = t[1];
#endif
        assertValid(r);
    }

};


#endif  /* MODREDC126_HPP */
