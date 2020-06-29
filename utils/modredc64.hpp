/* Some functions for modular arithmetic with residues and modulus in
   with up to 64 bits. Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining,
   this file must be included in the caller's source code with #include */

#ifndef MODREDC64_HPP
#define MODREDC64_HPP

/**********************************************************************/
#include <cstdlib>       // for size_t, llabs, NULL
#include <new>            // for operator new
#include <cstdint>
#include "macros.h"
#include "u64arith.h"
#include "modint.hpp"
#include "mod_stdop.hpp"
#include "macros.h"     // ASSERT MAYBE_UNUSED // IWYU pragma: keep


class ModulusREDC64 {
    /* Type definitions */
public:
    typedef Integer64 Integer;
    class Residue {
        friend class ModulusREDC64;
    protected:
        uint64_t r;
    public:
        typedef ModulusREDC64 Modulus;
        typedef Modulus::Integer Integer;
        typedef bool IsResidueType;
        Residue() = delete;
        Residue(const Modulus &m MAYBE_UNUSED) : r(0) {}
        Residue(const Modulus &m MAYBE_UNUSED, const Residue &s) : r(s.r) {}
        Residue(const Residue &&s) : r(s.r) {}
    protected:
        Residue &operator=(const Residue &s) {r = s.r; return *this;}
        Residue &operator=(const Integer &s) {r = 0; s.get(&r, 1); return *this;}
        Residue &operator=(const uint64_t s) {r = s; return *this;}
    };

    typedef ResidueStdOp<Residue> ResidueOp;

protected:
    /* Data members */
    uint64_t m;
    uint64_t invm;
    uint64_t mrecip;
    Residue one;

    /* Methods used internally */
    void assertValid(const Residue &a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a.r < m);}
    void assertValid(const uint64_t a MAYBE_UNUSED) const {ASSERT_EXPENSIVE (a < m);}

    /* Computes (a * 2^64) % m */
    void tomontgomery (Residue &r, const Residue &a) const {
        assertValid (a);
        int shift = u64arith_clz(m);
        uint64_t ml = m << shift, dummy;
        u64arith_divqr_2_1_1_recip_precomp(&dummy, &r.r, 0, a.r, ml, mrecip,
                                           shift);
    }

    /* Computes (a / 2^64) % m. Assumes a < m */

    void frommontgomery (uint64_t &r, const uint64_t a) const {
      uint64_t tlow, thigh;
      assertValid (a);
      tlow = a * invm;
      u64arith_mul_1_1_2 (&tlow, &thigh, tlow, m);
      r = thigh + ((a != 0) ? 1 : 0);
    }

    uint64_t get_u64 (const Residue &s) const {
        uint64_t r;
        assertValid (s);
        frommontgomery (r, s.r);
        return r;
    }

    /* Methods of the API */
public:
    static Integer getminmod() {return Integer(1);}
    static Integer getmaxmod() {return Integer(UINT64_MAX);}
    static void getminmod(Integer &r) {r = getminmod();}
    static void getmaxmod(Integer &r) {r = getmaxmod();}
    static bool valid(const Integer &m) {
        return getminmod() <= m && m <= getmaxmod() && m % 2 == 1;
    }
    /* Methods for the modulus */
    ModulusREDC64 (const uint64_t s) : m(s), invm(-u64arith_invmod (s)), one(*this) {
        int shift = u64arith_clz(m);
        uint64_t ml = m << shift, dummy;
        mrecip = u64arith_reciprocal_for_div(ml);
        u64arith_divqr_2_1_1_recip_precomp(&dummy, &one.r, 0, 1, ml, mrecip,
                                           shift);
    }
    ModulusREDC64(const ModulusREDC64 &s) : m(s.m), invm(s.invm), mrecip(s.mrecip), one(s) {one = s.one;}
    ModulusREDC64 (const Integer &s) : ModulusREDC64(s.getWord(0)) {}
    ~ModulusREDC64 () {}

    uint64_t getmod_u64 () const {return m;}
    void getmod (Integer &r) const {r = Integer(m);}

    /* Methods for residues */

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

    void set (Residue &r, const Residue &s) const {assertValid(s); r = s;}
    /* Puts in r the value of s * beta mod m, where beta is the word base.
       Note: s can be any uint64_t, in particular can be larger than m.
       When 0 <= s < m, use set_reduced for better efficiency. */
    void set (Residue &r, const uint64_t s) const {
      uint64_t plow, phigh;

      u64arith_mul_1_1_2 (&plow, &phigh, s, one.r);
      u64arith_redc (&r.r, plow, phigh, m, invm);
      tomontgomery (r, r);
    }
    void set (Residue &r, const Integer &s) const {set(r, s.getWord(0));}

    /* Sets the residue_t to the class represented by the integer s. Assumes that
       s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced (Residue &r, const uint64_t s) const {
      assertValid (s);
      r.r = s;
      tomontgomery (r, r);
    }
    void set_reduced (Residue &r, const Integer &s) const {set_reduced(r, s.getWord(0));}
    void set_int64 (Residue &r, const int64_t s) const {set(r, llabs(s)); if (s < 0) neg(r, r);}
    void set0 (Residue &r) const {r.r = 0;}
    void set1 (Residue &r) const {r = one;}
    void swap (Residue &a, Residue &b) const {uint64_t t = a.r; a.r = b.r; b.r = t;}

    void get (Integer &r, const Residue &s) const {
        assertValid (s);
        uint64_t t;
        frommontgomery (t, s.r);
        r = Integer(t);
    }

    bool equal (const Residue &a, const Residue &b) const {assertValid(a); assertValid(b); return (a.r == b.r);}
    bool is0 (const Residue &a) const {assertValid(a); return (a.r == 0);}
    bool is1 (const Residue &a) const {return equal(a, one);}
    void neg (Residue &r, const Residue &a) const {
        assertValid(a);
        if (a.r == 0)
            r.r = a.r;
        else
            r.r = m - a.r;
    }
    void add (Residue &r, const Residue &a, const Residue &b) const {u64arith_addmod_1_1(&r.r, a.r, b.r, m);}
    void add1 (Residue &r, const Residue &a) const {add(r, a, one);}
    void add (Residue &r, const Residue &a, const uint64_t b) const {
        Residue t(*this);

        assertValid (a);
        set (t, b);
        add (r, a, t);
    }
    void sub (Residue &r, const Residue &a, const Residue &b) const {
        u64arith_submod_1_1(&r.r, a.r, b.r, m);
    }
    void sub1 (Residue &r, const Residue &a) const {sub(r, a, one);}
    void sub (Residue &r, const Residue &a, const uint64_t b) const {
        Residue t(*this);

        assertValid (a);
        set (t, b);
        sub (r, a, t);
    }

    void mul (Residue &r, const Residue &a, const Residue &b) const {
        uint64_t plow, phigh;

        ASSERT_EXPENSIVE (m % 2 != 0);
        assertValid (a);
        assertValid (b);
        u64arith_mul_1_1_2 (&plow, &phigh, a.r, b.r);
        u64arith_redc (&r.r, plow, phigh, m, invm);
    }

    /* For a residue class a (mod m) and non-negative integer b, set r to
       the smallest non-negative integer in the residue class a*b (mod m). */
    void mul_u64_u64 (uint64_t &r, const Residue &a, const uint64_t b) const {
      uint64_t plow, phigh;

      ASSERT_EXPENSIVE (m % 2 != 0);
      assertValid (a);

      u64arith_mul_1_1_2 (&plow, &phigh, a.r, b);
      /* We have a <= m-1, b <= 2^64 - 1. Thus the product
         phigh:plow <= (m-1)*(2^64 - 1) = m*2^64 - 2^64 - m + 1,
         and with m >= 1,
         phigh:plow <= m*2^64 - 2^64, so phigh < m. */
      u64arith_redc (&r, plow, phigh, m, invm);
    }


    void sqr (Residue &r, const Residue &a) const {
        uint64_t plow, phigh;

        assertValid (a);

        u64arith_sqr_1_2 (&plow, &phigh, a.r);
        u64arith_redc (&r.r, plow, phigh, m, invm);
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

    bool next (Residue &r) const {return (++r.r == m);}
    bool finished (const Residue &r) const {return (r.r == m);}
    bool div2 (Residue &r, const Residue &a) const {r.r = u64arith_div2mod(a.r, m); return 1;}
    bool div3 (Residue &, const Residue &) const;
    bool div5 (Residue &, const Residue &) const;
    bool div7 (Residue &, const Residue &) const;
    bool div11 (Residue &, const Residue &) const;
    bool div13 (Residue &, const Residue &) const;
    void gcd (Integer &, const Residue &) const;
    void pow (Residue &, const Residue &, const uint64_t) const;
    void pow (Residue &, const Residue &, const uint64_t *, const size_t) const;
    void pow (Residue &, const Residue &, const Integer &) const;
    void pow2 (Residue &, const uint64_t) const;
    void pow2 (Residue &, const uint64_t *, const size_t) const;
    void pow2 (Residue &r, const Integer &) const;
    void pow3 (Residue &, const uint64_t) const;
    void V (Residue &r, Residue *rp1, const Residue &b,
            const uint64_t k) const;
    bool sprp (const Residue &) const;
    bool sprp2 () const;
    bool isprime () const;
    bool inv (Residue &, const Residue &) const;
    bool inv_odd (Residue &, const Residue &) const;
    bool intinv (Residue &, const Residue &) const;
    bool batchinv (Residue *, const Residue *, size_t, const Residue *) const;
    bool batchinv_u64 (uint64_t *, const uint64_t *, uint64_t, const size_t) const;
    bool batch_Q_to_Fp (uint64_t *, uint64_t, uint64_t, uint64_t, const uint64_t *, size_t) const;
    int jacobi (const Residue &) const;
protected:
    bool find_minus1 (Residue &r1, const Residue &minusone, const int po2) const;
};
#endif  /* MODREDC64_HPP */
