/* A class for modular arithmetic with residues and modulus of up to 64
 * bits. */

#ifndef MODMPZ_HPP
#define MODMPZ_HPP

/**********************************************************************/
#include <cstdint>
#include <gmp.h>          // for mp_limb_t, mpz_size, __gmpn_copyi, __gmpn_s...
#include <cstddef>        // for size_t, NULL
#include <new>            // for operator new
#include "macros.h"
#include "gmp_aux.h"      // for mpz_cmp_uint64, mpz_get_uint64, mpz_set_uint64
#include "cxx_mpz.hpp"
// #include "modint.hpp"
#include "mod_stdop.hpp"
#include "misc.h"

class ModulusMPZ {
    /* Type definitions */
public:
    typedef cxx_mpz Integer;
    class Residue {
        friend class ModulusMPZ;
    protected:
        mp_limb_t *r;
    public:
        typedef ModulusMPZ Modulus;
        typedef Modulus::Integer Integer;
        typedef bool IsResidueType;
        Residue() = delete;
        Residue(const Modulus &m) {
            r = new mp_limb_t[mpz_size(m.m) + 1];
            mpn_zero(r, mpz_size(m.m) + 1);
        }
        Residue(const Modulus &m, const Residue &s) {
            r = new mp_limb_t[mpz_size(m.m) + 1];
            mpn_copyi(r, s.r, mpz_size(m.m));
        }
        ~Residue() {
            delete[] r;
        }
        Residue(Residue &&s) : r(s.r) {
            s.r = NULL;
        }
        Residue(Residue const & s) = delete;
        Residue& operator=(Residue &&s) {
            delete[] r;
            r = s.r;
            s.r = NULL;
            return *this;
        }
        Residue& operator=(Residue const & s) = delete;
    };

    typedef ResidueStdOp<Residue> ResidueOp;

    /* Data members */
protected:
    mpz_t m;
    size_t bits;
    static const size_t limbsPerUint64 = iceildiv(64, GMP_NUMB_BITS);
    
    /** Convert a uint64_t to an array of mp_limb_t.
     * Always writes exactly limbsToWrite limbs. */
    static void uint64ToLimbs(mp_limb_t *r, const uint64_t s, const size_t limbsToWrite) {
        uint64_t t = s;
        for (size_t i = 0; i < limbsToWrite; i++) {
            r[i] = t & GMP_NUMB_MASK;
#if GMP_NUMB_BITS < 64
                t >>= GMP_NUMB_BITS;
#else
                t = 0;
#endif
        }
        ASSERT_ALWAYS(t == 0);
    }
    
    /* Methods used internally */
    /** Set a residue to m.
     *  This leaves the residue in a non-reduced state. */
    void setM(mp_limb_t *r) const {
        mpn_copyi(r, m->_mp_d, mpz_size(m));
    }
    /** Add M to a residue, return carry */
    mp_limb_t addM(mp_limb_t *r, const mp_limb_t *s) const {
        return mpn_add_n(r, s, m->_mp_d, mpz_size(m));
    }
    /** Subtract M from a residue, return borrow */
    mp_limb_t subM(mp_limb_t *r, const mp_limb_t *s) const {
        return mpn_sub_n(r, s, m->_mp_d, mpz_size(m));
    }
    /** Returns a negative value if s < m, 0 if s == m, and a positive value if s > m. */
    int cmpM(const mp_limb_t *s) const {
        return mpn_cmp(s, m->_mp_d, mpz_size(m));
    }
    uint64_t modM(const uint64_t s) const {
        if (bits > 64 || mpz_cmp_uint64(m, s) > 0)
            return s;
        /* m <= s, thus m fits in uint64_t */
        uint64_t m64 = mpz_get_uint64(m);
        /* We always make sure that the modulus in non-zero in the ctor
         */
        ASSERT_FOR_STATIC_ANALYZER(m64 != 0);
        return s % m64;
    }
    void set_residue_u64(Residue &r, const uint64_t s) const {
        ASSERT(mpz_cmp_uint64(m, s) > 0);
        const size_t limbsToWrite = MIN(limbsPerUint64, mpz_size(m));
        uint64ToLimbs(r.r, s, limbsToWrite);
        for (size_t i = limbsToWrite; i < mpz_size(m); i++)
            r.r[i] = 0;
    }
    void set_residue_mpz(Residue &r, const mpz_t s) const {
        ASSERT_ALWAYS(mpz_cmp(s, m) < 0);
        size_t written;
        mpz_export(r.r, &written, -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, s);
        ASSERT_ALWAYS(written <= mpz_size(m));
        for (size_t i = written; i < mpz_size(m); i++)
            r.r[i] = 0;
    }
    void set_mpz_residue(mpz_t r, const Residue &s) const {
        mpz_import(r, mpz_size(m), -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, s.r);
    }
    void assertValid(const Residue &a MAYBE_UNUSED) const {
        ASSERT_EXPENSIVE (cmpM(a.r) < 0);
    }
    void assertValid(const uint64_t a MAYBE_UNUSED) const {
        ASSERT_EXPENSIVE (mpz_cmp_uint64(m, a) > 0);
    }
    void assertValid(const mpz_t s MAYBE_UNUSED) const {
        ASSERT_EXPENSIVE (mpz_cmp(s, m) < 0);
    }
    void assertValid(const cxx_mpz &s MAYBE_UNUSED) const {
        assertValid((mpz_srcptr) s);
    }

    /* Methods of the API */
public:
    static void getminmod(Integer &r) {
        mpz_set_ui(r, 0);
    }
    static void getmaxmod(Integer &r) {
        mpz_set_ui(r, 0);
    }
    static bool valid(const Integer &m MAYBE_UNUSED) {
        return true;
    }

    ModulusMPZ(const uint64_t s) {
        ASSERT_ALWAYS(s > 0);
        mpz_init(m);
        mpz_set_uint64(m, s);
        bits = mpz_sizeinbase(m, 2);
    }
    ModulusMPZ(const ModulusMPZ &s) {
        mpz_init_set(m, s.m);
        bits = s.bits;
    }
    ModulusMPZ(const Integer &s) {
        ASSERT_ALWAYS(mpz_sgn(s) > 0);
        mpz_init(m);
        mpz_set(m, s);
        bits = mpz_sizeinbase(m, 2);
    }
    ~ModulusMPZ() {
        mpz_clear(m);
    }
    void getmod (Integer &r) const {
        mpz_set(r, m);
    }

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


    void set (Residue &r, const Residue &s) const {
        assertValid(s);
        mpn_copyi(r.r, s.r, mpz_size(m));
    }
    void set (Residue &r, const uint64_t s) const {
        const uint64_t sm = modM(s);
        set_residue_u64(r, sm);
    }
    void set (Residue &r, const Integer &s) const {
        cxx_mpz t;
        mpz_mod(t, s, m);
        set_residue_mpz(r, t);
    }
    /* Sets the Residue to the class represented by the integer s. Assumes that
       s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced (Residue &r, const uint64_t s) const {
        assertValid(s);
        set_residue_u64(r, s);
    }
    void set_reduced (Residue &r, const Integer &s) const {
        assertValid(r);
        set_residue_mpz(r, s);
    }
    void set_int64 (Residue &r, const int64_t s) const {
        const uint64_t u = modM(safe_abs64(s));
        set_residue_u64(r, u);
        if (s < 0)
            neg(r, r);
    }
    void set0 (Residue &r) const {
        mpn_zero(r.r, mpz_size(m));
    }
    void set1 (Residue &r) const {
        set0(r);
        if (mpz_cmp_ui(m, 1) > 0) {
            r.r[0] = 1;
        }
    }
    /* Exchanges the values of the two arguments */
    void swap (Residue &a, Residue &b) const {
        mp_limb_t *t = a.r;
        a.r = b.r;
        b.r = t;
    }
    void get (Integer &r, const Residue &s) const {
        assertValid(s);
        set_mpz_residue(r, s);
    }
    bool equal (const Residue &a, const Residue &b) const {
        assertValid(a);
        assertValid(b);
        return mpn_cmp(a.r, b.r, mpz_size(m)) == 0;
    }
    bool is0 (const Residue &a) const {
        assertValid(a);
        return mpn_zero_p(a.r, mpz_size(m));
    }
    bool is1 (const Residue &a) const {
        assertValid(a);
        if (mpz_cmp_ui(m, 1) == 0)
            return 1;
        return a.r[0] == 1 && (mpz_size(m) == 1 || mpn_zero_p(a.r + 1, mpz_size(m) - 1));
    }
    void neg (Residue &r, const Residue &a) const {
        assertValid(a);
        if (is0(a) != 0) {
            if (&r != &a)
                set0(r);
        } else {
            mp_limb_t bw = mpn_sub_n(r.r, m->_mp_d, r.r, mpz_size(m));
            ASSERT_ALWAYS(bw == 0);
        }
        assertValid(r);
    }
    void add (Residue &r, const Residue &a, const Residue &b) const {
        assertValid(a);
        assertValid(b);
        const mp_limb_t cy = mpn_add_n(r.r, a.r, b.r, mpz_size(m));
        if (cy || cmpM(r.r) >= 0) {
            const mp_limb_t bw = subM(r.r, r.r);
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void add1 (Residue &r, const Residue &a) const {
        assertValid(a);
        const mp_limb_t cy = mpn_add_1(r.r, a.r, 1, mpz_size(m));
        if (cy || cmpM(r.r) >= 0) {
            const mp_limb_t bw = subM(r.r, r.r);
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void add (Residue &r, const Residue &a, const uint64_t b) const {
        assertValid(a);
        const uint64_t bm = modM(b);
        mp_limb_t cy;
        
        if (limbsPerUint64 == 1 || bm < GMP_NUMB_MAX) {
            cy = mpn_add_1(r.r, a.r, (mp_limb_t) bm, mpz_size(m));
        } else {
            mp_limb_t t[limbsPerUint64];
            const size_t toWrite = MIN(limbsPerUint64, mpz_size(m));

            uint64ToLimbs(t, bm, limbsPerUint64);
            cy = mpn_add(r.r, a.r, mpz_size(m), t, toWrite);
        }
        if (cy || cmpM(r.r) >= 0) {
            mp_limb_t bw = subM(r.r, r.r);
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void sub (Residue &r, const Residue &a, const Residue &b) const {
        assertValid(a);
        assertValid(b);
        const mp_limb_t bw = mpn_sub_n(r.r, a.r, b.r, mpz_size(m));
        if (bw) {
            const mp_limb_t cy = addM(r.r, r.r);
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    void sub1 (Residue &r, const Residue &a) const {
        assertValid(a);
        const mp_limb_t bw = mpn_sub_1(r.r, a.r, 1, mpz_size(m));
        if (bw) {
            const mp_limb_t cy = addM(r.r, r.r);
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    void sub (Residue &r, const Residue &a, const uint64_t b) const {
        assertValid(a);
        const uint64_t bm = modM(b);
        mp_limb_t bw;
        
        if (limbsPerUint64 == 1 || bm < GMP_NUMB_MAX) {
            bw = mpn_sub_1(r.r, a.r, (mp_limb_t) bm, mpz_size(m));
        } else {
            mp_limb_t t[limbsPerUint64];
            const size_t toWrite = MIN(limbsPerUint64, mpz_size(m));

            uint64ToLimbs(t, bm, limbsPerUint64);
            bw = mpn_sub(r.r, a.r, mpz_size(m), t, toWrite);
        }
        if (bw) {
            mp_limb_t cy = addM(r.r, r.r);
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    void mul (Residue &r, const Residue &a, const Residue &b) const {
        const size_t nrWords = mpz_size(m);
        mp_limb_t Q[2];
        mp_limb_t *t;
        if (r.r == a.r || r.r == b.r) {
            t = new mp_limb_t[nrWords + 1];
        } else {
            t = r.r;
        }
        t[nrWords] = mpn_mul_1 (t, a.r, nrWords, b.r[nrWords - 1]);
        if (t[nrWords] != 0)
            mpn_tdiv_qr (Q, t, 0, t, nrWords + 1, m->_mp_d, nrWords);
        /* t <= (m-1) * beta */
        for (size_t iWord = nrWords - 1; iWord > 0; iWord--) {
            mpn_copyd(t + 1, t, nrWords);
            t[0] = 0;
            mp_limb_t msw = mpn_addmul_1 (t, a.r, nrWords, b.r[iWord - 1]);
            t[nrWords] += msw;
            mp_limb_t cy = t[nrWords] < msw;
            if (cy) {
                mp_limb_t bw = subM(t + 1, t + 1);
                ASSERT_ALWAYS(bw == cy);
            }
            if (t[nrWords] != 0)
                mpn_tdiv_qr (Q, t, 0, t, nrWords + 1, m->_mp_d, nrWords);
        }
        if (cmpM(t) >= 0) {
            mpn_tdiv_qr (Q, t, 0, t, nrWords, m->_mp_d, nrWords);
        }
        if (r.r == a.r || r.r == b.r) {
            mpn_copyi(r.r, t, nrWords);
            delete[] t;
        }
    }
    void sqr (Residue &r, const Residue &a) const {
        mul(r, a, a);
    }
    bool next (Residue &r) const {
        add1(r, r);
        return finished(r);
    }
    bool finished (const Residue &r) const {
        return is0(r);
    }
    bool div2 (Residue &r, const Residue &a) const {
        assertValid(a);
        if (mpz_even_p(m)) {
            return false;
        } else {
            if (a.r[0] % 2 == 0) {
                mp_limb_t lsb = mpn_rshift(r.r, a.r, mpz_size(m), 1);
                ASSERT_ALWAYS(lsb == 0);
            } else {
                mp_limb_t cy = addM(r.r, a.r);
                mp_limb_t lsb = mpn_rshift(r.r, r.r, mpz_size(m), 1);
                ASSERT_ALWAYS(lsb == 0);
                r.r[mpz_size(m) - 1] |= cy << (GMP_NUMB_BITS - 1);
            }
            assertValid(r);
            return true;
        }
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
    bool divn (Residue &, const Residue &, unsigned long) const;
    void gcd (Integer &, const Residue &) const;
    void pow (Residue &, const Residue &, const uint64_t) const;
    void pow (Residue &, const Residue &, const uint64_t *, const size_t) const;
    void pow (Residue &r, const Residue &b, const Integer &e) const;
    void pow2 (Residue &, const uint64_t) const;
    void pow2 (Residue &, const uint64_t *, const size_t) const;
    void pow2 (Residue &r, const Integer &e) const;
    void pow3 (Residue &, uint64_t) const;
    void V (Residue &, const Residue &, const uint64_t) const;
    void V (Residue &, const Residue &, const uint64_t *, const int) const;
    void V (Residue &r, const Residue &b, const Integer &e) const;
    void V (Residue &r, Residue *rp1, const Residue &b,
            const uint64_t k) const;
    bool sprp (const Residue &) const;
    bool sprp2 () const;
    bool isprime () const;
    bool inv (Residue &, const Residue &) const;
    bool inv_odd (Residue &, const Residue &) const;
    bool inv_powerof2 (Residue &, const Residue &) const;
    bool batchinv (Residue *, const Residue *, size_t, const Residue *) const;
    int jacobi (const Residue &) const;
protected:
    bool find_minus1 (Residue &r1, const Residue &minusone, const int po2) const;
    bool divn (Residue &, const Residue &, unsigned long, const mp_limb_t *, mp_limb_t) const;
};

#endif  /* MOD64_HPP */
