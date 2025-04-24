#ifndef CADO_UTILS_ARITHXX_MOD_MPZ_NEW_HPP
#define CADO_UTILS_ARITHXX_MOD_MPZ_NEW_HPP

/* A class for modular arithmetic with residues and modulus or arbitrary
 * size.
 */

#include <cstddef>
#include <cstdint>
#include <climits>

#include <memory>
#include <type_traits>
#include <utility>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "misc.h"
#include "arithxx_common.hpp"


struct arithxx_mod_mpz_new {
    class Modulus;
    class Residue;
    typedef cxx_mpz Integer;

    /* These two are not used for this layer, but for completeness... */
    typedef std::integral_constant<int, 4> mul_c_cutoff;
    typedef std::integral_constant<int, INT_MAX> overflow_bits;

    typedef std::false_type uses_montgomery_representation;

    typedef std::true_type even_moduli_allowed;
};

/* unlike what happens with the other types, here we're *not* deriving
 * from Residue_base, because our underling field is not an Integer
 * object but rather a raw pointer.
 */
class arithxx_mod_mpz_new::Residue
{
    typedef arithxx_mod_mpz_new layer;
    friend class layer::Modulus;
    friend struct arithxx_details::api<arithxx_mod_mpz_new>;

#if 0
    protected:
    std::unique_ptr<mp_limb_t[]> r_owned;
    mp_limb_t * r;

    public:
    /* These two are implemented as inlines at the end of this header file
     * (we must do so because the Modulus type isn't complete yet)
     */
    explicit Residue(Modulus const & m);
    Residue(Modulus const & m, Residue const & s);
    Residue(Residue const & o) = delete;
    Residue(Residue && o) noexcept : r(o.r) { o.r = nullptr; }
    Residue& operator=(Residue const & o) = delete;
    Residue& operator=(Residue && o) noexcept { std::swap(r, o.r); return *this; }
    ~Residue() { delete[] r; }
    Residue() = delete;
#else
    protected:
    std::unique_ptr<mp_limb_t[]> r;

    public:
    /* These two are implemented as inlines at the end of this header file
     * (we must do so because the Modulus type isn't complete yet)
     */
    explicit Residue(Modulus const & m);
    Residue(Modulus const & m, Residue const & s);
#endif
};

class arithxx_mod_mpz_new::Modulus
    : public arithxx_details::api<arithxx_mod_mpz_new>

{
    typedef arithxx_mod_mpz_new layer;
    friend class layer::Residue;

    friend struct arithxx_details::api<arithxx_mod_mpz_new>;

  protected:
    static constexpr mp_size_t limbsPerUint64 = iceildiv(64, GMP_NUMB_BITS);

    /* {{{ ctors, validity range, and asserts */
  public:
    static bool valid(Integer const & m MAYBE_UNUSED) {
        return m > 0 && mpz_odd_p(m);
    }

    explicit Modulus(Integer const & s)
        : api(s)
    {
        ASSERT_ALWAYS(mpz_sgn(s) > 0);
    }
    explicit Modulus(uint64_t const s)
        : Modulus(Integer(s))
    {
        ASSERT_ALWAYS(s > 0);
    }

  protected:
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(Residue const & a MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(cmpM(a.r) < 0);
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(uint64_t const a MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(mpz_cmp_uint64(m, a) > 0);
    }
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void assertValid(mpz_srcptr s MAYBE_UNUSED) const
    {
        ASSERT_EXPENSIVE(mpz_cmp(s, m) < 0);
    }
    void assertValid(cxx_mpz const & s MAYBE_UNUSED) const
    {
        assertValid(mpz_srcptr(s));
    }
    /* }}} */

  protected:

    /** Convert a uint64_t to an array of mp_limb_t.
     * Always writes exactly limbsToWrite limbs. */
    static void uint64ToLimbs(mp_limb_t * r, uint64_t const s,
                              size_t const limbsToWrite)
    {
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
    void setM(mp_limb_t * r) const { mpn_copyi(r, mpz_limbs_read(m), mpz_size(m)); }
    /** Add M to a residue, return carry */
    mp_limb_t addM(mp_limb_t * r, mp_limb_t const * s) const
    {
        return mpn_add_n(r, s, mpz_limbs_read(m), mpz_size(m));
    }
    /** Subtract M from a residue, return borrow */
    mp_limb_t subM(mp_limb_t * r, mp_limb_t const * s) const
    {
        return mpn_sub_n(r, s, mpz_limbs_read(m), mpz_size(m));
    }
    /** Returns a negative value if s < m, 0 if s == m, and a positive value if
     * s > m. */
    int cmpM(mp_limb_t const * s) const
    {
        return mpn_cmp(s, mpz_limbs_read(m), mpz_size(m));
    }
    uint64_t modM(uint64_t const s) const
    {
        if (mpz_sizeinbase(m, 2) > 64 || mpz_cmp_uint64(m, s) > 0)
            return s;
        /* m <= s, thus m fits in uint64_t */
        const uint64_t m64 = mpz_get_uint64(m);
        /* We always make sure that the modulus in non-zero in the ctor
         */
        ASSERT_FOR_STATIC_ANALYZER(m64 != 0);
        return s % m64;
    }
    void set_residue_u64(Residue & r, uint64_t const s) const
    {
        ASSERT(mpz_cmp_uint64(m, s) > 0);
        size_t const limbsToWrite = MIN(limbsPerUint64, mpz_size(m));
        uint64ToLimbs(r.r.get(), s, limbsToWrite);
        for (size_t i = limbsToWrite; i < mpz_size(m); i++)
            r.r[i] = 0;
    }
    void set_residue_mpz(Residue & r, mpz_srcptr s) const
    {
        ASSERT_ALWAYS(mpz_cmp(s, m) < 0);
        size_t written;
        mpz_export(r.r.get(), &written, -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, s);
        ASSERT_ALWAYS(written <= mpz_size(m));
        for (size_t i = written; i < mpz_size(m); i++)
            r.r[i] = 0;
    }
    void set_mpz_residue(mpz_ptr r, Residue const & s) const
    {
        mpz_import(r, mpz_size(m), -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS,
                   s.r.get());
    }

    /* Methods of the API */
  public:

    /* Methods for residues */

    /* {{{ set(*4), set_reduced(*2), set0, set1 */
    void set(Residue & r, Residue const & s) const
    {
        assertValid(s);
        mpn_copyi(r.r.get(), s.r.get(), mpz_size(m));
    }
    void set(Residue & r, uint64_t const s) const
    {
        uint64_t const sm = modM(s);
        set_residue_u64(r, sm);
    }
    void set(Residue & r, int64_t const s) const
    {
        uint64_t const u = modM(safe_abs64(s));
        set_residue_u64(r, u);
        if (s < 0)
            neg(r, r);
    }
    void set(Residue & r, Integer const & s) const
    {
        cxx_mpz t;
        mpz_mod(t, s, m);
        set_residue_mpz(r, t);
    }
    /* Sets the Residue to the class represented by the integer s. Assumes that
       s is reduced (mod m), i.e. 0 <= s < m */
    void set_reduced(Residue & r, uint64_t const s) const
    {
        assertValid(s);
        set_residue_u64(r, s);
    }
    void set_reduced(Residue & r, Integer const & s) const
    {
        assertValid(r);
        set_residue_mpz(r, s);
    }
    void set0(Residue & r) const { mpn_zero(r.r.get(), mpz_size(m)); }
    void set1(Residue & r) const
    {
        set0(r);
        if (mpz_cmp_ui(m, 1) > 0) {
            r.r[0] = 1;
        }
    }
    /* }}} */

    /* {{{ get equal is0 is1 */
    Integer get(Residue const & s) const
    {
        Integer r;
        assertValid(s);
        set_mpz_residue(r, s);
        return r;
    }
    bool equal(Residue const & a, Residue const & b) const
    {
        assertValid(a);
        assertValid(b);
        return mpn_cmp(a.r.get(), b.r.get(), mpz_size(m)) == 0;
    }
    bool is0(Residue const & a) const
    {
        assertValid(a);
        return mpn_zero_p(a.r.get(), mpz_size(m));
    }
    bool is1(Residue const & a) const
    {
        assertValid(a);
        if (mpz_cmp_ui(m, 1) == 0)
            return true;
        return a.r[0] == 1 &&
               (mpz_size(m) == 1 || mpn_zero_p(a.r.get() + 1, mpz_size(m) - 1));
    }
    /* }}} */

    /* {{{ neg add(*2) add1 sub(*2) sub1 div2 */
    void neg(Residue & r, Residue const & a) const
    {
        assertValid(a);
        if (is0(a) != 0) {
            if (&r != &a)
                set0(r);
        } else {
            const mp_limb_t bw = mpn_sub_n(r.r.get(), mpz_limbs_read(m), a.r.get(), mpz_size(m));
            ASSERT_ALWAYS(bw == 0);
        }
        assertValid(r);
    }
    void add(Residue & r, Residue const & a, Residue const & b) const
    {
        assertValid(a);
        assertValid(b);
        mp_limb_t const cy = mpn_add_n(r.r.get(), a.r.get(), b.r.get(), mpz_size(m));
        if (cy || cmpM(r.r.get()) >= 0) {
            mp_limb_t const bw = subM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void add1(Residue & r, Residue const & a) const
    {
        assertValid(a);
        mp_limb_t const cy = mpn_add_1(r.r.get(), a.r.get(), mpz_size(m), 1);
        if (cy || cmpM(r.r.get()) >= 0) {
            mp_limb_t const bw = subM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void add(Residue & r, Residue const & a, uint64_t const b) const
    {
        assertValid(a);
        uint64_t const bm = modM(b);
        mp_limb_t cy;

        if (limbsPerUint64 == 1 || bm < GMP_NUMB_MAX) {
            cy = mpn_add_1(r.r.get(), a.r.get(), mpz_size(m), bm);
        } else {
            mp_limb_t t[limbsPerUint64];
            mp_size_t const toWrite = MIN(limbsPerUint64, mpz_size(m));

            uint64ToLimbs(t, bm, limbsPerUint64);
            cy = mpn_add(r.r.get(), a.r.get(), mpz_size(m), t, toWrite);
        }
        if (cy || cmpM(r.r.get()) >= 0) {
            const mp_limb_t bw = subM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(bw == cy);
        }
        assertValid(r);
    }
    void sub(Residue & r, Residue const & a, Residue const & b) const
    {
        assertValid(a);
        assertValid(b);
        mp_limb_t const bw = mpn_sub_n(r.r.get(), a.r.get(), b.r.get(), mpz_size(m));
        if (bw) {
            mp_limb_t const cy = addM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    void sub1(Residue & r, Residue const & a) const
    {
        assertValid(a);
        mp_limb_t const bw = mpn_sub_1(r.r.get(), a.r.get(), mpz_size(m), 1);
        if (bw) {
            mp_limb_t const cy = addM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    void sub(Residue & r, Residue const & a, uint64_t const b) const
    {
        assertValid(a);
        uint64_t const bm = modM(b);
        mp_limb_t bw;

        if (limbsPerUint64 == 1 || bm < GMP_NUMB_MAX) {
            bw = mpn_sub_1(r.r.get(), a.r.get(), mpz_size(m), bm);
        } else {
            mp_limb_t t[limbsPerUint64];
            mp_size_t const toWrite = MIN(limbsPerUint64, mpz_size(m));

            uint64ToLimbs(t, bm, limbsPerUint64);
            bw = mpn_sub(r.r.get(), a.r.get(), mpz_size(m), t, toWrite);
        }
        if (bw) {
            const mp_limb_t cy = addM(r.r.get(), r.r.get());
            ASSERT_ALWAYS(cy == bw);
        }
        assertValid(r);
    }
    bool div2(Residue & r, Residue const & a) const
    {
        assertValid(a);
        if (mpz_even_p(m)) {
            return false;
        } else {
            if (a.r[0] % 2 == 0) {
                const mp_limb_t lsb = mpn_rshift(r.r.get(), a.r.get(), mpz_size(m), 1);
                ASSERT_ALWAYS(lsb == 0);
            } else {
                const mp_limb_t cy = addM(r.r.get(), a.r.get());
                const mp_limb_t lsb = mpn_rshift(r.r.get(), r.r.get(), mpz_size(m), 1);
                ASSERT_ALWAYS(lsb == 0);
                r.r[mpz_size(m) - 1] |= cy << (GMP_NUMB_BITS - 1);
            }
            assertValid(r);
            return true;
        }
    }
    /* }}} */

    /* {{{ mul sqr */
    void mul(Residue & r, Residue const & a, Residue const & b) const
    {
        mp_size_t const nrWords = mpz_size(m);
        mp_limb_t Q[2];
        mp_limb_t * t;
        if (r.r == a.r || r.r == b.r) {
            t = new mp_limb_t[nrWords + 1];
        } else {
            t = r.r.get();
        }
        t[nrWords] = mpn_mul_1(t, a.r.get(), nrWords, b.r[nrWords - 1]);
        if (t[nrWords] != 0)
            mpn_tdiv_qr(Q, t, 0, t, nrWords + 1, mpz_limbs_read(m), nrWords);
        /* t <= (m-1) * beta */
        for (mp_size_t iWord = nrWords - 1; iWord > 0; iWord--) {
            mpn_copyd(t + 1, t, nrWords);
            t[0] = 0;
            const mp_limb_t msw = mpn_addmul_1(t, a.r.get(), nrWords, b.r[iWord - 1]);
            t[nrWords] += msw;
            const mp_limb_t cy = t[nrWords] < msw;
            if (cy) {
                const mp_limb_t bw = subM(t + 1, t + 1);
                ASSERT_ALWAYS(bw == cy);
            }
            if (t[nrWords] != 0)
                mpn_tdiv_qr(Q, t, 0, t, nrWords + 1, mpz_limbs_read(m), nrWords);
        }
        if (cmpM(t) >= 0) {
            mpn_tdiv_qr(Q, t, 0, t, nrWords, mpz_limbs_read(m), nrWords);
        }
        if (r.r == a.r || r.r == b.r) {
            mpn_copyi(r.r.get(), t, nrWords);
            delete[] t;
        }
    }
    void sqr(Residue & r, Residue const & a) const
    {
        mul(r, a, a);
    }
    /* }}} */

  protected:
    bool divn(Residue &, Residue const &, unsigned long, mp_limb_t const *,
              mp_limb_t) const;
};

inline arithxx_mod_mpz_new::Residue::Residue(Modulus const & m)
    : r(new mp_limb_t[mpz_size(m.m) + 1])
{
    mpn_zero(r.get(), mpz_size(m.m) + 1);
}
inline arithxx_mod_mpz_new::Residue::Residue(Modulus const & m, Residue const & s)
    : r(new mp_limb_t[mpz_size(m.m) + 1])
{
    mpn_copyi(r.get(), s.r.get(), mpz_size(m.m));
}

#endif /* CADO_UTILS_ARITHXX_MOD_MPZ_NEW_HPP */
