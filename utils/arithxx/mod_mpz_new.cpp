#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>


#include <gmp.h>

#include "arith/ularith.h"    // IWYU pragma: keep
#include "arithxx/u64arith.h" // for u64arith_mul_1_1_2
#include "cxx_mpz.hpp"
#include "macros.h"
#include "mod_mpz_new.hpp"
#include "arithxx_common.hpp"

/* Only the .cpp source files that emit the non-inline symbols will
 * include this impl header file. So even though it does not look like
 * we're using it, in fact we are!  */
#include "arithxx_api_impl.hpp"      // IWYU pragma: keep

// scan-headers: stop here

/* {{{ pow */
/* r = b^e mod m. The base is read in place, and nothing is copied
 * except the result, whose allocation is all that is left.
 */
static void mod_mpz_new_powm(mp_limb_t * r, mpz_srcptr b, mpz_srcptr e,
        mpz_srcptr m)
{
    cxx_mpz R;
    mpz_powm(R, b, e, m);
    mp_size_t const n = mpz_size(m);
    mp_size_t const k = mpz_size(R);
    mpn_copyi(r, mpz_limbs_read(R), k);
    mpn_zero(r + k, n - k);
}

/* an mpz view of the exponent e, which is copied only if the limbs are
 * not 64-bit words */
static void mod_mpz_new_powm(mp_limb_t * r, mpz_srcptr b,
        uint64_t const * e, size_t const nrWords, mpz_srcptr m)
{
#if GMP_LIMB_BITS == 64
    mpz_t E;
    mpz_roinit_n(E, reinterpret_cast<mp_limb_t const *>(e), mp_size_t(nrWords));
    mod_mpz_new_powm(r, b, E, m);
#else
    mod_mpz_new_powm(r, b, cxx_mpz(e, nrWords), m);
#endif
}

template <>
void arithxx_details::api<arithxx_mod_mpz_new>::pow(
    Residue & r, Residue const & b, uint64_t const * e,
    size_t const nrWords) const
{
    auto const & me = downcast();
    mpz_t B;
    mpz_roinit_n(B, b.r.get(), mpz_size(me.m));
    mod_mpz_new_powm(r.r.get(), B, e, nrWords, me.m);
}

template <>
void arithxx_details::api<arithxx_mod_mpz_new>::pow2(
    Residue & r, uint64_t const * e,
    size_t const nrWords) const
{
    auto const & me = downcast();
    mp_limb_t const two = 2;
    mpz_t B;
    mpz_roinit_n(B, &two, 1);
    mod_mpz_new_powm(r.r.get(), B, e, nrWords, me.m);
}

/* The generic one-word version works with mul and sqr; mpz_powm is
 * faster. */
template <>
void arithxx_details::api<arithxx_mod_mpz_new>::pow(Residue & r, Residue const & b, uint64_t e) const
{
    pow(r, b, &e, 1);
}

template <>
void arithxx_details::api<arithxx_mod_mpz_new>::pow(Residue & r, Residue const & b, Integer const & e) const
{
    auto const & me = downcast();
    mpz_t B;
    mpz_roinit_n(B, b.r.get(), mpz_size(me.m));
    mod_mpz_new_powm(r.r.get(), B, e, me.m);
}

template <>
void arithxx_details::api<arithxx_mod_mpz_new>::pow2(Residue &r, const Integer &e) const
{
    auto const & me = downcast();
    mp_limb_t const two = 2;
    mpz_t B;
    mpz_roinit_n(B, &two, 1);
    mod_mpz_new_powm(r.r.get(), B, e, me.m);
}

/* }}} */

/* {{{ ModulusMPZ::divn */
bool arithxx_mod_mpz_new::Modulus::divn(Residue & r, Residue const & a, unsigned long const b,
                      mp_limb_t const * minvb,
                      mp_limb_t const binvw MAYBE_UNUSED) const
{
    mp_size_t const sz = mpz_size(m);
    unsigned long const mModB = mpz_tdiv_ui(m, b);

    if (mModB == 0) {
        return false;
    }

    mp_limb_t const aModB = mpn_mod_1(a.r.get(), sz, b);
    mp_limb_t const k = (minvb[mModB] * aModB) % b;

    if (r.r != a.r) {
        mpn_copyi(r.r.get(), a.r.get(), sz);
    }
    mp_limb_t cy MAYBE_UNUSED;
    cy = mpn_addmul_1(r.r.get(), m->_mp_d, sz, k);

#if defined(HAVE_MPN_DIVEXACT_1)
    mpn_divexact_1(r.r.get(), r.r.get(), sz, b);
#else
    const mp_limb_t binvlimb = binvw & GMP_NUMB_MASK;
    static_assert(GMP_LIMB_BITS <= 64,
                  "Currently can't handle GMP_LIMB_BITS > 64");

    mp_limb_t borrow = 0;
    mp_limb_t t0 = (r.r[0] * binvlimb) & GMP_NUMB_MASK;
    r.r[0] = t0;

    for (mp_size_t i = 1; i < sz; i++) {
        mp_limb_t t1, t2;
#if GMP_LIMB_BITS == 64
        // while sizeof(mp_limb_t) == sizeof(uint64_t), it does not mean
        // that mp_limb_t and uint64_t are typedef'ed to the same type.
        // One might be unsigned long, and the other unsigned long long,
        // in which case we need pointer casts to avoid a compiler
        // warning.
        u64arith_mul_1_1_2((uint64_t *)&t1, (uint64_t *)&t2, t0, b);
#elif GMP_LIMB_BITS == ULONG_BITS
        ularith_mul_ul_ul_2ul(&t1, &t2, t0, b);
#else
#error GMP limb is neither 64 bit nor same size as unsigned long
#endif
#if GMP_NAIL_BITS != 0
        t2 = (t2 << GMP_NAIL_BITS) | (t1 >> GMP_NUMB_BITS);
        t1 &= GMP_NUMB_MASK;
#endif
        ASSERT(t2 < GMP_NUMB_MASK); /* No overflow in next line */
        t2 += borrow;
        t0 = r.r[i] - t2;
        borrow = (t0 > r.r[i]);
        t0 = (t0 * binvlimb) & GMP_NUMB_MASK;
        r.r[i] = t0;
    }

    {
        mp_limb_t t1, t2;
#if GMP_LIMB_BITS == 64
        // see remark above.
        u64arith_mul_1_1_2((uint64_t *)&t1, (uint64_t *)&t2, t0, b);
#elif GMP_LIMB_BITS == ULONG_BITS
        ularith_mul_ul_ul_2ul(&t1, &t2, t0, b);
#else
#error GMP limb is neither 64 bit nor same size as unsigned long
#endif
        t2 += borrow;
        t0 = cy - t2;
        borrow = (t0 > t2);
        ASSERT_ALWAYS(borrow == 0);
        ASSERT_ALWAYS(t0 == 0);
    }

#endif

    return true;
}
/* }}} */

// {{{ div{3,5,7,11,13}
template <>
bool arithxx_details::api<arithxx_mod_mpz_new>::div3(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    Residue t(me);
    unsigned int const mMod3 = mpz_tdiv_ui(me.m, 3);
    unsigned int const aMod3 = mpn_mod_1(a.r.get(), mpz_size(me.m), 3);

    if (mMod3 == 0) {
        return false;
    }

    mp_limb_t cy = 0;
    if (aMod3 == 0) {
        me.set(t, a);
    } else if (aMod3 == mMod3) {
        cy = me.addM(t.r.get(), a.r.get());
        cy += me.addM(t.r.get(), t.r.get());
    } else {
        cy = me.addM(t.r.get(), a.r.get());
    }
    mp_limb_t const spill = mpn_divexact_by3(r.r.get(), t.r.get(), mpz_size(me.m));
    ASSERT_ALWAYS(cy != 0 || spill == 0);

    return true;
}

/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */

template <>
bool arithxx_details::api<arithxx_mod_mpz_new>::div5(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    /* inv_5[i] = -1/i (mod 5) */
    mp_limb_t const inv_5[5] = {0, 4, 2, 3, 1};
    auto const c = (mp_limb_t)UINT64_C(0xcccccccccccccccd); /* 1/5 (mod 2^64) */

    return me.divn(r, a, 5, inv_5, c);
}

/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */

template <>
bool arithxx_details::api<arithxx_mod_mpz_new>::div7(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    /* inv_7[i] = -1/i (mod 7) */
    mp_limb_t const inv_7[7] = {0, 6, 3, 2, 5, 4, 1};
    auto const c = (mp_limb_t)UINT64_C(0x6db6db6db6db6db7); /* 1/7 (mod 2^64) */
    return me.divn(r, a, 7, inv_7, c);
}

/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */

template <>
bool arithxx_details::api<arithxx_mod_mpz_new>::div11(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    /* inv_11[i] = -1/i (mod 11) */
    mp_limb_t const inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1};
    auto const c =
        (mp_limb_t)UINT64_C(0x2e8ba2e8ba2e8ba3); /* 1/11 (mod 2^64) */
    return me.divn(r, a, 11, inv_11, c);
}

/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */

template <>
bool arithxx_details::api<arithxx_mod_mpz_new>::div13(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    /* inv_13[i] = -1/i (mod 13) */
    mp_limb_t const inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1};
    auto const c =
        (mp_limb_t)UINT64_C(0x4ec4ec4ec4ec4ec5); /* 1/13 (mod 2^64) */
    return me.divn(r, a, 13, inv_13, c);
}

// }}}

template<>
void
arithxx_details::api<arithxx_mod_mpz_new>::gcd
(Integer & r, const Residue & a) const
{
    auto const & me = downcast();
    /* read the limbs of a in place, no copy */
    mpz_t A;
    mpz_roinit_n(A, a.r.get(), mpz_size(me.m));
    mpz_gcd(r, A, me.m);
}

template<>
bool arithxx_details::api<arithxx_mod_mpz_new>::inv(Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    /* read the limbs of a in place */
    mpz_t A;
    mpz_roinit_n(A, a.r.get(), mpz_size(me.m));
    cxx_mpz R;
    bool const exists = mpz_invert(R, A, me.m);
    if (exists) {
        mp_size_t const n = mpz_size(me.m);
        mp_size_t const k = mpz_size(R);
        mpn_copyi(r.r.get(), mpz_limbs_read(R), k);
        mpn_zero(r.r.get() + k, n - k);
    }
    return exists;
}

template<>
int arithxx_details::api<arithxx_mod_mpz_new>::jacobi(Residue const & a MAYBE_UNUSED) const
{
    auto const & me = downcast();
    cxx_mpz A;
    me.set_mpz_residue(A, a);
    return mpz_jacobi(A, me.m);
}

/* It's a bit awkward. We need to explicitly delete all the member
 * templates that make no sense *and* that we explicitly override anyway
 * in arithxx_mod_mpz_new::Modulus
 */
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set_reduced(Residue & r, Integer const & s) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set1(Residue & r) const = delete;
template<>
arithxx_mod_mpz_new::Integer arithxx_details::api<arithxx_mod_mpz_new>::get(Residue const & r) const = delete;
template<>
arithxx_mod_mpz_new::Integer arithxx_details::api<arithxx_mod_mpz_new>::get_raw(Residue const & r) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set(Residue & r, Residue const & s) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::neg(Residue & r, Residue const & s) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set(Residue & r, uint64_t) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set(Residue & r, Integer const &) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set_raw(Residue & r, uint64_t) const = delete;
template<>
void arithxx_details::api<arithxx_mod_mpz_new>::set_raw(Residue & r, Integer const &) const = delete;
template<>
bool arithxx_details::api<arithxx_mod_mpz_new>::is1(Residue const & r) const = delete;

template struct arithxx_details::api<arithxx_mod_mpz_new>;
