#include "cado.h" // IWYU pragma: keep
#include <gmp.h>
#include "u64arith.h"      // for u64arith_mul_1_1_2
#include "ularith.h"    // IWYU pragma: keep
#include "mod_mpz_new.hpp"

#define MOD_NO_SHARED_MOD_POW_UL 1
#define MOD_NO_SHARED_MOD_POW_MP 1
#define MOD_NO_SHARED_MOD_POW_INT 1

typedef ModulusMPZ Modulus;
#include "mod_common.cpp"
#include "macros.h"

void ModulusMPZ::pow (Residue &r, const Residue &b, const uint64_t e) const
{
    cxx_mpz R, B;
    set_mpz_residue(B, b);
    if (ULONG_BITS == 64) {
        mpz_powm_ui(R, B, (unsigned long) e, m);
    } else {
       cxx_mpz E(e);
       mpz_powm(R, B, E, m);
    }
    set_residue_mpz(r, R);
}

void ModulusMPZ::pow (Residue &r, const Residue &b, const Integer &e) const
{
    cxx_mpz R, B;
    set_mpz_residue(B, b);
    mpz_powm(R, B, e, m);
    set_residue_mpz(r, R);

}

void ModulusMPZ::pow (Residue &r, const Residue &b, const uint64_t *e, const size_t nrWords) const
{
    cxx_mpz E;
    E.set(e, nrWords);
    pow(r, b, E);
}

void ModulusMPZ::pow2 (Residue &r, const uint64_t e) const
{
    cxx_mpz R;
    cxx_mpz B(2);
    if (ULONG_BITS == 64) {
        mpz_powm_ui(R, B, e, m);
    } else {
       cxx_mpz E(e);
       mpz_powm(R, B, E, m);
    }
    set_residue_mpz(r, R);
}

void ModulusMPZ::pow2 (Residue &r, const Integer &e) const
{
    cxx_mpz R, B(2);
    mpz_powm(R, B, e, m);
    set_residue_mpz(r, R);

}

void ModulusMPZ::pow2 (Residue &r, const uint64_t *e, const size_t nrWords) const
{
    cxx_mpz E;
    E.set(e, nrWords);
    pow2(r, E);
}

bool
ModulusMPZ::div3(ModulusMPZ::Residue &r, ModulusMPZ::Residue const &a) const {
    Residue t(*this);
    const unsigned int mMod3 = mpz_tdiv_ui(m, 3);
    const unsigned int aMod3 = mpn_mod_1(a.r, mpz_size(m), 3);

    if (mMod3 == 0) {
        return false;
    }

    mp_limb_t cy = 0;
    if (aMod3 == 0) {
        set(t, a);
    } else if (aMod3 == mMod3) {
        cy = addM(t.r, a.r);
        cy += addM(t.r, t.r);
    } else {
        cy = addM(t.r, a.r);
    }
    mp_limb_t spill = mpn_divexact_by3(r.r, t.r, mpz_size(m));
    ASSERT_ALWAYS(cy != 0 || spill == 0);
    
    return true;
}

bool
ModulusMPZ::divn (Residue &r, const Residue &a, const unsigned long b,
                  const mp_limb_t *minvb, const mp_limb_t binvw MAYBE_UNUSED) const {
    const size_t sz = mpz_size(m);
    const unsigned long mModB = mpz_tdiv_ui(m, b);

    if (mModB == 0) {
        return false;
    }

    const mp_limb_t aModB = mpn_mod_1(a.r, sz, b);
    const mp_limb_t k = (minvb[mModB] * aModB) % b;

    if (r.r != a.r) {
        mpn_copyi(r.r, a.r, sz);
    }
    mp_limb_t cy MAYBE_UNUSED;
    cy = mpn_addmul_1 (r.r, m->_mp_d, sz, k);

#if defined(HAVE_MPN_DIVEXACT_1)
    mpn_divexact_1 (r.r, r.r, sz, b);
#else
    const mp_limb_t binvlimb = binvw & GMP_NUMB_MASK;
    static_assert(GMP_LIMB_BITS <= 64, "Currently can't handle GMP_LIMB_BITS > 64");

    mp_limb_t borrow = 0;
    mp_limb_t t0 = (r.r[0] * binvlimb) & GMP_NUMB_MASK;
    r.r[0] = t0;
    
    for (size_t i = 1; i < sz; i++) {
        mp_limb_t t1, t2;
#if GMP_LIMB_BITS == 64
        // while sizeof(mp_limb_t) == sizeof(uint64_t), it does not mean
        // that mp_limb_t and uint64_t are typedef'ed to the same type.
        // One might be unsigned long, and the other unsigned long long,
        // in which case we need pointer casts to avoid a compiler
        // warning.
        u64arith_mul_1_1_2((uint64_t *) &t1, (uint64_t *) &t2, t0, b);
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

    if (1) {
        mp_limb_t t1, t2;
#if GMP_LIMB_BITS == 64
        // see remark above.
        u64arith_mul_1_1_2((uint64_t *) &t1, (uint64_t *) &t2, t0, b);
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

/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */

bool
ModulusMPZ::div5 (Residue &r, const Residue &a) const
{
  /* inv_5[i] = -1/i (mod 5) */
  const mp_limb_t inv_5[5] = {0,4,2,3,1};
  const mp_limb_t c = (mp_limb_t) UINT64_C(0xcccccccccccccccd); /* 1/5 (mod 2^64) */
  
  return divn (r, a, 5, inv_5, c);
}


/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */

bool
ModulusMPZ::div7 (Residue &r, const Residue &a) const
{
  /* inv_7[i] = -1/i (mod 7) */
  const mp_limb_t inv_7[7] = {0,6,3,2,5,4,1};
  const mp_limb_t c = (mp_limb_t) UINT64_C(0x6db6db6db6db6db7); /* 1/7 (mod 2^64) */
  return divn (r, a, 7, inv_7, c);
}


/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */

bool
ModulusMPZ::div11 (Residue &r, const Residue &a) const
{
  /* inv_11[i] = -1/i (mod 11) */
  const mp_limb_t inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1}; 
  const mp_limb_t c = (mp_limb_t) UINT64_C(0x2e8ba2e8ba2e8ba3); /* 1/11 (mod 2^64) */
  return divn (r, a, 11, inv_11, c);
}


/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */

bool
ModulusMPZ::div13 (Residue &r, const Residue &a) const
{
  /* inv_13[i] = -1/i (mod 13) */
  const mp_limb_t inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1}; 
  const mp_limb_t c = (mp_limb_t) UINT64_C(0x4ec4ec4ec4ec4ec5); /* 1/13 (mod 2^64) */
  return divn (r, a, 13, inv_13, c);
}

void ModulusMPZ::gcd (Integer &r, const Residue &a) const {
    cxx_mpz A;
    
    set_mpz_residue(A, a);
    mpz_gcd(r, A, m);
}

bool ModulusMPZ::sprp(const Residue &a MAYBE_UNUSED) const
{
    return isprime();
}

bool ModulusMPZ::sprp2() const
{
    return isprime();
}

bool ModulusMPZ::isprime() const
{
    return mpz_probab_prime_p(m, 25);
}

bool ModulusMPZ::inv (Residue &r, const Residue &a) const {
    cxx_mpz A;
    set_mpz_residue(A, a);
    bool exists = mpz_invert(A, A, m);
    if (exists)
        set_residue_mpz(r, A);
    return exists;
}

int ModulusMPZ::jacobi(const Residue &a MAYBE_UNUSED) const {
    cxx_mpz A;
    set_mpz_residue(A, a);
    return mpz_jacobi(A, m);
}
