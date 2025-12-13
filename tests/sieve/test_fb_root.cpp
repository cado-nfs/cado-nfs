/* test fb_root_in_qlattice_31bits */
#include "cado.h" // IWYU pragma: keep

#include <cstdint>                 // for uint32_t, uint64_t
#include <cstring>                 // for strcmp
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <array>
#include <iostream>

#include <gmp.h>
#include "fmt/ostream.h"
#include "fmt/format.h"
#include "fmt/base.h"

#include "cado_poly.h"              // MAX_DEGREE
#include "cxx_mpz.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "gmp_aux.h"                // for mpz_add_int64, mpz_init_set_int64
#include "las-arith.hpp"            // for invmod_redc_32
#include "las-fbroot-qlattice.hpp"
#include "las-qlattice.hpp"
#include "macros.h"
#include "tests_common.h"
#include "timing.h"  // seconds

/*
 * R encodes (p:1) if R<p, or (1:R-p) if R >= p
 *
 * when R < p (affine root),
 *      return num/den%p with num=R*b1-a1 den=a0-R*b0
 * when R >= p (projective root), 
 *      return num/den%p with num=b1-(R-p)*a1 den=(R-p)*a0-b0
 *
 * It is assumed that R*b1-a1 and a0-R*b0 cannot be both multiples of p,
 * since the matrix (a0,b0,a1,b1) has determinant coprime to p, and
 * gcd(1,R) == 1
 *
 * if den in either expression is divisile by p, then the
 * resulting root is projective. Therefore, we return p + den/num
 * */

template<typename T>
static std::ostream& operator<<(std::ostream& os, fb_root_p1_t<T> const & R)
{
    if (R.proj)
        os << "(1:" << R.r << ")";
    else
        os << "(" << R.r << ":1)";
    return os;
}

/* declare these in the fmt namespace to work around a bug in g++-6
 * https://stackoverflow.com/questions/25594644/warning-specialization-of-template-in-different-namespace
 */
namespace fmt {
template <typename T> struct formatter<fb_root_p1_t<T>>: ostream_formatter {};
}


static fb_root_p1_t<cxx_mpz>
ref_fb_root_in_qlattice (fbprime_t p, fb_root_p1 R, qlattice_basis const & basis)
{
  cxx_mpz num, den;

  if (R.is_affine()) {
      den = -basis.b0;
      mpz_mul_ui (den, den, R.r);
      mpz_add_int64 (den, den, basis.a0);

      num = basis.b1;
      mpz_mul_ui (num, num, R.r);
      mpz_sub_int64 (num, num, basis.a1);
  } else {
      den = basis.a0;
      mpz_mul_ui (den, den, R.r);
      mpz_sub_int64 (den, den, basis.b0);

      num = -basis.a1;
      mpz_mul_ui (num, num, R.r);
      mpz_add_int64 (num, num, basis.b1);
  }

  if (mpz_gcd_ui(nullptr, den, p) != 1) {
      /* root is projective */
      mpz_invert (num, num, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
      return { num, true };
  } else {
      mpz_invert (den, den, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
      return { num, false };
  }
}

using basis_citer_t = std::vector<qlattice_basis>::const_iterator;

/* Make a random prime in the interval [2^(FBPRIME_BITS-1), 2^FBPRIME_BITS].
 * t is a temp variable that will get clobbered. */
static fbprime_t
make_random_fb_prime(mpz_ptr t)
{
    do {
        tests_common_urandomb (t, FBPRIME_BITS);
        mpz_nextprime (t, t);
    } while (mpz_sizeinbase (t, 2) != FBPRIME_BITS);
    return mpz_get_ui (t);
}

/* Make a random root mod p.
 * t is a temp variable that will get clobbered. */
static fbroot_t
make_random_fb_root(const fbprime_t p, mpz_ptr t)
{
    mpz_set_ui (t, p);
    tests_common_urandomm (t, t);
    fbroot_t const r = mpz_get_ui (t);
    ASSERT_ALWAYS (r < p);
    return r;
}

/* Compute -1/p[i] mod 2^32.
 * t and u are temp variables that will get clobbered. */
static uint32_t
make_invp32 (const fbprime_t p, mpz_ptr t, mpz_ptr u)
{
    mpz_ui_pow_ui (u, 2, 32);
    mpz_set_ui (t, p);
    mpz_neg (t, t);
    mpz_invert (t, t, u);
    return mpz_get_ui (t);
}

/* Compute -1/p[i] mod 2^64.
 * t and u are temp variables that will get clobbered. */
static uint64_t
make_invp64 (const fbprime_t p, mpz_t t, mpz_t u)
{
    mpz_ui_pow_ui (u, 2, 64);
    mpz_set_ui (t, p);
    mpz_neg (t, t);
    mpz_invert (t, t, u);
    return mpz_get_uint64 (t);
}

/* Make a random qlattice_basis where each coordinate is in [-2^31, 2^31-1]
 * if bits==31, or in [-2^63, 2^63-1] if bits == 63. */
static void
make_random_basis (qlattice_basis &basis, mpz_t t, const int bits)
{
    ASSERT_ALWAYS(bits == 31 || bits == 63);
    const int64_t offset = (bits == 31) ? INT32_MIN : INT64_MIN;
    tests_common_urandomb (t, bits + 1);
    basis.a0 = mpz_get_ui (t) + offset;
    tests_common_urandomb (t, bits + 1);
    basis.a1 = mpz_get_ui (t) + offset;
    tests_common_urandomb (t, bits + 1);
    basis.b0 = mpz_get_ui (t) + offset;
    tests_common_urandomb (t, bits + 1);
    basis.b1 = mpz_get_ui (t) + offset;
    if (bits == 31) {
        ASSERT_ALWAYS(basis.fits_31bits());
    }
}

static void
make_extremal_basis (qlattice_basis &basis, const int bits, const unsigned int count)
{
    const int64_t min = (bits == 31) ? INT32_MIN : INT64_MIN / 4;
    const int64_t max = (bits == 31) ? INT32_MAX : INT64_MAX / 4;

    basis.a0 = (count & 1) != 0 ? min : max;
    basis.a1 = (count & 2) != 0 ? min : max;
    basis.b0 = (count & 4) != 0 ? min : max; /* with bits==31, min-2 happens to pass, min-3 fails */
    basis.b1 = (count & 8) != 0 ? min : max; /* with bits==31, min-2 happens to pass, min-3 fails */
}

static void
print_error_and_exit(const fbprime_t p, fb_root_p1 const & Rab, fb_root_p1 const & rt, fb_root_p1_t<cxx_mpz> const & rref, const qlattice_basis &basis, const int bits)
{
    std::cerr
        << fmt::format("Error for p={}; R={}; {};\n", p, Rab, basis)
        << fmt::format("fb_root_in_qlattice_{}bits gives {}\n", bits, rt)
        << fmt::format("ref_fb_root_in_qlattice gives {}\n", rref);
    abort();
}

template<int Nr_roots>
static 
fb_entry_x_roots<Nr_roots> make_random_fb_entry_x_roots(
        mpz_ptr t, mpz_ptr u)
{
    fbprime_t const p = make_random_fb_prime(t);
    redc_invp_t const invq = make_invp32(p, t, u);
    fbroot_t rr[Nr_roots];
    for (int i_root = 0; i_root < Nr_roots; i_root++) {
        rr[i_root] = make_random_fb_root(p, t);
    }
    return fb_entry_x_roots<Nr_roots>(p, invq, rr);
}

template<int Nr_roots>
static void
test_chain_fb_root_in_qlattice_batch(basis_citer_t basis_begin,
                                basis_citer_t basis_end,
                                const unsigned long N, mpz_t t, mpz_t u,
                                const bool do_speed, const bool do_test,
                                const int bits)
{
    /* Test N random factor base entries with Nr_roots roots each against
     * each of the bases in [basis_begin, basis_end] */
    /* The reason why this is a recursively called function template is
     * that I'd like to call fb_entry_x_roots<Nr_roots>::transform() in
     * for the test eventually, but currently it is decided at compile
     * time (via SUPPORT_LARGE_Q macro) whether
     * fb_root_in_qlattice_31bits_batch() or
     * fb_root_in_qlattice_127bits_batch() gets called from there.
     * It should not be too hard to make that decision at run-time,
     * then this test code can be adapted accordingly.
    */

    std::vector<fb_entry_x_roots<Nr_roots>> fbv;
    fbv.reserve(N);

    /* Make N random factor base entries */
    for (unsigned long i = 0; i < N; i++)
        fbv.emplace_back(make_random_fb_entry_x_roots<Nr_roots>(t, u));

    if (do_speed) {
        double st = seconds ();
        fbroot_t fake_sum = 0; /* Fake sum to stop compiler from optimizing away
                                  everything due to unused results */
        for (auto basis_iter = basis_begin;
                basis_iter != basis_end; basis_iter++) {
            for (unsigned long i_fb = 0; i_fb < N; i_fb++) {
                typename fb_entry_x_roots<Nr_roots>::transformed_entry_t fbt;
                const fbprime_t p = fbv[i_fb].get_q();
                // fbv->transform_roots(fbt, *basis_iter);
                if (bits == 31) {
                    fb_root_in_qlattice_31bits_batch (fbt.roots.data(), p, fbv[i_fb].roots.data(), fbv[i_fb].invq, *basis_iter, Nr_roots);
                } else {
                    fb_root_in_qlattice_127bits_batch (fbt.roots.data(), p, fbv[i_fb].roots.data(), fbv[i_fb].invq, *basis_iter, Nr_roots);
                }
                for (unsigned long i_root = 0; i_root + 1 < Nr_roots + 1; i_root++)
                    fake_sum += fbt.get_r(i_root);
            }
        }
        st = seconds() - st;
        volatile fbroot_t const fake_sum_vol = fake_sum;
        if (fake_sum_vol) {}
        fmt::print ("fb_entry_x_roots<{}>::transform_roots with {}-bit basis: {} tests took {:.2f}s\n",
                Nr_roots, bits, N * N, st);
    }

    if (do_test) {
        /* Compare all the transformed roots with the reference implementation */
        for (auto basis_iter = basis_begin;
                basis_iter != basis_end; basis_iter++) {
            for (unsigned long i_fb = 0; i_fb < N; i_fb++) {
                typename fb_entry_x_roots<Nr_roots>::transformed_entry_t fbt;
                std::array<fb_root_p1_t<cxx_mpz>, Nr_roots> fbt_ref;
                const fbprime_t p = fbv[i_fb].get_q();

                /* Compute reference roots */
                bool had_any_proj = false;
                for (unsigned char i_root = 0; i_root + 1 < Nr_roots + 1; i_root++) {
                    const fbroot_t original_root = fbv[i_fb].get_r(i_root);
                    fbt_ref[i_root] = ref_fb_root_in_qlattice (p, original_root, *basis_iter);
                    had_any_proj |= fbt_ref[i_root].is_projective();
                }

                /* Compute batch transform of roots */
                fbt.p = p;
                const bool batch_worked = (bits == 31) ? fb_root_in_qlattice_31bits_batch (fbt.roots.data(), p, fbv[i_fb].roots.data(), fbv[i_fb].invq, *basis_iter, Nr_roots)
                                                       : fb_root_in_qlattice_127bits_batch (fbt.roots.data(), p, fbv[i_fb].roots.data(), fbv[i_fb].invq, *basis_iter, Nr_roots);

                /* If the batch transform did not work, verify that at least
                 * one transformed root is projective */
                if (!batch_worked) {
                    if (!had_any_proj) {
                        fmt::print(stderr, "Batch transform failed, but none of the "
                                "transformed roots are projective according to the "
                                "reference implementation\n");
                        abort();
                    }
                    /* If batch transform did not work and there actually was a
                     * projective transformed root, then all went as it should. */
                } else {
                    /* If batch transform worked, check that it agrees
                     * with reference roots
                    if (false) {
                        fmt::print(stderr, "p = {}, fbt.get_q() = {}\n",
                                       p, fbt.get_q());
                    }
                     * */
                    for (unsigned char i_root = 0; i_root + 1 < Nr_roots + 1; i_root++) {
                        if (fbt_ref[i_root].r != fbt.get_r(i_root) || fbt_ref[i_root].is_projective()) {
                            print_error_and_exit(p, fbv[i_fb].get_r(i_root), fbt.get_r(i_root),
                                    fbt_ref[i_root],
                                    *basis_iter, bits);
                        }
                    }
                }
            }
        }
    }

    /* Repeat the test with factor base entries that have one fewer root each */
    test_chain_fb_root_in_qlattice_batch<Nr_roots - 1>(basis_begin, basis_end, N, t, u, do_speed, do_test, bits);
}

/* Specialize the test for length 0 to terminate the recursion.
 *
 * (we used to specialize for length -1, but length 0 actually doesn't do
 * anything either)
 */
template<> void
test_chain_fb_root_in_qlattice_batch<0>(basis_citer_t basis_begin MAYBE_UNUSED,
                                         basis_citer_t basis_end MAYBE_UNUSED,
                                         const unsigned long N MAYBE_UNUSED,
                                         mpz_t t MAYBE_UNUSED, 
                                         mpz_t u MAYBE_UNUSED,
                                         const bool do_speed MAYBE_UNUSED,
                                         const bool do_test MAYBE_UNUSED,
                                         const int bits MAYBE_UNUSED)
{
}

static void test_one_root_31bits(const fbprime_t p, fb_root_p1 const & Rab, const uint32_t invp,
                     const qlattice_basis &basis)
{
    fb_root_p1 const r31 = fb_root_in_qlattice_31bits (p, Rab, invp, basis);
    auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
    if (rref != r31) {
        print_error_and_exit(p, Rab, r31, rref, basis, 31);
    }
}

static void
test_fb_root_in_qlattice_31bits (const bool test_timing,
    const bool test_correctness, const unsigned long N)
{
    std::vector<fbprime_t> p, R;
    std::vector<uint32_t> invp;
    fbprime_t r;
    cxx_mpz t, u;
    std::vector<qlattice_basis> basis;
    double st;

    p.assign(N, 0);
    R.assign(N, 0);
    invp.assign(N, 0);
    basis.assign(N, qlattice_basis());

    if (0 < N) {
        p[0] = 4294967291U;
        R[0] = p[0] - 1;
        invp[0] = make_invp32(p[0], t, u);
    }

    /* generate p[i], R[i], invp[i] for 0 <= i < N */
    for (unsigned long i = 0; i < N; i++) {
        if (i == 0) {
            p[0] = 4294967291U;
            R[0] = p[0] - 1;
        } else {
            p[i] = make_random_fb_prime(t);
            R[i] = make_random_fb_root(p[i], t);
        }
        /* invp[i] is -1/p[i] mod 2^32 */
        invp[i] = make_invp32(p[i], t, u);

        /* Fill the first up to 16 basis entries with extremal bases */
        if (i < 16)
            make_extremal_basis(basis[i], 31, i);
        else
            make_random_basis(basis[i], t, 31);

    }

    test_chain_fb_root_in_qlattice_batch<MAX_DEGREE>(basis.cbegin(), basis.cend(), N, t, u, test_timing, test_correctness, 31);

    /* Timing of fb_root_in_qlattice_31bits(), i.e., without batch inversion */
    if (test_timing) {
        /* efficiency test */
        st = seconds ();
        r = 0;
        for (unsigned long j = 0; j < N; j++)
            for (unsigned long i = 0; i < N; i++) {
                fb_root_p1 const Rab { R[i], false };
                auto Rij = fb_root_in_qlattice_31bits (p[i], Rab, invp[i], basis[j]);
                r += Rij.r;
            }
        st = seconds () - st;
        fmt::print ("fb_root_in_qlattice_31bits: {} tests took {:.2f}s (r={})\n",
                N * N, st, r);
    }

    /* Test of fb_root_in_qlattice_31bits(), i.e., without batch inversion */
    if (test_correctness) {
        for (unsigned long i = 0; i < N; i++)
            for (unsigned long j = 0; j < N; j++)
                test_one_root_31bits(p[i], R[i], invp[i], basis[j]);
    }
}

static void
test_fb_root_in_qlattice_127bits (const bool test_timing,
    const bool test_correctness, const unsigned long N)
{
    std::vector<fbprime_t> p, R;
    std::vector<uint32_t> invp32;
    std::vector<uint64_t> invp64;
    fbprime_t r;
    cxx_mpz t, u;
    std::vector<qlattice_basis> basis;
    double st;

    p.assign(N, 0);
    R.assign(N, 0);
    invp32.assign(N, 0);
    invp64.assign(N, 0);
    basis.assign(N, qlattice_basis());

    /* generate p[i], R[i], invp32[i], invp64[i] for 0 <= i < N */
    for (unsigned long i = 0; i < N; i++) {
        p[i] = make_random_fb_prime(t);
        R[i] = make_random_fb_root(p[i], t);
        invp32[i] = make_invp32(p[i], t, u);
        invp64[i] = make_invp64(p[i], t, u);

        /* Fill the first up to 16 basis entries with extremal bases */
        if (i < 16)
            make_extremal_basis(basis[i], 63, i);
        else
            make_random_basis(basis[i], t, 63);
    }

    test_chain_fb_root_in_qlattice_batch<MAX_DEGREE>(basis.cbegin(), basis.cend(), N, t, u, test_timing, test_correctness, 127);

    if (test_timing) {
        /* efficiency test */
        st = seconds ();
        r = 0;
        for (unsigned long j = 0; j < N; j++)
            for (unsigned long i = 0; i < N; i++) {
                fb_root_p1 const Rab = R[i];
                auto Rij = fb_root_in_qlattice_127bits (p[i], Rab, invp64[i], basis[j]);
                r += Rij.r;
            }
        st = seconds () - st;
        fmt::print ("fb_root_in_qlattice_127bits: {} tests took {:.2f}s (r={})\n",
                N * N, st, r);
    }

    if (test_correctness) {
        /* correctness test */
        for (unsigned long i = 0; i < N; i++) {
            for (unsigned long j = 0; j < N; j++) {
                fb_root_p1 const Rab = R[i];
                fb_root_p1 const r127 = fb_root_in_qlattice_127bits (p[i], R[i], invp64[i], basis[j]);
                auto rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
                if (rref != r127)
                    print_error_and_exit(p[i], Rab, r127, rref, basis[j], 127);
            }
        }
    }

}

/* exercise bugs in fb_root_in_qlattice_127bits */
static void
bug20200225 ()
{
    {
        /* exercises bug in assembly part of invmod_redc_32 (starting
         * around line 326 with 2nd #ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM) */
        unsigned long const a = 8088625;
        unsigned long const b = 2163105767;
        uint64_t const expected = 2062858318;
        uint64_t const got = invmod_redc_32 (a, b, -invmod_po2(b));
        if (expected != got) {
            fmt::print (stderr, "Error in invmod_redc_32 for a={} b={}\n", a, b);
            fmt::print (stderr, "Expected {}\n", expected);
            fmt::print (stderr, "Got      {}\n", got);
            exit (1);
        }
    }

    {
        unsigned long const a = 76285;
        unsigned long const b = 2353808591;
        uint64_t const expected = 2102979166;
        uint64_t const got = invmod_redc_32(a, b, -invmod_po2(b));
        if (expected != got) {
            fmt::print (stderr, "Error in invmod_redc_32 for a:={}; b:={};\n",
                    a, b);
            fmt::print (stderr, "Expected {}\n", expected);
            fmt::print (stderr, "Got      {}\n", got);
            exit (1);
        }
    }

    {
        fbprime_t const p = 3628762957;
        fb_root_p1 const Rab { 1702941053 };
        /* This is a redc64 value. Our code only works with redc32
         * values, though */
        uint64_t const invp = 5839589727713490555UL;
        qlattice_basis const basis { 
            -2503835703516628395L, 238650852,
                -3992552824749287692L, 766395543
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
        fb_root_p1 const r127 = fb_root_in_qlattice_127bits (
                /* the cast to redc_invp_t narrows invp to a redc32
                 * value, which is what the routine expects! */
                p, Rab, (redc_invp_t) invp, basis);
        if (rref != r127)
            print_error_and_exit(p, Rab, r127, rref, basis, 127);
    }

    {
        /* exercises bug in invmod_redc_32, already tested directly above */
        fbprime_t const p = 2163105767;
        fb_root_p1 const Rab = 1743312141;
        uint64_t const invp = 3235101737;
        qlattice_basis const basis {
            -30118114923155082L, -749622022,
                -2851499432479966615L, -443074848,
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);

        fb_root_p1 const r31 = fb_root_in_qlattice_31bits (p, Rab, invp, basis);
        if (rref != r31)
            print_error_and_exit(p, Rab, r31, rref, basis, 31);
    }

    {
        fbprime_t const p = 3725310689;
        fb_root_p1 const Rab = 2661839516;
        /* This is a redc64 value. Our code only works with redc32
         * values, though */
        uint64_t const invp = 1066179678986106591UL;
        qlattice_basis const basis {
            3008222006914909739L, 877054135,
            3170231873717741170L, 932375769,
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
        fb_root_p1 const r127 = fb_root_in_qlattice_127bits (
                /* the cast to redc_invp_t narrows invp to a redc32
                 * value, which is what the routine expects! */
                p, Rab, (redc_invp_t) invp, basis);
        if (rref != r127)
            print_error_and_exit(p, Rab, r127, rref, basis, 127);
    }

    {
        /* This is playing with the 32-bit limit */
        fbprime_t const p = 3486784401;       // 3^20
        fb_root_p1 const Rab = 2009510725;
        uint64_t const invp = 898235023;
        qlattice_basis const basis {
            -1353180941,
                -5,
                -223660881,
                8,
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
        fb_root_p1 const r127 = fb_root_in_qlattice_127bits (p, Rab, invp, basis);
        if (rref != r127)
            print_error_and_exit(p, Rab, r127, rref, basis, 127);
    }
    {
        fbprime_t const p = 3;
        fb_root_p1 const Rab = 2;
        uint64_t const invp = 1431655765;
        qlattice_basis const basis {
            -14730287151, 11, -6528529, -2,
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
        fb_root_p1 const r127 = fb_root_in_qlattice_127bits (p, Rab, invp, basis);
        if (rref != r127)
            print_error_and_exit(p, Rab, r127, rref, basis, 127);
    }
    {
        /* This is with
         * f1=-201985095924438464843258385780-2092848637636605284780511*x
         *    +200523270267920133662*x^2+383823019593599*x^3-6458886018*x^4
         *    +10080*x^5
         * and side-1 q=598766043859183; rho=150546891162130 (within a
         * descent), leading to
         *      a0=-1285560395691; b0=15859721; a1=-23568698; b1=-175;
         *      I=1024; J=32768
         */
        fbprime_t const p = 3;
        fb_root_p1 const Rab = 2;
        uint64_t const invp = 1431655765;
        qlattice_basis const basis {
            -1285560395691, 15859721, -23568698, -175,
        };
        auto rref = ref_fb_root_in_qlattice (p, Rab, basis);
        fb_root_p1 const r127 = fb_root_in_qlattice_127bits (p, Rab, invp, basis);
        if (rref != r127)
            print_error_and_exit(p, Rab, r127, rref, basis, 127);
    }

    {
        fbprime_t const p = 4273;
        constexpr size_t nroots = 3;
        std::array<fbroot_t, nroots> const Rab { 3698, 3454, 827 };
        uint64_t const invp = 3379283887;
        qlattice_basis const basis {
            -22223178967276, 82882403, -2164021285451, -27018774,
        };
        std::array<fbroot_t, nroots> Rt {};
        const bool t = fb_root_in_qlattice_127bits_batch(Rt.data(), p, Rab.data(), invp, basis, nroots);
        ASSERT_ALWAYS(t);
        for(size_t i = 0 ; i < nroots ; i++) {
            auto rref = ref_fb_root_in_qlattice (p, Rab[i], basis);
            if (rref.proj || rref.r != Rt[i])
                print_error_and_exit(p, Rab[i], Rt[i], rref, basis, 127);
        }
    }
}

/* The usual tests command line parameters "-seed" and "-iter" are accepted.
   Giving the "-check" parameter does only correctness tests but not timing.
   Giving only "-time" does only timing. Giving neither or both does both. */

// coverity[root_function]
int main(int argc, char const * argv[])
{
    unsigned long N = 100;
    int test_correctness = 1, test_timing = 1;

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_CHECK | PARSE_TIME);
    tests_common_get_iter(&N);
    tests_common_get_check_and_time(&test_correctness, &test_timing);

    bug20200225 ();
    test_fb_root_in_qlattice_31bits (test_timing != 0, test_correctness != 0, N);
    test_fb_root_in_qlattice_127bits (test_timing != 0, test_correctness != 0, N);

    return 0;
}
