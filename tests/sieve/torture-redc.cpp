#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include "macros.h"
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cinttypes>          // for PRIu32, PRId64, PRIu64
#include <sys/time.h>
#include <gmp.h>               // for gmp_randstate_t, gmp_urandomb_ui, mpz_...
#include "gmp_aux.h"           // for mpz_get_uint64, mpz_fits_sint64_p, mpz...
#include "cxx_mpz.hpp"
#include "las-arith.hpp"

/* smaller p bits are good to spot some corner cases that happen only with
 * probability 1/p
 * */
int minimum_p_bits = 10;

/*
 *
 * g++ -O3 -DNDEBUG -W -Wall -I build/localhost -I utils -I . torture-redc.cpp  -lgmp && ./a.out
 *
 */

// [[THEORY]]  Signed redc_32 based on 64-bit arithmetic
// [[THEORY]]  Assume:
// [[THEORY]]    * p is an odd number < 2^32.
// [[THEORY]]    * invp is -1/p mod 2^32.
// [[THEORY]]    * x is some signed integer in ]-2^32*p, 2^32*p[
// [[THEORY]]  Compute:
// [[THEORY]]    * x/2^32 mod p as an integer in [0, p[

int redc_32_preconditions(const int64_t x, const uint32_t p, const uint32_t invp)
{
    cxx_mpz xx(x);
    cxx_mpz pp(p);
    cxx_mpz ii(invp);
    cxx_mpz z = ii*pp+1;

    cxx_mpz B = 1;
    B <<= 32;

    if (!(p & 1)) return 0;
    // if (!(p < B)) return 0;
    // if (!(ii < B)) return 0;
    if (!((z % B) == 0)) return 0;
    if (x >= 0 && !(xx-B*pp < 0)) return 0;
    if (x <= 0 && !(xx+B*pp > 0)) return 0;

    return 1;
}

template<typename T>
int redc_32_reference(const T x, const uint32_t p)
{
    cxx_mpz xx(x);
    cxx_mpz pp(p);
    cxx_mpz B = 1;
    B <<= 32;
    cxx_mpz iB;
    mpz_invert(iB, B, pp);
    xx = xx % pp;
    xx = xx * iB;
    /* we want a nonnegative representative */
    mpz_fdiv_r(xx, xx, pp);
    return mpz_get_uint64(xx);
}


/* x is at most B*p */
int redc_u32_preconditions(const uint64_t x, const uint32_t p, const uint32_t invp)
{
    cxx_mpz xx(x);
    cxx_mpz pp(p);
    cxx_mpz ii(invp);
    cxx_mpz z = ii*pp+1;

    cxx_mpz B = 1;
    B <<= 32;

    if (!(p & 1)) return 0;
    // if (!(p < B)) return 0;
    // if (!(ii < B)) return 0;
    if (!((z % B) == 0)) return 0;
    if (!(xx-B*pp < 0)) return 0;

    return 1;
}

int redc_32_postconditions(uint32_t u, const int64_t x, const uint32_t p, const uint32_t invp MAYBE_UNUSED)
{
    cxx_mpz uu(u);
    cxx_mpz xx(x);
    cxx_mpz pp(p);

    cxx_mpz B = 1;
    B <<= 32;

    if (!(u < p)) return 0;
    if (!((uu*B-xx)%pp == 0)) return 0;

    return 1;
}

int redc_u32_postconditions(uint32_t u, const uint64_t x, const uint32_t p, const uint32_t invp MAYBE_UNUSED)
{
    cxx_mpz uu(u);
    cxx_mpz xx(x);
    cxx_mpz pp(p);

    cxx_mpz B = 1;
    B <<= 32;

    if (!(u < p)) return 0;
    if (!((uu*B-xx)%pp == 0)) return 0;

    return 1;
}

/* This version seems to be ok with 31-bit p */
uint32_t
ok31_redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t)x * invp;
  int32_t u = (x + (uint64_t)t * (uint64_t)p) >> 32;
  // if x > 0, u might be too large by p,
  // if x < 0, u might be too small by p.
  t = u;
  if (x < 0 && (int32_t) t < 0) u = t + p;
  t -= p;
  if (x > 0 && (int32_t) t >= 0) u = t;
  if (UNLIKELY((uint32_t) u >= p))
    return 42;
  return u;
}

inline uint32_t
oldbuggy_redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
  uint32_t t = (uint32_t) x * invp;                            /* t = x * invp mod 2^32 */
  uint32_t u = (x + (uint64_t)t * (uint64_t)p) >> 32;
  /* x + t*p is bounded by 2^32*p-1+(2^32-1)*p < 2*2^32*p:
     we might want p < 2^31 so that there is no overflow */
  t = u - p;
  if ((int32_t) t >= 0) u = t;
  return u;
}

#if defined(GENUINE_GNUC)
#pragma GCC optimize ("no-tree-loop-vectorize")
#endif
template <bool CARRY>
int test_redc_32(gmp_randstate_t rstate, size_t N, bool check, bool signed_x = true)
{
    constexpr unsigned int loops = 1024;
    constexpr int maximum_p_bits = CARRY ? 32 : 31;
    std::vector<uint32_t> ps;
    std::vector<int64_t> xs;
    std::vector<uint32_t> ips;
    std::vector<uint32_t> us;
    std::vector<uint32_t> rs;
    ps.reserve(N);
    xs.reserve(N);
    ips.reserve(N);
    us.reserve(N);

    for(size_t i = 0 ; i < N ; i++) {
        /* number of bits in [minimum_p_bits..32] */
        size_t pbits = minimum_p_bits + gmp_urandomb_ui(rstate, minimum_p_bits) % (maximum_p_bits+1-minimum_p_bits);
        cxx_mpz p;
        mpz_rrandomb(p, rstate, pbits);
        uint64_t pi = mpz_get_uint64(p);
        if (!pi)
            pi++;
        else if (!(pi & 1))
            pi--;
        p = pi;
        ps.push_back(pi);

        cxx_mpz B = 1;
        B <<= 32;

        cxx_mpz ip;
        mpz_neg(ip, p);
        mpz_invert(ip, ip, B);
        ips.push_back(mpz_get_uint64(ip));

        cxx_mpz xmax = B * p;
        if (signed_x && (xmax >> 63 != 0)) {
            /* XXX XXX XXX This is quite ugly. The constraints written in
             * redc_32 for the x interval are obviously quite loose when
             * 2^32*p exceeds 2^63. In practice, we expect that we'll
             * never reach this size if p is a large factor base prime.
             */
            xmax = 1;
            xmax <<= 63;
        }
        cxx_mpz x;
        mpz_urandomm(x, rstate, xmax);
        if (signed_x) {
            if (gmp_urandomb_ui(rstate, 1) & 1)
                mpz_neg(x, x);
            ASSERT_ALWAYS(mpz_fits_sint64_p(x));
            xs.push_back(mpz_get_int64(x));
            ASSERT_ALWAYS(redc_32_preconditions(xs[i], ps[i], ips[i]));
        } else {
            ASSERT_ALWAYS(mpz_fits_uint64_p(x));
            xs.push_back(mpz_get_uint64(x));
            ASSERT_ALWAYS(redc_u32_preconditions(xs[i], ps[i], ips[i]));
        }
    }

    clock_t clk0 = clock();

    if (check) {
        if (signed_x) {
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
            for(size_t i = 0 ; i < N ; i++)
                us.push_back(redc_32<CARRY>(xs[i], ps[i], ips[i]));

            for(size_t i = 0 ; i < N ; i++)
                rs.push_back(redc_32_reference(xs[i], ps[i]));
        } else {
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
            for(size_t i = 0 ; i < N ; i++)
                us.push_back(redc_u32<CARRY>(xs[i], ps[i], ips[i]));

            for(size_t i = 0 ; i < N ; i++)
                rs.push_back(redc_32_reference((uint64_t) xs[i], ps[i]));
        }
    } else {
        uint32_t fake_sum = 0;
        for (unsigned int loop = 0; loop < loops; loop++) {
            if (signed_x) {
#if defined(ALIGN_LOOP_32)
                __asm__ volatile (".p2align 5");
#endif
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
                for(size_t i = 0 ; i < N ; i++)
                    fake_sum += redc_32<CARRY>(xs[i], ps[i], ips[i]);
            } else {
#if defined(ALIGN_LOOP_32)
                __asm__ volatile (".p2align 5");
#endif
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
                for(size_t i = 0 ; i < N ; i++)
                    fake_sum += redc_u32<CARRY>(xs[i], ps[i], ips[i]);
            }
        }
        volatile uint32_t fake_sum_vol = fake_sum;
        if (fake_sum_vol) {}
    }

    clock_t clk1 = clock();

    if (check) {
        if (signed_x) {
            for(size_t i = 0 ; i < N ; i++) {
                if (us[i] != rs[i] || !redc_32_postconditions(us[i], xs[i], ps[i], ips[i])) {
                    fprintf(stderr, "ERROR: redc_32<%s>("
                            "%" PRId64 ", "
                            "%" PRIu32 ", "
                            "%" PRIu32 ") "
                            "returns "
                            "%" PRIu32
                            " instead of expected %" PRIu32 "\n",
                            CARRY ? "true" : "false",
                            xs[i], ps[i], ips[i], us[i],
                            rs[i]
                           );
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            for(size_t i = 0 ; i < N ; i++) {
                if (us[i] != rs[i] || !redc_u32_postconditions(us[i], xs[i], ps[i], ips[i])) {
                    fprintf(stderr, "ERROR: redc_u32<%s>("
                            "%" PRIu64 ", "
                            "%" PRIu32 ", "
                            "%" PRIu32 ") "
                            "returns "
                            "%" PRIu32
                            " instead of expected %" PRIu32 "\n",
                            CARRY ? "true" : "false",
                            (uint64_t) xs[i], ps[i], ips[i], us[i],
                            rs[i]
                           );
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    const char * fname[2] = { "redc_u32", "redc_32" };
    printf("%s<%s>: %zu tests in %.4fs\n",
            fname[signed_x], CARRY ? "true" : "false", (check) ? N : N*(size_t)loops, ((double)(clk1-clk0))/CLOCKS_PER_SEC);

    return 0;
}

template <bool CARRY>
int test_redc_u32(gmp_randstate_t rstate, size_t N, bool check)
{
    return test_redc_32<CARRY>(rstate, N, check, false);
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    size_t Nmax = 1e5;
    bool check = true;
    for( ; argc > 1 ; argv++,argc--) {
        if (strcmp(argv[1], "--minimum-p-bits") == 0) {
            argv++,argc--;
            minimum_p_bits = atoi(argv[1]);
        } else if (strcmp(argv[1], "-t") == 0) {
            check = false;
        } else {
            // coverity[tainted_argument]
            Nmax = atol(argv[1]);
        }
    }
    if (minimum_p_bits > 32) {
        fprintf(stderr, "--minimum_p_bits accepts at most 32\n");
        exit(EXIT_FAILURE);
    }

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    for(size_t N = 1 ; N < Nmax ; N *= 2) {
        test_redc_32<false>(rstate, N, check);
        test_redc_u32<false>(rstate, N, check);
        test_redc_32<true>(rstate, N, check);
        test_redc_u32<true>(rstate, N, check);
    }
    gmp_randclear(rstate);
}

