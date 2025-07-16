#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <exception>
#include <map>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <sys/types.h>
#include <unistd.h>
#if defined(HAVE_AVX512F) || defined(HAVE_AVX2) || defined(HAVE_AVX) ||        \
    defined(HAVE_SSE41)
#include <x86intrin.h>
#endif

#include "fmt/base.h"
#include "fmt/format.h"
#include <gmp.h>

#include "getprime.h"
#include "las-plattice.hpp"
#include "macros.h"

// see bug #30052
#if GNUC_VERSION_ATLEAST(12, 0, 0) && !GNUC_VERSION_ATLEAST(13, 0, 0)
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wuninitialized"
#pragma GCC diagnostic ignored "-Winit-self"
#endif

// scan-headers: stop here

/* see plattice.sage */

/* see https://gitlab.inria.fr/cado-nfs/cado-nfs/-/merge_requests/43 for
 * the different versions as well as the code performance.
 */

/* glossary of the different versions:
 *
 * simplistic --> two_legs (partial unrolling)
 * simplistic --> swapping loop (maintain a flip counter, and use a single
 * code path)
 * swapping_loop2 is functionally similar to swapping_loop, but is written
 * in a way that allows for easier simd translation.
 * two_legs --> using_64bit_mul; but it's not done the right way.
 * two_legs --> reduce_plattice_asm ; production code which combines
 * partial unrolling and 64-bit.
 * mimick_production_noasm: a C version of the asm code.
 */

/* Each implementation is now in a separate header file in
 * tests/sieve/reduce-plattice, and the two that are used in production
 * are in sieve/
 *
 * The plattice_proxy is used to allow the old-interface functions to
 * peek/poke into the plattice protected data members.
 *
 * The different implementations are also written as members as
 * plattice_proxy (and included from its definition), in order to mimick
 * correctly what could happen if we put them inside plattice_info.
 */

#include "reduce-plattice/plattice-proxy.hpp"

#include "reduce-plattice/las-reduce-plattice-reference.hpp"
#include "reduce-plattice/las-reduce-plattice-reference2.hpp"

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) &&                                \
    !(defined(__APPLE_CC__) && defined(__llvm__) && __APPLE_CC__ == 5621) &&   \
    !defined(__INTEL_COMPILER)
/* icc can't compile the avx512 code if we happen to enable this code.
 * https://community.intel.com/t5/Intel-C-Compiler/asm-callq-and-kand-mask8-intrinsic-generate-vkmovb-quot-no-such/m-p/1140906
 */
#define TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
#endif

#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
#include "reduce-plattice/las-reduce-plattice-reference2_asm.hpp"
#endif

#include "reduce-plattice/las-reduce-plattice-simd.hpp"

struct post_condition_error : public std::exception {
};
struct disagreement : public std::exception {
};

struct a_test {
    unsigned long p;
    unsigned long q;
    unsigned long r;
    int k;
};

struct test_wrap {
    uint32_t I = 512;
    bool nocheck = false;
    bool timing = false;
    bool failed = false;
    std::vector<a_test> tests;
};

struct call_production {
    static constexpr size_t const batch_count = 1;
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static constexpr char const * what = "production";
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        asm volatile("");
        L.reduce(I);
        asm volatile("");
    }
};

struct call_production_noasm {
    static constexpr size_t const batch_count = 1;
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static constexpr char const * what = "mimick_production_noasm";
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        asm volatile("");
        L.mimick_production_noasm(I);
        asm volatile("");
    }
};

struct call_two_legs {
    static constexpr size_t const batch_count = 1;
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static constexpr char const * what = "two_legs";
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        L.reduce_plattice_two_legs(I);
    }
};

struct call_simplistic {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "simplistic";
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        L.reduce_plattice_simplistic(I);
    }
};

struct call_using_64bit_mul {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "using_64bit_mul";
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        L.using_64bit_mul(I);
    }
};

struct call_swapping_loop {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "swapping_loop";
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        L.swapping_loop(I);
    }
};

struct call_swapping_loop2 {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "swapping_loop2";
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static void call(plattice_proxy & L, uint32_t I, a_test = a_test())
    {
        L.swapping_loop2(I);
    }
};

struct call_old_reference {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "old_reference";
    static constexpr bool const has_known_bugs = true;
    static constexpr bool const old_interface = true;
    static void call(plattice_proxy & L, uint32_t I, a_test const & aa)
    {
        reference(&L, aa.q, aa.r, I);
    }
};

struct call_reference2 {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "reference2";
    static constexpr bool const has_known_bugs = true;
    static constexpr bool const old_interface = true;
    static void call(plattice_proxy & L, uint32_t I, a_test const & aa)
    {
        reference2(&L, aa.q, aa.r, I);
    }
};

#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
struct call_reference2_asm {
    static constexpr size_t const batch_count = 1;
    static constexpr char const * what = "reference2_asm";
    static constexpr bool const has_known_bugs = true;
    static constexpr bool const old_interface = true;
    static void call(plattice_proxy & L, uint32_t I, a_test const & aa)
    {
        reference2_asm(&L, aa.q, aa.r, I);
    }
};
#endif

template <size_t N> struct call_simd_base {
    static constexpr size_t const batch_count = N;
    static constexpr bool const has_known_bugs = false;
    static constexpr bool const old_interface = false;
    static void call(plattice_proxy * pL, uint32_t I) { simd<N>(pL, I); }
};

struct call_simd_avx512 : public call_simd_base<16> {
    static constexpr char const * what = "simd-avx512";
};

struct call_simd_avx2 : public call_simd_base<8> {
    static constexpr char const * what = "simd-avx2";
};

template <size_t N> struct call_simd : public call_simd_base<N> {
};
template <> struct call_simd<4> : public call_simd_base<4> {
    static constexpr char const * what = "simd-synthetic<4> (sse-4.1, maybe)";
};
template <> struct call_simd<2> : public call_simd_base<2> {
    static constexpr char const * what = "simd-synthetic<2>";
};

template <typename T>
static inline unsigned long test_inner(plattice_proxy * L, test_wrap const & tw,
                                       a_test const * aa)
    requires T::old_interface
{
    T::call(*L, tw.I, *aa);
    return L->mi0;
}

template <typename T>
static inline unsigned long test_inner(plattice_proxy * L, test_wrap const & tw,
                                       a_test const * aa, size_t const M)
    requires T::old_interface
{
    unsigned long dummy_local = 0;
    for (size_t i = 0; i < M; i++)
        dummy_local += test_inner<T>(L + i, tw, aa + i);
    return dummy_local;
}

template <typename T>
static inline unsigned long test_inner(plattice_proxy * L, test_wrap const & tw,
                                       a_test const * aa)
    requires(!T::old_interface && T::batch_count == 1)
{
    bool const proj = aa->r >= aa->q;
    L->initial_basis(aa->q, proj ? (aa->r - aa->q) : aa->r, proj);
    T::call(*L, tw.I);
    return L->mi0;
}

template <typename T>
static inline unsigned long test_inner(plattice_proxy * L, test_wrap const & tw,
                                       a_test const * aa, size_t const M)
    requires(!T::old_interface && T::batch_count == 1)
{
    unsigned long dummy_local = 0;
    for (size_t i = 0; i < M; i++)
        dummy_local += test_inner<T>(L + i, tw, aa + i);
    return dummy_local;
}

template <typename T>
static inline unsigned long test_inner(plattice_proxy * L, test_wrap const & tw,
                                       a_test const * aa,
                                       size_t N = T::batch_count)
    requires(!T::old_interface && T::batch_count != 1)
{
    size_t j = 0;
    bool normal = true;
    /* prepare the lattices for all calls. */
    for (j = 0; j < N; j++) {
        // unsigned long p = tests[j].p;
        unsigned long const q = aa[j].q;
        unsigned long const r = aa[j].r;
        bool const proj = r >= q;
        L[j].initial_basis(q, proj ? (r - q) : r, proj);
        if (L[j].needs_special_treatment(tw.I))
            normal = false;
    }
    if (normal) {
        T::call(L, tw.I);
    } else {
        for (j = 0; j < N; j++)
            L[j].reduce(tw.I);
    }
    return L[0].mi0;
}

template <typename T> static void test_correctness(test_wrap & tw)
{
    if (tw.nocheck)
        return;
    constexpr size_t N = T::batch_count;
    int nfailed = 0;
    std::vector<a_test> const & tests = tw.tests;
    std::string thiscode = T::what;
    if (T::has_known_bugs)
        thiscode += " (known buggy)";

    for (size_t i = 0; i < tests.size(); i += N) {
        plattice_proxy Lref[N];
        plattice_proxy L[N];
        size_t j = 0;
        char const * when = "pre";
        try {
            if (i + N <= tests.size()) {
                test_inner<call_simplistic>(Lref, tw, &(tests[i]), N);
                test_inner<T>(L, tw, &(tests[i]), N);
            } else {
                size_t const M = tests.size() - i;
                test_inner<call_simplistic>(Lref, tw, &(tests[i]), M);
                test_inner<T>(L, tw, &(tests[i]), M);
            }
            when = "post";
            for (j = 0; j < N && i + j < tests.size(); j++) {
                if (!L[j].check_post_conditions(tw.I))
                    throw post_condition_error();
                if (L[j] != Lref[j])
                    throw disagreement();
            }
        } catch (plattice_proxy::error const & e) {
            a_test const & aa = tests[i + j];
            fmt::print(stderr,
                       "{}: failed in-algorithm check"
                       " for p^k={}^{} r={}\n",
                       thiscode, aa.p, aa.k, aa.r);
            if (!T::has_known_bugs)
                tw.failed = true;
        } catch (post_condition_error const & e) {
            if (nfailed < 16) {
                a_test const & aa = tests[i + j];
                fmt::print(stderr,
                           "{}: failed check ({}) for "
                           "p^k={}^{} r={}\n",
                           thiscode, when, aa.p, aa.k, aa.r);
                if (!T::has_known_bugs)
                    tw.failed = true;
                if (++nfailed >= 16)
                    fmt::print(
                        stderr,
                        "{}:"
                        " stopped reporting errors, go fix your program\n",
                        thiscode);
            }
        } catch (disagreement const & e) {
            if (nfailed < 16) {
                a_test const & aa = tests[i + j];
                std::string msg = fmt::format(
                    "{}: disagreement with simplistic on p^k={}^{} r={}\n",
                    thiscode, aa.p, aa.k, aa.r);
                msg += fmt::format("simplistic: [({}, {}), ({}, {})]\n",
                                   Lref[j].get_i0(), Lref[j].get_j0(),
                                   Lref[j].get_i1(), Lref[j].get_j1());
                msg += fmt::format("{}: [({}, {}), ({}, {})]\n", thiscode,
                                   L[j].get_i0(), L[j].get_j0(), L[j].get_i1(),
                                   L[j].get_j1());
                fputs(msg.c_str(), stderr);
                if (!T::has_known_bugs)
                    tw.failed = true;
                if (++nfailed >= 16)
                    fmt::print(
                        stderr,
                        "{}:"
                        " stopped reporting errors, go fix your program\n",
                        thiscode);
            }
        }
    }
}

template <typename T>
static inline unsigned long test_speed(test_wrap & tw)
    requires(T::batch_count != 1)
{
    test_correctness<T>(tw);
    if (!tw.timing)
        return 0;
    clock_t const clk0 = clock();
    unsigned long dummy_local = 0;
    size_t i;
    for (i = 0; i + T::batch_count <= tw.tests.size(); i += T::batch_count) {
        plattice_proxy L[T::batch_count];
        dummy_local += test_inner<T>(L, tw, &tw.tests[i]);
    }
    if (i < tw.tests.size()) {
        plattice_proxy L[T::batch_count];
        dummy_local += test_inner<T>(L, tw, &tw.tests[i], tw.tests.size() - i);
    }
    clock_t const clk1 = clock();
    if (tw.timing) {
        fmt::print("# {}: {} tests in {:.4f}s\n", T::what, tw.tests.size(),
                   ((double)(clk1 - clk0)) / CLOCKS_PER_SEC);
    }
    return dummy_local;
}

template <typename T>
static inline unsigned long test_speed(test_wrap & tw)
    requires(T::batch_count == 1)
{
    test_correctness<T>(tw);
    if (!tw.timing)
        return 0;
    clock_t const clk0 = clock();
    unsigned long dummy_local = 0;
    for (a_test const & aa: tw.tests) {
        plattice_proxy L[T::batch_count];
        dummy_local += test_inner<T>(L, tw, &aa);
    }
    clock_t const clk1 = clock();
    if (tw.timing) {
        fmt::print("# {}: {} tests in {:.4f}s\n", T::what, tw.tests.size(),
                   ((double)(clk1 - clk0)) / CLOCKS_PER_SEC);
    }
    return dummy_local;
}

int main(int argc, char const * argv[])
{
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    test_wrap tw;

    uint32_t B = 1e6;
    size_t ntests = 1000;
    unsigned long seed = 0;
    bool quiet = false;

    for (; argc > 1; argv++, argc--) {
        if (strcmp(argv[1], "-B") == 0) {
            argv++, argc--;
            B = atoi(argv[1]);
        } else if (strcmp(argv[1], "-I") == 0) {
            argv++, argc--;
            tw.I = atoi(argv[1]);
        } else if (strcmp(argv[1], "-ntests") == 0) {
            argv++, argc--;
            ntests = atoi(argv[1]);
        } else if (strcmp(argv[1], "-seed") == 0) {
            argv++, argc--;
            seed = atol(argv[1]);
        } else if (strcmp(argv[1], "-q") == 0) {
            quiet = true;
        } else if (strcmp(argv[1], "-N") == 0) {
            tw.nocheck = true;
        } else if (strcmp(argv[1], "-T") == 0) {
            tw.timing = true;
        }
    }

    std::vector<std::tuple<unsigned long, int>> prime_powers;
    /* count primes up to B, and add prime powers as well */
    prime_info pi;
    prime_info_init(pi);
    for (unsigned long p = 2; p < B; p = getprime_mt(pi)) {
        int k = 1;
        unsigned long q = p;
        for (; q < B; q *= p, k++) {
            if (q < tw.I)
                continue;
            prime_powers.emplace_back(p, k);
        }
    }
    prime_info_clear(pi);

    std::map<int, std::map<std::string, int>> stats;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    if (!seed) {
        seed = getpid();
        if (!quiet)
            fmt::print("seeding with pid = {}\n", seed);
    }
    gmp_randseed_ui(rstate, seed);

    std::vector<a_test> & tests = tw.tests;

    /* interesting with I=0x20000 */
    tests.emplace_back(a_test {5, 390625, 771795, 8});

    /* we have sporadic failures with this one */
    tests.emplace_back(a_test {1579, 1579, 1579, 1});

    for (; tw.tests.size() < ntests;) {
        unsigned long const j = gmp_urandomm_ui(rstate, prime_powers.size());
        unsigned long p;
        int k;
        std::tie(p, k) = prime_powers[j];
        unsigned long q = 1;
        for (int s = k; s--; q *= p)
            ;
        unsigned long r = gmp_urandomm_ui(rstate, q + q / p);
        bool const proj = r >= q;
        if (proj)
            r = q + p * (r - q);

        a_test const aa {p, q, r, k};

        tests.emplace_back(aa);

        if (!quiet) {
            std::string desc = aa.r > aa.q ? "proj" : "affine";
            if (aa.p == 2)
                desc += "+even";
            stats[aa.k][desc]++;
        }
    }

    {
        /* verify that all tests are in range */
        auto jt = tests.begin();
        for (auto const & aa: tests) {
            plattice_proxy L;
            bool const proj = aa.r >= aa.q;
            L.initial_basis(aa.q, proj ? (aa.r - aa.q) : aa.r, proj);
            if (L.check_pre_conditions(tw.I)) {
                *jt++ = aa;
            } else {
                fmt::print(stderr,
                           "skipping out-of-range test for p^k={}^{} r={}\n",
                           aa.p, aa.k, aa.r);
            }
        }
        tests.erase(jt, tests.end());
    }

    /* declare it as volatile so that the compiler doesn't outsmart us.
     */
    unsigned long volatile dummy MAYBE_UNUSED = 0;

    dummy += test_speed<call_production>(tw);
    dummy += test_speed<call_production_noasm>(tw);

    dummy += test_speed<call_simplistic>(tw);
    dummy += test_speed<call_two_legs>(tw);
    dummy += test_speed<call_using_64bit_mul>(tw);
    dummy += test_speed<call_swapping_loop>(tw);
    dummy += test_speed<call_swapping_loop2>(tw);

    if (tw.timing) {
        /* as far as correctness is concerned, we're not interested in
         * this code. We *know* that it is buggy. But when we report
         * timings, we want to know how it fares !
         */
        dummy += test_speed<call_old_reference>(tw);
        dummy += test_speed<call_reference2>(tw);
#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
        dummy += test_speed<call_reference2_asm>(tw);
#endif
    }

#ifdef HAVE_AVX512F
    dummy += test_speed<call_simd_avx512>(tw);
#endif

#ifdef HAVE_AVX2
    dummy += test_speed<call_simd_avx2>(tw);
#endif

    dummy += test_speed<call_simd<4>>(tw);
    dummy += test_speed<call_simd<2>>(tw);

    if (!quiet) {
        /* Finish with that one, but here we don't care about the
         * correctness, we only want to do the frequency table of quotients.
         */
        std::map<int, unsigned long> T;
        for (auto const & aa: tests) {
            bool const proj = aa.r >= aa.q;
            plattice_proxy L;
            L.initial_basis(aa.q, proj ? (aa.r - aa.q) : aa.r, proj);
            if (L.needs_special_treatment(tw.I)) {
                T[-1]++;
                L.reduce_with_vertical_vector(tw.I);
            } else {
                L.instrumented_two_legs(tw.I, T);
            }
        }
        for (auto & kv: stats) {
            fmt::print("{}: {}\n", kv.first,
                       join(kv.second.begin(), kv.second.end(), " ",
                            [](decltype(kv.second)::value_type const & v) {
                                return fmt::format("{}:{}", v.first, v.second);
                            }));
        }
        unsigned long sumT = 0;
        for (auto const & x: T) {
            sumT += x.second;
        }
        unsigned long cumulative = 0;
        for (auto const & x: T) {
            fmt::print("{}: {} ({:.1f}%%) ({:.1f}%%)\n", x.first, x.second,
                       100.0 * x.second / sumT,
                       100.0 * (cumulative += x.second) / sumT);
            if (x.first >= 16)
                break;
        }
    }
    gmp_randclear(rstate);
    return tw.failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

// see bug #30052
#if GNUC_VERSION_ATLEAST(12, 0, 0) && !GNUC_VERSION_ATLEAST(13, 0, 0)
#pragma GCC diagnostic pop
#endif
