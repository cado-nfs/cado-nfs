#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <memory>
#include <vector>
#include <string>

#include "fmt/base.h"
#include <gmp.h>

#include "gf2x-fft.h"
#include "gmp_aux.h"
#include "lingen_mat_types.hpp"
#include "macros.h"
#include "misc.h"
#include "portability.h"
#include "utils_cxx.hpp"

// scan-headers: stop here

static void usage()
{
    fmt::print(stderr, "Usage: ./bench_polmatmul [--nrep <k>] <N> <m> <n>\n");
    exit(1);
}

/* Counted in unsigned longs -- must exceed all caches by large, so as to
 * avoid cheating on benches */
#define DATA_POOL_SIZE (1 << 23)

/* handy globals */
static unsigned int nrep = 10;

#if 0
void bench_polmul(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    fmt::print("gf2x_mul ({} x {}):", d1, d2);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        gf2x_mul(h, f, nw1, g, nw2);
    }
    fmt::print(" {}", (double) (cputime()-tt)/nrep);
    fmt::print("\n");

    delete[] f;
    delete[] g;
    delete[] h;
}
#endif

#define NREPS_MAX 1000
#define REPEAT_TIME_MAX 0.1

template <typename F>
static double one_bench(F const & f)
{
    clock_t const tt = clock();
    unsigned int r;
    clock_t const clocklim = tt + lround(REPEAT_TIME_MAX * CLOCKS_PER_SEC);
    for (r = 0; r < NREPS_MAX && clock() < clocklim; r++) {
        f();
    }
    if (!r) return 0;
    return double(clock() - tt) / r / CLOCKS_PER_SEC;
}

#if 0
void bench_c128(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt0;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    c128_info_t o;

    c128_init(o, d1, d2);

    c128_t * tf = c128_alloc(o, 1);
    c128_t * tg = c128_alloc(o, 1);
    c128_t * th = c128_alloc(o, 1);

    fmt::print("c128 ({} x {}):", d1, d2);

    tt0 = cputime();
    double v;

    double v = one_bench(c128_dft, o, tf, f, d1);
    double v = one_bench(c128_dft, o, tg, g, d2);
    double v = one_bench(c128_compose, o, th, tf, tg);
    double v = one_bench(c128_ift, o, h, d1 + d2, th);

    c128_clear(o);
    fmt::print(" tot {}", (double) (cputime()-tt0)/nrep);
    fmt::print("\n");
    delete[] f;
    delete[] g;
    delete[] h;
}
#endif

/* Since strassen is a recursive algorithm that falls back on matrices
 * whose sizes are halved, and since for our applications, the matrices
 * do have a large power of 2 in their dimension, we consider different
 * possible dimensions of the form (d1*2^i,d2*2^i) times (d2*2^i,d3*2^i),
 * for all values of i from 0 to some bound, and for all value of d from
 * 1 to some power of 2 (of course 0 for d makes no sense).
 */

#define BITS_IN_DIM_D 3
#define BITS_IN_DIM_I 4

/* This is in case 1 << BITS_IN_DIM_I is unacceptably large. */
#define MAX_I_FOR_TUNING 8 /* < (1 << BITS_IN_DIM_I) */

// dimension must be <= d * 2^i
// d <= (1 << BITS_IN_DIM_D)
// i < (1 << BITS_IN_DIM_I)
//
// one tuning table entry contains
// 3*BITS_IN_DIM_D+BITS_IN_DIM_I bits
//
// representing a max dimension (1 << BITS_IN_DIM_D) * (1 << ((1 <<
// BITS_IN_DIM_I) - 1)) ; however not all dimensions can be represented
// up to this size, for obvious reasons.
//
// 1,3 means at most 256        ;  6 bits per entry index
// 2,3 means at most 512        ;  9 bits per entry index
// 1,4 means at most 65536      ;  7 bits per entry index
// 2,4 means at most 131072     ;  10 bits per entry index
// 3,4 means at most 262144     ;  13 bits per entry index

// contains the value n0 such that for nbits >= n0, we'd better use strassen

#define SELECTOR_NB_INDICES (1 << (3 * BITS_IN_DIM_D + BITS_IN_DIM_I))

struct my_strassen_selector {

    unsigned int strassen_threshold[SELECTOR_NB_INDICES];

    /* give a number-of-words threshold value for multiplying two
     * matrices of sizes m*n and n*p, depending on BITS_IN_DIM_D and
     * BITS_IN_DIM_I */
    static unsigned int threshold_index(unsigned int m, unsigned int n,
                                        unsigned int p)
    {
        unsigned int const combined = m | n | p;
        unsigned int b = cado_ctzl(combined);
        if (b >= (1UL << BITS_IN_DIM_I)) {
            // arrange so that we
            b = (1UL << BITS_IN_DIM_I) - 1;
        }
        m >>= b;
        m -= 1;
        ASSERT(m < (1 << BITS_IN_DIM_D));
        n >>= b;
        n -= 1;
        ASSERT(n < (1 << BITS_IN_DIM_D));
        p >>= b;
        p -= 1;
        ASSERT(p < (1 << BITS_IN_DIM_D));

        unsigned int index = 0;
        index = m;
        index = (index << BITS_IN_DIM_D) | n;
        index = (index << BITS_IN_DIM_D) | p;
        index = (index << BITS_IN_DIM_I) | b;
        return index;
    }

    unsigned int const & threshold(unsigned int m, unsigned int n,
                                   unsigned int p) const
    {
        return strassen_threshold[threshold_index(m, n, p)];
    }

    unsigned int & threshold(unsigned int m, unsigned int n, unsigned int p)
    {
        return strassen_threshold[threshold_index(m, n, p)];
    }

    int operator()(unsigned int m, unsigned int n, unsigned int p,
                   unsigned int nbits) const
    {
        int const answer = nbits >= threshold(m, n, p);
        // fmt::print("([{},{},{},{}]=>{})",m,n,p,nbits,answer);
        return answer;
    }
    void dump(std::string const & name) const
    {
        fmt::print("#ifdef  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS\n");
        fmt::print("#define {}_STRASSEN_THRESHOLDS_D {}\n", name,
                   BITS_IN_DIM_D);
        fmt::print("#define {}_STRASSEN_THRESHOLDS_I {}\n", name,
                   BITS_IN_DIM_I);
        fmt::print("#define {}_STRASSEN_THRESHOLDS {{\t\\\n\t", name);
        unsigned int disp = 0;
        for (unsigned int i = 0; i < SELECTOR_NB_INDICES; i++) {
            if (strassen_threshold[i]) {
                fmt::print(" {{{}, {}}},", i, (int)strassen_threshold[i]);
                if (++disp % 4 == 0) {
                    fmt::print("\t\\\n\t");
                }
            }
        }
        fmt::print("}}\n");
        fmt::print("#else   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */\n");
        auto s = name;
        std::ranges::transform(s, s.begin(),
                       [](char c) { return std::tolower(c); });
        fmt::print(
            "template<> unsigned int foo<{}>::default_selector_data[] = {{\n",
            s);
        fmt::print("\t/* D = {} */\n", BITS_IN_DIM_D);
        fmt::print("\t/* I = {} */\n", BITS_IN_DIM_I);
        fmt::print("\t");
        disp = 0;
        for (unsigned int i = 0; i < SELECTOR_NB_INDICES; i++) {
            if (strassen_threshold[i]) {
                fmt::print(" [{}] = {},", i, (int)strassen_threshold[i]);
                if (++disp % 4 == 0) {
                    fmt::print("\t\\\n\t");
                }
            }
        }
        fmt::print("}}\n");
        fmt::print("#endif   /*  STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS */\n");
    }
};

struct logbook {
    double t1;
    double t2;
    double c;
    double i;
};

template <typename T> struct foo {
    static logbook l;
    static my_strassen_selector s;
    static unsigned int default_selector_data[];
};

template <> logbook foo<gf2x_fake_fft_info>::l = logbook();
template <> logbook foo<gf2x_cantor_fft_info>::l = logbook();
template <> logbook foo<gf2x_ternary_fft_info>::l = logbook();
template <>
my_strassen_selector foo<gf2x_fake_fft_info>::s = my_strassen_selector();
template <>
my_strassen_selector foo<gf2x_cantor_fft_info>::s = my_strassen_selector();
template <>
my_strassen_selector foo<gf2x_ternary_fft_info>::s = my_strassen_selector();

#define STRASSEN_THRESHOLDS_AS_CPP_CONSTANTS
#include "strassen-thresholds.hpp" // IWYU pragma: keep
/*
template<> unsigned int foo<gf2x_fake_fft_info>::default_selector_data[] = {};
template<> unsigned int foo<gf2x_cantor_fft_info>::default_selector_data[] = {};
template<> unsigned int foo<gf2x_ternary_fft_info>::default_selector_data[] =
{};
*/

template <typename T> struct pointed_type {
};
template <typename T> struct pointed_type<T *> {
    using type = T;
};
template <typename T> struct pointed_type<T const *> {
    using type = T;
};

template <typename fft_type>
static int depth(fft_type const & o, size_t d1, size_t d2, size_t d3)
{
    int k;
    my_strassen_selector const & s(foo<fft_type>::s);
    const size_t nbits = o.size0_bytes() * CHAR_BIT;
    for (k = 0; (((d1 | d2 | d3) >> k) & 1UL) == 0; k++) {

        if (!s(d1 >> k, d2 >> k, d3 >> k, nbits))
            break;
    }
    return k;
}

template <typename fft_type>
static unsigned long nmults(fft_type const & o, size_t d1, size_t d2, size_t d3)
{
    const int d = depth(o, d1, d2, d3);
    unsigned long nmuls = d1 * d2 * d3;
    for (int k = 0; k < d; k++) {
        nmuls /= 8;
        nmuls *= 7;
    }
    return nmuls;
}

/* One step:
 *
 * pi matrix (m+n) * (m+n), length ceil(d/2)*m/(m+n) ; two transforms.
 * E matrix m * (m+n), length ceil(d/2)*(1+m/(m+n)) = d-ceil(d/2)*n/(m+n) ; one
 * transform
 *
 * length of E*pi d
 * length of pi*pi d*m/(m+n)
 *
 */

static inline size_t W(size_t x)
{
    return (x + ULONG_BITS - 1) / ULONG_BITS;
}
static inline size_t I(size_t x)
{
    return x / ULONG_BITS;
}
static inline size_t R(size_t x)
{
    return x % ULONG_BITS;
}
static inline unsigned long MASK(size_t x)
{
    return (1UL << R(x)) - 1UL;
}

static std::vector<unsigned long> random_bitstring(size_t n,
                                                   cxx_gmp_randstate & state)
{
    std::vector<unsigned long> data(W(n));
    memfill_random(data.data(), W(n) * sizeof(unsigned long), state);
    if (R(n))
        data[I(n)] &= MASK(n);
    return data;
}

template <typename T>
static void fft_times(double & dft1, double & dft2, double & compose,
                      double & ift, T & o, unsigned long n1, unsigned long n2,
                      cxx_gmp_randstate & state)
{
    typename T::ptr f = o.alloc(1);
    typename T::ptr g = o.alloc(1);
    /* Make sure we're timing for at least a second. */
    unsigned long const n3 = n1 + n2 - 1;

    using elt = typename pointed_type<typename T::ptr>::type;

    auto temp1 = std::make_unique<elt[]>(o.size1_bytes() / sizeof(elt));
    auto temp2 = std::make_unique<elt[]>(o.size2_bytes() / sizeof(elt));

    {
        auto data1 = random_bitstring(n1, state);
        dft1 = one_bench([&](){o.dft(f, data1.data(), n1, temp1.get());});
        fmt::print(" {:.2f}", dft1 * 1.0e6);
    }
    {
        auto data2 = random_bitstring(n2, state);
        dft2 = one_bench([&](){o.dft(g, data2.data(), n2, temp1.get());});
        fmt::print(" {:.2f}", dft2 * 1.0e6);
    }
    {
        compose = one_bench([&](){o.compose(g, g, f, temp2.get());});
        fmt::print(" {:.2f}", compose * 1.0e6);
    }
    {
        auto data3 = random_bitstring(n3, state);
        ift = one_bench([&](){o.ift(data3.data(), n3, g, temp1.get());});
        fmt::print(" {:.2f}", ift * 1.0e6);
    }

    o.free(f, 1);
    o.free(g, 1);
}

#if 0
/* One step:
 *
 * pi matrix (m+n) * (m+n), length d/2*m/(m+n) ; two transforms.
 * E matrix m * (m+n), length d/2*(1+m/(m+n)) = d-d/2*n/(m+n) ; one transform
 *
 * length of E*pi d
 * length of pi*pi d*m/(m+n)
 *
 */

    template<typename T>
void bench_polmatmul(const char * s,
        unsigned int n1, unsigned int n2,
        unsigned int d1, unsigned int d2)
{
    unsigned int nc1 = d1 + 1;
    unsigned int nc2 = d2 + 1;
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);

    fmt::print("polmatmul<{}> (dim {} x {} deg {} x {})", s, n1, n2, d1, d2);

    polmat f(n1, n2, nc1);
    polmat g(n2, n2, nc2);
    polmat h(n1, n2, nc1 + nc2 - 1);

    f.randomize(state);
    g.randomize(state);

    T o(d1, d2);

    tpolmat<T> tf(n1,n2,o);
    tpolmat<T> tg(n2,n2,o);
    tpolmat<T> th(n1,n2,o);

    double foo<T>::l.t1 = one_bench(transform, tf, f, o, d1);
    double foo<T>::l.t2 = one_bench(transform, tg, g, o, d2);
    double foo<T>::l.c = one_bench(compose, th, tf, tg, o);
    double foo<T>::l.i = one_bench(itransform, h, th, o, d1 + d2);

    double t = foo<T>::l.t1 + foo<T>::l.t2 + foo<T>::l.c + foo<T>::l.i;

    fmt::print(" tot {}", t);
    fmt::print("\n");
}
#endif

/* The fft_type const& argument was meant as a "base" fft type with which
 * we would strive to be compatible. Unfortunately this functionality got
 * killed in gf2x in 2019
 */
template <typename fft_type>
static void tune_strassen1(fft_type const &, unsigned int d1, unsigned int d2,
                           unsigned int d3, size_t maxlen,
                           cxx_gmp_randstate & state)
{
    my_strassen_selector & s(foo<fft_type>::s);

    unsigned int const bp = cado_ctzl(d1 | d2 | d3);
    unsigned int dd1 = d1 >> bp;
    unsigned int dd2 = d2 >> bp;
    unsigned int dd3 = d3 >> bp;
    fmt::print(" -- tuning strassen for {}x{} * {}x{}\n", d1, d2, d2, d3);
    fmt::print(" -- max {} levels above {}x{} * {}x{} --\n", bp, dd1, dd2, dd2,
               dd3);
    fmt::print(" -- maximum length is {} bits\n",
               iceildiv(maxlen, ULONG_BITS) * ULONG_BITS);
    ASSERT((dd1 - 1) < (1 << BITS_IN_DIM_D));
    ASSERT((dd2 - 1) < (1 << BITS_IN_DIM_D));
    ASSERT((dd3 - 1) < (1 << BITS_IN_DIM_D));
    // for the very small matrices, we just _cannot_ use strassen because
    // there is an odd dimension hanging around.
    s.threshold(dd1, dd2, dd3) = UINT_MAX;
    fmt::print("Strassen threshold for {}x{} * {}x{} is {} bits [index {}]\n",
               dd1, dd2, dd2, dd3, UINT_MAX,
               my_strassen_selector::threshold_index(dd1, dd2, dd3));
    dd1 <<= 1;
    dd2 <<= 1;
    dd3 <<= 1;
    // If we know that for smaller matrices, strassen was paying off from
    // nbits on, then of course it is going to pay off at worst at this
    // size.
    // old_wt stands for old threshold in words.
    unsigned int const max_wt = iceildiv(maxlen, ULONG_BITS);
    unsigned int const min_wt = 0;
    unsigned int old_wt = max_wt;
    size_t earliest_good_strassen = UINT_MAX;
    for (; dd1 <= d1; dd1 <<= 1, dd2 <<= 1, dd3 <<= 1) {
        fmt::print("Tuning for {}x{} * {}x{}\n", dd1, dd2, dd2, dd3);
        if (old_wt == min_wt) {
            fmt::print("[keeping low threshold]\n");
            s.threshold(dd1, dd2, dd3) = 0;
            continue;
        }
        // in order to avoid accumulating errors, we introduce some
        // looseness here. Setting wt1 to twice the old wt value will
        // force re-checking at the previous cutoff value. If ever this
        // cutoff value was wrong, we'll settle for one which sits above
        // in the interval.
        unsigned int wt1 = 2 * old_wt;
        unsigned int wt0 = min_wt;
        for (; wt1 - wt0 > 1 && wt0 < max_wt;) {
            // assuming cubic is better for 64 * t0, ans strassen better
            // for 64 * t1
            unsigned int const wt = (wt1 + wt0) / 2;
            // unsigned int t1 = wt1 * ULONG_BITS;
            // unsigned int t0 = wt0 * ULONG_BITS;
            unsigned int t = wt * ULONG_BITS;
            if (t == 0) {
                t = 2;
            }
            fft_type o(t / 2, t / 2);
            unsigned int tr = o.size0_bytes() * CHAR_BIT;
            tpolmat<fft_type> tf(dd1, dd2, o);
            tpolmat<fft_type> tg(dd2, dd3, o);
            tpolmat<fft_type> th(dd1, dd3, o);

            tf.randomize(state);
            tg.randomize(state);
            fmt::print("{} [{}]", t, tr);
            double t_cubic;
            double t_strassen;
            // force cubic
            s.threshold(dd1, dd2, dd3) = UINT_MAX;

            {
                auto F = [&](){compose_inner(th, tf, tg, o, s);};
                t_cubic = one_bench(F);
                fmt::print(" {:.2g}", t_cubic);
                // force strassen
                s.threshold(dd1, dd2, dd3) = 0;
                t_strassen = one_bench(F);
                fmt::print(" {:.2g}", t_strassen);
            }
            if (t_cubic <= t_strassen) {
                fmt::print(" cubic [{:.2f}]\n", t_cubic / t_strassen);
                wt0 = wt;
            } else {
                fmt::print(" strassen [{:.2f}]\n", t_cubic / t_strassen);
                wt1 = wt;
                earliest_good_strassen = tr;
            }
        }
        wt1 = std::min(wt1, max_wt);
        old_wt = wt1;
        if (wt1 == max_wt)
            wt1 = UINT_MAX;
        if (wt1 == min_wt)
            wt1 = 0;
        if (wt1 == 0)
            earliest_good_strassen = 0;
        if (wt1 == UINT_MAX)
            earliest_good_strassen = UINT_MAX;
        s.threshold(dd1, dd2, dd3) = earliest_good_strassen;
        fmt::print(
            "Strassen threshold for {}x{} * {}x{} is {} bits [index {}]\n", dd1,
            dd2, dd2, dd3, s.threshold(dd1, dd2, dd3),
            my_strassen_selector::threshold_index(dd1, dd2, dd3));
    }
    fmt::print("\n");
}

template <typename fft_type>
static void tune_strassen(fft_type const & base, size_t maxlen,
                          cxx_gmp_randstate & state)
{
    my_strassen_selector & s(foo<fft_type>::s);
    for (unsigned int i = 1; i <= (1 << BITS_IN_DIM_D); i++) {
        for (unsigned int j = 1; j <= (1 << BITS_IN_DIM_D); j++) {
            for (unsigned int k = 1; k <= (1 << BITS_IN_DIM_D); k++) {
                tune_strassen1(base, i << MAX_I_FOR_TUNING,
                               j << MAX_I_FOR_TUNING, k << MAX_I_FOR_TUNING,
                               maxlen, state);
            }
        }
    }
    for (int i = 0; i < SELECTOR_NB_INDICES; i++) {
        unsigned int t = s.strassen_threshold[i];
        if (t)
            fmt::print("{} {}\n", i, t);
    }
}

template <typename fft_type>
static void plot_compose(char const * name MAYBE_UNUSED, unsigned int n1,
                         unsigned int n2, unsigned int n3, unsigned long wt,
                         cxx_gmp_randstate & state)
{
    my_strassen_selector & s(foo<fft_type>::s);

    unsigned int nx1 = n1;
    unsigned int nx2 = n2;
    unsigned int nx3 = n3;
    unsigned long nmul = 1;
    for (; !((nx1 | nx2 | nx3) & 1UL);) {
        if (!s(nx1, nx2, nx3, wt * ULONG_BITS))
            break;
        nx1 >>= 1;
        nx2 >>= 1;
        nx3 >>= 1;
        nmul *= 7;
    }
    nmul *= nx1 * nx2 * nx3;
    fmt::print("For {}x{} * {}x{}, strassen begins at {}x{} * {}x{} ; roughly "
               "{} mults\n",
               n1, n2, n2, n3, nx1, nx2, nx2, nx3, nmul);

    for (; !((nx1 | nx2 | nx3) & 1UL); nx1 >>= 1, nx2 >>= 1, nx3 >>= 1)
        ;

    /* Start by a measurement of the unit time */
    fft_type o(wt / 2, wt / 2);
    double single_compose;
    {
        typename fft_type::ptr pf = o.alloc(1);
        typename fft_type::ptr pg = o.alloc(1);
        typename fft_type::ptr ph = o.alloc(1);
        single_compose = one_bench([&](){o.compose(ph, pf, pg);});
        o.clear(pf, 1);
        o.clear(pg, 1);
        o.clear(ph, 1);
    }

    /* Now bench all possible matrix timings. */
    nmul = nx1 * nx2 * nx3;
    for (; nx1 < n1;) {
        tpolmat<fft_type> tf(nx1, nx2, o);
        tpolmat<fft_type> tg(nx2, nx3, o);
        tpolmat<fft_type> th(nx1, nx3, o);

        tf.randomize(state);
        tg.randomize(state);
        double v = one_bench(compose_inner, th, tf, tg, o, s);
        fmt::print("{}x{} * {}x{} : {:.2e} = {:.2f}%% of {} * {:.2e}\n",
                nx1, nx2, nx2, nx3,
                v, 100.0 * double_ratio(v, nmul, single_compose),
                nmul, single_compose);
        nx1 <<= 1;
        nx2 <<= 1;
        nx3 <<= 1;
        nmul *= s(nx1, nx2, nx3, wt * ULONG_BITS) ? 7 : 8;
    }
}

typedef struct {
    char const * engine;
    int nstrassen;
    double dft;
    double ift;
    double compose;
} level_info;

static bool operator<(level_info const & a, level_info const & b)
{
    return (a.dft + a.compose + a.ift) < (b.dft + b.compose + b.ift);
}

template <typename fft_type>
static level_info
bench_one_polmm_projected_sub(fft_type & o, unsigned long d1, unsigned long d2,
                              unsigned long d3, unsigned long n1,
                              unsigned long n2, cxx_gmp_randstate & state)
{
    level_info res;
    double dft1, dft2, compose, ift;
    fft_times(dft1, dft2, compose, ift, o, n1, n2, state);
    double total_compose_time = compose * nmults(o, d1, d2, d3);
    int const d = depth(o, d1, d2, d3);

    res.nstrassen = d;
    res.dft = dft1 * double(d1 * d2) + dft2 * double(d2 * d3);
    res.ift = ift * double(d1 * d3);
    res.compose = total_compose_time;
    res.engine = fft_type::name;

    fmt::print(" --> dft: {:.2f} ; mul(d{}): {:.2f} ; ift: {:.2f}\n",
               (dft1 * double(d1 * d2) + dft2 * double(d2 * d3)), d, total_compose_time,
               ift * double(d1 * d3));
    return res;
}

static level_info bench_one_polmm_projected(unsigned long d1, unsigned long d2,
                                            unsigned long d3, unsigned long n1,
                                            unsigned long n2,
                                            cxx_gmp_randstate & state)
{
    fmt::print("Timings {}x{} ({}-bit entries)"
               " times {}x{} ({}-bit entries) [projected timings]\n",
               d1, d2, n1, d2, d3, n2);

    gf2x_cantor_fft_info oc(n1, n2);
    fmt::print("c128:");
    level_info lc = bench_one_polmm_projected_sub(oc, d1, d2, d3, n1, n2, state);
    lc.engine = "c128";

    gf2x_fake_fft_info of(n1, n2);
    fmt::print("fake:");
    level_info lf = bench_one_polmm_projected_sub(of, d1, d2, d3, n1, n2, state);
    lf.engine = "fake";

    gf2x_ternary_fft_info os(n1, n2);
    fmt::print("tfft({}):", 81);
    level_info ls = bench_one_polmm_projected_sub(os, d1, d2, d3, n1, n2, state);
    os.adjust(GF2X_FFT_ADJUST_DEPTH, 81);
    ls.engine = "tfft(81)";

    return std::min(std::min(lc, lf), ls);
}

template <typename fft_type>
static level_info
bench_one_polmm_complete_sub(fft_type & o, unsigned long d1, unsigned long d2,
                             unsigned long d3, unsigned long n1,
                             unsigned long n2, cxx_gmp_randstate & state)
{
    my_strassen_selector & s(foo<fft_type>::s);

    polmat f(d1, d2, n1);
    polmat g(d2, d2, n2);
    polmat h(d1, d3, n1 + n2 - 1);
    tpolmat<fft_type> tf(d1, d2, o);
    tpolmat<fft_type> tg(d2, d3, o);
    tpolmat<fft_type> th(d1, d3, o);

    f.randomize(state);
    g.randomize(state);
    tf.randomize(state);
    tg.randomize(state);

    logbook & l(foo<fft_type>::l);

    level_info res;

    res.engine = fft_type::name;

    {
        l.t1 = one_bench([&](){transform(tf, f, o, n1);});
        fmt::print(" dft1: {:.2f}", l.t1 * 1);
    }
    {
        l.t2 = one_bench([&](){transform(tg, g, o, n2);});
        fmt::print(" dft2: {:.2f}", l.t2 * 1);
    }
    {
        l.c = one_bench([&](){compose_inner(th, tf, tg, o, s);});
        fmt::print(" compose: {:.2f}", l.c * 1);
    }

#if 0
    fmt::print(" [");
    /* Do composition by checking how deep strassen might pay off. */
    for(unsigned int k = 0 ; (((d1 | d2 | d3) >> k) & 1UL) == 0 ; k++) {
        // Force naive
        s.threshold(d1>>k,d2>>k,d3>>k) = UINT_MAX;
        double naive = one_bench(compose_inner, th, tf, tg, o, s);
        // Force strassen
        s.threshold(d1>>k,d2>>k,d3>>k) = 0;
        double strassen = one_bench(compose_inner, th, tf, tg, o, s);
        if (naive < strassen) {
            l.c = naive;
            break;
        }
        l.c = strassen;
    }
    fmt::print("]");
#endif

    {
        l.i = one_bench([&](){itransform(h, th, o, n1 + n2 - 1);});
        fmt::print(" ift: {:.2f}", l.i * 1);
    }
    fmt::print("\n");

    res.nstrassen = depth(o, d1, d2, d3);
    res.dft = l.t1 + l.t2;
    res.compose = l.c;
    res.ift = l.i;

    return res;
}

static level_info bench_one_polmm_complete(unsigned long d1, unsigned long d2,
                                           unsigned long d3, unsigned long n1,
                                           unsigned long n2, cxx_gmp_randstate & state)
{
    fmt::print("Timings {}x{} ({}-bit entries)"
               " times {}x{} ({}-bit entries) [complete timings]\n",
               d1, d2, n1, d2, d3, n2);

    gf2x_cantor_fft_info oc(n1, n2);
    fmt::print("c128:");
    level_info lc = bench_one_polmm_complete_sub(oc, d1, d2, d3, n1, n2, state);
    lc.engine = "c128";

    gf2x_fake_fft_info of(n1, n2);
    fmt::print("fake:");
    level_info lf = bench_one_polmm_complete_sub(of, d1, d2, d3, n1, n2, state);
    lf.engine = "fake";

    gf2x_ternary_fft_info os(n1, n2);
    fmt::print("tfft({}):", 81);
    level_info ls = bench_one_polmm_complete_sub(os, d1, d2, d3, n1, n2, state);
    os.adjust(GF2X_FFT_ADJUST_DEPTH, 81);
    ls.engine = "tfft(81)";

    return std::min(std::min(lc, lf), ls);
}

static void tune_strassen_global(unsigned long m, unsigned long n,
                                 unsigned long N, cxx_gmp_randstate & state)
{
    unsigned long const b = m + n;
#if 1
    /* Tune for E * pi */
    /* disabling, as in reality it's a middle product, and we're not
     * really doing a middle product, which is annoying.
     */
    if (false) {
        size_t const d = (N * b / m / n);

        unsigned long const dl = d - d / 2;
        unsigned long const pi_l_len = dl * m / b; // always <= dl
        /* it's a middle product. So we have something like
         * d * pi_l_len -> d/2, chopping off all dl first coeffs. So
         * that dl-pi_l_len coeffs of E are unused, which is the (chop)
         * value below */
        unsigned long const chop = dl - pi_l_len;

        size_t const n1 = d - chop;
        size_t const n2 = pi_l_len;
        fmt::print("Top-level multiplications E*pi: len {}, {}x{} * len "
                   "{}, {}x{} (alpha={:.2f})\n",
                   n1, m, b, n2, b, b, double_ratio(n1, n2));

        const gf2x_cantor_fft_info oc(n1, n2);
        fmt::print("=== c128 ===\n");
        tune_strassen1(oc, m, b, b, n1 + n2 - 1, state);

        const gf2x_fake_fft_info of(n1, n2);
        fmt::print("=== fake ===\n");
        tune_strassen1(of, m, b, b, n1 + n2 - 1, state);

        gf2x_ternary_fft_info os(n1, n2);
        os.adjust(GF2X_FFT_ADJUST_DEPTH, 81);
        fmt::print("=== tfft({}) ===\n", 81);
        tune_strassen1(os, m, b, b, n1 + n2 - 1, state);

        fmt::print("Options for composition at top level E*pi\n");
        fmt::print("c128: {} levels of strassen, {} pol.muls\n",
                   depth(oc, m, b, b), nmults(oc, m, b, b));
        fmt::print("fake: {} levels of strassen, {} pol.muls\n",
                   depth(of, m, b, b), nmults(of, m, b, b));
        fmt::print("gf2x: {} levels of strassen, {} pol.muls\n",
                   depth(os, m, b, b), nmults(os, m, b, b));
    }
    /* Tune for pi * pi */
    {
        ASSERT_ALWAYS(b && m && n);
        size_t const d = (N * b / m / n);

        unsigned long const dl = d - d / 2;
        unsigned long const pi_l_len = dl * m / b; // always <= dl
        size_t const n1 = pi_l_len;
        size_t const n2 = pi_l_len;

        const gf2x_cantor_fft_info oc(n1, n2);
        fmt::print("=== c128 ===\n");
        tune_strassen1(oc, b, b, b, n1 + n2 - 1, state);

        const gf2x_fake_fft_info of(n1, n2);
        fmt::print("=== fake ===\n");
        tune_strassen1(of, b, b, b, n1 + n2 - 1, state);

        gf2x_ternary_fft_info os(n1, n2);
        os.adjust(GF2X_FFT_ADJUST_DEPTH, 81);
        fmt::print("=== tfft({}) ===\n", 81);
        tune_strassen1(os, b, b, b, n1 + n2 - 1, state);

        fmt::print("Options for composition at top level pi*pi\n");
        fmt::print("c128: {} levels of strassen, {} pol.muls\n",
                   depth(oc, b, b, b), nmults(oc, b, b, b));
        fmt::print("fake: {} levels of strassen, {} pol.muls\n",
                   depth(of, b, b, b), nmults(of, b, b, b));
        fmt::print("gf2x: {} levels of strassen, {} pol.muls\n",
                   depth(os, b, b, b), nmults(os, b, b, b));
    }
    foo<gf2x_cantor_fft_info>::s.dump("C128");
    foo<gf2x_fake_fft_info>::s.dump("FAKE");
    foo<gf2x_ternary_fft_info>::s.dump("gf2x_ternary_fft");
#endif
}

static void do_polmm_timings(unsigned long m, unsigned long n, unsigned long N, cxx_gmp_randstate & state)
{
    unsigned long const b = m + n;

    fmt::print("Some timings are in microseconds\n");
    fmt::print("Matrix size {}\n", N);
    // unsigned long c = m - n;

    std::vector<level_info> results;

    for (unsigned long level = 0;; level++) {
        unsigned long const d = (N * b / m / n) >> level;
        if (d <= 64) {
            fmt::print("[{}] is below threshold (length {})\n", level, d);
            break;
        }

        if (DATA_POOL_SIZE * ULONG_BITS / 2 < d) {
            fmt::print(stderr, "Please increa DATA_POOL_SIZE to at least {}\n",
                       W(2 * d));
            exit(1);
        }

        fmt::print("[{}] length of E is {} ({} MB)\n", level, d,
                   d * m * b >> 23);
        unsigned long const dl = d - d / 2;
        unsigned long const pi_l_len = dl * m / b; // always <= dl
        fmt::print("[{}] Top-level length of pi_left is {} ({} MB)\n", level,
                   pi_l_len, pi_l_len * b * b >> 23);
        unsigned long const chop = dl - pi_l_len;
        fmt::print("[{}] Number of chopped of bits at top level is {}\n",
                   level, chop);
        fmt::print("[{}] Top-level degree of truncated E is {}\n", level,
                   d - chop);
        fmt::print("[{}] Degree of product E'*pi_left is {}\n", level,
                   d - chop + pi_l_len);
        fmt::print("[{}] Degree of product pi_left*pi_right is {}\n", level,
                   d * m / b);

        // bench_one_polmm_projected(m, b, b, d-chop, pi_l_len, state);
        const level_info r =
            bench_one_polmm_projected(b, b, b, pi_l_len, pi_l_len, state);
        results.push_back(r);
        // bench_one_polmm_complete(m, b, b, d-chop, pi_l_len);
        bench_one_polmm_complete(b, b, b, pi_l_len, pi_l_len, state);

#if 0
        /* Start by benching degree d polynomial multiplication */
        bench_polmul(d, d);

        /* Then time for cantor, separating the different steps */
        bench_c128(d, d);

        /* Now the matrix versions */
        bench_polmatmul<gf2x_cantor_fft_info>("gf2x_cantor_fft", n, n, d, d);
        bench_polmatmul<gf2x_fake_fft_info>("gf2x_fake_fft", n, n, d, d);
        bench_polmatmul<gf2x_ternary_fft_info>("gf2x_ternary_fft", n, n, d, d);

        /* Now do these again, but for unbalanced computations */

        bench_polmul(d, d/4);
        bench_c128(d, d/4);
        bench_polmatmul<gf2x_cantor_fft_info>("gf2x_cantor_fft", n/2, n, d, d/4);
        bench_polmatmul<gf2x_fake_fft_info>("gf2x_fake_fft", n/2, n, d, d/4);
        bench_polmatmul<gf2x_ternary_fft_info>("gf2x_ternary_fft", n/2, n, d, d/4);
#endif
    }

    fmt::print("Summary for {} levels\n", results.size());
    double total = 0;
    for (unsigned int i = 0; i < results.size(); i++) {
        double loc =
            (results[i].dft + results[i].compose + results[i].ift) * (1 << i);
        total += loc;
        fmt::print(
            "[{}] {} (d{}) {:.2f}+{:.2f}+{:.2f} --> {:.2f}, total {:.2f}\n", i,
            results[i].engine, results[i].nstrassen, results[i].dft,
            results[i].compose, results[i].ift, loc, total);
    }
}

int main(int argc, char const * argv[])
{
    /* This is a bench program. We want the output quick.  */
    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    unsigned long m = 0, n = 0, N = 0;
    argc--, argv++;

    for (; argc;) {
        if (strcmp(argv[0], "--nrep") == 0) {
            argc--, argv++;
            if (!argc)
                usage();
            char * p;
            nrep = strtol(argv[0], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argc--, argv++;
            continue;
        }
        if (N == 0) {
            char * p;
            N = strtol(argv[0], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argc--, argv++;
            continue;
        }
        if (m == 0) {
            char * p;
            m = strtol(argv[0], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argc--, argv++;
            continue;
        }
        if (n == 0) {
            char * p;
            n = strtol(argv[0], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argc--, argv++;
            continue;
        }
        usage();
    }
    if (!N || !m || !n) usage();

    unsigned long const b = m + n;

#if 0
    /* tune strassen for all sizes... Takes a looong time */
#define MAX_CONSIDERED_THRESHOLD 1000000
    {
        size_t d = (N * b / m / n);
        size_t n1 = d;
        size_t n2 = d * m / b;
        gf2x_cantor_fft_info oc(n1, n2); fmt::print("c128");
        tune_strassen(oc, MAX_CONSIDERED_THRESHOLD, state);
        gf2x_fake_fft_info of(n1, n2); fmt::print("fake");
        tune_strassen(of, MAX_CONSIDERED_THRESHOLD, state);
        gf2x_ternary_fft_info os(n1, n2, 81); fmt::print("tfft({})", 81);
        tune_strassen(of, MAX_CONSIDERED_THRESHOLD, state);
    }
#endif

    foo<gf2x_cantor_fft_info>::s.threshold(b, b, b) = UINT_MAX;
    foo<gf2x_fake_fft_info>::s.threshold(b, b, b) = UINT_MAX;
    foo<gf2x_ternary_fft_info>::s.threshold(b, b, b) = UINT_MAX;

    cxx_gmp_randstate state;
    gmp_randseed_ui(state, time(nullptr));

    do_polmm_timings(m, n, N, state);

    tune_strassen_global(m, n, N, state);
    do_polmm_timings(m, n, N, state);

    return 0;
}
