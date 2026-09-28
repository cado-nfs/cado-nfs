#include "cado.h" // IWYU pragma: keep

/* Tuner for the compile-time choice search_survivors_uses_patterns, which is
 * made in sieve/las-unsieve-sse2.cpp.
 *
 * Besides the generic variant, the survivor search has two variants that use
 * an SSE bound pattern to kill in bulk the positions whose real abscissa is a
 * multiple of 3 (resp. of 5), on the rows where that matters. They are a loss
 * on both the microarchitectures we have measured, but the margin is a
 * property of the target, and the break-even survivor density is not far
 * enough away to take on trust.
 *
 * This program times the very functions that las runs, with the choice forced
 * both ways, and compares the outcome with the compile-time default. It exits
 * with a non-zero status when that default is off by more than -margin, so
 * that a port to a new microarchitecture says so rather than leaving us to
 * guess. It also checks that the two variants report the same survivors.
 *
 * The stake is small: the survivor search is well under one percent of a las
 * run. This is a tuner, not a benchmark of anything that matters on its own.
 */

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <vector>

#include "las-sublat.hpp"
#include "las-unsieve.hpp"
#include "macros.h"
#include "params.hpp"

/* A bucket region is 2^LOG_BUCKET_REGION bytes per side, and that is the
 * working set the survivor search really sees. We keep the buffer to that
 * size, and get the variety of row indices we need by sweeping j across calls
 * rather than across lines in memory: rows are recycled modulo the number of
 * lines the buffer holds.
 */
static constexpr int log_region = 16;

/* How many distinct rows the sweep visits. Any multiple of 15 gives the
 * pattern-3 and pattern-5 rows their right share.
 */
static constexpr unsigned int nrows = 240;

struct rng {
    uint64_t s;
    explicit rng(uint64_t seed) : s(seed ? seed : 1) {}
    uint64_t operator()() {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        return s;
    }
};

struct bench {
    unsigned int nlines;        /* lines held in the buffer */
    size_t linelen;             /* bytes per line */
    int i0, i1;
    unsigned char bound[2] = { 100, 100 };
    unsigned int td_max = 100;
    std::vector<unsigned char> S[2];
    /* offsets of the planted low values, one list per line of the buffer */
    std::vector<std::vector<size_t>> planted;
    j_divisibility_helper j_div;

    static constexpr unsigned char high = 200;  /* above the bound */
    static constexpr unsigned char low = 50;    /* below the bound */

    bench(int logI, double density, uint64_t seed, uint32_t jmax)
        : nlines(1U << (log_region > logI ? log_region - logI : 0))
        , linelen(size_t(1) << (logI < log_region ? logI : log_region))
        , i0(-(1 << (logI - 1)))
        , i1(i0 + int(linelen))
        , planted(nlines)
        , j_div(jmax)
    {
        /* Both sides are filled well above the bound, then positions that
         * pass the bound test are planted at the requested density. That
         * density is what decides how often the scan leaves its inner loop,
         * so it matters much more than the values themselves; 1e-4 is what
         * las reports on the inputs we have looked at (survivor ratios from
         * 3e-5 to 2.5e-4 across F9 and RSA-768).
         */
        rng r(seed);
        size_t const n = size_t(nlines) * linelen;
        for (auto & side : S) side.assign(n, high);
        for (size_t k = size_t(double(n) * density); k--; ) {
            size_t const x = r() % n;
            S[0][x] = S[1][x] = low;
            planted[x / linelen].push_back(x);
        }
    }

    /* The search raises to 255 the positions that it rejects, among them the
     * ones that trial division finds non-coprime -- and those are exactly the
     * ones the bound patterns exist to dispose of cheaply. Scanning the same
     * line twice would therefore measure the patterns' cost without their
     * benefit, so we put the planted values back before each scan. There are
     * only a handful of them, so this costs nothing worth subtracting. The
     * positions that were merely raised from `high` to 255 can be left alone:
     * they fail the bound test either way.
     */
    void refresh(unsigned int j) {
        for (size_t const x : planted[j % nlines])
            S[0][x] = S[1][x] = low;
    }

    unsigned char * line(int side, unsigned int j) {
        return S[side].data() + (j % nlines) * linelen;
    }
};

template<bool P, bool oneside>
static void one_line(bench & B, unsigned int j, sublat_runtime_t sublat,
        std::vector<uint32_t> & sv)
{
    if constexpr (oneside) {
        search_survivors_in_line_sse2_oneside_choice<P>(B.line(0, j),
                B.bound[0], j, B.i0, B.i1, 0, B.j_div, B.td_max, sv, sublat);
    } else {
        unsigned char * const SS[2] = { B.line(0, j), B.line(1, j) };
        search_survivors_in_line_sse2_choice<P>(SS, B.bound, j, B.i0, B.i1,
                0, B.j_div, B.td_max, sv, sublat);
    }
}

/* Returns nanoseconds per byte scanned, best of nrounds. */
template<bool P, bool oneside>
static double measure(bench & B, sublat_runtime_t sublat, int nrounds,
        double target)
{
    std::vector<uint32_t> sv;
    sv.reserve(1 << 16);

    auto pass = [&](int nrep) {
        clock_t const tt = clock();
        for (int i = 0; i < nrep; i++) {
            sv.clear();
            for (unsigned int j = 1; j <= nrows; j++) {
                B.refresh(j);
                one_line<P, oneside>(B, j, sublat, sv);
            }
        }
        return double(clock() - tt) / CLOCKS_PER_SEC;
    };

    /* warm up, and pick a repeat count that brings one round to about
     * -time seconds */
    int nrep = 1;
    double t;
    while ((t = pass(nrep)) < target / 8 && nrep < (1 << 20))
        nrep *= 8;
    nrep = int(double(nrep) * target / (t > 0 ? t : target)) + 1;

    double best = -1;
    for (int round = 0; round < nrounds; round++) {
        t = pass(nrep) / (double(nrep) * double(nrows) * double(B.linelen));
        if (best < 0 || t < best) best = t;
    }
    return best * 1e9;
}

/* The two variants must agree on the survivors: the positions that the
 * pattern kills are exactly the ones the generic variant rejects by trial
 * division. The planted values are put back around each scan, so a divergence
 * on one row does not cascade into the next.
 */
template<bool oneside>
static bool check_same(bench & B, sublat_runtime_t sublat)
{
    std::vector<uint32_t> a, b;
    for (unsigned int j = 1; j <= nrows; j++) {
        B.refresh(j);
        size_t const na = a.size();
        one_line<false, oneside>(B, j, sublat, a);
        B.refresh(j);
        size_t const nb = b.size();
        one_line<true, oneside>(B, j, sublat, b);
        if (a.size() - na == b.size() - nb
                && memcmp(a.data() + na, b.data() + nb,
                    (a.size() - na) * sizeof(uint32_t)) == 0)
            continue;
        fprintf(stderr, "# the two variants disagree on row j=%u"
                " (sublat m=%u, i0=%u, j0=%u, %s):"
                " %zu survivors without the patterns, %zu with\n",
                j, sublat.m, sublat.i0, sublat.j0,
                oneside ? "one side" : "two sides",
                a.size() - na, b.size() - nb);
        return false;
    }
    return true;
}

struct outcome {
    char const * what;
    double t_off, t_on;
};

int main(int argc, char const * argv[])
{
    int logI = 15;
    double density = 1e-4;
    double margin = 0.05;
    double target = 0.1;
    int nrounds = 5;
    int report_only = 0;
    unsigned long seed = 42;

    cxx_param_list pl;
    pl.declare_usage("I", "log of the line length (default 15)");
    pl.declare_usage("density", "survivor density (default 1e-4)");
    pl.declare_usage("margin", "tolerated gap to the best choice (default .05)");
    pl.declare_usage("rounds", "measurement rounds, best is kept (default 5)");
    pl.declare_usage("time", "seconds per measurement round (default .1)");
    pl.declare_usage("seed", "random seed");
    pl.declare_usage("report-only", "report, but exit successfully anyway");
    pl.configure_switch("-report-only");
    pl.process_command_line(argc, argv, false);
    pl.parse("I", logI);
    pl.parse("density", density);
    pl.parse("margin", margin);
    pl.parse("rounds", nrounds);
    pl.parse("time", target);
    pl.parse("seed", seed);
    pl.parse("-report-only", report_only);

    ASSERT_ALWAYS(logI >= 5 && logI <= log_region);

    /* j_div is indexed by the real row jj = m*j + j0, and its size must be a
     * power of two. */
    uint32_t jmax = 2;
    while (jmax <= 6 * nrows) jmax *= 2;

    bench B(logI, density, seed, jmax);

    printf("# survivor search: bound patterns for the primes 3 and 5\n");
    printf("# logI=%d, %u lines of %zu bytes, %u rows, density %.1e,"
            " best of %d rounds\n",
            logI, B.nlines, B.linelen, nrows, density, nrounds);

    static constexpr sublat_runtime_t plain { 1, 0, 0 };
    static constexpr sublat_runtime_t sublat2 { 2, 1, 0 };

    int rc = 0;
    if (!check_same<false>(B, plain)) rc = 1;
    if (!check_same<false>(B, sublat2)) rc = 1;
    if (!check_same<true>(B, plain)) rc = 1;
    if (rc) {
        fprintf(stderr, "# the two variants are not equivalent,"
                " there is no point in timing them\n");
        return rc;
    }

    outcome const res[] = {
        { "two sides, no sublattice",
            measure<false, false>(B, plain, nrounds, target),
            measure<true, false>(B, plain, nrounds, target) },
        { "two sides, sublat=2",
            measure<false, false>(B, sublat2, nrounds, target),
            measure<true, false>(B, sublat2, nrounds, target) },
        { "one side, no sublattice",
            measure<false, true>(B, plain, nrounds, target),
            measure<true, true>(B, plain, nrounds, target) },
    };

    printf("#\n# %-24s %10s %10s   %s\n", "configuration",
            "generic", "pattern", "best");
    for (auto const & o : res) {
        bool const on_is_best = o.t_on < o.t_off;
        double const gap = on_is_best
            ? o.t_off / o.t_on - 1 : o.t_on / o.t_off - 1;
        printf("# %-24s %10.4f %10.4f   %-7s by %5.1f%%%s\n", o.what,
                o.t_off, o.t_on, on_is_best ? "pattern" : "generic",
                100 * gap,
                on_is_best == search_survivors_uses_patterns ? ""
                : (gap > margin ? "   <== mistuned" : "   (within margin)"));
        if (on_is_best != search_survivors_uses_patterns && gap > margin
                && &o != &res[2])
            rc = 1;
    }
    printf("# (nanoseconds per byte scanned, lower is better)\n");

    /* A sanity check on the measurement itself. The one-sided search does
     * strictly less work than the two-sided one, so it cannot be slower; when
     * it is, this binary's code placement is deciding the outcome rather than
     * the loops, and nothing here is worth reading. That is what an unaligned
     * build did on Skylake-SP, which is why the tuner asks for
     * -falign-loops=64.
     */
    if (res[2].t_off > res[0].t_off) {
        printf("#\n# The one-sided search came out slower than the two-sided"
                " one, which it cannot\n# be. This binary is measuring its own"
                " code placement. Refusing to conclude.\n");
        return report_only ? 0 : 1;
    }

    /* Two sides is what las runs for ordinary factoring, so it decides; the
     * one-sided figure is there to be looked at, not to vote. */
    bool const want = res[0].t_on < res[0].t_off;
    if ((res[1].t_on < res[1].t_off) != want) {
        printf("#\n# The two-sided configurations disagree with each other,"
                " so the effect is too\n# small to act on here. Leaving the"
                " default alone.\n");
        printf("UNSIEVE_PATTERNS=%d\n", int(search_survivors_uses_patterns));
        return 0;
    }
    if (rc) {
        printf("# search_survivors_uses_patterns is %s in"
                " sieve/las-unsieve-sse2.cpp, and the\n"
                "# measurements above say it should be %s on this"
                " microarchitecture. Rebuild\n"
                "# with -DUNSIEVE_PATTERNS=%d to act on that, and please"
                " report the outcome:\n"
                "# the compile-time default wants a case for this target.\n",
                search_survivors_uses_patterns ? "true" : "false",
                search_survivors_uses_patterns ? "false" : "true",
                int(want));
    }
    printf("UNSIEVE_PATTERNS=%d\n", int(want));
    return report_only ? 0 : rc;
}
