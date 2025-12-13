#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <algorithm>
#include <iostream>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "ecm.h"
#include "macros.h"
#include "smooth_detect.hpp"


// For a bug in ecm ?
static double default_B1done;

#define EA_THRESHOLD 0.8 // Heuristic for early-abort: low = keep many

// Expected effort to extract of prime of p bits.
// This has been computed for B1min = 200.0, but this is rather stable
// for interesting values of p.
#define MAX_PBIT 100
static const double expected_effort[MAX_PBIT] = {
    200,      200,      200,       200,      200,
    200,      200,      200,       200,      200, // 0 - 9
    200,      200,      200,       200,      200,
    200,      200,      211,       219,      243, // 10 - 19
    296,      280,      351,       380,      482,
    588,      653,      719,       834,      985, // 20 - 29
    1280,     1489,     2352,      2324,     2714,
    3813,     5141,     5891,      7986,     8816, // 30 - 39
    10087,    12264,    14191,     18827,    25633,
    25491,    37803,    39392,     44290,    51925,    // 40-49
    87203,    79943,    110007,    121644,   147602,   // 50 - 54
    174995,   199245,   257190,    279228,   345960,   // 55 - 59
    351682,   530703,   640140,    759310,   775311,   // 60 - 64
    960249,   1267879,  1174122,   1272107,  1589907,  // 65 - 69
    2258437,  2004235,  2903959,   3002629,  3888904,  // 70 - 74
    4373729,  4899345,  5218152,   6269843,  7063446,  // 75 - 79
    9553542,  9891138,  10623352,  13795248, 17574109, // 80 - 84
    18790448, 23529670, 24757303,  30897420, 31188626, // 85 - 89
    36647830, 41007546, 45692079,  53881307, 69709984, // 90 - 94
    81857570, 84194044, 107900749, 94609433, 136660173 // 95 - 99
};

static double
get_time()
{
    return (double)(clock()) / CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////
// Candidate: structure to store a pair of numbers currently being
// factored.
///////////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& os, descent_init_candidate const& c)
{
    cxx_mpz const & u0 = c.u0;
    cxx_mpz const & v0 = c.v0;
    cxx_mpz const & u = c.u;
    cxx_mpz const & v = c.v;
    return os
           << fmt::format("Candidate e = {}\n", c.e)
           << fmt::format("u0={}\nv0={}\n", u0, v0)
           << fmt::format(
                "u={} ({} bits) largest prime so far of {} bits)\n",
                u,
                mpz_sizeinbase(c.u, 2),
                c.lpu)
           << fmt::format(
                "v={} ({} bits) largest prime so far of {} bits)\n",
                v,
                mpz_sizeinbase(c.v, 2),
                c.lpv)
           << fmt::format("effort={:.0f}\n", c.effort);
}

// if effort is such that all primes up to b-bits have been removed,
// and if an unfactored part has less than 3*b bits, then there are only
// two factors. If size is more than 2*bound, it can not be smooth.
// Heuristic: we consider that if the effort is twice the average for
// detecting p bits, then all p bits prime have been removed.
bool
descent_init_candidate::is_probably_not_smooth(unsigned int bound) const
{
    double const eff = effort;
    unsigned int bits = 0;
    while ((bits < MAX_PBIT) && (2 * expected_effort[bits] < eff)) {
        bits++;
    }
    if (bits == MAX_PBIT) {
        return false; // can not conclude; no data for this effort
    }
    // bits -= 3; // take some margin

    for (unsigned int k = 2; k < 4; ++k) {
        unsigned long const bu = mpz_sizeinbase(u, 2);
        if ((bu < (k + 1) * bits) && (bu > k * bound)) {
            //      fmt::print("Probably not smooth, level {}: {} bits!\n", k, bu);
            return true;
        }
        unsigned long const bv = mpz_sizeinbase(v, 2);
        if ((bv < (k + 1) * bits) && (bv > k * bound)) {
            //      fmt::print("Probably not smooth, level {}: {} bits!\n", k, bv);
            return true;
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////
// Pool: contains a list of candidates that are still interesting to
// try to factor
/////////////////////////////////////////////////////////////////////////////

// The list of candidates must always be sorted

// Remove candidates that have a prime with bitsize more than lmax,
// and those that
// are fully factored. If there are still more than max_size elements,
// keep those with smallest cost.

static void
purge(std::vector<descent_init_candidate>& P,
      unsigned int max_size,
      unsigned int lmax)
{
    ASSERT(std::ranges::is_sorted(P));
    auto w = P.begin();

    for (auto r = P.begin(); r != P.end(); ++r) {
        if (!r->is_factored()
                && !r->is_probably_not_smooth(lmax)
                && r->maxlp() <= lmax)
        {
            std::swap(*r, *w++);
            if ((w - P.begin()) > max_size)
                break;
        }
    }
}

static std::ostream&
operator<<(std::ostream& os, std::vector<descent_init_candidate> const& P)
{
    size_t i = 0;
    for (auto const& C : P) {
        os << fmt::format("{}: {} {} ({})\n",
                          i,
                          mpz_sizeinbase(C.u, 2),
                          mpz_sizeinbase(C.v, 2),
                          C.maxlp());
        i++;
    }
    return os;
}

///////////////////////////////////////////////////////////////////////////
// Stats: keep track of statistics about the exepected number of bits
// obtained after running x curves.
///////////////////////////////////////////////////////////////////////////////

class ecm_stats
{
#define MAX_CPT 2048
    double aver_gain[MAX_CPT] {
        0,
    }; // average number of bits removed after i curves
    unsigned long nb_test[MAX_CPT] {
        0,
    }; // size of the sample on which this avearge was done
    public:
    void update(double gain, unsigned int i)
    {
        if (i >= MAX_CPT) {
            abort();
        }
        double newav = aver_gain[i] * double(nb_test[i]) + gain;
        nb_test[i]++;
        newav /= double(nb_test[i]);
        aver_gain[i] = newav;
    }

    bool is_below_average(double gain, unsigned int i) const
    {
        if (nb_test[i] <= 100)
            return false; // not enough data to conclude
        else
            return gain < EA_THRESHOLD * aver_gain[i];
    }

    friend std::ostream& operator<<(std::ostream& os, ecm_stats const& s);
};

std::ostream&
operator<<(std::ostream& os, ecm_stats const& s)
{
    os << ("[");
    for (int i = 1; i < MIN(500, MAX_CPT); ++i) {
        if (s.nb_test[i] <= 100)
            break;
        os << fmt::format(" {:.0f}", s.aver_gain[i]);
    }
    os << " ]\n";
    return os;
}

///////////////////////////////////////////////////////////////////////////
// General data structure for an instance of the search of smooth numbers
///////////////////////////////////////////////////////////////////////////

struct context
{
    std::vector<descent_init_candidate> pool;
    ecm_stats stats;
    const void* param_next_cand = nullptr;
    int (*next_cand)(descent_init_candidate&,
                     const void*); // new candidate put in first arg.
    unsigned long target = 0;      // smoothness bound (in bits)
    double current_effort = 0;     // current effort per candidate.
    double max_effort = 0;
    unsigned long max_pool_size = 0;
    double minB1 = 0;

    void increase_effort()
    {
        current_effort += sqrt(current_effort);
        current_effort = MIN(current_effort, max_effort);
    }
};

static double
remove_small_factors(mpz_t z)
{
    double gain = 0.0;
    while (mpz_divisible_ui_p(z, 2)) {
        mpz_divexact_ui(z, z, 2);
        gain += 1.0;
    }
    while (mpz_divisible_ui_p(z, 3)) {
        mpz_divexact_ui(z, z, 3);
        gain += log2(3.0);
    }
    return gain;
}

// get a B1, so that we can quickly cover the target effort
static double
get_B1_from_effort(double effort, double minB1)
{
    double B1 = minB1;
    double S = B1;
    while (S < effort) {
        B1 += sqrt(B1);
        S += B1;
    }
    return B1;
}

static void
my_ecm_factor(cxx_mpz& f, cxx_mpz& z, double B1)
{
    ecm_params ecm_par;
    ecm_init(ecm_par);
    long const sig = random();
    mpz_set_ui(ecm_par->sigma, sig);
    ecm_par->B1done = default_B1done; /* issue with ECM 6.4.x */
    ecm_factor(f, z, B1, ecm_par);
    ecm_clear(ecm_par);
}

// One step of smoothness detection: get a new candidate, run a bunch of
// ECMs, update the pool, and the stats.
// Return value:
//   1: smooth
//   0: non-smooth
//   -1: early stop, no more candidate to test.

static int
smooth_detect_one_step(descent_init_candidate& winner, context& ctx)
{
    descent_init_candidate C;
    // Get a new candidate, not obviously non-smooth
    double gain_u = 0;
    double gain_v = 0;
    do {
        int const ret = ctx.next_cand(C, ctx.param_next_cand);
        if (ret == 0) {
            return -1; // early stop, no more candidates.
        }
        gain_u += remove_small_factors(C.u);
        gain_v += remove_small_factors(C.v);
        C.check_prime();
    } while (C.maxlp() >= ctx.target);

    // Start a loop of ECM
    double B1 = ctx.minB1; // initial B1
    int cpt = 0;           // number of curves tried on this number
    cxx_mpz f;             // output of ecm
    while (!C.is_probably_not_smooth(ctx.target) && !C.is_factored() &&
           C.effort < ctx.current_effort) {
        cpt++;
        // u-side
        if (mpz_cmp_ui(C.u, 1) != 0) {
            my_ecm_factor(f, C.u, B1);
            gain_u += log2(mpz_get_d(f));
            ctx.stats.update(gain_u, cpt);
            if (mpz_cmp_ui(f, 1) > 0) {
                mpz_divexact(C.u, C.u, f);
                C.lpu = MAX(C.lpu, mpz_sizeinbase(f, 2));
                C.check_prime();
                if (C.lpu >= ctx.target)
                    break;
            }
        }
        // v-side
        if (mpz_cmp_ui(C.v, 1) != 0) {
            my_ecm_factor(f, C.v, B1);
            gain_v += log2(mpz_get_d(f));
            ctx.stats.update(gain_v, cpt);
            if (mpz_cmp_ui(f, 1) > 0) {
                mpz_divexact(C.v, C.v, f);
                C.lpv = MAX(C.lpv, mpz_sizeinbase(f, 2));
                C.check_prime();
                if (C.lpv >= ctx.target)
                    break;
            }
        }

        // if both gain are below average, then abort this candidate
        if (ctx.stats.is_below_average(gain_u, cpt) &&
            ctx.stats.is_below_average(gain_v, cpt)) {
            // mark it as unsmooth, by increasing lpu
            C.lpu = UINT_MAX;
            break;
        }
        // remember current effort for this number, and increase B1.
        C.effort += B1;
        B1 += sqrt(B1);
    }

    // If we had a prime larger than expected, then abort.
    if (C.maxlp() > ctx.target) {
        ctx.increase_effort();
        return 0;
    }

    // Did we factor the candidate completely?
    // If so, print result, otherwise, insert in pool
    if (C.is_factored()) {
        std::cout << C;
        std::swap(winner, C);
        return 1;
    } else {
        ctx.pool.push_back(C);
    }

    B1 = get_B1_from_effort(ctx.current_effort, ctx.minB1);

    // Run ECM on the numbers in pool, so that they all have received more
    // or less the same effort.
    //
    // at this point we want to use our understanding of which candidates
    // are best. we'll run a full sort on the remaining candidates
    // eventually, so it seems okay-ish to do it twice.

    std::ranges::sort(ctx.pool);
    size_t i = 0;
    for (auto& c : ctx.pool) {
        double effort = ctx.current_effort;
        // more effort for the most promising candidates!
        // the correcting factors below are completely heuristic...
        if (ctx.pool.size() > 5) {
            if (i == 0) {
                effort *= 2;
            }
            if (i == 1) {
                effort *= 1.6;
            }
            if (i == 2) {
                effort *= 1.3;
            }
            if (i == 3) {
                effort *= 1.1;
            }
        }
        while (!c.is_factored() && c.effort < effort) {
            if (mpz_cmp_ui(c.u, 1) != 0) {
                my_ecm_factor(f, c.u, B1);
                if (mpz_cmp_ui(f, 1) > 0) {
                    mpz_divexact(c.u, c.u, f);
                    c.lpu = MAX(c.lpu, mpz_sizeinbase(f, 2));
                }
            }
            if (mpz_cmp_ui(c.v, 1) != 0) {
                my_ecm_factor(f, c.v, B1);
                if (mpz_cmp_ui(f, 1) > 0) {
                    mpz_divexact(c.v, c.v, f);
                    c.lpv = MAX(c.lpv, mpz_sizeinbase(f, 2));
                }
            }
            c.check_prime();
            c.effort += B1;
            if (c.maxlp() > ctx.target)
                break;
        }
        i++;
    }

    // Go through pool: detect winners, purge bad candidates, keep best.
    for (auto& c : ctx.pool) {
        if (c.is_factored() && c.maxlp() <= ctx.target) {
            std::cout << c;
            std::swap(winner, c);
            return 1;
        }
    }
    std::ranges::sort(ctx.pool);
    purge(ctx.pool, ctx.max_pool_size, ctx.target);

    ctx.increase_effort();
    return 0;
}

descent_init_candidate
smooth_detect(int (*next_cand)(descent_init_candidate&, const void*),
              const void* param_next_cand,
              unsigned long target,
              smooth_detect_params const& param)
{
    /* fix issue with ECM 6.4.x */
    {
        ecm_params params;
        ecm_init(params);
        default_B1done = params->B1done;
        ecm_clear(params);
    }

    /*
    const smooth_detect_param_s default_param = {2000.0, 1e20, 10, 1, 100.0};
    if (param == NULL) {
      param = &default_param;
    }
    */

    // Create a context.
    context ctx;
    ctx.param_next_cand = param_next_cand;
    ctx.next_cand = next_cand;
    ctx.target = target;
    ctx.current_effort = param.min_effort;
    ctx.max_effort = param.max_effort;
    ctx.max_pool_size = param.max_pool_size;
    ctx.minB1 = param.minB1;

    descent_init_candidate C;

    double const tm = get_time();
    int cpt = 0;
    int found = 0;
    while (found == 0) {
        found = smooth_detect_one_step(C, ctx);
        cpt++;
        if (param.verbose && (cpt % 20 == 0)) {
            fmt::print("***** Pool status after {} candidates in {:.1f}s\n",
                   cpt,
                   get_time() - tm);
            fmt::print("current_effort = {:.0f}\n", ctx.current_effort);
            fmt::print("current max B1 = {:.0f}\n",
                   get_B1_from_effort(ctx.current_effort, ctx.minB1));
            fmt::print("current stats:\n");
            std::cout << ctx.stats;
            std::cout << ctx.pool;
        }
    }
    return C;
}
