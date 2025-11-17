/* program to determine the size of the L1 cache */
#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>

#include <array>
#include <limits>
#include <vector>
#include <memory>
#include <algorithm>

#include "timing.h"     // microseconds

/* we assume the L1 cache has size 2^k */
static constexpr size_t spin_max = 100000000;

extern "C" {
size_t cachesize_guess(int verbose);
}

size_t
cachesize_guess (int verbose)
{
    std::vector<std::pair<int, uint64_t>> results;

    uint64_t mintime = UINT64_MAX;

    for (int logsize = 10 ; logsize <= 20 ; logsize++) {
        const size_t n = 1 << logsize;
        const size_t mask = n - 1;
        auto s0 = std::make_unique<char[]>(2 * n);
        /* align s on a multiple of n */
        char * s = s0.get() + (n - ((uintptr_t) s0.get() & mask));
        for (size_t j = 0; j < n; j++)
            s[j] = 0;
        /* We compute k(j) = 2*j + 3*j^2 mod n.
           We have k(j+1) - k(j) = 6*j + 5. */
        for (size_t j = 0, k = 0, q = 5; j < spin_max / 10 ; j++) {
            /* invariant: q = 6j+5 */
            s[k] ++;
            k = (k + q) & mask;
            q += 6;
        }
        uint64_t t = microseconds ();
        for (size_t j = 0, k = 0, q = 5; j < spin_max ; j++) {
            /* invariant: q = 6j+5 */
            s[k] ++;
            k = (k + q) & mask;
            q += 6;
        }
        t = microseconds () - t;
        results.emplace_back(logsize, t);
        mintime = std::min(t, mintime);
        fprintf (stderr, "size=%zu time=%" PRIu64 "\n", n, t);
    }

    /* for all x, record y = time(x)/mintime as well as x^2 and x*y */
    std::vector<std::array<double, 4>> results_scaled;
    for(auto [ i, t ] : results) {
        const double x = i;
        const double y = double(t) / double(mintime);
        const double xx = x * x;
        const double xy = x * y;
        results_scaled.push_back({x, y, xx, xy});

        /*
        const size_t n = 1 << i;
        fprintf (stderr, "size=%lu time=%" PRIu64 "(%1.2f)\n", n, t,
                (double) t / (double) mintime);
                */
    }

    size_t best_i = 0;
    double best_distance = std::numeric_limits<double>::max();

    /*
    for (size_t i = 0 ; i < results.size() ; i++) {
        auto [ logsize, t ] = results[i];
        const auto [ x, y, xx, xy ] = results_scaled[i];
        printf("%d %.3f\n", logsize, y);
    }
    */

    std::vector<double> distances;

    for (size_t i = 0 ; i < results.size() ; i++) {
        /* compute a least-squares-fit of the values below i and above i
         * separately. This is often off by a factor of two, but it isn't
         * as bad as the previous heuristics that went easily off-track
         * with a 16x error. */
        double distance = 0;
        if (i > 0) {
            /* left */
            const size_t N = i + 1;
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for(size_t j = 0 ; j <= i ; j++) {
                const auto [ x, y, xx, xy ] = results_scaled[j];
                sx += x;
                sy += y;
                sxx += xx;
                sxy += xy;
            }
            const double m = (N * sxy - sx * sy) / (N * sxx - sx * sx);
            const double b = (sy - m * sx) / N;
            for(size_t j = 0 ; j <= i ; j++) {
                const auto [ x, y, xx, xy ] = results_scaled[j];
                const double eps = m * x + b - y;
                distance += eps * eps;
            }
        }
        if (i + 1 < results.size()) {
            /* right */
            const size_t N = results.size() - i;
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for(size_t j = i ; j < results.size() ; j++) {
                const auto [ x, y, xx, xy ] = results_scaled[j];
                sx += x;
                sy += y;
                sxx += xx;
                sxy += xy;
            }
            const double m = (N * sxy - sx * sy) / (N * sxx - sx * sx);
            const double b = (sy - m * sx) / N;
            for(size_t j = i ; j < results.size() ; j++) {
                const auto [ x, y, xx, xy ] = results_scaled[j];
                const double eps = m * x + b - y;
                distance += eps * eps;
            }
        }
        // auto [ logsize, t ] = results[i];
        // printf("%d %.4f\n", logsize, distance);
        if (distance < best_distance) {
            best_i = i;
            best_distance = distance;
        }
        distances.push_back(distance);
    }

#if 0
    /* see what is the best time we have *around* best_i. It is quite
     * easy for this test to end up being off by one, but then the
     * closest min is probably a good guess.
     * (well, results are a bit disappointing)
     */
    size_t i0 = std::max(best_i, size_t(1)) - 1;
    size_t i1 = std::min(best_i + 1, results.size() - 1);

    mintime = UINT64_MAX;
    for(size_t i = i0 ; i <= i1 ; i++) {
        if (distances[i] >= best_distance + best_distance / 2)
            continue;
        auto [ logsize, t ] = results[i];
        // fprintf (stderr, "size=%zu time=%" PRIu64 "\n", (size_t) (1 << logsize) , t);
        if (t < mintime) {
            best_i = i;
            mintime = t;
        }
    }
#endif

    const auto [ logsize, t ] = results[best_i];
    const size_t n = 1 << logsize;
    if (verbose) {
        printf ("#define L1_CACHE_SIZE %zu\n", n);
    }
    return n;
}
