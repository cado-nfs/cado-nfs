#include "cado.h" // IWYU pragma: keep

/* The bucket_hash<> class below is a c++ implementation of the exact
 * same algorithm as in polyselect_shash.c. (see comments there for the
 * description of the algorithm)
 *
 * Speed is like 5-10% slower in general compared to the legacy C code,
 * and the slowdown seems pretty well spread between the push and the
 * find_collision steps. I don't have an obvious explanation for this
 * offhand.
 *
 * The prefetching mechanism that is used in the original C code seems to
 * have only very small effect.
 *
 * Performance varies quite dramatically from a machine to another.
 */

#define EMIT_ADDRESSABLE_shash_add
#define EXPOSE_DEPRECATED_polyselect_shash_find_collision

#include <ctime>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "polyselect_shash.h"
#include "misc.h"
#include "utils_cxx.hpp"
#include "bucket_hash.hpp"

static void onetest(int64_t umax, size_t pushed_entries, unsigned long seed, int nruns)
{
    const size_t expected_entries = pushed_entries;
    cxx_gmp_randstate rstate;

    fmt::print("test with umax={} pushed_entries={} seed={} nruns={}\n",
            umax, pushed_entries, seed, nruns);

    std::vector<size_t> collisions_c_code;
    {
        gmp_randseed_ui(rstate, seed);

        polyselect_shash_t H;
        polyselect_shash_init(H, expected_entries);

        std::vector<size_t> found;
        clock_t st = clock();
        clock_t tfill = 0, tsearch = 0;
        for(int i = 0 ; i < nruns ; i++) {
            polyselect_shash_reset(H);
            for(size_t i = 0 ; i < pushed_entries ; i++) {
                auto const v = i64_random(rstate) % (2*umax) - umax;
                polyselect_shash_add(H, v);
            }
            tfill += clock() - st; st = clock();
            if (polyselect_shash_find_collision(H))
                found.push_back(i);
            tsearch += clock() - st; st = clock();
        }
        fmt::print("{} {} {:.2f} {:.2f} {:.2f} {}\n",
                nruns, found.size(),
                (double) (tfill + tsearch) / CLOCKS_PER_SEC,
                (double) (tfill) / CLOCKS_PER_SEC,
                (double) (tsearch) / CLOCKS_PER_SEC,
                join(found, " "));

        polyselect_shash_clear(H);

        collisions_c_code = found;
    }

    std::vector<size_t> collisions_cxx_code;
    {
        gmp_randseed_ui(rstate, seed);
        cado::bucket_hash<cado::polyselect_shash_config> H(expected_entries);

        std::vector<size_t> found;
        clock_t st = clock();
        clock_t tfill = 0, tsearch = 0;
        for(int i = 0 ; i < nruns ; i++) {
            H.reset();
            for(size_t i = 0 ; i < pushed_entries ; i++) {
                try {
                    auto const v = i64_random(rstate) % (2*umax) - umax;
                    H.push(v);
                } catch(std::runtime_error const& e) {
                    fprintf(stderr, "overflow after %zu-th push\n", i);
                    break;
                }
            }
            tfill += clock() - st; st = clock();
            if (H.has_collision())
                found.push_back(i);
            tsearch += clock() - st; st = clock();
        }
        fmt::print("{} {} {:.2f} {:.2f} {:.2f} {}\n",
                nruns, found.size(),
                (double) (tfill + tsearch) / CLOCKS_PER_SEC,
                (double) (tfill) / CLOCKS_PER_SEC,
                (double) (tsearch) / CLOCKS_PER_SEC,
                join(found, " "));

        collisions_cxx_code = found;
    }
    if (collisions_c_code != collisions_cxx_code) {
        fprintf(stderr, "The two implementations don't give matching results\n");
        exit(EXIT_FAILURE);
    }

    /* This code is functionally equivalent, but 6 to 7 times slower on
     * my laptop. */
    if (false) {
        gmp_randseed_ui(rstate, seed);
        std::vector<int64_t> H;

        std::vector<size_t> found;
        clock_t st = clock();
        clock_t tfill = 0, tsearch = 0;
        for(int i = 0 ; i < nruns ; i++) {
            H.clear();
            for(size_t i = 0 ; i < pushed_entries ; i++) {
                auto const v = i64_random(rstate) % (2*umax) - umax;
                H.push_back(v);
            }
            tfill += clock() - st; st = clock();
            std::ranges::sort(H);
            for(auto it = ++H.begin() ; it != H.end() ; ++it) {
                if (it[0] == it[-1]) {
                    found.push_back(i);
                    break;
                }
            }
            tsearch += clock() - st; st = clock();
            // printf("%d %d\n", i, found);
        }
        fmt::print("{} {} {:.2f} {:.2f} {:.2f} {}\n",
                nruns, found.size(),
                (double) (tfill + tsearch) / CLOCKS_PER_SEC,
                (double) (tfill) / CLOCKS_PER_SEC,
                (double) (tsearch) / CLOCKS_PER_SEC,
                join(found, " "));

        if (found != collisions_c_code) {
            fmt::print(stderr, "slow reference code disagrees with other implementations\n");
            // exit(EXIT_FAILURE);
        }
    }

}

int main()
{
    /* the number of entries that go in the hash table is pushed_entries,
     * and the value range is [-umax, umax].
     * So the expected number of collisions is pushed_entries^2 / 2 / 2umax.
     * With pushed_entries=1e6 and umax=2e12, this means an
     * expected value of 1/8. When repeated 100 times, this should mean
     * about 12 test with collisions.
     */
    {
        polyselect_shash_t H;
        polyselect_shash_init(H, 1e6);

        /* this triggers a collision, although it is a fake one. Because
         * the ordering of the inserts in the legacy C code is not
         * consistent, what we're seeing this is actually the third value
         * being inserted first, and the first one being inserted last,
         * yielding a collision.
         *
         * Getting a collision when the value space exceeds 2^40 (which
         * is the case here) isn't much of a surprise, though!
         */
        polyselect_shash_add(H, INT64_C(0x47087285be));
        polyselect_shash_add(H, INT64_C(-0xb9f78d7942));
        polyselect_shash_add(H, INT64_C(-0x4dde60d7a42));
        if (polyselect_shash_find_collision(H))
            fmt::print("corner case check: collision found\n");
        else
            fmt::print("corner case check: no collision found\n");
        polyselect_shash_clear(H);
    }

    onetest(2e6, 2e3, 0xdeadbeef, 10);
    onetest(2e7, 1e4, 0xdeadbeef, 10);
    onetest(1e8, 1e4, 0xdeadbeef, 10);
    onetest(2e4, 1e2, 1, 10);
    onetest(2e6, 1e3, 1, 10);
    onetest(2e8, 1e4, 1, 10);
    onetest(2e10, 1e5, 1, 10);
    onetest(2e12, 1e6, 1, 100);

    return 0;
}

