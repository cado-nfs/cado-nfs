#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <string>
#include <utility>  // for pair, swap
#include <vector>
#include <functional>
#include <algorithm>
#include <sys/time.h>
#include "macros.h"
#include "iqsort.h"

static inline void
renumber_sort_ul (unsigned long *r, size_t n)
{
    unsigned long rmin;

    if (UNLIKELY (n < 2))
        return;

    if (UNLIKELY (n == 2)) {
        if (r[0] < r[1]) {
            rmin = r[0];
            r[0] = r[1];
            r[1] = rmin;
        }
        return;
    }

    for (size_t i = n; --i;) {
        size_t min = i;
        rmin = r[min];
        for (size_t j = i; j--;) {
            unsigned long rj = r[j];
            if (UNLIKELY (rj < rmin)) {
                min = j;
                rmin = rj;
            }
        }
        if (LIKELY (min != i)) {
            r[min] = r[i];
            r[i] = rmin;
        }
    }
}

static inline void
hard_mergesort (unsigned long *r, size_t n)
{
    if (n > 256)
        return std::sort(r, r + n, std::greater<unsigned long>());
    unsigned long temp[256];
    for(size_t i = 0 ; i + 1 < n ; i += 2) {
        if (r[i] < r[i+1])
            std::swap(r[i], r[i+1]);
    }
    /* all sequences of length d are sorted */
    for(size_t d = 2; d < n ; d<<=1) {
        for(size_t i = 0 ; i < n ; i += (d << 1)) {
            size_t i0 = 0;
            size_t i1 = 0;
            size_t j = 0;
            for( ; i0 < d && (i + i0 < n) && i1 < d && (i + d + i1 < n) ; ) {
                if (r[i + i0] < r[i + d + i1]) {
                    temp[j++] = r[i + d + i1++];
                } else {
                    temp[j++] = r[i + i0++];
                }
            }
            for( ; i0 < d && (i + i0 < n) ; ) {
                temp[j++] = r[i + i0++];
            }
            for( ; i1 < d && (i + d + i1 < n) ; ) {
                temp[j++] = r[i + d + i1++];
            }
            std::copy(temp, temp + j, r + i);
        }
    }
}

static int gcmp(const void * pa, const void * pb)
{
    unsigned long a = * (const unsigned long *) pa;
    unsigned long b = * (const unsigned long *) pb;
    return (a < b) - (b < a);
}

int main(int argc, char * argv[])
{
    double exp_growth = 0.1;
    bool verbose = false;
    size_t logsize = 18;
    size_t maxspan = 1 << 9;
    gmp_randstate_t state;

    argv++, argc--;
    for( ; argc ; argc--,argv++) {
        std::string s = *argv;
        if (s == "--verbose" || s == "-v") {
            verbose = true;
            continue;
        }
        if (s == "--exp-growth") {
            exp_growth = atof(argv[1]);
            argc--,argv++;
            continue;
        }
        if (s == "--log-bigsize") {
            logsize = atol(argv[1]);
            argc--,argv++;
            continue;
        }
        if (s == "--max-span") {
            maxspan = atol(argv[1]);
            ASSERT_ALWAYS(maxspan > 0);
#ifdef __COVERITY__
            __coverity_mark_pointee_as_sanitized__(&maxspan, LOOP_BOUND);
#endif
            argc--,argv++;
            continue;
        }
        fprintf(stderr, "Unexpected argument %s\n", s.c_str());
        exit(EXIT_FAILURE);
    }
    gmp_randinit_default (state);

    std::vector<std::pair<unsigned int, std::string>> best;
    std::string last_best = "none";
    std::vector<unsigned long> blah(1 << logsize);
    for(unsigned int v = 2 ; v < maxspan ; v += 1 + v * exp_growth) {
        clock_t t0;
        clock_t d;
        std::vector<std::pair<double, std::string>> record;
        std::string s;

        if (verbose) printf("%u", v);

        s = "std::sort";
        for(auto & x : blah)
            x = gmp_urandomm_ui(state, 1000);
        t0 = clock();

        auto beg = blah.begin();

        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            std::sort(beg + i, beg + i+v, std::greater<unsigned long>());
        d = clock() - t0;
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            if (!std::is_sorted(beg + i, beg + i+v, std::greater<unsigned long>()))
                abort();
        record.emplace_back(d, s);
        if (verbose) printf(" %s %.3g", s.c_str(), d / (double) CLOCKS_PER_SEC * v / blah.size());

        s = "qsort";
        t0 = clock();
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            qsort(&(blah[i]), v, sizeof(unsigned long), gcmp);
        d = clock() - t0;
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            if (!std::is_sorted(beg+i, beg+i+v, std::greater<unsigned long>()))
                abort();
        record.emplace_back(d, s);
        if (verbose) printf(" %s %.3g", s.c_str(), d / (double) CLOCKS_PER_SEC * v / blah.size());

        s = "iqsort";
        t0 = clock();
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v) {
#define islt(a,b) (*(a) > *(b))
            // coverity[escape]
            QSORT(unsigned long, &(blah[i]), v, islt);
#undef islt  
        }
        d = clock() - t0;
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            if (!std::is_sorted(beg+i, beg+i+v, std::greater<unsigned long>()))
                abort();
        record.emplace_back(d, s);
        if (verbose) printf(" %s %.3g", s.c_str(), d / (double) CLOCKS_PER_SEC * v / blah.size());

        s = "custom-insert";
        for(auto & x : blah)
            x = gmp_urandomm_ui(state, 1000);
        t0 = clock();
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            renumber_sort_ul(&(blah[i]), v);
        d = clock() - t0;
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            if (!std::is_sorted(beg+i, beg+i+v, std::greater<unsigned long>()))
                abort();
        record.emplace_back(d, s);
        if (verbose) printf(" %s %.3g", s.c_str(), d / (double) CLOCKS_PER_SEC * v / blah.size());

        s = "custom-merge";
        for(auto & x : blah)
            x = gmp_urandomm_ui(state, 1000);
        t0 = clock();
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            hard_mergesort(&(blah[i]), v);
        d = clock() - t0;
        for(unsigned int i = 0 ; i + v <= blah.size() ; i += v)
            if (!std::is_sorted(beg+i, beg+i+v, std::greater<unsigned long>()))
                abort();
        record.emplace_back(d, s);
        if (verbose) printf(" %s %.3g", s.c_str(), d / (double) CLOCKS_PER_SEC * v / blah.size());

        if (verbose) printf("\n");

        auto winner = *std::min_element(record.begin(), record.end());
        if (winner.second != last_best) {
            last_best = winner.second;
            best.emplace_back(v, winner.second);
        }
    }
    for(auto x : best)
        printf("n >= %u : %s\n", x.first, x.second.c_str());
    gmp_randclear (state);
    return 0;
}
