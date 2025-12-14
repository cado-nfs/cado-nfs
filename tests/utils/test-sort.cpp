#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <string>

#include <utility>
#include <vector>
#include <functional>
#include <algorithm>
#include <ranges>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "macros.h"
#include "iqsort.h"
#include "utils_cxx.hpp"

struct test_env {
    std::vector<unsigned long> blah;
    cxx_gmp_randstate state;
    std::vector<std::pair<double, std::string>> record;
    explicit test_env(size_t N)
        : blah(N)
    {
        fill();
    }
    void fill() {
        for(auto & x : blah)
            x = gmp_urandomm_ui(state, 1000);
    }
};

struct custom_insert {
    static constexpr const char * name = "custom-insert";
    void operator() (unsigned long *r, size_t n) const
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
                unsigned long const rj = r[j];
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
};

struct custom_mergesort {
    static constexpr const char * name = "custom-merge";
    void operator() (unsigned long *r, size_t n) const
    {
        if (n > 256) {
            std::sort(r, r + n, std::greater<unsigned long>());
            return;
        }
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
};

static double get_speed(clock_t d, size_t nitems)
{
    return double_ratio(d, CLOCKS_PER_SEC, nitems);
}

struct iqsort {
    static constexpr const char * name = "iqsort";
    void operator()(unsigned long * data, size_t v) const
    {
#define islt(a,b) (*(a) > *(b))
        // coverity[escape]
        QSORT(unsigned long, data, v, islt);
#undef islt  
    }
};

struct c_qsort {
    static constexpr const char * name = "qsort";
    static int gcmp(const void * pa, const void * pb)
    {
        auto const a = * static_cast<const unsigned long *>(pa);
        auto const b = * static_cast<const unsigned long *>(pb);
        return (a < b) - (b < a);
    }

    void operator()(unsigned long * data, size_t v) const {
        qsort(data, v, sizeof(unsigned long), gcmp);
    }
};

struct stdsort {
    static constexpr const char * name = "std::sort";
    void operator()(unsigned long * data, size_t v) const
    {
        std::sort(data, data + v, std::greater<unsigned long>());
    }
};

template<typename T>
static void test_one(T const & F, test_env & E, size_t v, bool verbose)
{
    E.fill();
    const clock_t t0 = clock();
    for(unsigned int i = 0 ; i + v <= E.blah.size() ; i += v)
        F(&(E.blah[i]), v);
    const clock_t d = clock() - t0;
    for(size_t i = 0 ; i + v <= E.blah.size() ; i += v)
        if (!std::ranges::is_sorted(
                    E.blah | std::views::take(v),
                    std::greater<unsigned long>()))
            abort();
    if (verbose)
        fmt::print(" {} {:.3g}", T::name, get_speed(d, E.blah.size()));
    E.record.emplace_back(d, T::name);
}


int main(int argc, char const * argv[])
{
    double exp_growth = 0.1;
    size_t logsize = 18;
    size_t maxspan = 1 << 9;
    bool verbose = false;

    argv++, argc--;
    for( ; argc ; argc--,argv++) {
        std::string const s = *argv;
        if (s == "--verbose" || s == "-v") {
            verbose = true;
            continue;
        }
        if (s == "--exp-growth") {
            char * p;
            exp_growth = strtod(argv[1], &p);
            ASSERT_ALWAYS(*p == '\0');
            argc--,argv++;
            continue;
        }
        if (s == "--log-bigsize") {
            char * p;
            logsize = strtol(argv[1], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argc--,argv++;
            continue;
        }
        if (s == "--max-span") {
            char * p;
            maxspan = strtol(argv[1], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
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

    std::vector<std::pair<unsigned int, std::string>> best;
    std::string last_best = "none";
    test_env E(1 << logsize);
    for(size_t v = 2 ; v <= maxspan ; v += 1 + size_t(double(v) * exp_growth)) {

        E.record.clear();

        if (verbose) printf("%zu", v);

        test_one(stdsort(), E, v, verbose);
        test_one(c_qsort(), E, v, verbose);
        test_one(iqsort(), E, v, verbose);
        test_one(custom_insert(), E, v, verbose);
        test_one(custom_mergesort(), E, v, verbose);
        if (verbose) printf("\n");

        auto winner = *std::ranges::min_element(E.record);
        if (winner.second != last_best) {
            last_best = winner.second;
            best.emplace_back(v, winner.second);
        }
    }
    for(auto const & x : best)
        printf("n >= %u : %s\n", x.first, x.second.c_str());
    return 0;
}
