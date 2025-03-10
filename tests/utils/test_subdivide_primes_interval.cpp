#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <climits>

#include <algorithm>
#include <vector>

#include "fmt/base.h"

#include "misc.h"
#include "getprime.h"
#include "params.h"

static void decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "p0", "prime lower bound");
    param_list_decl_usage(pl, "p1", "prime upper bound");
    param_list_decl_usage(pl, "n", "number of intervals");
}

int main(int argc, char const * argv[])
{
    const char * progname = argv[0];

    size_t n = 10;
    size_t p0 = 0;
    size_t p1 = 1000000;

    cxx_param_list pl;

    decl_usage(pl);

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        fmt::print(stderr, "Unknown option: {}\n", argv[0]);
        param_list_print_usage(pl, progname, stderr);
        return EXIT_FAILURE;
    }

    param_list_parse(pl, "p0", p0);
    param_list_parse(pl, "p1", p1);
    param_list_parse(pl, "n", n);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, progname, stderr);
        return EXIT_FAILURE;
    }

    std::vector<unsigned long> splits = subdivide_primes_interval(p0, p1, n);

    prime_info pi;
    prime_info_init(pi);
    size_t max = 0, min = SIZE_MAX;
    for(size_t i = 1 ; i < n ; i++) {
        size_t j = 0;
        const unsigned long q0 = splits[i-1];
        const unsigned long q1 = splits[i];
        prime_info_seek(pi, q0);
        for( ; getprime_mt(pi) < q1 ; j++);
        const double est = nprimes_interval(double(q0), double(q1));
        const size_t exact = j;
        fmt::print("adjusted interval {}/{}: [{}, {}) ({:.1f}%%), estimate {:.0f} primes, exact {}\n",
                i-1, n, q0, q1,
                100.0 * double(q1 - q0) / double(p1 - p0) * double(n),
                est, exact);
        max = std::max(max, exact);
        min = std::min(min, exact);
    }
    fmt::print("adjusted max/min = {:.3f}\n", double(max) / double(min));
    prime_info_clear(pi);

    fmt::print("Here follows what simplistic equal-size intervals would give\n");
    max = 0, min = SIZE_MAX;
    prime_info_init(pi);
    for(size_t i = 1 ; i < n ; i++) {
        size_t j = 0;
        const unsigned long q0 = p0 + (i-1) * (p1 - p0) / n;
        const unsigned long q1 = p0 + (i) * (p1 - p0) / n;
        prime_info_seek(pi, q0);
        for( ; getprime_mt(pi) < q1 ; j++);
        double est = nprimes_interval(double(q0), double(q1));
        const size_t exact = j;
        fmt::print("arithmetic interval {}/{}: [{}, {}), estimate {:.0f} primes, exact {}\n",
                i-1, n, q0, q1,
                est, exact);
        max = std::max(max, exact);
        min = std::min(min, exact);
    }
    fmt::print("arithmetic max/min = {:.3f}\n", double(max) / double(min));
    prime_info_clear(pi);

    return 0;
}

