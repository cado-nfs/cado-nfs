#include "cado.h"
#include "misc.h"
#include "getprime.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

int main(int argc, char * argv[])
{
    size_t n = 10;
    size_t p0 = 0;
    size_t p1 = 1000000;

    for(argv++,argc-- ; argc ; argv++,argc--) {
        if (argc > 1 && strcmp(argv[0], "-p0") == 0) {
            p0 = atol(argv[1]);
            argv++,argc--;
            continue;
        } else if (argc > 1 && strcmp(argv[0], "-p1") == 0) {
            p1 = atol(argv[1]);
            argv++,argc--;
            continue;
        } else if (argc > 1 && strcmp(argv[0], "-n") == 0) {
            n = atol(argv[1]);
            argv++,argc--;
            continue;
        } else {
            fprintf(stderr, "Unexpected arg: %s\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    std::vector<unsigned long> splits = subdivide_primes_interval(p0, p1, n);

    prime_info pi;
    prime_info_init(pi);
    size_t max = 0, min = SIZE_MAX;
    for(size_t i = 1 ; i < n ; i++) {
        size_t j = 0;
        unsigned long q0 = splits[i-1];
        unsigned long q1 = splits[i];
        prime_info_seek(pi, q0);
        for(unsigned long p ; (p = getprime_mt(pi)) < q1 ; j++);
        double est = nprimes_interval(q0, q1);
        size_t exact = j;
        printf("adjusted interval %zu/%zu: [%lu, %lu) (%.1f%%), estimate %.0f primes, exact %zu\n",
                i-1, n, q0, q1,
                100.0 * (q1 - q0) / (p1 - p0) * n,
                est, exact);
        if (exact > max) max = exact;
        if (exact < min) min = exact;
    }
    printf("adjusted max/min = %.3f\n", (double) max/min);
    prime_info_clear(pi);

    printf("Here follows what simplistic equal-size intervals would give\n");
    max = 0, min = SIZE_MAX;
    prime_info_init(pi);
    for(size_t i = 1 ; i < n ; i++) {
        size_t j = 0;
        unsigned long q0 = p0 + (i-1) * (p1 - p0) / n;
        unsigned long q1 = p0 + (i) * (p1 - p0) / n;
        prime_info_seek(pi, q0);
        for(unsigned long p ; (p = getprime_mt(pi)) < q1 ; j++);
        double est = nprimes_interval(q0, q1);
        size_t exact = j;
        printf("arithmetic interval %zu/%zu: [%lu, %lu), estimate %.0f primes, exact %zu\n",
                i-1, n, q0, q1,
                est, exact);
        if (exact > max) max = exact;
        if (exact < min) min = exact;
    }
    printf("arithmetic max/min = %.3f\n", (double) max/min);
    prime_info_clear(pi);

    return 0;
}

