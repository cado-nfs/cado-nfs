#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <memory>
// iwyu wants it for allocator_traits<>::value_type, which seems weird
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "mpz_poly.h"

bool operator<(cxx_mpz_poly const& a, cxx_mpz_poly const& b) {
    return mpz_poly_cmp(a, b) < 0;
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    if (argc > 1) {
        srand(atoi(argv[1]));
    }

    std::vector<cxx_mpz_poly> v;

    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_poly x;
        int jmax = rand() % 16;
        for(int j = 0 ; j < jmax ; j++) {
            mpz_poly_setcoeff_si(x, j, (rand() - (RAND_MAX / 2)));
        }
        v.push_back(x);
    }
    sort(v.begin(), v.end());

    for(size_t i = 0 ; i < v.size() ; i++) {
        mpz_poly_fprintf(stdout, v[i]);
    }
}

