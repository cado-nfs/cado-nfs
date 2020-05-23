#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "mpz_poly.h"
#include "tests_common.h"

using namespace std;

bool operator<(cxx_mpz_poly const& a, cxx_mpz_poly const& b) {
    return mpz_poly_cmp(a, b) < 0;
}

int main(int argc, char * argv[])
{
    if (argc > 1) {
        srand(atoi(argv[1]));
    }

    vector<cxx_mpz_poly> v;

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

