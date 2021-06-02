#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <iosfwd>     // for std
#include <utility>    // for pair, make_pair
#include <map>
#include <algorithm>
#include <gmp.h>
#include "mpz_mat.h"

using namespace std;

// coverity[root_function]
int main(int argc, char * argv[])
{
    if (argc > 1) {
        srand(atoi(argv[1]));
    }

    map<unsigned long, cxx_mpz_mat> v;

    mpz_t det;
    mpz_init(det);
    unsigned long p = 1009;

    /* This generates many matrices at random, and puts in the std::map
     * above the relationship with the determinant of the reduction mod
     * 1009 of the leading submatrix of their HNF
     */
    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_mat M;

        mpz_mat_realloc(M, rand() % 16 + 2, rand() % 16 + 2);
        for(unsigned int i = 0 ; i < M->m ; i++) {
            for(unsigned int j = 0 ; j < M->n ; j++) {
                mpz_set_si(mpz_mat_entry(M, i, j), (rand() - (RAND_MAX / 2)));
            }
        }
        mpz_mat_mod_ui(M, M, p);
        cxx_mpz_mat M1 = M, T, M2;
        mpz_mat_hnf_backend(M1, T);
        // mpz_mat_fprint(stdout, M1);
        unsigned int d = std::min(M1->m, M1->n);
        mpz_mat_realloc(M2, d, d);
        mpz_mat_submat_swap(M2, 0, 0, M1, 0, 0, d, d);
        mpz_mat_determinant_triangular(det, M2);
        mpz_mod_ui(det, det, p);
        v.insert(make_pair(mpz_get_ui(det), M));
    }

    for(map<unsigned long, cxx_mpz_mat>::const_iterator it = v.begin() ; it != v.end() ; it++) {
        printf("[det=%ld] ", it->first);
        mpz_mat_fprint(stdout, it->second);
    }

    mpz_clear(det);
}
