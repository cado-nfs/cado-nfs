#include "cado.h" // IWYU pragma: keep
                  //
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <map>
#include <vector>

#include "fmt/base.h"
#include <gmp.h>

#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "timing.h"

static void silly_hnf_test(unsigned long seed)
{
    if (seed)
        srand(seed);

    std::map<unsigned long, cxx_mpz_mat> v;

    cxx_mpz det;
    unsigned long const p = 1009;

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
        cxx_mpz_mat M1 = M, M2;
        mpz_mat_hermite_form(M1);
        // mpz_mat_fprint(stdout, M1);
        unsigned int const d = std::min(M1->m, M1->n);
        mpz_mat_realloc(M2, d, d);
        mpz_mat_submat_swap(M2, 0, 0, M1, 0, 0, d, d);
        mpz_mat_determinant_triangular(det, M2);
        mpz_mod_ui(det, det, p);
        v.emplace(mpz_get_ui(det), M);
    }

    for(auto const & [d, M] : v) {
        printf("[det=%ld] ", d);
        mpz_mat_fprint(stdout, M);
    }
}

static void hnf_timing_test(unsigned long seed, bool check_performance = false)
{
    cxx_gmp_randstate state;
    if (seed)
        gmp_randseed_ui(state, seed);

    std::vector<int> bits_choices = { 8, 32, 64 };
    std::vector<int> dims_choices = { 4, 8, 16 };

    if (check_performance) {
        bits_choices = { 8, 32, 64, 128 };
        dims_choices = { 4, 8, 16, 32, 64, 128 /*, 256 */ };
    }

    for(int const bits : bits_choices) {
        cxx_mpz p;
        for(p = 1 ; !mpz_probab_prime_p(p, 2) ; mpz_urandomb(p, state, bits));
        for(unsigned int const n : dims_choices) {
            cxx_mpz_mat H(2 * n, n);
            for(unsigned int i = 0 ; i < H->m ; i++)
                for(unsigned int j = 0 ; j < H->n ; j++)
                    mpz_urandomm(H(i, j), state, p);
            // fmt::print("A="); mpz_mat_fprint(stdout, H); fmt::print("\n");
            double t = -wct_seconds();
            mpz_mat_hermite_form(H);
            t += wct_seconds();
            fmt::print("{} {} {} {}\n",
                    bits, H->m, H->n, t);
            // fmt::print("H="); mpz_mat_fprint(stdout, H); fmt::print("\n");
            // fmt::print("T="); mpz_mat_fprint(stdout, T); fmt::print("\n");
        }
    }

    /* magma comparison code
        for bits in [8,32,64,128] do
            p:=RandomPrime(bits);
            for n in [8,16,32,64,128,256] do
                A:=Matrix(2*n, n, [Integers()!Random(GF(p)):i in [1..2*n^2]]);
                t:=Cputime();
                H:=HermiteForm(A);
                t:=Cputime()-t;
                printf "%o %o %o %o\n", bits, 2*n, n, t;
            end for;
        end for;
     */
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    unsigned int seed = 0;
    if (argc > 1) {
        seed = atoi(argv[1]);
    }

    silly_hnf_test(seed);

    hnf_timing_test(seed);
}
