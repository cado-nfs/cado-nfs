#include "cado.h" // IWYU pragma: keep
                  //
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <algorithm>
#include <map>
#include <vector>

#include "fmt/base.h"
#include <gmp.h>

#include "mpz_mat.h"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "timing.h"
#include "macros.h"
#include "getprime.h"

static void silly_hnf_test(unsigned long seed)
{
    cxx_gmp_randstate state;
    if (seed)
        gmp_randseed_ui(state, seed);

    std::map<unsigned long, cxx_mpz_mat> v;

    cxx_mpz det;
    unsigned long const p = 1009;

    /* This generates many matrices at random, and puts in the std::map
     * above the relationship with the determinant of the reduction mod
     * 1009 of the leading submatrix of their HNF
     */
    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_mat M;
        unsigned int const m = gmp_urandomm_ui(state, p) % 16 + 2;
        unsigned int const n = gmp_urandomm_ui(state, p) % 16 + 2;
        mpz_mat_realloc(M, m, n);
        mpz_mat_urandomm(M, state, cxx_mpz(1009));

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

static void more_tests(unsigned long seed)
{
    cxx_gmp_randstate state;
    if (seed)
        gmp_randseed_ui(state, seed);

    for(int i = 0 ; i < 1000 ; i++) {
        cxx_mpz p;
        for(p = 1 ; !mpz_probab_prime_p(p, 2) ; mpz_urandomb(p, state, 16));

        cxx_mpz_mat A(2, 2);
        mpz_mat_urandomm(A, state, p);
        cxx_mpz detA, tr;
        mpz_mat_trace(tr, A);
        mpz_mat_determinant(detA, A);

        cxx_mpz_poly chi { detA, -tr, 1 };
        {
            cxx_mpz_mat Z;
            mpz_poly_eval_mpz_mat(Z, chi, A);
            ASSERT_ALWAYS(mpz_mat_is_zero(Z));
            mpz_poly_eval_mpz_mat_mod_ui(Z, chi, A, 11);
            ASSERT_ALWAYS(mpz_mat_is_zero(Z));

            cxx_mpz q = 1009;
            cxx_mpz_poly chi_q, psi_q;
            mpz_poly_mod_mpz(chi_q, chi, q, nullptr);
            mpz_poly_sub(psi_q, chi, chi_q);
            mpz_poly_divexact_mpz(psi_q, psi_q, q);
            mpz_poly_eval_mpz_mat_mod_mpz(Z, chi_q, A, q);
            ASSERT_ALWAYS(mpz_mat_is_zero(Z));
            mpz_poly_eval_mpz_mat(Z, chi_q, A);
            mpz_mat_mul_si(Z, Z, -1);
            mpz_mat_divexact_ui(Z, Z, mpz_get_ui(q));
            cxx_mpz_mat psi_A;
            mpz_poly_eval_mpz_mat(psi_A, psi_q, A);
            ASSERT_ALWAYS(mpz_mat_cmp(psi_A, Z) == 0);
        }

        cxx_mpz_mat B = A;
        mpz_mat_add_ui(B, 1);   // B = A + I

        /* now make sure that det(B) = chi(-1) */
        cxx_mpz detB;
        mpz_mat_determinant(detB, B);

        cxx_mpz chi_m1;
        mpz_poly_eval(chi_m1, chi, cxx_mpz(-1));

        ASSERT_ALWAYS(mpz_cmp(detB, chi_m1) == 0);

        {
            /* reversing rows of a 2x2 matrix should change the sign of
             * the determinant */
            cxx_mpz_mat Ar;

            auto unchanged = [&](auto B, cxx_mpz const & s = 1) {
                cxx_mpz d, ref;
                mpz_mat_determinant(d, B);
                mpz_mul(ref, detA, s);
                return mpz_cmp(d, ref) == 0;
            };

            mpz_mat_reverse_rows(Ar, A);
            ASSERT_ALWAYS(unchanged(Ar, -1));

            mpz_mat_reverse_columns(Ar, A);
            ASSERT_ALWAYS(unchanged(Ar, -1));

            Ar = A;
            mpz_mat_addmulrow(Ar, 0, 1, cxx_mpz(3));
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_addmulrow(Ar, 1, 0, cxx_mpz(-17));
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_submulrow(Ar, 0, 1, cxx_mpz(3));
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_submulrow(Ar, 1, 0, cxx_mpz(-17));
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_addrow(Ar, 1, 0);
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_subrow(Ar, 0, 1);
            ASSERT_ALWAYS(unchanged(Ar));
            mpz_mat_mulrow(Ar, 0, cxx_mpz(11));
            ASSERT_ALWAYS(unchanged(Ar, 11));

            mpz_mat_mul_mpz(Ar, A, cxx_mpz(7));
            ASSERT_ALWAYS(unchanged(Ar, 7*7));
            mpz_mat_mul_ui(Ar, A, 7);
            ASSERT_ALWAYS(unchanged(Ar, 7*7));
            mpz_mat_mul_si(Ar, A, -7);
            ASSERT_ALWAYS(unchanged(Ar, 7*7));
        }


        {
            /* mpq_mat_determinant_triangular needs specific testing */
            cxx_mpq_mat Xq = A;
            mpq_mat_div_mpz(Xq, Xq, cxx_mpz(17));
            mpq_set_ui(Xq(1, 0), 0, 1);
            cxx_mpq detXq;
            mpq_mat_determinant_triangular(detXq, Xq);
            cxx_mpq diag;
            mpz_mul(mpq_numref(diag), A(0,0), A(1,1));
            mpz_ui_pow_ui(mpq_denref(diag), 17, 2);
            mpq_canonicalize(diag);
            ASSERT_ALWAYS(mpq_cmp(diag, detXq) == 0);
        }

        {
            cxx_mpq_mat Aq = A;
            mpq_mat_div_mpz(Aq, Aq, cxx_mpz(17));

            cxx_mpq detAq;
            mpq_mat_determinant(detAq, Aq);

            auto unchanged = [&](auto B, cxx_mpq const & s = 1) {
                cxx_mpq d, ref;
                mpq_mat_determinant(d, B);
                mpq_mul(ref, detAq, s);
                return mpq_cmp(d, ref) == 0;
            };

            mpq_mat_transpose(Aq, Aq);
            ASSERT_ALWAYS(unchanged(Aq));

            cxx_mpq_mat Aqr = Aq;
            mpq_mat_addmulrow(Aqr, 0, 1, cxx_mpq(3));
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_addmulrow(Aqr, 1, 0, cxx_mpq(-17));
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_submulrow(Aqr, 0, 1, cxx_mpq(3));
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_submulrow(Aqr, 1, 0, cxx_mpq(-17));
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_addrow(Aqr, 1, 0);
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_subrow(Aqr, 0, 1);
            ASSERT_ALWAYS(unchanged(Aqr));
            mpq_mat_mulrow(Aqr, 0, cxx_mpq(11));
            ASSERT_ALWAYS(unchanged(Aqr, 11));

            mpq_mat_mul_mpq(Aqr, Aq, cxx_mpq(17, 42));
            ASSERT_ALWAYS(unchanged(Aqr, cxx_mpq(17*17, 42*42)));
            mpq_mat_mul_ui(Aqr, Aq, 17);
            ASSERT_ALWAYS(unchanged(Aqr, 17*17));
            mpq_mat_mul_si(Aqr, Aq, -17);
            ASSERT_ALWAYS(unchanged(Aqr, 17*17));
            mpq_mat_div_mpq(Aqr, Aq, cxx_mpq(17, 42));
            ASSERT_ALWAYS(unchanged(Aqr, cxx_mpq(42*42, 17*17)));
            mpq_mat_div_ui(Aqr, Aq, 17);
            ASSERT_ALWAYS(unchanged(Aqr, cxx_mpq(1, 17*17)));
            mpq_mat_div_si(Aqr, Aq, -17);
            ASSERT_ALWAYS(unchanged(Aqr, cxx_mpq(1, 17*17)));
        }
            

        if (mpz_cmp_ui(detA, 0) == 0) continue;

        {
            /* now make a block matrix A,B,C,D. We know that B and A commute,
             * so the determinant should be det(D*A-C*B)
             */
            cxx_mpz_mat C(2, 2);
            cxx_mpz_mat D(2, 2);
            mpz_mat_urandomm(C, state, p);
            mpz_mat_urandomm(D, state, p);
            cxx_mpz_mat ABCD, L = C;
            mpz_mat_horizontal_join(ABCD, A, B);

            /* A has same number of rows, but fewer columns than ABCD, so
             * it must compare below */
            ASSERT_ALWAYS(mpz_mat_cmp(A, ABCD) < 0);

            mpz_mat_horizontal_join(L, L, D);
            mpz_mat_vertical_join(ABCD, ABCD, L);

            /* L has fewer rows than ABCD, so it must compare below */
            ASSERT_ALWAYS(mpz_mat_cmp(L, ABCD) < 0);

            {
                cxx_mpz_mat X;
                mpz_mat_pow_ui(X, ABCD, 0);
                mpz_mat_mul(X, X, ABCD);
                ASSERT_ALWAYS(mpz_mat_cmp(X, ABCD) == 0);
            }

            cxx_mpz detABCD;
            mpz_mat_determinant(detABCD, ABCD);
            cxx_mpz_mat P, Q;
            mpz_mat_mul(P, D, A);
            mpz_mat_mul(Q, C, B);
            mpz_mat_sub(P, P, Q);
            cxx_mpz detP;
            mpz_mat_determinant(detP, P);
            ASSERT_ALWAYS(mpz_cmp(detP, detABCD) == 0);

            /* it's a very weak test of add/sub, but at least it's
             * something.
             */
            mpz_mat_add(P, P, Q);
            mpz_mat_sub(P, P, Q);
            mpz_mat_determinant(detP, P);
            ASSERT_ALWAYS(mpz_cmp(detP, detABCD) == 0);

            {
                /* reversing rows of a 2x2 matrix should change the sign of
                 * the determinant */
                cxx_mpz_mat ABCDr;
                cxx_mpz detABCDr;

                mpz_mat_reverse_rows(ABCDr, ABCD);
                mpz_mat_determinant(detABCDr, ABCDr);
                ASSERT_ALWAYS(mpz_cmp(detABCD, detABCDr) == 0);

                mpz_mat_reverse_columns(ABCDr, ABCD);
                mpz_mat_determinant(detABCDr, ABCDr);
                ASSERT_ALWAYS(mpz_cmp(detABCD, detABCDr) == 0);
            }

        }

        {
            /* do the same with mpq_mat */
            cxx_mpq_mat Aq = A;
            cxx_mpq_mat Bq = B;
            mpq_mat_div_mpz(Aq, Aq, cxx_mpz(17));
            mpq_mat_div_mpz(Bq, Bq, cxx_mpz(42));
            cxx_mpq_mat Cq(2, 2);
            cxx_mpq_mat Dq(2, 2);
            mpq_mat_urandomm(Cq, state, p);
            mpq_mat_urandomm(Dq, state, p);
            cxx_mpq_mat ABCDq, Lq = Cq;
            mpq_mat_horizontal_join(ABCDq, Aq, Bq);

            ASSERT_ALWAYS(mpq_mat_cmp(Aq, ABCDq) < 0);
            mpq_mat_horizontal_join(Lq, Lq, Dq);
            mpq_mat_vertical_join(ABCDq, ABCDq, Lq);
            ASSERT_ALWAYS(mpq_mat_cmp(Lq, ABCDq) < 0);
            cxx_mpq detABCDq;
            mpq_mat_determinant(detABCDq, ABCDq);

            cxx_mpq_mat Pq, Qq;
            mpq_mat_mul(Pq, Dq, Aq);
            mpq_mat_mul(Qq, Cq, Bq);
            mpq_mat_sub(Pq, Pq, Qq);
            cxx_mpq detPq;
            mpq_mat_determinant(detPq, Pq);

            ASSERT_ALWAYS(mpq_cmp(detPq, detABCDq) == 0);

            /* test swaprows */
            unsigned int i0 = gmp_urandomm_ui(state, 4);
            unsigned int i1 = gmp_urandomm_ui(state, 4);
            mpq_mat_swaprows(ABCDq, i0, i1);
            mpq_mat_determinant(detABCDq, ABCDq);
            if (i0 != i1)
                mpq_neg(detABCDq, detABCDq);
            ASSERT_ALWAYS(mpq_cmp(detPq, detABCDq) == 0);

            /* test transpose of rectangular matrices */
            cxx_mpq_mat U1(4, 2), U2(4, 2);
            mpq_mat_submat_set(U1, 0, 0, ABCDq, 0, 0, 4, 2);
            mpq_mat_submat_set(U2, 0, 0, ABCDq, 0, 2, 4, 2);
            mpq_mat_transpose(U1, U1);
            mpq_mat_transpose(U2, U2);
            mpq_mat_vertical_join(U1, U1, U2);
            mpq_mat_transpose(U1, U1);
            ASSERT_ALWAYS(mpq_mat_cmp(ABCDq, U1) == 0);
        }

        {
            cxx_mpz_poly T { tr, -2 };
            cxx_mpz res;
            mpz_poly_resultant(res, T, chi);
            prime_info pi;
            prime_info_init(pi);
            for (unsigned long p = 2; p < 100; p = getprime_mt (pi)) {
                if (p == 2)
                    continue;
                if (!mpz_divisible_ui_p(res, p))
                    continue;
                /* trace(A)-2*x = trace(A-x) and chi = det(x-A) have a
                 * common root mod p (which is trace(A)/2). So A-x has
                 * zero determinant and trace, thus it is nilpotent mod p */
                cxx_mpz minus_x = tr * ((p - 1) / 2);
                cxx_mpz_mat Y = A;
                mpz_mat_add_mpz(Y, minus_x);
                cxx_mpz dY, tY;
                mpz_mat_determinant(dY, Y);
                mpz_mat_trace(tY, Y);
                ASSERT_ALWAYS(mpz_divisible_p(dY, cxx_mpz(p)));
                ASSERT_ALWAYS(mpz_divisible_p(tY, cxx_mpz(p)));
                cxx_mpz_mat Z;
                Z = Y;
                mpz_mat_pow_ui_mod_ui(Z, Z, 2, p*p);
                ASSERT_ALWAYS(mpz_mat_p_valuation_ui(Z, p) > 0);
                Z = Y;
                mpz_mat_pow_ui(Z, Z, 2);
                ASSERT_ALWAYS(mpz_mat_p_valuation_ui(Z, p) > 0);
                Z = Y;
                mpz_mat_mul_mod_ui(Z, Z, Z, p*p);
                ASSERT_ALWAYS(mpz_mat_p_valuation_ui(Z, p) > 0);
            }
            prime_info_clear(pi);
        }
    }
}


// coverity[root_function]
int main(int argc, char const * argv[])
{
    unsigned int seed = 0;
    if (argc > 1) {
        seed = atoi(argv[1]);
    }

    silly_hnf_test(seed);

    hnf_timing_test(seed, false);

    more_tests(seed);
}
