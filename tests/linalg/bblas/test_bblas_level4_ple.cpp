#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include "bblas_mat8.hpp"
#include <cstdio>                // for printf
#include <cstdint>
#include <algorithm>
#include <string>                 // for basic_string, string
#include <utility>                // for pair
#include <vector>                 // for vector

#include <gmp.h>                  // for gmp_urandomm_ui, mpz_clear, mpz_div...

#include "bblas_bitmat.hpp"       // for bitmat<>::vector_type, bitmat
#include "bblas_level4_ple_internal.hpp"
#include "bpack.hpp"              // for bpack, bpack_view, bpack_const_view
#include "cxx_mpz.hpp"
#include "gmp_aux.h"              // for mpz_get_uint64, memfill_random
#include "macros.h"               // for ASSERT_ALWAYS
#include "test_bblas_base.hpp"    // for test_bblas_base, test_bblas_base::t...
#include "time_bblas_common.hpp"
#include "test_bblas_level4.hpp"

template<typename T>
int test_bblas_level4::test_PLE_find_pivot(unsigned int m, unsigned int n)/*{{{*/
{
    using matrix = bitmat<T>;
    constexpr const unsigned int B = matrix::width;
    bpack<T> A(m, n);
    PLE<T> ple(A.view());

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        std::vector<unsigned int> f;
        std::vector<unsigned int> u;
        for(unsigned int j = 0 ; j < n ; j++) {
            f.push_back(gmp_urandomb_ui(rstate, 30));
            u.push_back(1 + gmp_urandomm_ui(rstate, 3*m-1));
        }
        A = 0;
        for(unsigned int j = 0 ; j < n ; j++) {
            for(unsigned int i = 0 ; i < m ; i++) {
                T b = (f[j]+i) % u[j] == 0;
                A.X[(i/B)*(n/B)+(j/B)][i%B] ^= (b << (j%B));
            }
        }
        for(unsigned int ell = 0 ; ell < 100 ; ell++) {
            unsigned int const jj = gmp_urandomm_ui(rstate, n);
            unsigned int const ii = gmp_urandomm_ui(rstate, m);
            int const p = ple.find_pivot(ii / B, jj / B, ii % B, jj % B);
            unsigned int const x = f[jj]+ii;
            unsigned int const d = u[jj];
            /* first r>=0 such that x+r is a multiple of d */
            unsigned int const r = (d - 1) - ((x-1) % d);
            if ((ii + r) >= m) {
                ASSERT_ALWAYS(p == -1);
            } else {
                ASSERT_ALWAYS(p >= 0);
                ASSERT_ALWAYS((unsigned int) p == (ii+r));
            }
        }
    }

    return 0;
}/*}}}*/

template<typename T>
int test_bblas_level4::test_PLE_propagate_pivot(unsigned int m, unsigned int n)/*{{{*/
{
    using matrix = bitmat<T>;
    constexpr const unsigned int B = matrix::width;
    bpack<T> A(m, n);
    PLE<T> ple(A.view());

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        unsigned int const kjj0 = gmp_urandomm_ui(rstate, n - 1) + 1;
        unsigned int const kjj1 = gmp_urandomm_ui(rstate, n + 1 - kjj0) + kjj0;
        unsigned int const pjj = gmp_urandomm_ui(rstate, kjj0);
        unsigned int const pii = gmp_urandomm_ui(rstate, m);
        A = 0;
        T examples[n/B][4];
        for(unsigned int bj = 0 ; bj < n/B ; bj++) {
            T a = 0;
            if (kjj0 / B < bj)
                a ^= ~T(0);
            else if (kjj0 / B == bj)
                a ^= -(T(1) << (kjj0 % B));
            if (kjj1 / B < bj)
                a ^= ~T(0);
            else if (kjj1 / B == bj)
                a ^= -(T(1) << (kjj1 % B));
            examples[bj][0] = a;
            examples[bj][1] = a;
            examples[bj][2 + (pii & 1)] = 0;
            examples[bj][2 + (1 ^ (pii & 1))] = a;
            if (pjj / B == bj) {
                examples[bj][pii & 1] ^= T(1) << (pjj % B);
                examples[bj][2 + (pii & 1)] ^= T(1) << (pjj % B);
            }
        }

        for(unsigned int bj = 0 ; bj < n/B ; bj++) {
            for(unsigned int ii = 0 ; ii < m ; ii++) {
                A.X[(ii/B)*(n/B)+bj][ii%B] = examples[bj][ii & 1];
            }
        }
        ple.propagate_pivot(pii / B, pjj / B, pii % B, pjj % B);
        unsigned int const bj = pjj / B;
        for(unsigned int ii = 0 ; ii < m ; ii++) {
            T a = examples[bj][2 * (ii > pii) + (ii & 1)];
            ASSERT_ALWAYS(A.X[(ii/B)*(n/B)+bj][ii%B] == a);
        }
    }

    return 0;
}/*}}}*/

template<typename T>
int test_bblas_level4::test_PLE_propagate_row_permutations(unsigned int m, unsigned int n)/*{{{*/
{
    using matrix = bitmat<T>;
    constexpr const unsigned int B = matrix::width;
    bpack<T> A(m, n);
    PLE<T> ple(A.view());

    for(unsigned int ii = 0, kk = 0 ; ii < m ; ii++) {
        for(unsigned int bj = 0 ; bj < n/B ; bj++, kk++) {
            unsigned int  const bi =  ii / B;
            unsigned int   const i =  ii % B;
            A.X[bi * n/B + bj][i] = kk;
        }
    }

    for(unsigned int k = 1 ; k <= 1000 ; k++) {
        unsigned int const ii0 = gmp_urandomm_ui(rstate, m);
        unsigned int const ii1 = std::min((mp_limb_t) k, gmp_urandomm_ui(rstate, m + 1 - ii0)) + ii0;
        unsigned int const bj0 = gmp_urandomm_ui(rstate, n/B);

        std::vector<unsigned int> q;

        for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
            unsigned int const pii = gmp_urandomm_ui(rstate, m - ii) + ii;
            q.push_back(pii);
            if (ii == pii) continue;
            unsigned int  const bi =  ii / B;
            unsigned int   const i =  ii % B;
            unsigned int const pbi = pii / B;
            unsigned int  const pi = pii % B;
            bitmat<T> & Y = A.X[bi * n/B + bj0];
            bitmat<T> & pY = A.X[pbi * n/B + bj0];
            T c = Y[i] ^ pY[pi];
            Y[i] ^= c;
            pY[pi] ^= c;
        }

        ASSERT_ALWAYS(q.size() == (unsigned int) (ii1 - ii0));

        ple.propagate_row_permutations(ii1, bj0, q.begin(), q.end());

        for(unsigned int ii = 0 ; ii < m ; ii++) {
            for(unsigned int bj = 0 ; bj < n/B ; bj++) {
                unsigned int  const bi =  ii / B;
                unsigned int   const i =  ii % B;
                ASSERT_ALWAYS(A.X[bi * n/B + bj][i] == T(A.X[bi * n/B][i] + bj));
            }
        }
    }

    return 0;
}/*}}}*/

template<typename T>
int test_bblas_level4::test_PLE_move_L_fragments(unsigned int m, unsigned int n)/*{{{*/
{
    using matrix = bitmat<T>;
    constexpr const unsigned int B = matrix::width;
    bpack<T> A(m, n);
    PLE<T> ple(A.view());

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        /* first generate the transpose of the bitmat<T> that we will
         * reduce.
         * We'll generate matrices with the following shape. The space is
         * only here as an extra information indicating the diagonal
         *
         * 1 010101... 
         * 01 010101...
         * 10/0 0000
         * 01/10101010...
         * 101/10101010...
         *
         * where the right part of some rows is intentionally shifted
         * down. The number of pivot rows generated is r, which is less
         * than min(m,n). We check that the lower triangular part of the
         * output bitmat<T> has the expected shape, which is:
         *
         * --------..
         * 0-------..
         * 10------..
         * and so on, with m rows and r columns.
         *
         */

        unsigned int const r = gmp_urandomm_ui(rstate, std::min(m, n) + 1);
        std::vector<unsigned int> pivs;
        /* generating binomial(n, r) is like generating with repetition r
         * things among n-(r-1) = n-r+1.  */
        pivs.reserve(r);
        for(unsigned int k = 0 ; k < r ; k++)
            pivs.push_back(gmp_urandomm_ui(rstate, n-r+1));
        std::ranges::sort(pivs);
        for(unsigned int k = 0 ; k < r ; k++)
            pivs[k] += k;
        {
            bpack<T> tA(n, m);
            tA = 0;
            auto ppiv = pivs.begin();
            unsigned int rr = 0;
            for(unsigned int jj = 0 ; jj < n ; jj++) {
                bool is_piv= false;
                unsigned int const bj = jj / B;
                if (ppiv < pivs.end() && jj == *ppiv) {
                    is_piv = true;
                    ppiv++;
                }
                /* generate rr bits of 1010/0101 pattern, depending on jj.
                 * and then:
                 * if is_piv, m-rr bits of 1010/0101 pattern, depending on rr, 
                 * if !is_piv, m-rr bits of 0's
                 */
                mpz_t NN,PP;
                mpz_init(NN);
                mpz_init(PP);
                /* jj even: rr bits of 101010 (starting with 1) is NN, with:
                 *      if rr is even, 3NN+1=(1<<rr)
                 *      if rr is odd, 3NN+1=(2<<rr)
                 * jj odd: rr bits of 010101 (starting with 0) is NN, with:
                 *      if rr is even, 3NN+2=(2<<rr)
                 *      if rr is odd, 3NN+2=(1<<rr)
                 *
                 * so that it's 3NN+1+(jj&1) = 1 << (rr + ((rr^jj)&1))
                 */
                if (jj < r) {
                    mpz_ui_pow_ui(NN, 2, rr + ((rr ^ jj) & 1));
                    mpz_sub_ui(NN,NN,1 + (jj & 1));
                    ASSERT_ALWAYS(mpz_divisible_ui_p(NN, 3));
                    mpz_divexact_ui(NN,NN,3);
                } else {
                    mpz_set_ui(NN, 0);
                }
                if (is_piv) {
                    /* rest: we'll make an m-bit mask of 1010..., i.e. 3NN+1=(1<<m)
                    */
                    mpz_ui_pow_ui(PP, 2, m);
                    mpz_sub_ui(PP, PP, 1);
                    ASSERT_ALWAYS(mpz_divisible_ui_p(PP, 3));
                    mpz_divexact_ui(PP,PP,3);
                    mpz_mul_2exp(PP, PP, rr);
                    mpz_add(NN, NN, PP);
                    mpz_clear(PP);
                    rr++;
                }
                static_assert(B <= 64, "we need B<=64 since we use mpz_get_uint64");
                for(unsigned int bi = 0 ; bi < m/B ; bi++) {
                    tA.X[bi + bj * (m/B)][jj % B] = mpz_get_uint64(NN);
                    mpz_fdiv_q_2exp(NN, NN, B);
                }
                mpz_clear(NN);
            }
            for(unsigned int bi = 0 ; bi < m/B ; bi++) {
                for(unsigned int bj = 0 ; bj < n/B ; bj++) {
                    bitmat<T>::transpose(A.X[bi * (n/B) + bj], tA.X[bi + bj * (m/B)]);
                }
            }
        }
        
        /* Now pivs is the set of _columns_ that we have to move. Let's
         * do that. As per the specification of move_L_fragments, we'll
         * do it in chunks where the source blocks and destination blocks
         * are constant.
         */

        for(unsigned int rr = 0, jj = 0, rr1 = 0 ; rr < r ; rr = rr1) {
            for( ; rr < r && pivs[rr] == jj ; rr++, jj++);
            if (rr == r) break;
            rr1 = rr + 1;
            unsigned int const br = rr / B;
            unsigned int const bp = pivs[rr] / B;
            for( ; rr1 < r && rr1 / B == br && pivs[rr1] / B == bp ; rr1++);
            /* move rr1 - rr columns together */
            std::vector<unsigned int> const Q(pivs.begin() + rr, pivs.begin() + rr1);

            ple.move_L_fragments(rr, Q);
        }
        
        /* Now, check that we have what we should have.
         */
        /* ii even: z=min(ii, r) bits of 101010 (starting with 1) is NN, with:
         *      if z is even, 3NN+1=(1<<z)
         *      if z is odd, 3NN+1=(2<<z)
         * ii odd: z=min(ii, r) bits of 010101 (starting with 0) is NN, with:
         *      if z is even, 3NN+2=(2<<z)
         *      if z is odd, 3NN+2=(1<<z)
         * so that it's 3NN+1+(ii&1) = 1 << (z + ((z^ii)&1))
         */
        for(unsigned int ii = 0 ; ii < m ; ii++) {
            cxx_mpz NN;
            unsigned int z = std::min(ii, r);
            mpz_ui_pow_ui(NN, 2, z + ((z ^ ii) & 1));
            mpz_sub_ui(NN,NN,1 + (ii & 1));
            ASSERT_ALWAYS(mpz_divisible_ui_p(NN, 3));
            mpz_divexact_ui(NN,NN,3);
            for(unsigned int bj = 0 ; bj < n/B ; bj++, z -= B) {
                T c = A.X[(ii/B)*(n/B)+bj][ii%B];
                c ^= mpz_get_uint64(NN);
                if (z < B) c &= (T(1) << z) - 1;
                mpz_fdiv_q_2exp(NN, NN, B);
                ASSERT_ALWAYS(c == 0);
                if (z < B) break;
            }
        }
    }

    return 0;
}/*}}}*/

template<typename T>
int test_bblas_level4::test_PLE(unsigned int m, unsigned int n)
{
    for(unsigned int k = 0 ; k < 100 ; k++) {
        bpack<T> A(m, n);
        A.fill_random(rstate);

        /* The main important thing is the fact of enabling the
         * debug_stuff object, which does invariant checks throughout the
         * PLE computation. As a matter of fact, the final test that
         * we do here is not even needed, since it's already one of those
         * checks triggered by debug_stuff.
         */
        PLE<T> ple(A.view());
        typename PLE<T>::debug_stuff D(ple);
        std::vector<unsigned int> pivs = ple(&D);

        D.check(ple.const_view(), pivs.begin(), pivs.size());
    }

    return 0;
}

template<typename T>
void test_bblas_level4::meta_ple()
{
    using matrix = bitmat<T>;
    constexpr const unsigned int B = matrix::width;
    std::vector<std::pair<unsigned int, unsigned int>> const mns
    {{
         {B,B},
         {B,2*B},
         {B,3*B},
         {2*B,B},
         {2*B,2*B},
         {2*B,3*B},
         {2*B,10*B},
         {10*B,2*B},
         {24*B,16*B},
         {2*B,100*B},
         {100*B,2*B},
         {2*B,1000*B},
         {1000*B,2*B},
         {32*B,32*B},
     }};

    for(auto x : mns) {
        unsigned int m = x.first;
        unsigned int n = x.second;
        if (m + n < 10*B) {
            test_PLE_find_pivot<T>(m, n);
            test_PLE_propagate_pivot<T>(m, n);
            test_PLE_propagate_row_permutations<T>(m, n);
            test_PLE_move_L_fragments<T>(m, n);
            test_PLE<T>(m, n);
        }

        typename matrix::vector_type X ((m/B)*(n/B), 0);

        auto randomize = [&]() { memfill_random(&X[0], m/B*n/B*sizeof(bitmat<T>), rstate); };
        auto do_ple = [&](bitmat<T> * X, unsigned int mm, unsigned int nn) {
            return bpack_view<T>(X, mm, nn).ple().size();
        };

        printf(" -- PLE(m=%u, n=%u)\n", m, n);

        // TIME1N_SPINS(randomize(), 2, do_ple, &X[0], m/B, n/B);
        std::string const what = "PLE";       // complete me

        if (m + n >= 10 * B) {
            randomize();
            TIME1(2, do_ple, X.data(), m/B, n/B);
        } else {
            bblas_timer(4, what).time1n_classify(n, randomize, do_ple, X.data(), m/B, n/B);
        }
#ifdef TIME_PLE
        PLE<T>::print_and_flush_stats();
#endif
    }

}

test_bblas_base::tags_t test_bblas_level4::ple_tags { "ple", "l4" };
void test_bblas_level4::ple()
{
    meta_ple<uint64_t>();
    meta_ple<uint8_t>();
}

