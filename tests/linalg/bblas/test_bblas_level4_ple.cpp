#include "cado.h"
#include "test_bblas_level4.hpp"
#include "level4_ple_internal.hpp"
#include <cstring>
#include <algorithm>

int test_bblas_level4::test_PLE_find_pivot(unsigned int m, unsigned int n)/*{{{*/
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        std::vector<unsigned int> f;
        std::vector<unsigned int> u;
        for(unsigned int j = 0 ; j < n ; j++) {
            f.push_back(gmp_urandomb_ui(rstate, 30));
            u.push_back(1 + gmp_urandomm_ui(rstate, 3*m-1));
        }
        memset((void *) X, 0, (m/B) * (n/B) * sizeof(mat64));
        for(unsigned int j = 0 ; j < n ; j++) {
            for(unsigned int i = 0 ; i < m ; i++) {
                uint64_t b = (f[j]+i) % u[j] == 0;
                X[(i/B)*(n/B)+(j/B)][i%B] ^= (b << (j%B));
            }
        }
        for(unsigned int ell = 0 ; ell < 100 ; ell++) {
            unsigned int jj = gmp_urandomm_ui(rstate, n);
            unsigned int ii = gmp_urandomm_ui(rstate, m);
            int p = ple.find_pivot(ii / B, jj / B, ii % B, jj % B);
            unsigned int x = f[jj]+ii;
            unsigned int d = u[jj];
            /* first r>=0 such that x+r is a multiple of d */
            unsigned int r = (d - 1) - ((x-1) % d);
            if ((ii + r) >= m) {
                ASSERT_ALWAYS(p == -1);
            } else {
                ASSERT_ALWAYS(p >= 0);
                ASSERT_ALWAYS((unsigned int) p == (ii+r));
            }
        }
    }
    delete[] X;

    return 0;
}/*}}}*/

int test_bblas_level4::test_PLE_propagate_pivot(unsigned int m, unsigned int n)/*{{{*/
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        unsigned int kjj0 = gmp_urandomm_ui(rstate, n - 1) + 1;
        unsigned int kjj1 = gmp_urandomm_ui(rstate, n + 1 - kjj0) + kjj0;
        unsigned int pjj = gmp_urandomm_ui(rstate, kjj0);
        unsigned int pii = gmp_urandomm_ui(rstate, m);
        memset((void *) X, 0, (m/B) * (n/B) * sizeof(mat64));
        uint64_t examples[n/B][4];
        for(unsigned int bj = 0 ; bj < n/B ; bj++) {
            uint64_t a = 0;
            if (kjj0 / B < bj)
                a ^= ~UINT64_C(0);
            else if (kjj0 / B == bj)
                a ^= -(UINT64_C(1) << (kjj0 % B));
            if (kjj1 / B < bj)
                a ^= ~UINT64_C(0);
            else if (kjj1 / B == bj)
                a ^= -(UINT64_C(1) << (kjj1 % B));
            examples[bj][0] = a;
            examples[bj][1] = a;
            examples[bj][2 + (pii & 1)] = 0;
            examples[bj][2 + (1 ^ (pii & 1))] = a;
            if (pjj / B == bj) {
                examples[bj][pii & 1] ^= UINT64_C(1) << (pjj % B);
                examples[bj][2 + (pii & 1)] ^= UINT64_C(1) << (pjj % B);
            }
        }

        for(unsigned int bj = 0 ; bj < n/B ; bj++) {
            for(unsigned int ii = 0 ; ii < m ; ii++) {
                X[(ii/B)*(n/B)+bj][ii%B] = examples[bj][ii & 1];
            }
        }
        ple.propagate_pivot(pii / B, pjj / B, pii % B, pjj % B);
        unsigned int bj = pjj / B;
        for(unsigned int ii = 0 ; ii < m ; ii++) {
            uint64_t a = examples[bj][2 * (ii > pii) + (ii & 1)];
            ASSERT_ALWAYS(X[(ii/B)*(n/B)+bj][ii%B] == a);
        }
    }
    delete[] X;

    return 0;
}/*}}}*/

int test_bblas_level4::test_PLE_propagate_permutations(unsigned int m, unsigned int n)/*{{{*/
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int ii = 0, kk = 0 ; ii < m ; ii++) {
        for(unsigned int bj = 0 ; bj < n/B ; bj++, kk++) {
            unsigned int  bi =  ii / B;
            unsigned int   i =  ii % B;
            X[bi * n/B + bj][i] = kk;
        }
    }

    for(unsigned int k = 1 ; k <= 1000 ; k++) {
        unsigned int ii0 = gmp_urandomm_ui(rstate, m);
        unsigned int ii1 = std::min((mp_limb_t) k, gmp_urandomm_ui(rstate, m + 1 - ii0)) + ii0;
        unsigned int bj0 = gmp_urandomm_ui(rstate, n/B);

        std::vector<unsigned int> q;

        for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
            unsigned int pii = gmp_urandomm_ui(rstate, m - ii) + ii;
            q.push_back(pii);
            if (ii == pii) continue;
            unsigned int  bi =  ii / B;
            unsigned int   i =  ii % B;
            unsigned int pbi = pii / B;
            unsigned int  pi = pii % B;
            mat64 & Y = X[bi * n/B + bj0];
            mat64 & pY = X[pbi * n/B + bj0];
            uint64_t c = Y[i] ^ pY[pi];
            Y[i] ^= c;
            pY[pi] ^= c;
        }

        ASSERT_ALWAYS(q.size() == (unsigned int) (ii1 - ii0));

        ple.propagate_permutations(ii1, bj0, q.begin(), q.end());

        for(unsigned int ii = 0 ; ii < m ; ii++) {
            for(unsigned int bj = 0 ; bj < n/B ; bj++) {
                unsigned int  bi =  ii / B;
                unsigned int   i =  ii % B;
                ASSERT_ALWAYS(X[bi * n/B + bj][i] - X[bi * n/B][i] == bj);
            }
        }
    }

    delete[] X;

    return 0;
}/*}}}*/

int test_bblas_level4::test_PLE_move_L_fragments(unsigned int m, unsigned int n)/*{{{*/
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        /* first generate the transpose of the matrix that we will
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
         * output matrix has the expected shape, which is:
         *
         * --------..
         * 0-------..
         * 10------..
         * and so on, with m rows and r columns.
         *
         */

        unsigned int r = gmp_urandomm_ui(rstate, std::min(m, n) + 1);
        std::vector<unsigned int> pivs;
        /* generating binomial(n, r) is like generating with repetition r
         * things among n-(r-1) = n-r+1.  */
        for(unsigned int k = 0 ; k < r ; k++)
            pivs.push_back(gmp_urandomm_ui(rstate, n-r+1));
        std::sort(pivs.begin(), pivs.end());
        for(unsigned int k = 0 ; k < r ; k++)
            pivs[k] += k;
        mat64 * tX = new mat64[(n/B)*(m/B)];
        std::fill_n(tX, (n/B)*(m/B), 0);
        auto ppiv = pivs.begin();
        unsigned int rr = 0;
        for(unsigned int jj = 0 ; jj < n ; jj++) {
            bool is_piv= false;
            unsigned int bj = jj / B;
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
            for(unsigned int bi = 0 ; bi < m/B ; bi++) {
                tX[bi + bj * (m/B)][jj % B] = mpz_get_uint64(NN);
                mpz_fdiv_q_2exp(NN, NN, 64);
            }
            mpz_clear(NN);
        }
        for(unsigned int bi = 0 ; bi < m/B ; bi++) {
            for(unsigned int bj = 0 ; bj < n/B ; bj++) {
                mat64_transpose(X[bi * (n/B) + bj], tX[bi + bj * (m/B)]);
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
            unsigned int br = rr / B;
            unsigned int bp = pivs[rr] / B;
            for( ; rr1 < r && rr1 / B == br && pivs[rr1] / B == bp ; rr1++);
            /* move rr1 - rr columns together */
            std::vector<unsigned int> Q(pivs.begin() + rr, pivs.begin() + rr1);

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
            mpz_t NN;
            mpz_init(NN);
            unsigned int z = ii;
            if (z >= r) z = r;
            mpz_ui_pow_ui(NN, 2, z + ((z ^ ii) & 1));
            mpz_sub_ui(NN,NN,1 + (ii & 1));
            ASSERT_ALWAYS(mpz_divisible_ui_p(NN, 3));
            mpz_divexact_ui(NN,NN,3);
            for(unsigned int bj = 0 ; bj < n/B ; bj++, z -= B) {
                uint64_t c = X[(ii/B)*(n/B)+bj][ii%B];
                c ^= mpz_get_uint64(NN);
                if (z < B) c &= (UINT64_C(1) << z) - 1;
                mpz_fdiv_q_2exp(NN, NN, 64);
                ASSERT_ALWAYS(c == 0);
                if (z < B) break;
            }
            mpz_clear(NN);
        }
    }
    delete[] X;

    return 0;
}/*}}}*/

int test_bblas_level4::test_PLE(unsigned int m, unsigned int n)
{
    const unsigned int B = 64;
    for(unsigned int k = 0 ; k < 100 ; k++) {
        std::vector<mat64> X ((m/B)*(n/B), 0);
        for(unsigned int bi = 0 ; bi < m/B ; bi++) {
            for(unsigned int bj = 0 ; bj < n/B ; bj++) {
                mat64_fill_random(X[bi*(n/B)+bj], rstate);
            }
        }
        std::vector<mat64> Xc = X;

        /* The main importaant thing is the fact of enabling the
         * debug_stuff object, which does invariant checks throughout the
         * PLE computation. As a matter of fact, the final test that
         * we do here is not even needed, since it's already one of those
         * checks triggered by debug_stuff.
         */
        PLE ple(&X[0], m/B, n/B);
        PLE::debug_stuff D(ple);
        std::vector<unsigned int> pivs = ple(&D);

        D.start_check(X);
        D.apply_permutations(pivs);
        unsigned int r = pivs.size();
        auto LL = D.get_LL(r);
        auto UU = D.get_UU(r);
        ASSERT_ALWAYS(D.complete_check(LL, UU));
    }

    return 0;
}



test_bblas_base::tags_t test_bblas_level4::ple_tags { "ple", "l4" };
void test_bblas_level4::ple()
{
    std::vector<std::pair<unsigned int, unsigned int>> mns
    {{
         {64,128},
         {64,64},
         {64,192},
         {128,64},
         {128,128},
         {128,192},
     }};

    for(auto x : mns) {
        test_PLE_find_pivot(x.first, x.second);
        test_PLE_propagate_pivot(x.first, x.second);
        test_PLE_propagate_permutations(x.first, x.second);
        test_PLE_move_L_fragments(x.first, x.second);
        test_PLE(x.first, x.second);   // pass
    }
}
