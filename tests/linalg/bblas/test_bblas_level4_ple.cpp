#include "cado.h"
#include "test_bblas_level4.hpp"
#include "level4_ple_internal.hpp"
#include <cstring>

int test_bblas_level4::test_PLE_find_pivot(unsigned int m, unsigned int n)
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int k = 0 ; k < 10 ; k++) {
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
}

int test_bblas_level4::test_PLE_propagate_pivot(unsigned int m, unsigned int n)
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
}

int test_bblas_level4::test_PLE_propagate_permutations(unsigned int m, unsigned int n)
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple(X, m/B, n/B);

    for(unsigned int k = 1 ; k <= 100 ; k++) {
    for(unsigned int ii = 0, kk = 0 ; ii < m ; ii++) {
        for(unsigned int bj = 0 ; bj < n/B ; bj++, kk++) {
            unsigned int  bi =  ii / B;
            unsigned int   i =  ii % B;
            X[bi * n/B + bj][i] = kk;
        }
    }

    unsigned int ii0 = gmp_urandomm_ui(rstate, m);
    unsigned int ii1 = std::min((mp_limb_t) k, gmp_urandomm_ui(rstate, m + 1 - ii0)) + ii0;
    unsigned int bj0 = gmp_urandomm_ui(rstate, n/B);
    
    unsigned int q[ii1 - ii0];

    unsigned int * q0 = q;
    unsigned int * q1 = q0 + ii1 - ii0;

    for(unsigned int ii = ii0 ; ii < ii1 ; ii++) {
        unsigned int pii = gmp_urandomm_ui(rstate, m - ii) + ii;
        q[ii - ii0] = pii;
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

    ple.propagate_permutations(ii1, bj0, q0, q1);

    for(unsigned int ii = 0 ; ii < m ; ii++) {
        for(unsigned int bj = 0 ; bj < n/B ; bj++) {
            unsigned int  bi =  ii / B;
            unsigned int   i =  ii % B;
            ASSERT_ALWAYS(X[bi * n/B + bj][i] - X[bi * n/B][i] == bj);
        }
    }
    }

    return 0;
}

test_bblas_base::tags_t test_bblas_level4::ple_tags { "ple", "l4" };
void test_bblas_level4::ple()
{
    test_PLE_find_pivot(64, 64);
    test_PLE_propagate_pivot(64, 64);
    test_PLE_propagate_permutations(64, 64);
    test_PLE_find_pivot(128, 192);
    test_PLE_propagate_pivot(128, 192);
    test_PLE_propagate_permutations(128, 192);
}
