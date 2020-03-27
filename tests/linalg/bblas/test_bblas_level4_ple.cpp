#include "cado.h"
#include "test_bblas_level4.hpp"
#include "level4_ple_internal.hpp"
#include <cstring>

int test_bblas_level4::test_PLE_find_pivot(unsigned int m, unsigned int n)
{
    const unsigned int B = 64;
    mat64 * X = new mat64[(m/B)*(n/B)];
    PLE ple((mat64_ptr) X, m/B, n/B);

    for(unsigned int k = 0 ; k < 10 ; k++) {
        std::vector<unsigned int> f;
        std::vector<unsigned int> u;
        for(unsigned int j = 0 ; j < n ; j++) {
            f.push_back(gmp_urandomb_ui(rstate, 30));
            u.push_back(1 + gmp_urandomm_ui(rstate, 3*m-1));
        }
        memset(X, 0, (m/B) * (n/B) * sizeof(mat64));
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

test_bblas_base::tags_t test_bblas_level4::ple_tags { "ple", "l4" };
void test_bblas_level4::ple()
{
    test_PLE_find_pivot(64, 64);
    test_PLE_find_pivot(128, 192);
}
