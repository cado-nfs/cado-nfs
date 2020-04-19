#include "cado.h"
#include "test_bpack.hpp"
#include <cstring>
#include <algorithm>

template<typename matrix>
int test_bpack::test_invert_triangular(unsigned int m, unsigned int n)/*{{{*/
{
    constexpr const unsigned int B = matrix::width;
    // typedef typename matrix::datatype U;
    matrix * X = matrix::alloc((m/B)*(n/B));
    bpack<matrix> P(X, m/B, n/B);

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        memfill_random((void *) X, (m/B) * (n/B) * sizeof(matrix), rstate);
        for(unsigned int j = 0 ; j < n ; j++) {
            for(unsigned int i = 0 ; i < m ; i++) {
                if (j > i)
                    X[(i/B)*(n/B)+(j/B)] = 0;
                else if (j == i) {
                    // X[(i/B)*(n/B)+(j/B)].make_unit_lower_triangular();
                }
            }
        }
    }
    matrix::free(X);

    return 0;
}/*}}}*/

template<typename matrix>
void test_bpack::meta_bpack()
{
    constexpr const unsigned int B = matrix::width;
    std::vector<std::pair<unsigned int, unsigned int>> mns
    {{
         {B,B},
         {B,2*B},
         {B,3*B},
         {2*B,B},
         {2*B,2*B},
         {2*B,3*B},
         {2*B,10*B},
         {10*B,2*B},
         {2*B,100*B},
         {100*B,2*B},
         {2*B,1000*B},
         {1000*B,2*B},
         {32*B,32*B},
     }};

    for(auto x : mns) {
        unsigned int m = x.first;
        unsigned int n = x.second;
        test_invert_triangular<matrix>(m, n);
    }

}

test_bblas_base::tags_t test_bpack::do_bpack_tags { "bpack", "l4" };
void test_bpack::do_bpack()
{
    meta_bpack<mat64>();
    meta_bpack<mat8>();
}

