#include "cado.h"
#include "test_bpack.hpp"
#include <cstring>
#include <algorithm>

template<typename matrix>
int test_bpack::test_invert_triangular(unsigned int m, unsigned int n)/*{{{*/
{
    bpack<matrix> P(m, n);

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        P.fill_random(rstate);
        P.make_unit_lowertriangular();
    }

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

