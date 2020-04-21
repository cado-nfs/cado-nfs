#include "cado.h"
#include "test_bpack.hpp"
#include "time_bblas_common.hpp"
#include <cstring>
#include <algorithm>

template<typename matrix>
int test_bpack::test_invert_triangular(unsigned int m, unsigned int n)/*{{{*/
{
    if (n > m) return 0;

    bpack<matrix> P(m, n);
    bpack<matrix> U(m, n);

        P.fill_random(rstate);
        P.make_unit_lowertriangular();

        U = P;
        U.invert_lower_triangular();

        ASSERT_ALWAYS(U.is_lowertriangular());
        ASSERT_ALWAYS(U.triangular_is_unit());

        /* must now do some sort of mul_lt_ge U*P ; we'll do P*U, in
         * fact, so that we can work on U and not on P.
         */
        for(unsigned int bi = U.nrowblocks() ; bi-- ; ) {
            for(unsigned int bj = U.ncolblocks() ; bj-- ; ) {
                /* (PU)_{i,j} is the sum for j<=k<=i of P_{i,k}U_{k,j}
                 * except for bi>=ncolblocks, where it is:
                 *
                 * U_{i,j}+sum for j<=k<ncolblocks P_{i,k}U_{k,j}
                 *
                 */
                matrix S;
                if (bi >= U.ncolblocks())
                    S = U.cell(bi, bj);
                else
                    matrix::mul(S, P.cell(bi, bi), U.cell(bi, bj));
                for(unsigned int bk = bj ; bk < bi && bk < U.ncolblocks() ; bk++) {
                    matrix::addmul(S, P.cell(bi, bk), U.cell(bk, bj));
                }
                ASSERT_ALWAYS(S == (bi == bj));
            }
        }

        auto do_ilp = [](bpack<matrix> & U, bpack<matrix> const & P) {
            U = P;
            U.invert_lower_triangular();
        };
        TIME1(2, do_ilp, U, P);

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
        if (n > m) continue;
        printf(" -- ILP(m=%u, n=%u)\n", m, n);
        test_invert_triangular<matrix>(m, n);
    }

}

test_bblas_base::tags_t test_bpack::do_bpack_tags { "bpack", "l4" };
void test_bpack::do_bpack()
{
    meta_bpack<mat64>();
    meta_bpack<mat8>();
}

