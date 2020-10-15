#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include "bblas_mat8.hpp"
// IWYU pragma: no_include "bblas_mat64.hpp"
#include <cstdint>                // for uint64_t, uint8_t
#include <utility>                // for pair
#include "bblas_bitmat.hpp"
#include "bpack.hpp"              // for bpack
#include "macros.h"               // for ASSERT_ALWAYS
#include "test_bblas_base.hpp"
#include "test_bpack.hpp"
#include "time_bblas_common.hpp"

template<typename T>
int test_bpack::test_invert_triangular(unsigned int m, unsigned int n)/*{{{*/
{
    if (n > m) return 0;

    bpack<T> P(m, n);
    bpack<T> U(m, n);

    P.fill_random(rstate);
    P.make_unit_lowertriangular();

    U = P;
    U.invert_lower_triangular();

    ASSERT_ALWAYS(U.is_lowertriangular());
    ASSERT_ALWAYS(U.triangular_is_unit());

    /* must now do some sort of mul_lt_ge U*P ; we'll do P*U, in
     * fact, so that we can work on U and not on P.
     */
    bpack<T>::mul_lt_ge(P, U);

    for(unsigned int bi = U.nrowblocks() ; bi-- ; ) {
        for(unsigned int bj = U.ncolblocks() ; bj-- ; ) {
            ASSERT_ALWAYS(U.cell(bi, bj) == (bi == bj));
        }
    }

    auto do_mul_lt_ge = [](bpack<T> & U, bpack<T> const & P) {
        bpack<T>::mul_lt_ge(P, U);
    };
    auto do_ilt = [](bpack<T> & U, bpack<T> const & P) {
        U = P;
        U.invert_lower_triangular();
    };
    TIME1(2, do_ilt, U, P);
    TIME1(2, do_mul_lt_ge, U, P);

    return 0;
}/*}}}*/

template<typename T>
void test_bpack::meta_bpack()
{
    constexpr const unsigned int B = bitmat<T>::width;
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
        if (n > m) continue;
        printf(" -- ILT(m=%u, n=%u)\n", m, n);
        test_invert_triangular<T>(m, n);
    }

}

test_bblas_base::tags_t test_bpack::do_bpack_tags { "bpack", "l4" };
void test_bpack::do_bpack()
{
    meta_bpack<uint64_t>();
    meta_bpack<uint8_t>();
}

