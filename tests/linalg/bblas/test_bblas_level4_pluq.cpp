#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdio>                // for printf
#include <algorithm>              // for fill_n
#include <memory>                 // for allocator_traits<>::value_type
#include <string>                 // for basic_string
#include <vector>                 // for vector
#include "bblas_bitmat.hpp"       // for bitmat<>::vector_type, bitmat
#include "bblas_level3b.hpp"      // for mul_6464_6464
#include "bblas_mat64.hpp"
#include "bblas_level4.hpp"
#include "gmp_aux.h"              // for memfill_random
#include "macros.h"               // for ASSERT_ALWAYS
#include "test_bblas_base.hpp"    // for test_bblas_base, test_bblas_base::t...
#include "test_bblas_level4.hpp"
#include "bblas_perm_matrix.hpp"

/* PLUQ helpers -- well we're not computing exactly PLUQ 
 *
 * PLUQ says: Any m*n matrix A with rank r , can be written A = P*L*U*Q
 * where P and Q are two permutation matrices, of dimension respectively
 * m*m and n*n, L is m*r unit lower triangular and U is r*n upper
 * triangular.
 *
 * Here we compute p,l,u,q such that p*l*a*transpose(q) = an upper
 * triangular matrix (and u is l*a).
 *
 * p*l is lower triangular.
 */
void check_pluq(perm_matrix_ptr p, mat64 * l, mat64 * u, perm_matrix_ptr q, mat64 * m, int n) /*{{{*/
{
    constexpr const unsigned int B = mat64::width;
    mat64::vector_type pm((n/B)*(n/B));
    perm_matrix_get_matrix(pm.data(), p);

    perm_matrix qt;
    perm_matrix_init(qt, n);
    perm_matrix_transpose(qt, q);

    mat64::vector_type qmt((n/B)*(n/B));
    perm_matrix_get_matrix(qmt.data(), qt);

    /* compute p*u*transpose(q) */
    mat64::vector_type pu((n/B)*(n/B));
    std::fill_n(pu.begin(), (n/B)*(n/B), 0);

    for(unsigned int i = 0 ; i < (n/B) ; i++ )
    for(unsigned int j = 0 ; j < (n/B) ; j++ )
    for(unsigned int k = 0 ; k < (n/B) ; k++ ) {
        mat64::addmul(pu[i*(n/B)+j], pm[i*(n/B)+k], u[k*(n/B)+j]);
    }

    mat64::vector_type puq((n/B)*(n/B));
    std::fill_n(puq.data(), (n/B)*(n/B), 0);

    for(unsigned int i = 0 ; i < (n/B) ; i++ )
    for(unsigned int j = 0 ; j < (n/B) ; j++ )
    for(unsigned int k = 0 ; k < (n/B) ; k++ ) {
        mul_6464_6464(puq[i*(n/B)+j], pu[i*(n/B)+k], qmt[k*(n/B)+j]);
    }
    
    /* at this point puq = p*u*transpose(q) should be a upper triangular,
     * with normalized diagonal. */
    for(unsigned int i = 0 ; i < (n/B) ; i++ ) {
        ASSERT_ALWAYS(puq[i*(n/B)+i].is_uppertriangular());
    }

    mat64::vector_type lm((n/B)*(n/B));
    std::fill_n(lm.data(), (n/B)*(n/B), 0);

    for(unsigned int i = 0 ; i < (n/B) ; i++ )
    for(unsigned int j = 0 ; j < (n/B) ; j++ )
    for(unsigned int k = 0 ; k <= i ; k++ ) {
        mat64::addmul(lm[i*(n/B)+j], l[i*(n/B)+k], m[k*(n/B)+j]);
    }

    for(unsigned int i = 0 ; i < (n/B) ; i++ ) {
        ASSERT_ALWAYS(l[i*(n/B)+i].is_lowertriangular());
        ASSERT_ALWAYS(l[i*(n/B)+i].triangular_is_unit());
        for(unsigned int j = 0 ; j < (n/B) ; j++ ) {
            ASSERT_ALWAYS(lm[i*(n/B)+j] == u[i*(n/B)+j]);
        }
    }
    perm_matrix_clear(qt);
}
/*}}}*/

test_bblas_base::tags_t test_bblas_level4::pluq_tags { "pluq", "l4" };/*{{{*/
void test_bblas_level4::pluq() {
    mat64 m[4],l[4],u[4];
    {
        perm_matrix p,q;
        perm_matrix_init(p, 64);
        perm_matrix_init(q, 64);
        memfill_random(m, sizeof(mat64), rstate);
        PLUQ64(p,l[0],u[0],q,m[0]);
        check_pluq(p,l,u,q,m,64);
        perm_matrix_clear(p);
        perm_matrix_clear(q);
    }

    {
        perm_matrix p,q;
        perm_matrix_init(p, 128);
        perm_matrix_init(q, 128);
        memfill_random(m, 4 * sizeof(mat64), rstate);
        PLUQ128(p,l,u,q,m);
        check_pluq(p,l,u,q,m,128);
        perm_matrix_clear(p);
        perm_matrix_clear(q);
    }

    printf("PLUQ128 stub executed ok\n");

    // perm_matrix_6464(m[0]); printf("\n");
    // perm_matrix_6464(u); printf("\n");
    // perm_matrix_6464(l); printf("\n");
    // perm_matrix_6464(p); printf("\n");
    // perm_matrix_6464(t); printf("\n");


    // mul_N64_T6464(t,u,p,64);
    // perm_matrix_6464(t); printf("\n");

#if 0
    memset(e,0,sizeof(e));
    memfill_random(m[0], (64) * sizeof(uint32_t), rstate);
    memfill_random(m[1], (64) * sizeof(uint32_t), rstate);
    memfill_random(m[2], (64) * sizeof(uint32_t), rstate);
    memfill_random(m[3], (64) * sizeof(uint32_t), rstate);
    // perm_matrix_6464(m[0]); printf("\n");
    full_echelon_6464_imm(mm,e[0],m[0]);

    /* 
       mul_6464

       perm_matrix_6464(mm); printf("\n");
       perm_matrix_6464(e); printf("\n"); printf("\n");
       */
#endif
}/*}}}*/

