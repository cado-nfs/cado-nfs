#include "cado.h"
#include "test_bblas_level4.hpp"
#include "perm_matrix.hpp"
#include <cstring>

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
    mat64 pm[(n/64)*(n/64)];
    perm_matrix_get_matrix(pm, p);

    perm_matrix qt;
    perm_matrix_init(qt, n);
    perm_matrix_transpose(qt, q);

    mat64 qmt[(n/64)*(n/64)];
    perm_matrix_get_matrix(qmt, qt);

    /* compute p*u*transpose(q) */
    mat64 pu[(n/64)*(n/64)];
    memset(pu, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pm[i*(n/64)+k], u[k*(n/64)+j]);
        mat64_add(pu[i*(n/64)+j], pu[i*(n/64)+j], tmp);
    }

    mat64 puq[(n/64)*(n/64)];
    memset(puq, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pu[i*(n/64)+k], qmt[k*(n/64)+j]);
        mat64_add(puq[i*(n/64)+j], puq[i*(n/64)+j], tmp);
    }
    
    /* at this point puq = p*u*transpose(q) should be a upper triangular,
     * with normalized diagonal. */
    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_uppertriangular(puq[i*(n/64)+i]));
    }

    mat64 lm[(n/64)*(n/64)];
    memset(lm, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k <= i ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, l[i*(n/64)+k], m[k*(n/64)+j]);
        mat64_add(lm[i*(n/64)+j], lm[i*(n/64)+j], tmp);
    }

    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_lowertriangular(l[i*(n/64)+i]));
        ASSERT_ALWAYS(mat64_triangular_is_unit(l[i*(n/64)+i]));
        for(int j = 0 ; j < (n/64) ; j++ ) {
            ASSERT_ALWAYS(lm[i*(n/64)+j] == u[i*(n/64)+j]);
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

