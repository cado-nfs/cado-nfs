#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <limits.h>     /* for INT_MAX */
#include <ostream>
#include "mpz_mat.h"
#include "portability.h"
#include "utils.h"

/*{{{ entry access*/
mpz_ptr mpz_mat_entry(mpz_mat_ptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpz_srcptr mpz_mat_entry_const(mpz_mat_srcptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpq_ptr mpq_mat_entry(mpq_mat_ptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}

mpq_srcptr mpq_mat_entry_const(mpq_mat_srcptr M, unsigned int i, unsigned int j)
{
    return M->x[i * M->n + j];
}
/*}}}*/
/*{{{ init/clear/realloc*/
void mpz_mat_init(mpz_mat_ptr M, unsigned int m, unsigned int n)
{
    M->x = (mpz_t*) ((m && n) ? malloc(m * n * sizeof(mpz_t)) : NULL);

    M->m = m;
    M->n = n;
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_init(mpz_mat_entry(M, i, j));
}

void mpz_mat_clear(mpz_mat_ptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_clear(mpz_mat_entry(M, i, j));
    free(M->x);
}

void mpz_mat_realloc(mpz_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n) return;
    mpz_mat_clear(M);
    mpz_mat_init(M, m, n);
}

void mpq_mat_init(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    M->x = (mpq_t*) ((m && n) ? malloc(m * n * sizeof(mpq_t)) : NULL);
    M->m = m;
    M->n = n;
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_init(mpq_mat_entry(M, i, j));
}

void mpq_mat_clear(mpq_mat_ptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_clear(mpq_mat_entry(M, i, j));
    free(M->x);
}

void mpq_mat_realloc(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n) return;
    mpq_mat_clear(M);
    mpq_mat_init(M, m, n);
}

/*}}}*/
/*{{{ operations on submatrices, and swaps*/
void mpz_mat_submat_swap(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpz_swap(
                    mpz_mat_entry(A0, i0 + i, j0 + j),
                    mpz_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpz_mat_submat_set(
        mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpz_mat_srcptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpz_set(
                    mpz_mat_entry(A0, i0 + i, j0 + j),
                    mpz_mat_entry_const(A1, i1 + i, j1 + j));
        }
    }
}

void mpq_mat_swap(mpz_mat_ptr A, mpz_mat_ptr B)
{
    mpz_mat C;
    memcpy(C, A, sizeof(mpz_mat));
    memcpy(A, B, sizeof(mpz_mat));
    memcpy(B, C, sizeof(mpz_mat));
}

void mpq_mat_submat_swap(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_ptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpq_swap(
                    mpq_mat_entry(A0, i0 + i, j0 + j),
                    mpq_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpq_mat_submat_set(
        mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
        mpq_mat_srcptr A1, unsigned int i1, unsigned int j1,
        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for(unsigned int i = 0 ; i < dm ; i++) {
        for(unsigned int j = 0 ; j < dn ; j++) {
            mpq_set(
                    mpq_mat_entry(A0, i0 + i, j0 + j),
                    mpq_mat_entry_const(A1, i1 + i, j1 + j));
        }
    }
}

void mpz_mat_swap(mpz_mat_ptr A, mpz_mat_ptr B)
{
    mpz_mat C;
    memcpy(C, A, sizeof(mpz_mat));
    memcpy(A, B, sizeof(mpz_mat));
    memcpy(B, C, sizeof(mpz_mat));
}

void mpq_mat_swap(mpq_mat_ptr A, mpq_mat_ptr B)
{
    mpq_mat C;
    memcpy(C, A, sizeof(mpq_mat));
    memcpy(A, B, sizeof(mpq_mat));
    memcpy(B, C, sizeof(mpq_mat));
}

/*}}}*/
/*{{{ set/set_ui*/
void mpz_mat_set(mpz_mat_ptr dst, mpz_mat_srcptr src)
{
    mpz_mat_realloc(dst, src->m, src->n);
    for(unsigned int i = 0 ; i < src->m ; i++)
        for(unsigned int j = 0 ; j < src->n ; j++)
            mpz_set(mpz_mat_entry(dst, i, j), mpz_mat_entry_const(src, i, j));
}

void mpz_mat_set_ui(mpz_mat_ptr M, unsigned long a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_set_ui(mpz_mat_entry(M, i, j), i == j ? a : 0);
}

void mpz_mat_set_mpz(mpz_mat_ptr M, mpz_srcptr a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        mpz_set(mpz_mat_entry(M, i, i), a);
    }
}

void mpz_mat_add_ui(mpz_mat_ptr M, unsigned long a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        mpz_ptr mii = mpz_mat_entry(M, i, i);
        mpz_add_ui(mii, mii, a);
    }
}

void mpz_mat_add_mpz(mpz_mat_ptr M, mpz_srcptr a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        mpz_ptr mii = mpz_mat_entry(M, i, i);
        mpz_add(mii, mii, a);
    }
}

void mpq_mat_set_ui(mpq_mat_ptr M, unsigned long a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpq_set_ui(mpq_mat_entry(M, i, j), i == j ? a : 0, 1);
}

void mpq_mat_set(mpq_mat_ptr dst, mpq_mat_srcptr src)
{
    mpq_mat_realloc(dst, src->m, src->n);
    for(unsigned int i = 0 ; i < src->m ; i++)
        for(unsigned int j = 0 ; j < src->n ; j++)
            mpq_set(mpq_mat_entry(dst, i, j), mpq_mat_entry_const(src, i, j));
}

void mpz_mat_urandomm(mpz_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++)
            mpz_urandomm(mpz_mat_entry(M, i, j), state, p);
}

void mpq_mat_urandomm(mpq_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for(unsigned int i = 0 ; i < M->m ; i++)
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_urandomm(mpq_numref(mpq_mat_entry(M, i, j)), state, p);
            mpz_set_ui(mpq_denref(mpq_mat_entry(M, i, j)), 1);
        }
}



/*}}}*/
/*{{{ Joining matrices */
void mpz_mat_vertical_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2)
{
    if (N == M1 || N == M2) {
        mpz_mat Nx;
        mpz_mat_init(Nx, 0, 0);
        mpz_mat_vertical_join(Nx, M1, M2);
        mpz_mat_swap(Nx, N);
        mpz_mat_clear(Nx);
        return;
    }
    ASSERT_ALWAYS(M1->n == M2->n);
    mpz_mat_realloc(N,M1->m+M2->m,M1->n);
    mpz_mat_submat_set(N,0,0,M1,0,0,M1->m,M1->n);
    mpz_mat_submat_set(N,M1->m,0,M2,0,0,M2->m,M2->n);
}
void mpq_mat_vertical_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2)
{
    if (N == M1 || N == M2) {
        mpq_mat Nx;
        mpq_mat_init(Nx, 0, 0);
        mpq_mat_vertical_join(Nx, M1, M2);
        mpq_mat_swap(Nx, N);
        mpq_mat_clear(Nx);
        return;
    }
    ASSERT_ALWAYS(M1->n == M2->n);
    mpq_mat_realloc(N,M1->m+M2->m,M1->n);
    mpq_mat_submat_set(N,0,0,M1,0,0,M1->m,M1->n);
    mpq_mat_submat_set(N,M1->m,0,M2,0,0,M2->m,M2->n);
}
void mpz_mat_horizontal_join(mpz_mat_ptr N, mpz_mat_srcptr M1, mpz_mat_srcptr M2)
{
    if (N == M1 || N == M2) {
        mpz_mat Nx;
        mpz_mat_init(Nx, 0, 0);
        mpz_mat_horizontal_join(Nx, M1, M2);
        mpz_mat_swap(Nx, N);
        mpz_mat_clear(Nx);
        return;
    }
    ASSERT_ALWAYS(M1->m == M2->m);
    mpz_mat_realloc(N,M1->m,M1->n + M2->n);
    mpz_mat_submat_set(N,0,0,M1,0,0,M1->m,M1->n);
    mpz_mat_submat_set(N,0,M1->n,M2,0,0,M2->m,M2->n);
}
void mpq_mat_horizontal_join(mpq_mat_ptr N, mpq_mat_srcptr M1, mpq_mat_srcptr M2)
{
    if (N == M1 || N == M2) {
        mpq_mat Nx;
        mpq_mat_init(Nx, 0, 0);
        mpq_mat_horizontal_join(Nx, M1, M2);
        mpq_mat_swap(Nx, N);
        mpq_mat_clear(Nx);
        return;
    }
    ASSERT_ALWAYS(M1->m == M2->m);
    mpq_mat_realloc(N,M1->m,M1->n + M2->n);
    mpq_mat_submat_set(N,0,0,M1,0,0,M1->m,M1->n);
    mpq_mat_submat_set(N,0,M1->n,M2,0,0,M2->m,M2->n);
}
/*}}}*/
/*{{{ determinant, trace and transposition */
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M)/*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(t, 0);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_add(t, t, mpz_mat_entry_const(M,i,i));
}
/*}}}*/
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M)/*{{{*/
{
    // We assume that M is triangular
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(d, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpz_mul(d, d, mpz_mat_entry_const(M,i,i));
}
/*}}}*/
void mpz_mat_transpose(mpz_mat_ptr D, mpz_mat_srcptr M)/*{{{*/
{
    if (D != M) {
        mpz_mat_realloc(D,M->n,M->m);
        for(unsigned int i = 0 ; i < M->m ; i++){
            for(unsigned int j = 0 ; j < M->n ; j++){
                mpz_srcptr mij = mpz_mat_entry_const(M,i,j);
                mpz_ptr dji = mpz_mat_entry(D,j,i);
                mpz_set(dji, mij);
            }
        }
    } else if (M->m != M->n) {
        /* transpose a rectangular matrix in place. Rather annoying to do
         * with real swaps, right ? */
        mpz_mat Mc;
        mpz_mat_init(Mc, 0, 0);
        mpz_mat_set(Mc, M);
        mpz_mat_transpose(D, Mc);
        mpz_mat_clear(Mc);
    } else {
        for(unsigned int i = 0 ; i < M->m ; i++){
            for(unsigned int j = i + 1 ; j < M->n ; j++){
                mpz_ptr dij = mpz_mat_entry(D,i,j);
                mpz_ptr dji = mpz_mat_entry(D,j,i);
                mpz_swap(dji, dij);
            }
        }
    }
}/*}}}*/

void mpz_mat_reverse_rows(mpz_mat_ptr B, mpz_mat_srcptr A)
{
    if (A != B) {
        mpz_mat_realloc(B, A->m,A->n);
        for(unsigned int i = 0 ; i < A->m ; i++){
            mpz_mat_submat_set(B, i, 0, A, A->m - 1 - i, 0, 1, A->n);
        }
    } else {
        for(unsigned int i = 0 ; i < A->m - 1 - i ; i++){
            mpz_mat_submat_swap(B, i, 0, B, A->m - 1 - i, 0, 1, A->n);
        }
    }
}

void mpz_mat_reverse_columns(mpz_mat_ptr B, mpz_mat_srcptr A)
{
    if (A != B) {
        mpz_mat_realloc(B, A->m,A->n);
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpz_mat_submat_set(B, 0, j, A, 0, A->n - 1 - j, A->m, 1);
        }
    } else {
        for(unsigned int j = 0 ; j < A->n - 1 - j ; j++){
            mpz_mat_submat_swap(B, 0, j, B, 0, A->n - 1 - j, A->m, 1);
        }
    }
}

void mpq_mat_trace(mpq_ptr t, mpq_mat_srcptr M)/*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(t, 0, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpq_add(t, t, mpq_mat_entry_const(M,i,i));
}
/*}}}*/
void mpq_mat_determinant_triangular(mpq_ptr d, mpq_mat_srcptr M)/*{{{*/
{
    // We assume that M is triangular
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(d, 1, 1);
    for(unsigned int i = 0 ; i < M->n ; i++)
        mpq_mul(d, d, mpq_mat_entry_const(M,i,i));
}
/*}}}*/
void mpq_mat_transpose(mpq_mat_ptr D, mpq_mat_srcptr M)/*{{{*/
{
    if (D != M) {
        mpq_mat_realloc(D,M->n,M->m);
        for(unsigned int i = 0 ; i < M->m ; i++){
            for(unsigned int j = 0 ; j < M->n ; j++){
                mpq_srcptr mij = mpq_mat_entry_const(M,i,j);
                mpq_ptr dji = mpq_mat_entry(D,j,i);
                mpq_set(dji, mij);
            }
        }
    } else if (M->m != M->n) {
        /* transpose a rectangular matrix in place. Rather annoying to do
         * with real swaps, right ? */
        mpq_mat Mc;
        mpq_mat_init(Mc, 0, 0);
        mpq_mat_set(Mc, M);
        mpq_mat_transpose(D, Mc);
        mpq_mat_clear(Mc);
    } else {
        for(unsigned int i = 0 ; i < M->m ; i++){
            for(unsigned int j = i + 1 ; j < M->n ; j++){
                mpq_ptr dij = mpq_mat_entry(D,i,j);
                mpq_ptr dji = mpq_mat_entry(D,j,i);
                mpq_swap(dji, dij);
            }
        }
    }
}/*}}}*/
/*}}}*/
/*{{{ miscellaneous */

/* convert to integer matrix divided by lcm of denominator.
 * return 1
 * if den==NULL, return 0 if denominator happens to not be 1 (in which
 * case the matrix returned is undefined).
 */

int mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M)/*{{{*/
{
    int ret = 1;
    if (den) mpz_set_ui(den, 1);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_srcptr Mij =  mpq_denref(mpq_mat_entry_const(M, i, j));
            if (den) {
                mpz_lcm(den, den, Mij);
            } else if (mpz_cmp_ui(Mij, 1) != 0) {
                ret = 0;
            }
        }
    }
    if (!ret) return ret;
    mpz_mat_realloc(num, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_ptr dst = mpz_mat_entry(num, i, j);
            mpq_srcptr src = mpq_mat_entry_const(M, i, j);
            if (den) {
                mpz_divexact(dst, den, mpq_denref(src));
                mpz_mul(dst, dst, mpq_numref(src));
            } else {
                mpz_set(dst, mpq_numref(src));
            }
        }
    }
    return ret;
}
/*}}}*/
void mpq_mat_set_mpz_mat(mpq_mat_ptr N, mpz_mat_srcptr M)/*{{{*/
{
    mpq_mat_realloc(N,M->m,M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpq_set_z(mpq_mat_entry(N,i,j),mpz_mat_entry_const(M,i,j));
        }
    }
}
/*}}}*/
void mpq_mat_set_mpz_mat_denom(mpq_mat_ptr N, mpz_mat_srcptr M, mpz_srcptr d)/*{{{*/
{
    mpq_mat_realloc(N,M->m,M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpq_ptr nij = mpq_mat_entry(N,i,j);
            mpz_srcptr mij = mpz_mat_entry_const(M,i,j);
            mpz_set(mpq_numref(nij), mij);
            mpz_set(mpq_denref(nij), d);
            mpq_canonicalize(nij);
        }
    }
}
/*}}}*/
void mpz_mat_mod_ui(mpz_mat_ptr dst, mpz_mat_srcptr src, unsigned long p)/*{{{*/
{
    mpz_mat_realloc(dst, src->m, src->n); /* no-op if dst == src */
    for (unsigned int i = 0 ; i < src->m ; i++){
        for (unsigned int j = 0 ; j < src->n ; j++){
            mpz_fdiv_r_ui(mpz_mat_entry(dst,i,j),mpz_mat_entry_const(src,i,j),p);
        }
    }
}/*}}}*/
void mpz_mat_mod_mpz(mpz_mat_ptr dst, mpz_mat_srcptr src, mpz_srcptr p)/*{{{*/
{
    mpz_mat_realloc(dst, src->m, src->n); /* no-op if dst == src */
    for (unsigned int i = 0 ; i < src->m ; i++){
        for (unsigned int j = 0 ; j < src->n ; j++){
            mpz_fdiv_r(mpz_mat_entry(dst,i,j),mpz_mat_entry_const(src,i,j),p);
        }
    }
}/*}}}*/

int mpz_mat_p_valuation(mpz_mat_srcptr A, mpz_srcptr p)
{
    int val = INT_MAX;
    mpz_t c;
    mpz_init(c);
    for (unsigned int i = 0 ; val && i < A->m ; i++){
        for (unsigned int j = 0 ; val && j < A->n ; j++){
            int v = 0;
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            if (mpz_size(aij) == 0) continue;
            mpz_set(c, aij);
            for( ; v < val && mpz_divisible_p(c, p) ; v++)
                mpz_fdiv_q(c, c, p);
            val = v;
        }
    }
    mpz_clear(c);
    return val;
}

int mpz_mat_p_valuation_ui(mpz_mat_srcptr A, unsigned long p)
{
    int val = INT_MAX;
    mpz_t c;
    mpz_init(c);
    for (unsigned int i = 0 ; val && i < A->m ; i++){
        for (unsigned int j = 0 ; val && j < A->n ; j++){
            int v = 0;
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            if (mpz_size(aij) == 0) continue;
            mpz_set(c, aij);
            for( ; v < val && mpz_divisible_ui_p(c, p) ; v++)
                mpz_fdiv_q_ui(c, c, p);
            val = v;
        }
    }
    mpz_clear(c);
    return val;
}

/*}}}*/

/* {{{ row-level operations (for Gauss and friends, mostly) */
// Return 1 if the k-th line of M is null, 0 else
int mpz_mat_isnull_row(mpz_mat_srcptr M, unsigned int k){/*{{{*/
    unsigned int j = 0;
    ASSERT_ALWAYS(k < M->m);
    while((j < M->n) && !mpz_cmp_si(mpz_mat_entry_const(M,k,j),0)){
        j++;
    }
    if(j == M->n){
        return 1;
    }
    return 0;
}
/*}}}*/
void mpz_mat_swaprows(mpz_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_swap(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
void mpq_mat_swaprows(mpq_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_swap(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/
void mpz_mat_rotatedownrows(mpz_mat_ptr M, unsigned int i0, unsigned int k)/*{{{*/
{
    /* apply a circular shift on rows [i0...i0+k[ */
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i0 + k <= M->m);
    if (k <= 1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        for(unsigned int s = k - 1 ; s-- ; ) {
            mpz_swap(mpz_mat_entry(M, i0 + s, j), mpz_mat_entry(M, i0 + s + 1, j));
        }
    }
}
/*}}}*/
void mpq_mat_rotatedownrows(mpq_mat_ptr M, unsigned int i0, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i0 + k <= M->m);
    if (k <= 1) return;
    for(unsigned int j = 0 ; j < M->n ; j++) {
        for(unsigned int s = k - 1 ; s-- ; ) {
            mpq_swap(mpq_mat_entry(M, i0 + s, j), mpq_mat_entry(M, i0 + s + 1, j));
        }
    }
}
/*}}}*/
void mpz_mat_permuterows(mpz_mat_ptr M, unsigned int * perm)/*{{{*/
{
    /* put row perm[k] in row k */
    mpz_mat Mx;
    mpz_mat_init(Mx, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpz_swap(mpz_mat_entry(Mx, i, j),
                    mpz_mat_entry(M, perm[i], j));
        }
    }
    mpz_mat_swap(M, Mx);
    mpz_mat_clear(Mx);
}
/*}}}*/
void mpq_mat_permuterows(mpq_mat_ptr M, unsigned int * perm)/*{{{*/
{
    mpq_mat Mx;
    mpq_mat_init(Mx, M->m, M->n);
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            mpq_swap(mpq_mat_entry(Mx, i, j),
                    mpq_mat_entry(M, perm[i], j));
        }
    }
    mpq_mat_swap(M, Mx);
    mpq_mat_clear(Mx);
}
/*}}}*/
static int permutation_signature(const unsigned int * perm, unsigned int n)/*{{{*/
{
    /* can this be done in O(n) only ?? */
    int sign = 1;
    for(unsigned int i = 0 ; i < n ; i++)
        for(unsigned int j = i ; j < n ; j++)
            if (perm[j] < perm[i]) sign*=-1;
    return sign;
}
/*}}}*/

/* add lambda times row i1 to row i0 */
void mpz_mat_addmulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_addmul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
void mpz_mat_addmulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_addmul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}
/*}}}*/
/* add lambda times row i1 to row i0 */
void mpq_mat_addmulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        /* we have no mpq_addmul, unfortunately */
        mpq_mul(tmp,
                lambda,
                mpq_mat_entry(M, i1, j));
        mpq_add(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                tmp);
    }
    mpq_clear(tmp);
}
/*}}}*/
/* subtract lambda times row i1 to row i0 */
void mpz_mat_submulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
void mpz_mat_submulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1, mpz_srcptr lambda, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j),
                lambda,
                mpz_mat_entry(M, i1, j));
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}
/*}}}*/
/* subtract lambda times row i1 to row i0 */
void mpq_mat_submulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1, mpq_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        /* we have no mpq_submul, unfortunately */
        mpq_mul(tmp,
                lambda,
                mpq_mat_entry(M, i1, j));
        mpq_sub(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                tmp);
    }
    mpq_clear(tmp);
}
/*}}}*/

/* add row i1 to row i0 */
void mpz_mat_addrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_add(mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* add row i1 to row i0 */
void mpq_mat_addrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_add(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/

/* subtract row i1 to row i0 */
void mpz_mat_subrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_sub(mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* subtract row i1 to row i0 */
void mpq_mat_subrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_sub(mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* multiply row i0 by lambda */
void mpz_mat_mulrow(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
    }
}
/*}}}*/
void mpz_mat_mulrow_mod(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}
/*}}}*/
/* multiply row i0 by lambda */
void mpq_mat_mulrow(mpq_mat_ptr M, unsigned int i0, mpq_srcptr lambda)/*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for(unsigned int j = 0 ; j < M->n ; j++) {
        mpq_mul(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j), lambda);
    }
}
/*}}}*/
/* this computes an additive combination of n rows into row [didx] of the
 * initial matrix. We require that this destination row be cleared
 * initially.
 */
void mpz_mat_combinerows(mpz_mat_ptr M, unsigned int didx, unsigned int sidx,/*{{{*/
        mpz_srcptr * lambda, unsigned int n)
{
    for(unsigned int j = 0 ; j < M->n ; j++) {
        ASSERT_ALWAYS(mpz_cmp_ui(mpz_mat_entry(M, didx, j), 0) == 0);
    }
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_mat_addmulrow(M, didx, sidx + i, lambda[i]);
    }
}
/*}}}*/
/* }}} */
/*{{{ I/O*/
void mpz_mat_fprint(FILE * stream, mpz_mat_srcptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            gmp_fprintf(stream, " %Zd", mpz_mat_entry_const(M, i, j));
        }
        fprintf(stream, "\n");
    }
}

void mpq_mat_fprint(FILE * stream, mpq_mat_srcptr M)
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            gmp_fprintf(stream, " %Qd", mpq_mat_entry_const(M, i, j));
        }
        fprintf(stream, "\n");
    }
}

void mpq_mat_fprint_as_mpz(FILE* f, mpq_mat_srcptr M)
{
    mpz_t denom;
    mpz_mat N;

    mpz_init(denom);
    mpz_mat_init(N,0,0);

    mpq_mat_numden(N, denom, M);
    
    mpz_mat_fprint(f,N);
    gmp_fprintf(f, "Denominator is : %Zd\n",denom);
    
    mpz_mat_clear(N);
    mpz_clear(denom);
}
/*}}}*/
/*{{{ comparison */
int mpz_mat_cmp(mpz_mat_srcptr M, mpz_mat_srcptr N)/*{{{*/
{
    /* shall we abort or return something for incompatible matrices ?? */
    if (M->m != N->m) return M->m - N->m;
    if (M->n != N->n) return M->n - N->n;
    // ASSERT_ALWAYS((M->m == N->m) && (M->n == N->n));
    unsigned int m = M->m;
    unsigned int n = M->n;
    unsigned int i,j;
    for (i = 0 ; i < m ; i++){
        for (j = 0 ; j < n ; j++){
            int k = mpz_cmp(mpz_mat_entry_const(M,i,j),mpz_mat_entry_const(N,i,j));
            if(k!=0) return k;
        }
    }
    return 0;
}
/*}}}*/
int mpq_mat_cmp(mpq_mat_srcptr M, mpq_mat_srcptr N)/*{{{*/
{
    /* shall we abort or return something for incompatible matrices ?? */
    if (M->m != N->m) return M->m - N->m;
    if (M->n != N->n) return M->n - N->n;
    // ASSERT_ALWAYS((M->m == N->m) && (M->n == N->n));
    unsigned int m = M->m;
    unsigned int n = M->n;
    unsigned int i,j;
    for (i = 0 ; i < m ; i++){
        for (j = 0 ; j < n ; j++){
            int k = mpq_cmp(mpq_mat_entry_const(M,i,j),mpq_mat_entry_const(N,i,j));
            if(k!=0) return k;
        }
    }
    return 0;
}
/*}}}*/

/* TODO (perhaps) :
 * mp[qz]_mat_is_{upper,lower}_triangular
 * mp[qz]_mat_is_diagonal
 */
/*}}}*/
/*{{{ add and subtract */
void mpz_mat_add(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpz_mat_realloc(C,A->m,A->n);       /* no-op if C == A or C == B */
    
    for(unsigned int i = 0 ; i < A->m ; i++){
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            mpz_srcptr bij = mpz_mat_entry_const(B, i, j);
            mpz_ptr cij = mpz_mat_entry(C, i, j);
            mpz_add(cij,aij,bij);
        }
    }
}
/* }}} */
void mpq_mat_add(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpq_mat_realloc(C,A->m,A->n);       /* no-op if C == A or C == B */
    
    for(unsigned int i = 0 ; i < A->m ; i++){
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_srcptr bij = mpq_mat_entry_const(B, i, j);
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_add(cij,aij,bij);
        }
    }
}
/*}}}*/
void mpz_mat_sub(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpz_mat_realloc(C,A->m,A->n);       /* no-op if C == A or C == B */
    
    for(unsigned int i = 0 ; i < A->m ; i++){
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            mpz_srcptr bij = mpz_mat_entry_const(B, i, j);
            mpz_ptr cij = mpz_mat_entry(C, i, j);
            mpz_sub(cij,aij,bij);
        }
    }
}
/* }}} */
void mpq_mat_sub(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpq_mat_realloc(C,A->m,A->n);       /* no-op if C == A or C == B */
    
    for(unsigned int i = 0 ; i < A->m ; i++){
        for(unsigned int j = 0 ; j < A->n ; j++){
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_srcptr bij = mpq_mat_entry_const(B, i, j);
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_sub(cij,aij,bij);
        }
    }
}
/*}}}*/
/*}}}*/
/*{{{ multiplication */
void mpz_mat_mul(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->n == B->m);
    if (C == A || C == B) {
        mpz_mat D;
        mpz_mat_init(D, A->m, B->n);
        mpz_mat_mul(D, A, B);
        mpz_mat_swap(D, C);
        mpz_mat_clear(D);
        return;
    }
    mpz_mat_realloc(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j),
                        mpz_mat_entry_const(A, i, k),
                        mpz_mat_entry_const(B, k, j));
            }
        }
    }
}/*}}}*/

void mpz_mat_mul_mod_ui(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B, unsigned long p)/*{{{*/
{
    ASSERT_ALWAYS(A->n == B->m);
    if (C == A || C == B) {
        mpz_mat D;
        mpz_mat_init(D, A->m, B->n);
        mpz_mat_mul_mod_ui(D, A, B, p);
        mpz_mat_swap(D, C);
        mpz_mat_clear(D);
        return;
    }
    mpz_mat_realloc(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j),
                        mpz_mat_entry_const(A, i, k),
                        mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_mod_ui(C,C,p);
}/*}}}*/

void mpz_mat_mul_mod_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(A->n == B->m);
    if (C == A || C == B) {
        mpz_mat D;
        mpz_mat_init(D, A->m, B->n);
        mpz_mat_mul_mod_mpz(D, A, B, p);
        mpz_mat_swap(D, C);
        mpz_mat_clear(D);
        return;
    }
    mpz_mat_realloc(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j),
                        mpz_mat_entry_const(A, i, k),
                        mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_mod_mpz(C,C,p);
}/*}}}*/

void mpz_mat_mul_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, mpz_ptr k)/*{{{*/
{
    mpz_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_mul(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
}/*}}}*/

void mpz_mat_divexact_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, mpz_srcptr k)/*{{{*/
{
    mpz_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_divexact(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
}/*}}}*/

void mpz_mat_divexact_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long k)/*{{{*/
{
    mpz_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_divexact_ui(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
}/*}}}*/

void mpz_mat_mul_si(mpz_mat_ptr B, mpz_mat_srcptr A, long k)/*{{{*/
{
    mpz_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_mul_si(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
}/*}}}*/

void mpz_mat_mul_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long k)/*{{{*/
{
    mpz_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpz_mul_ui(mpz_mat_entry(B, i, j),mpz_mat_entry_const(A,i,j),k);
        }
    }
}/*}}}*/


void mpq_mat_mul(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B)/*{{{*/
{
    ASSERT_ALWAYS(A->n == B->m);
    if (C == A || C == B) {
        mpq_mat D;
        mpq_mat_init(D, A->m, B->n);
        mpq_mat_mul(D, A, B);
        mpq_mat_swap(D, C);
        mpq_mat_clear(D);
        return;
    }
    mpq_t z;
    mpq_init(z);
    mpq_mat_realloc(C, A->m, B->n);
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_set_ui(cij, 0, 1);
            for(unsigned int k = 0 ; k < B->m ; k++) {
                mpq_srcptr aik = mpq_mat_entry_const(A, i, k);
                mpq_srcptr bkj = mpq_mat_entry_const(B, k, j);
                mpq_mul(z, aik, bkj);
                mpq_add(cij, cij, z);
            }
            mpq_canonicalize(cij);
        }
    }
    mpq_clear(z);
}/*}}}*/

void mpq_mat_mul_mpz(mpq_mat_ptr B, mpq_mat_srcptr A, mpz_srcptr k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij), mpq_denref(aij));
            mpz_mul(mpq_numref(bij), mpq_numref(aij), k);
            mpq_canonicalize(bij);
        }
    }
}/*}}}*/
void mpq_mat_mul_mpq(mpq_mat_ptr B, mpq_mat_srcptr A, mpq_srcptr k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpq_mul(bij, aij, k);
        }
    }
}
/*}}}*/
void mpq_mat_mul_si(mpq_mat_ptr B, mpq_mat_srcptr A, long k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij),mpq_denref(aij));
            mpz_mul_si(mpq_numref(bij),mpq_numref(aij),k);
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_mul_ui(mpq_mat_ptr B, mpq_mat_srcptr A, unsigned long k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij),mpq_denref(aij));
            mpz_mul_ui(mpq_numref(bij),mpq_numref(aij),k);
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_mpz(mpq_mat_ptr B, mpq_mat_srcptr A, mpz_srcptr k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul(mpq_denref(bij), mpq_denref(aij), k);
            mpz_set(mpq_numref(bij), mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
}/*}}}*/
void mpq_mat_div_mpq(mpq_mat_ptr B, mpq_mat_srcptr A, mpq_srcptr k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul(mpq_numref(bij), mpq_numref(aij), mpq_denref(k));
            mpz_mul(mpq_denref(bij), mpq_denref(aij), mpq_numref(k));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_si(mpq_mat_ptr B, mpq_mat_srcptr A, long k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul_si(mpq_denref(bij),mpq_denref(aij), k);
            mpz_set(mpq_numref(bij),mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_ui(mpq_mat_ptr B, mpq_mat_srcptr A, unsigned long k)/*{{{*/
{
    mpq_mat_realloc(B,A->m,A->n);       /* no-op if A == B */
    for(unsigned int i = 0 ; i < B->m ; i++) {
        for(unsigned int j = 0 ; j < B->n ; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul_ui(mpq_denref(bij),mpq_denref(aij), k);
            mpz_set(mpq_numref(bij),mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
/* }}} */
/*{{{ powering */
// Returns A^n for n >= 0, A being a square matrix
void mpz_mat_pow_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n)/*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_set_ui(B, 1);
        return;
    }
    if (B == A) {
        mpz_mat C;
        mpz_mat_init(C, A->m, A->n);
        mpz_mat_pow_ui(C, A, n);
        mpz_mat_swap(B, C);
        mpz_mat_clear(C);
        return;
    }
    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > n ; k >>= 1);
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for( ; k >>= 1 ; ) {
        mpz_mat_mul(Q, B, B);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul(Q, B, A);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
}/*}}}*/

void mpz_mat_pow_ui_mod_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n, unsigned long p)/*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_set_ui(B, 1);
        return;
    }
    if (B == A) {
        mpz_mat C;
        mpz_mat_init(C, A->m, A->n);
        mpz_mat_pow_ui_mod_ui(C, A, n, p);
        mpz_mat_swap(B, C);
        mpz_mat_clear(C);
        return;
    }
    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > n ; k >>= 1);
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for( ; k >>= 1 ; ) {
        mpz_mat_mul_mod_ui(Q, B, B, p);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul_mod_ui(Q, B, A, p);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
}/*}}}*/

void mpz_mat_pow_ui_mod_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_set_ui(B, 1);
        return;
    }
    if (B == A) {
        mpz_mat C;
        mpz_mat_init(C, A->m, A->n);
        mpz_mat_pow_ui_mod_mpz(C, A, n, p);
        mpz_mat_swap(B, C);
        mpz_mat_clear(C);
        return;
    }
    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > n ; k >>= 1);
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for( ; k >>= 1 ; ) {
        mpz_mat_mul_mod_mpz(Q, B, B, p);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul_mod_mpz(Q, B, A, p);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
}/*}}}*/
/*}}}*/
/* {{{ polynomial evaluation */

void mpz_poly_eval_mpz_mat(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f)/*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int n = M->n;
    int d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat(X, M, f);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    mpz_mat_set_mpz(D, f->coeff[f->deg]);
    for(int i = f->deg - 1 ; i >= 0 ; i--) {
        mpz_mat_mul(D, M, D);
        mpz_mat_add_mpz(D, f->coeff[i]);
    }
}
/*}}}*/
void mpz_poly_eval_mpz_mat_mod_ui(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f, unsigned long p)/*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int n = M->n;
    int d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat_mod_ui(X, M, f, p);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    mpz_mat_set_ui(D, mpz_fdiv_ui(f->coeff[f->deg], p));
    for(int i = f->deg - 1 ; i >= 0 ; i--) {
        mpz_mat_mul_mod_ui(D, M, D, p);
        mpz_mat_add_ui(D, mpz_fdiv_ui(f->coeff[i], p));
        mpz_mat_mod_ui(D, D, p);
    }
}
/*}}}*/
void mpz_poly_eval_mpz_mat_mod_mpz(mpz_mat_ptr D, mpz_mat_srcptr M, mpz_poly_srcptr f, mpz_srcptr p)/*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int n = M->n;
    int d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat_mod_mpz(X, M, f, p);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_fdiv_r(tmp, f->coeff[f->deg], p);
    mpz_mat_set_mpz(D, tmp);
    for(int i = f->deg - 1 ; i >= 0 ; i--) {
        mpz_mat_mul_mod_mpz(D, M, D, p);
        mpz_fdiv_r(tmp, f->coeff[i], p);
        mpz_mat_add_mpz(D, tmp);
        mpz_mat_mod_mpz(D, D, p);
    }
    mpz_clear(tmp);
}
/*}}}*/
/*}}}*/

/* {{{ gaussian reduction over the rationals
 * this is a backend for row gaussian reduction. T receives the
 * transformation matrix. M is modified in place.
 * this includes reordering of the rows.
 *
 * this works reasonably well over Fp, but for rationals we suffer from
 * coefficient blowup.
 */
void mpq_mat_gauss_backend(mpq_mat_ptr M, mpq_mat_ptr T)
{
    unsigned int m = M->m;
    unsigned int n = M->n;
    ASSERT_ALWAYS(M != T);
    mpq_mat_realloc(T, m, m);
    mpq_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpq_t tmp1, tmp2;
    mpq_init(tmp1);
    mpq_init(tmp2);
    for(unsigned int j = 0 ; j < n ; j++) {
        unsigned int i1;
        for(i1 = rank ; i1 < M->m ; i1++) {
            if (mpq_cmp_ui(mpq_mat_entry(M, i1, j), 0, 1) != 0) break;
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        mpq_mat_swaprows(M, rank, i1);
        mpq_mat_swaprows(T, rank, i1);
        i1 = rank;
        /* canonicalize this row */
        mpq_inv(tmp1, mpq_mat_entry(M, rank, j));
        mpq_mat_mulrow(M, rank, tmp1);
        mpq_mat_mulrow(T, rank, tmp1);
        for(unsigned int i0 = 0 ; i0 < m ; i0++) {
            if (i0 == rank) continue;
            /* we've normalized row rank, so all it takes is to multiply by
             * M[i0, j]
             */
            mpq_set(tmp2, mpq_mat_entry(M, i0, j));
            mpq_mat_submulrow(M, i0, rank, tmp2);
            mpq_mat_submulrow(T, i0, rank, tmp2);
        }
        rank++;
    }
    mpq_clear(tmp1);
    mpq_clear(tmp2);
}
/* }}} */
/* {{{ Gaussian reduction over Z/pZ
 *
 * We sort of return something even in the case where p is not prime.
 *
 * M is the reduced matrix (modification in place of the input), T the
 * transformation matrix (so that T * input-M = output-M)
 *
 * Note that if M is rectangular, this the function we want to call in
 * order to get a row-reduced echelon form of M.
 *
 * T may be NULL in case we don't give a penny about the transformation
 * matrix.
 */
void mpz_mat_gauss_backend_mod_mpz(mpz_mat_ptr M, mpz_mat_ptr T, mpz_srcptr p)
{
    unsigned int m = M->m;
    unsigned int n = M->n;
    ASSERT_ALWAYS(M != T);
    if (T) mpz_mat_realloc(T, m, m);
    if (T) mpz_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpz_t gcd, tmp2, tmp3;
    mpz_init(gcd);
    mpz_init(tmp2);
    mpz_init(tmp3);
    for(unsigned int j = 0 ; j < n ; j++) {
        unsigned int i1 = M->m;
        mpz_set(gcd, p);
        for(unsigned int i = rank ; i < M->m ; i++) {
            mpz_gcd(tmp2, mpz_mat_entry(M, i, j), p);
            if (mpz_cmp(tmp2, gcd) < 0) {
                i1 = i;
                mpz_set(gcd, tmp2);
                if (mpz_cmp_ui(gcd, 1) == 0) break;
            }
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        mpz_mat_swaprows(M, rank, i1);
        if (T) mpz_mat_swaprows(T, rank, i1);
        i1 = rank;
        /* canonicalize this row */
        /* gcd is gcd(M[rank, j], p) */
        mpz_divexact(tmp2, mpz_mat_entry(M, rank, j), gcd);
        mpz_divexact(tmp3, p, gcd);
        mpz_invert(tmp2, tmp2, tmp3);
        mpz_mat_mulrow_mod(M, rank, tmp2, p);
        if (T) mpz_mat_mulrow_mod(T, rank, tmp2, p);
        ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry(M, rank, j), gcd) == 0);
        for(unsigned int i0 = 0 ; i0 < m ; i0++) {
            if (i0 == rank) continue;
            /* we've normalized row rank, so all it takes is to multiply by
             * M[i0, j]
             */
            mpz_fdiv_q(tmp2, mpz_mat_entry(M, i0, j), gcd);
            /* for i0 > rank, we should have exact quotient */
            mpz_mat_submulrow_mod(M, i0, rank, tmp2, p);
            if (T) mpz_mat_submulrow_mod(T, i0, rank, tmp2, p);
            ASSERT_ALWAYS(mpz_cmp_ui(mpz_mat_entry(M, i0, j), 0) == 0);
        }
        rank++;
    }
    mpz_clear(gcd);
    mpz_clear(tmp2);
    mpz_clear(tmp3);
}
/* }}} */

void mpz_mat_gauss_backend_mod_ui(mpz_mat_ptr M, mpz_mat_ptr T, unsigned long p)
{
    mpz_t pz;
    mpz_init_set_ui(pz, p);
    mpz_mat_gauss_backend_mod_mpz(M, T, pz);
    mpz_clear(pz);
}


/* {{{ integer heap routines */
/* input: heap[0..n-1[ a heap, and heap[n-1] unused.
 * output: heap[0..n[ a heap, includes x
 */
typedef int (*cmp_t)(void *, unsigned int, unsigned int);

static void push_heap(unsigned int * heap, unsigned int x, unsigned int n, cmp_t cmp, void * arg)
{
    unsigned int i = n - 1;
    for(unsigned int i0; i; i = i0) {
        i0 = (i-1) / 2;
        if (cmp(arg, x, heap[i0]) > 0)
            heap[i] = heap[i0];
        else
            break;
    }
    heap[i] = x;
}

/* put heap[0] in x, and keep the heap property on the resulting
 * heap[0..n-1[ ; heap[n-1] is unused on output.
 */
static void pop_heap(unsigned int * x, unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    *x = heap[0];
    unsigned int y = heap[n-1];
    /* make heap[0..n-1[ a heap, now */
    /* couldn't there be various strategies here ? */
    unsigned int i = 0;
    for(unsigned int i1 ; 2*i+2 <  n - 1 ; i = i1) {
        unsigned int left = 2*i + 1;
        unsigned int right = 2*i + 2;
        if (cmp(arg, heap[left], heap[right]) > 0)
            i1 = left;
        else
            i1 = right;
        if (cmp(arg, y, heap[i1]) > 0)
            break;
        heap[i] = heap[i1];
    }
    /* there we know that i has no right child, maybe even no left child.
     */
    if (2*i + 1 < n - 1 && cmp(arg, y, heap[2*i+1]) <= 0) {
        heap[i] = heap[2*i+1];
        i=2*i+1;
    }
    heap[i] = y;
}

/* make heap[0..n[ a heap. */
static void make_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = 1 ; i <= n ; i++) {
        unsigned int t = heap[i-1];
        push_heap(heap, t, i, cmp, arg);
    }
}

/* must be a heap already */
static void sort_heap(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    for(unsigned int i = n ; i ; i--) {
        unsigned int t;
        pop_heap(&t, heap, i, cmp, arg);
        heap[i-1] = t;
    }
}

static void sort_heap_inverse(unsigned int * heap, unsigned int n, cmp_t cmp, void * arg)
{
    sort_heap(heap, n, cmp, arg);
    for(unsigned int i = 0, j = n-1 ; i < j ; i++, j--) {
        int a = heap[i];
        int b = heap[j];
        heap[i] = b;
        heap[j] = a;
    }
}
/* }}} */
/* {{{ Hermite normal form over the integers
 * The algorithm here is randomly cooked from a stupid idea I had. The
 * transformation matrices produced seem to be reasonably good, somewhat
 * better than what magma gives (for overdetermined matrices, of course
 * -- otherwise they're unique). However we do so in an embarrassingly
 * larger amount of time (more than 1000* for 100x100 matrices).
 * Fortunately, our intent is not to run this code on anything larger
 * than 30*30, in which case it's sufficient.
 */
/*{{{ hnf helper algorithm (for column matrices) */
struct mpz_mat_hnf_helper_heap_aux {
    double * dd;
    mpz_mat_ptr A;
};

static int mpz_mat_hnf_helper_heap_aux_cmp(struct mpz_mat_hnf_helper_heap_aux * data, unsigned int a, unsigned int b)
{
    /* The heap will have the largest elements on top according to this
     * measure. Thus the comparison on the ->a field must be as per
     * mpz_cmp
     */
    mpz_srcptr xa = mpz_mat_entry_const(data->A, a, 0);
    mpz_srcptr xb = mpz_mat_entry_const(data->A, b, 0);

    assert(mpz_sgn(xa) >= 0);
    assert(mpz_sgn(xb) >= 0);
    int r = mpz_cmp(xa, xb);
    if (r) return r;
    /* among combinations giving the same result, we want the "best" to
     * be the one with *smallest* norm. So instead of the sign of da-db,
     * we want the sign of db-da, here.
     */
    double da = data->dd[a];
    double db = data->dd[b];
    return (db > da) - (da > db);
}

/* this is quite the same as computing the hnf on a column matrix, with
 * the added functionality that entries in a[0..n0[ are reduced with
 * respect to the entries in a[n0..n[. This is equivalent to saying that
 * we apply a transformation matrix which is:
 *       I_n0 R
 * T ==  0    S
 * where S is unimodular. The resulting vector T*A has entries [n0..n[
 * equal to (g,0,...,0), while entries [0..n0[ are reduced modulo g.
 *
 * return +1 or -1 which is the determinant of T.
 */
static int mpz_mat_hnf_helper(mpz_mat_ptr T, mpz_mat_ptr dT, mpz_ptr * a, unsigned int n0, unsigned int n)
{
    int signdet=1;
    ASSERT_ALWAYS(dT != T);
    mpz_mat_realloc(dT, n, n);
    mpz_mat_set_ui(dT, 1);

    if (n == n0) return signdet;

    mpz_t q, r2;
    mpz_mat A;
    mpz_mat_init(A, n, 1);
    mpz_init(q);
    mpz_init(r2);
    unsigned int * heap = (unsigned int *) malloc(n * sizeof(unsigned int));
    struct mpz_mat_hnf_helper_heap_aux cmpdata[1];
    cmpdata->A = A;
    cmpdata->dd = (double*) malloc(n * sizeof(double));
    cmp_t cmp = (cmp_t) &mpz_mat_hnf_helper_heap_aux_cmp;
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_ptr Ai = mpz_mat_entry(A, i, 0);
        mpz_set(Ai, a[i]);
        cmpdata->dd[i] = 0;
        for(unsigned int j = 0 ; j < n ; j++) {
            double d = mpz_get_d(mpz_mat_entry(T, i, j));
            cmpdata->dd[i] += d * d;
        }
        heap[i] = i;
        if (i < n0)
            continue;
        if (mpz_sgn(Ai) == -1) {
            mpz_neg(Ai, Ai);
            signdet = -signdet;
            mpz_set_si(mpz_mat_entry(dT, i, i), -1);
            for(unsigned int j = 0 ; j < n ; j++) {
                mpz_neg(mpz_mat_entry(T, i, j), mpz_mat_entry(T, i, j));
            }
        }
    }
    make_heap(heap + n0, n - n0, cmp, cmpdata);
    unsigned int head;
    for( ; ; ) {
        /* extract the head */
        pop_heap(&head, heap + n0, n - n0, cmp, cmpdata);
        mpz_ptr r0 = mpz_mat_entry(A, head, 0);
        if (mpz_cmp_ui(r0, 0) == 0) {
            /* we're done ! */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        /* reduce A[0..n0[ wrt r0 == A[head] */
        for(unsigned int i = 0 ; i < n0 ; i++) {
            mpz_fdiv_qr(q, r2, mpz_mat_entry(A, i, 0), r0);
            mpz_swap(r2, mpz_mat_entry(A, i, 0));
            mpz_mat_submulrow(T, i, head, q);
            mpz_mat_submulrow(dT, i, head, q);
        }
        if (n0 == n - 1) {
            /* we're actually computing the gcd of one single integer.
             * That makes no sense of course.
             */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        mpz_ptr r1 = mpz_mat_entry(A, heap[n0], 0);
        if (mpz_cmp_ui(r1, 0) == 0) {
            /* we're done ! */
            push_heap(heap + n0, head, n - n0, cmp, cmpdata);
            break;
        }
        mpz_fdiv_qr(q, r2, r0, r1);
        mpz_swap(r0, r2);
        mpz_mat_submulrow(T, head, heap[n0], q);
        mpz_mat_submulrow(dT, head, heap[n0], q);
        /* adjust the norm of the row in T */
        cmpdata->dd[head] = 0;
        for(unsigned int j = 0 ; j < n ; j++) {
            double d = mpz_get_d(mpz_mat_entry(T, head, j));
            cmpdata->dd[head] += d * d;
        }
        push_heap(heap + n0, head, n - n0, cmp, cmpdata);
    }
    sort_heap_inverse(heap + n0, n - n0, cmp, cmpdata);
    mpz_mat_permuterows(T, heap);
    mpz_mat_permuterows(dT, heap);
    mpz_mat_permuterows(A, heap);
    signdet*= permutation_signature(heap, n);
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_ptr Ai = mpz_mat_entry(A, i, 0);
        mpz_swap(a[i], Ai);
    }
    free(cmpdata->dd);
    mpz_mat_clear(A);
    mpz_clear(q);
    mpz_clear(r2);
    free(heap);

    return signdet;
}
/* {{{ this computes the row hnf on a column C.
 * The result is always of the form C'=(gcd(C),0,0,...,0). The transform
 * matrix T such that C'=T*C is most interesting.  a is modified in
 * place.
 */
void mpz_gcd_many(mpz_mat_ptr dT, mpz_ptr * a, unsigned int n)
{
    mpz_mat T;
    mpz_mat_init(T, n, n);
    mpz_mat_set_ui(T, 1);
    mpz_mat_hnf_helper(T, dT, a, 0, n);
    mpz_mat_clear(T);
}
/*}}}*/
/*}}}*/

/* return +1 or -1, which is the determinant of the transformation matrix
 * T */
int mpz_mat_hnf_backend(mpz_mat_ptr M, mpz_mat_ptr T)
{
    ASSERT_ALWAYS(M != T);
    if (T == NULL) {
        mpz_mat xT;
        mpz_mat_init(xT,0,0);
        int r = mpz_mat_hnf_backend(M, xT);
        mpz_mat_clear(xT);
        return r;
    }
    int signdet = 1;
    unsigned int m = M->m;
    unsigned int n = M->n;
    mpz_mat_realloc(T, m, m);
    mpz_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpz_t quo;
    mpz_init(quo);
    mpz_ptr * colm = (mpz_ptr *) malloc(m * sizeof(mpz_t));
    mpz_mat dT, Mx, My;
    mpz_mat_init(dT, 0, 0);
    mpz_mat_init(Mx, 0, 0);
    mpz_mat_init(My, 0, 0);
    for(unsigned int j = 0 ; j < n && rank < m; j++) {
        for(unsigned int i = 0 ; i < m ; i++) {
            colm[i] = mpz_mat_entry(M, i, j);
        }
        signdet *= mpz_mat_hnf_helper(T, dT, colm, rank, m);

        /* apply dT to the submatrix of size m * (n - 1 - j) of M */
        mpz_mat_realloc(Mx, m, n - 1 - j);
        mpz_mat_submat_swap(Mx, 0, 0, M, 0, j + 1, m, n - 1 - j);
        mpz_mat_mul(My, dT, Mx);
        mpz_mat_submat_swap(M, 0, j + 1, My, 0, 0, m, n - 1 - j);

        mpz_srcptr gcd = colm[rank];
        if (mpz_cmp_ui(gcd, 0) == 0) /* this column is zero -- nothing to do */
            continue;
        rank++;
    }
    free(colm);
    mpz_clear(quo);
    mpz_mat_clear(dT);
    mpz_mat_clear(Mx);
    mpz_mat_clear(My);

    return signdet;
}
/* }}} */
/*{{{ kernel*/
// This is supposed to compute the Kernel of M mod p and to store it in
// the matrix K. If r is the rank of M, and M is a square matrix n*n, K
// is a n*(n-r) matrix
// 
// self-assignment is ok.
void mpz_mat_kernel_mod_mpz(mpz_mat_ptr K, mpz_mat_srcptr M, mpz_srcptr p)
{
    mpz_mat T, H;
    unsigned int r;

    mpz_mat_init(T,M->m,M->m);
    mpz_mat_init(H,M->m,M->n);
    mpz_mat_set(H,M);

    // Storing a reduced matrix of M in H with gauss
    mpz_mat_gauss_backend_mod_mpz(H,T,p);

    // Finding the rank of M and H
    r = H->m;
    while((r>0) && mpz_mat_isnull_row(H,r-1)){
        r--;
    }
    // Kernel is of dimension n-r, and a basis of the kernel is given in the n-r last rows of T
    // We shall keep the convention of magma
    
    // Reallocating K with n-r columns
    mpz_mat_realloc(K, H->m-r, H->m);
    mpz_mat_submat_swap(    K, 0, 0,
                            T, r, 0, 
                            H->m - r, H->m);

    mpz_mat_clear(T);
    mpz_mat_clear(H);
}
void mpz_mat_kernel_mod_ui(mpz_mat_ptr K, mpz_mat_srcptr M, unsigned long p)
{
    mpz_t p2;
    mpz_init_set_ui(p2,p);
    mpz_mat_kernel_mod_mpz(K,M,p2);
    mpz_clear(p2);
}
/* }}} */
/*{{{ inversion*/
void mpq_mat_inv(mpq_mat_ptr dst, mpq_mat_srcptr src)
{
    ASSERT_ALWAYS(src->m == src->n);

    mpq_mat aux;
    mpq_mat_init(aux,src->m,src->n);
    mpq_mat_set(aux,src);

    /* This is self-assignment compatible (done reading src before
     * dst is touched). */
    mpq_mat_gauss_backend(aux,dst);

    mpq_mat_clear(aux);
}
/* }}} */

/*
 * For now, it is just a wrapper to use utils/lll.c.
 */
void mpz_mat_LLL(mpz_ptr det, mpz_mat_ptr M, mpz_mat_ptr U, mpz_srcptr a,
    mpz_srcptr b)
{
  mat_Z M_tmp;
  LLL_init(&M_tmp, M->m, M->n);
  for (unsigned int row = 1; row < M->m + 1; row++) {
    for (unsigned int col = 1; col < M->n + 1; col++) {
      mpz_set(M_tmp.coeff[row][col], mpz_mat_entry_const(M, row - 1, col - 1));
    }
  }

  if (U) {
    mat_Z U_tmp;
    LLL_init(&U_tmp, U->m, U->m);

    LLL(det, M_tmp, &U_tmp, a, b);

    for (unsigned int row = 1; row < U->m + 1; row++) {
      for (unsigned int col = 1; col < U->n + 1; col++) {
        mpz_set(mpz_mat_entry(U, row - 1, col - 1), U_tmp.coeff[row][col]);
      }
    }
    LLL_clear(&U_tmp);
  } else {
    LLL(det, M_tmp, NULL, a, b);
  }

  for (unsigned int row = 1; row < M->m + 1; row++) {
    for (unsigned int col = 1; col < M->n + 1; col++) {
      mpz_set(mpz_mat_entry(M, row - 1, col - 1), M_tmp.coeff[row][col]);
    }
  }

  LLL_clear(&M_tmp);
}

using namespace std;

ostream& operator<<(ostream& os, cxx_mpz_mat const& M)/*{{{*/
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            if (j) os << " ";
            os << mpz_mat_entry_const(M, i, j);
        }
        if (i < M->m - 1) os << "\n";
    }
    return os;
}
/*}}}*/
ostream& operator<<(ostream& os, cxx_mpq_mat const& M)/*{{{*/
{
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            if (j) os << " ";
            os << mpq_mat_entry_const(M, i, j);
        }
        if (i < M->m - 1) os << "\n";
    }
    return os;
}
/*}}}*/
