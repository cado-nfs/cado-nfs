#include "cado.h" // IWYU pragma: keep

#include <climits> /* for INT_MAX */
#include <cstdio>  // FILE // IWYU pragma: keep
#include <cstdlib>
#include <cstring>

#include <ostream>
#include <vector>
#include <string>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"
#include "fmt/ostream.h"

#include "cxx_mpz.hpp"
#include "fmt_helper_sagemath.hpp"
#include "lll.h" // mat_Z, LLL
#include "macros.h"
#include "mpz_mat.h"
#include "mpz_poly.h"

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
    M->x = (mpz_t *)((m && n) ? malloc(m * n * sizeof(mpz_t)) : NULL);

    M->m = m;
    M->n = n;
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpz_init(mpz_mat_entry(M, i, j));
}

void mpz_mat_clear(mpz_mat_ptr M)
{
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpz_clear(mpz_mat_entry(M, i, j));
    free(M->x);
}

void mpz_mat_realloc(mpz_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n)
        return;
    mpz_mat_clear(M);
    mpz_mat_init(M, m, n);
}

void mpq_mat_init(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    M->x = (mpq_t *)((m && n) ? malloc(m * n * sizeof(mpq_t)) : NULL);
    M->m = m;
    M->n = n;
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpq_init(mpq_mat_entry(M, i, j));
}

void mpq_mat_clear(mpq_mat_ptr M)
{
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpq_clear(mpq_mat_entry(M, i, j));
    free(M->x);
}

void mpq_mat_realloc(mpq_mat_ptr M, unsigned int m, unsigned int n)
{
    if (M->m == m && M->n == n)
        return;
    mpq_mat_clear(M);
    mpq_mat_init(M, m, n);
}

/*}}}*/
/*{{{ operations on submatrices, and swaps*/
void mpz_mat_submat_swap(mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
                         mpz_mat_ptr A1, unsigned int i1, unsigned int j1,
                         unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for (unsigned int i = 0; i < dm; i++) {
        for (unsigned int j = 0; j < dn; j++) {
            mpz_swap(mpz_mat_entry(A0, i0 + i, j0 + j),
                     mpz_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpz_mat_submat_set(mpz_mat_ptr A0, unsigned int i0, unsigned int j0,
                        mpz_mat_srcptr A1, unsigned int i1, unsigned int j1,
                        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for (unsigned int i = 0; i < dm; i++) {
        for (unsigned int j = 0; j < dn; j++) {
            mpz_set(mpz_mat_entry(A0, i0 + i, j0 + j),
                    mpz_mat_entry_const(A1, i1 + i, j1 + j));
        }
    }
}

void mpq_mat_swap(mpq_mat_ptr A, mpq_mat_ptr B)
{
    mpz_mat C;
    memcpy(C, A, sizeof(mpz_mat));
    memcpy(A, B, sizeof(mpz_mat));
    memcpy(B, C, sizeof(mpz_mat));
}

void mpq_mat_submat_swap(mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
                         mpq_mat_ptr A1, unsigned int i1, unsigned int j1,
                         unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for (unsigned int i = 0; i < dm; i++) {
        for (unsigned int j = 0; j < dn; j++) {
            mpq_swap(mpq_mat_entry(A0, i0 + i, j0 + j),
                     mpq_mat_entry(A1, i1 + i, j1 + j));
        }
    }
}

void mpq_mat_submat_set(mpq_mat_ptr A0, unsigned int i0, unsigned int j0,
                        mpq_mat_srcptr A1, unsigned int i1, unsigned int j1,
                        unsigned int dm, unsigned int dn)
{
    ASSERT_ALWAYS(i0 + dm <= A0->m);
    ASSERT_ALWAYS(i1 + dm <= A1->m);
    ASSERT_ALWAYS(j0 + dn <= A0->n);
    ASSERT_ALWAYS(j1 + dn <= A1->n);
    for (unsigned int i = 0; i < dm; i++) {
        for (unsigned int j = 0; j < dn; j++) {
            mpq_set(mpq_mat_entry(A0, i0 + i, j0 + j),
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

/*}}}*/
/*{{{ set/set_ui*/
void mpz_mat_set(mpz_mat_ptr dst, mpz_mat_srcptr src)
{
    mpz_mat_realloc(dst, src->m, src->n);
    for (unsigned int i = 0; i < src->m; i++)
        for (unsigned int j = 0; j < src->n; j++)
            mpz_set(mpz_mat_entry(dst, i, j), mpz_mat_entry_const(src, i, j));
}

void mpz_mat_set_ui(mpz_mat_ptr M, unsigned long a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpz_set_ui(mpz_mat_entry(M, i, j), i == j ? a : 0);
}

void mpz_mat_set_mpz(mpz_mat_ptr M, mpz_srcptr a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        mpz_set(mpz_mat_entry(M, i, i), a);
    }
}

void mpz_mat_add_ui(mpz_mat_ptr M, unsigned long a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        mpz_ptr mii = mpz_mat_entry(M, i, i);
        mpz_add_ui(mii, mii, a);
    }
}

void mpz_mat_add_mpz(mpz_mat_ptr M, mpz_srcptr a)
{
    ASSERT_ALWAYS(M->m == M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        mpz_ptr mii = mpz_mat_entry(M, i, i);
        mpz_add(mii, mii, a);
    }
}

void mpq_mat_set_ui(mpq_mat_ptr M, unsigned long a)
{
    if (a) {
        ASSERT_ALWAYS(M->m == M->n);
    }
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpq_set_ui(mpq_mat_entry(M, i, j), i == j ? a : 0, 1);
}

void mpq_mat_set(mpq_mat_ptr dst, mpq_mat_srcptr src)
{
    mpq_mat_realloc(dst, src->m, src->n);
    for (unsigned int i = 0; i < src->m; i++)
        for (unsigned int j = 0; j < src->n; j++)
            mpq_set(mpq_mat_entry(dst, i, j), mpq_mat_entry_const(src, i, j));
}

void mpz_mat_urandomm(mpz_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            mpz_urandomm(mpz_mat_entry(M, i, j), state, p);
}

void mpq_mat_urandomm(mpq_mat_ptr M, gmp_randstate_t state, mpz_srcptr p)
{
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++) {
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
    mpz_mat_realloc(N, M1->m + M2->m, M1->n);
    mpz_mat_submat_set(N, 0, 0, M1, 0, 0, M1->m, M1->n);
    mpz_mat_submat_set(N, M1->m, 0, M2, 0, 0, M2->m, M2->n);
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
    mpq_mat_realloc(N, M1->m + M2->m, M1->n);
    mpq_mat_submat_set(N, 0, 0, M1, 0, 0, M1->m, M1->n);
    mpq_mat_submat_set(N, M1->m, 0, M2, 0, 0, M2->m, M2->n);
}
void mpz_mat_horizontal_join(mpz_mat_ptr N, mpz_mat_srcptr M1,
                             mpz_mat_srcptr M2)
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
    mpz_mat_realloc(N, M1->m, M1->n + M2->n);
    mpz_mat_submat_set(N, 0, 0, M1, 0, 0, M1->m, M1->n);
    mpz_mat_submat_set(N, 0, M1->n, M2, 0, 0, M2->m, M2->n);
}
void mpq_mat_horizontal_join(mpq_mat_ptr N, mpq_mat_srcptr M1,
                             mpq_mat_srcptr M2)
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
    mpq_mat_realloc(N, M1->m, M1->n + M2->n);
    mpq_mat_submat_set(N, 0, 0, M1, 0, 0, M1->m, M1->n);
    mpq_mat_submat_set(N, 0, M1->n, M2, 0, 0, M2->m, M2->n);
}
/*}}}*/
/*{{{ determinant, trace and transposition */
void mpz_mat_trace(mpz_ptr t, mpz_mat_srcptr M) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(t, 0);
    for (unsigned int i = 0; i < M->n; i++)
        mpz_add(t, t, mpz_mat_entry_const(M, i, i));
}
/*}}}*/
void mpz_mat_determinant_triangular(mpz_ptr d, mpz_mat_srcptr M) /*{{{*/
{
    // We assume that M is triangular
    ASSERT_ALWAYS(M->m == M->n);
    mpz_set_ui(d, 1);
    for (unsigned int i = 0; i < M->n; i++)
        mpz_mul(d, d, mpz_mat_entry_const(M, i, i));
}
/*}}}*/
void mpz_mat_determinant(mpz_ptr d, mpz_mat_srcptr A)
{
    cxx_mpz_mat H;
    mpz_mat_set(H, A);
    int s = mpz_mat_hermite_form(H);
    mpz_mat_determinant_triangular(d, H);
    mpz_mul_si(d, d, s);
}

void mpz_mat_transpose(mpz_mat_ptr D, mpz_mat_srcptr M) /*{{{*/
{
    if (D != M) {
        mpz_mat_realloc(D, M->n, M->m);
        for (unsigned int i = 0; i < M->m; i++) {
            for (unsigned int j = 0; j < M->n; j++) {
                mpz_srcptr mij = mpz_mat_entry_const(M, i, j);
                mpz_ptr dji = mpz_mat_entry(D, j, i);
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
        for (unsigned int i = 0; i < M->m; i++) {
            for (unsigned int j = i + 1; j < M->n; j++) {
                mpz_ptr dij = mpz_mat_entry(D, i, j);
                mpz_ptr dji = mpz_mat_entry(D, j, i);
                mpz_swap(dji, dij);
            }
        }
    }
} /*}}}*/

void mpz_mat_reverse_rows(mpz_mat_ptr B, mpz_mat_srcptr A)
{
    if (A != B) {
        mpz_mat_realloc(B, A->m, A->n);
        for (unsigned int i = 0; i < A->m; i++) {
            mpz_mat_submat_set(B, i, 0, A, A->m - 1 - i, 0, 1, A->n);
        }
    } else {
        for (unsigned int i = 0; i < A->m - 1 - i; i++) {
            mpz_mat_submat_swap(B, i, 0, B, A->m - 1 - i, 0, 1, A->n);
        }
    }
}

void mpz_mat_reverse_columns(mpz_mat_ptr B, mpz_mat_srcptr A)
{
    if (A != B) {
        mpz_mat_realloc(B, A->m, A->n);
        for (unsigned int j = 0; j < A->n; j++) {
            mpz_mat_submat_set(B, 0, j, A, 0, A->n - 1 - j, A->m, 1);
        }
    } else {
        for (unsigned int j = 0; j < A->n - 1 - j; j++) {
            mpz_mat_submat_swap(B, 0, j, B, 0, A->n - 1 - j, A->m, 1);
        }
    }
}

void mpq_mat_trace(mpq_ptr t, mpq_mat_srcptr M) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(t, 0, 1);
    for (unsigned int i = 0; i < M->n; i++)
        mpq_add(t, t, mpq_mat_entry_const(M, i, i));
}
/*}}}*/
void mpq_mat_determinant_triangular(mpq_ptr d, mpq_mat_srcptr M) /*{{{*/
{
    // We assume that M is triangular
    ASSERT_ALWAYS(M->m == M->n);
    mpq_set_ui(d, 1, 1);
    for (unsigned int i = 0; i < M->n; i++)
        mpq_mul(d, d, mpq_mat_entry_const(M, i, i));
}
/*}}}*/
void mpq_mat_determinant(mpq_ptr d, mpq_mat_srcptr M) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    cxx_mpz_mat num;
    cxx_mpz den;
    mpq_mat_numden(num, den, M);
    mpz_mat_determinant(mpq_numref(d), num);
    mpz_pow_ui(mpq_denref(d), den, M->n);
    mpq_canonicalize(d);
}
/*}}}*/
void mpq_mat_transpose(mpq_mat_ptr D, mpq_mat_srcptr M) /*{{{*/
{
    if (D != M) {
        mpq_mat_realloc(D, M->n, M->m);
        for (unsigned int i = 0; i < M->m; i++) {
            for (unsigned int j = 0; j < M->n; j++) {
                mpq_srcptr mij = mpq_mat_entry_const(M, i, j);
                mpq_ptr dji = mpq_mat_entry(D, j, i);
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
        for (unsigned int i = 0; i < M->m; i++) {
            for (unsigned int j = i + 1; j < M->n; j++) {
                mpq_ptr dij = mpq_mat_entry(D, i, j);
                mpq_ptr dji = mpq_mat_entry(D, j, i);
                mpq_swap(dji, dij);
            }
        }
    }
} /*}}}*/
/*}}}*/
/*{{{ miscellaneous */

/* convert to integer matrix divided by lcm of denominator.
 *
 * returns 1 if the denominator is 1, 0 otherwise.
 *
 * if num==NULL, only the denominator is returned.
 *
 * if den==NULL, only the return value is meaningful (num is not changed)
 */

int mpq_mat_numden(mpz_mat_ptr num, mpz_ptr den, mpq_mat_srcptr M) /*{{{*/
{
    int ret = 1;
    if (den)
        mpz_set_ui(den, 1);
    for (unsigned int i = 0; i < M->m; i++) {
        for (unsigned int j = 0; j < M->n; j++) {
            mpz_srcptr Mij = mpq_denref(mpq_mat_entry_const(M, i, j));
            if (den) {
                mpz_lcm(den, den, Mij);
            } else if (mpz_cmp_ui(Mij, 1) != 0) {
                ret = 0;
            }
        }
    }
    if (!ret || !num)
        return ret;
    mpz_mat_realloc(num, M->m, M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        for (unsigned int j = 0; j < M->n; j++) {
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
void mpq_mat_set_mpz_mat(mpq_mat_ptr N, mpz_mat_srcptr M) /*{{{*/
{
    mpq_mat_realloc(N, M->m, M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        for (unsigned int j = 0; j < M->n; j++) {
            mpq_set_z(mpq_mat_entry(N, i, j), mpz_mat_entry_const(M, i, j));
        }
    }
}
/*}}}*/
void mpq_mat_set_mpz_mat_denom(mpq_mat_ptr N, mpz_mat_srcptr M,
                               mpz_srcptr d) /*{{{*/
{
    mpq_mat_realloc(N, M->m, M->n);
    for (unsigned int i = 0; i < M->m; i++) {
        for (unsigned int j = 0; j < M->n; j++) {
            mpq_ptr nij = mpq_mat_entry(N, i, j);
            mpz_srcptr mij = mpz_mat_entry_const(M, i, j);
            mpz_set(mpq_numref(nij), mij);
            mpz_set(mpq_denref(nij), d);
            mpq_canonicalize(nij);
        }
    }
}
/*}}}*/
void mpz_mat_mod_ui(mpz_mat_ptr dst, mpz_mat_srcptr src,
                    unsigned long p) /*{{{*/
{
    mpz_mat_realloc(dst, src->m, src->n); /* no-op if dst == src */
    for (unsigned int i = 0; i < src->m; i++) {
        for (unsigned int j = 0; j < src->n; j++) {
            mpz_fdiv_r_ui(mpz_mat_entry(dst, i, j),
                          mpz_mat_entry_const(src, i, j), p);
        }
    }
} /*}}}*/
void mpz_mat_mod_mpz(mpz_mat_ptr dst, mpz_mat_srcptr src, mpz_srcptr p) /*{{{*/
{
    mpz_mat_realloc(dst, src->m, src->n); /* no-op if dst == src */
    for (unsigned int i = 0; i < src->m; i++) {
        for (unsigned int j = 0; j < src->n; j++) {
            mpz_fdiv_r(mpz_mat_entry(dst, i, j), mpz_mat_entry_const(src, i, j),
                       p);
        }
    }
} /*}}}*/

int mpz_mat_p_valuation(mpz_mat_srcptr A, mpz_srcptr p)
{
    int val = INT_MAX;
    cxx_mpz c;
    for (unsigned int i = 0; val && i < A->m; i++) {
        for (unsigned int j = 0; val && j < A->n; j++) {
            int v = 0;
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            if (mpz_size(aij) == 0)
                continue;
            mpz_set(c, aij);
            for (; v < val && mpz_divisible_p(c, p); v++)
                mpz_fdiv_q(c, c, p);
            val = v;
        }
    }
    return val;
}

int mpz_mat_p_valuation_ui(mpz_mat_srcptr A, unsigned long p)
{
    int val = INT_MAX;
    cxx_mpz c;
    for (unsigned int i = 0; val && i < A->m; i++) {
        for (unsigned int j = 0; val && j < A->n; j++) {
            int v = 0;
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            if (mpz_size(aij) == 0)
                continue;
            mpz_set(c, aij);
            for (; v < val && mpz_divisible_ui_p(c, p); v++)
                mpz_fdiv_q_ui(c, c, p);
            val = v;
        }
    }
    return val;
}

/*}}}*/

/* {{{ row-level operations (for Gauss and friends, mostly) */
// Return 1 if the k-th line of M is null, 0 else
int mpz_mat_isnull_row(mpz_mat_srcptr M, unsigned int k)
{ /*{{{*/
    unsigned int j = 0;
    ASSERT_ALWAYS(k < M->m);
    while ((j < M->n) && !mpz_cmp_si(mpz_mat_entry_const(M, k, j), 0)) {
        j++;
    }
    if (j == M->n) {
        return 1;
    }
    return 0;
}
/*}}}*/
void mpz_mat_swaprows(mpz_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1)
        return;
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_swap(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
void mpq_mat_swaprows(mpq_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    if (i0 == i1)
        return;
    for (unsigned int j = 0; j < M->n; j++) {
        mpq_swap(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/
/*}}}*/

/* add lambda times row i1 to row i0 */
void mpz_mat_addmulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1,
                       mpz_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_addmul(mpz_mat_entry(M, i0, j), lambda, mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* add lambda times row i1 to row i0 */
void mpq_mat_addmulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1,
                       mpq_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for (unsigned int j = 0; j < M->n; j++) {
        /* we have no mpq_addmul, unfortunately */
        mpq_mul(tmp, lambda, mpq_mat_entry(M, i1, j));
        mpq_add(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j), tmp);
    }
    mpq_clear(tmp);
}
/*}}}*/
/* subtract lambda times row i1 to row i0 */
void mpz_mat_submulrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1,
                       mpz_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j), lambda, mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
static void mpz_mat_submulrow_mod(mpz_mat_ptr M, unsigned int i0, unsigned int i1,
                           mpz_srcptr lambda, mpz_srcptr p) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_submul(mpz_mat_entry(M, i0, j), lambda, mpz_mat_entry(M, i1, j));
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}
/*}}}*/
/* subtract lambda times row i1 to row i0 */
void mpq_mat_submulrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1,
                       mpq_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    mpq_t tmp;
    mpq_init(tmp);
    for (unsigned int j = 0; j < M->n; j++) {
        /* we have no mpq_submul, unfortunately */
        mpq_mul(tmp, lambda, mpq_mat_entry(M, i1, j));
        mpq_sub(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j), tmp);
    }
    mpq_clear(tmp);
}
/*}}}*/

/* add row i1 to row i0 */
void mpz_mat_addrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_add(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* add row i1 to row i0 */
void mpq_mat_addrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpq_add(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/

/* subtract row i1 to row i0 */
void mpz_mat_subrow(mpz_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_sub(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j),
                mpz_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* subtract row i1 to row i0 */
void mpq_mat_subrow(mpq_mat_ptr M, unsigned int i0, unsigned int i1) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    ASSERT_ALWAYS(i1 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpq_sub(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j),
                mpq_mat_entry(M, i1, j));
    }
}
/*}}}*/
/* multiply row i0 by lambda */
void mpz_mat_mulrow(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
    }
}
/*}}}*/
void mpz_mat_mulrow_mod(mpz_mat_ptr M, unsigned int i0, mpz_srcptr lambda,
                        mpz_srcptr p) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpz_mul(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), lambda);
        mpz_mod(mpz_mat_entry(M, i0, j), mpz_mat_entry(M, i0, j), p);
    }
}
/*}}}*/
/* multiply row i0 by lambda */
void mpq_mat_mulrow(mpq_mat_ptr M, unsigned int i0, mpq_srcptr lambda) /*{{{*/
{
    ASSERT_ALWAYS(i0 < M->m);
    for (unsigned int j = 0; j < M->n; j++) {
        mpq_mul(mpq_mat_entry(M, i0, j), mpq_mat_entry(M, i0, j), lambda);
    }
}
/*}}}*/

/* }}} */
/*{{{ comparison */
int mpz_mat_cmp(mpz_mat_srcptr M, mpz_mat_srcptr N) /*{{{*/
{
    /* shall we abort or return something for incompatible matrices ?? */
    if (M->m != N->m)
        return M->m - N->m;
    if (M->n != N->n)
        return M->n - N->n;
    // ASSERT_ALWAYS((M->m == N->m) && (M->n == N->n));
    unsigned int const m = M->m;
    unsigned int const n = M->n;
    unsigned int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            int const k = mpz_cmp(mpz_mat_entry_const(M, i, j),
                                  mpz_mat_entry_const(N, i, j));
            if (k != 0)
                return k;
        }
    }
    return 0;
}
/*}}}*/
int mpq_mat_cmp(mpq_mat_srcptr M, mpq_mat_srcptr N) /*{{{*/
{
    /* shall we abort or return something for incompatible matrices ?? */
    if (M->m != N->m)
        return M->m - N->m;
    if (M->n != N->n)
        return M->n - N->n;
    // ASSERT_ALWAYS((M->m == N->m) && (M->n == N->n));
    unsigned int const m = M->m;
    unsigned int const n = M->n;
    unsigned int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            int const k = mpq_cmp(mpq_mat_entry_const(M, i, j),
                                  mpq_mat_entry_const(N, i, j));
            if (k != 0)
                return k;
        }
    }
    return 0;
}
/*}}}*/

int mpz_mat_is_zero(mpz_mat_srcptr M)
{
    for (unsigned int i = 0; i < M->m; i++)
        for (unsigned int j = 0; j < M->n; j++)
            if (mpz_cmp_ui(mpz_mat_entry_const(M, i, j), 0) != 0)
                return 0;
    return 1;
}

/* TODO (perhaps) :
 * mp[qz]_mat_is_{upper,lower}_triangular
 * mp[qz]_mat_is_diagonal
 */
/*}}}*/
/*{{{ add and subtract */
void mpz_mat_add(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B) /*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpz_mat_realloc(C, A->m, A->n); /* no-op if C == A or C == B */

    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < A->n; j++) {
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            mpz_srcptr bij = mpz_mat_entry_const(B, i, j);
            mpz_ptr cij = mpz_mat_entry(C, i, j);
            mpz_add(cij, aij, bij);
        }
    }
}
/* }}} */
void mpq_mat_add(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B) /*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpq_mat_realloc(C, A->m, A->n); /* no-op if C == A or C == B */

    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < A->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_srcptr bij = mpq_mat_entry_const(B, i, j);
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_add(cij, aij, bij);
        }
    }
}
/*}}}*/
void mpz_mat_sub(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B) /*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpz_mat_realloc(C, A->m, A->n); /* no-op if C == A or C == B */

    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < A->n; j++) {
            mpz_srcptr aij = mpz_mat_entry_const(A, i, j);
            mpz_srcptr bij = mpz_mat_entry_const(B, i, j);
            mpz_ptr cij = mpz_mat_entry(C, i, j);
            mpz_sub(cij, aij, bij);
        }
    }
}
/* }}} */
void mpq_mat_sub(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B) /*{{{*/
{
    ASSERT_ALWAYS(A->m == B->m);
    ASSERT_ALWAYS(A->n == B->n);
    mpq_mat_realloc(C, A->m, A->n); /* no-op if C == A or C == B */

    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < A->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_srcptr bij = mpq_mat_entry_const(B, i, j);
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_sub(cij, aij, bij);
        }
    }
}
/*}}}*/
/*}}}*/
/*{{{ multiplication */
void mpz_mat_mul(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B) /*{{{*/
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
    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for (unsigned int k = 0; k < B->m; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j), mpz_mat_entry_const(A, i, k),
                           mpz_mat_entry_const(B, k, j));
            }
        }
    }
} /*}}}*/

void mpz_mat_mul_mod_ui(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B,
                        unsigned long p) /*{{{*/
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
    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for (unsigned int k = 0; k < B->m; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j), mpz_mat_entry_const(A, i, k),
                           mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_mod_ui(C, C, p);
} /*}}}*/

void mpz_mat_mul_mod_mpz(mpz_mat_ptr C, mpz_mat_srcptr A, mpz_mat_srcptr B,
                         mpz_srcptr p) /*{{{*/
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
    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_set_ui(mpz_mat_entry(C, i, j), 0);
            for (unsigned int k = 0; k < B->m; k++) {
                mpz_addmul(mpz_mat_entry(C, i, j), mpz_mat_entry_const(A, i, k),
                           mpz_mat_entry_const(B, k, j));
            }
        }
    }
    mpz_mat_mod_mpz(C, C, p);
} /*}}}*/

void mpz_mat_mul_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, mpz_ptr k) /*{{{*/
{
    mpz_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_mul(mpz_mat_entry(B, i, j), mpz_mat_entry_const(A, i, j), k);
        }
    }
} /*}}}*/

void mpz_mat_divexact_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, mpz_srcptr k) /*{{{*/
{
    mpz_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_divexact(mpz_mat_entry(B, i, j), mpz_mat_entry_const(A, i, j),
                         k);
        }
    }
} /*}}}*/

void mpz_mat_divexact_ui(mpz_mat_ptr B, mpz_mat_srcptr A,
                         unsigned long k) /*{{{*/
{
    mpz_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_divexact_ui(mpz_mat_entry(B, i, j),
                            mpz_mat_entry_const(A, i, j), k);
        }
    }
} /*}}}*/

void mpz_mat_mul_si(mpz_mat_ptr B, mpz_mat_srcptr A, long k) /*{{{*/
{
    mpz_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_mul_si(mpz_mat_entry(B, i, j), mpz_mat_entry_const(A, i, j), k);
        }
    }
} /*}}}*/

void mpz_mat_mul_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long k) /*{{{*/
{
    mpz_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpz_mul_ui(mpz_mat_entry(B, i, j), mpz_mat_entry_const(A, i, j), k);
        }
    }
} /*}}}*/

void mpq_mat_mul(mpq_mat_ptr C, mpq_mat_srcptr A, mpq_mat_srcptr B) /*{{{*/
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
    for (unsigned int i = 0; i < A->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_ptr cij = mpq_mat_entry(C, i, j);
            mpq_set_ui(cij, 0, 1);
            for (unsigned int k = 0; k < B->m; k++) {
                mpq_srcptr aik = mpq_mat_entry_const(A, i, k);
                mpq_srcptr bkj = mpq_mat_entry_const(B, k, j);
                mpq_mul(z, aik, bkj);
                mpq_add(cij, cij, z);
            }
            mpq_canonicalize(cij);
        }
    }
    mpq_clear(z);
} /*}}}*/

void mpq_mat_mul_mpz(mpq_mat_ptr B, mpq_mat_srcptr A, mpz_srcptr k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij), mpq_denref(aij));
            mpz_mul(mpq_numref(bij), mpq_numref(aij), k);
            mpq_canonicalize(bij);
        }
    }
} /*}}}*/
void mpq_mat_mul_mpq(mpq_mat_ptr B, mpq_mat_srcptr A, mpq_srcptr k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpq_mul(bij, aij, k);
        }
    }
}
/*}}}*/
void mpq_mat_mul_si(mpq_mat_ptr B, mpq_mat_srcptr A, long k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij), mpq_denref(aij));
            mpz_mul_si(mpq_numref(bij), mpq_numref(aij), k);
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_mul_ui(mpq_mat_ptr B, mpq_mat_srcptr A, unsigned long k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_set(mpq_denref(bij), mpq_denref(aij));
            mpz_mul_ui(mpq_numref(bij), mpq_numref(aij), k);
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_mpz(mpq_mat_ptr B, mpq_mat_srcptr A, mpz_srcptr k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul(mpq_denref(bij), mpq_denref(aij), k);
            mpz_set(mpq_numref(bij), mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
} /*}}}*/
void mpq_mat_div_mpq(mpq_mat_ptr B, mpq_mat_srcptr A, mpq_srcptr k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul(mpq_numref(bij), mpq_numref(aij), mpq_denref(k));
            mpz_mul(mpq_denref(bij), mpq_denref(aij), mpq_numref(k));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_si(mpq_mat_ptr B, mpq_mat_srcptr A, long k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul_si(mpq_denref(bij), mpq_denref(aij), k);
            mpz_set(mpq_numref(bij), mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
void mpq_mat_div_ui(mpq_mat_ptr B, mpq_mat_srcptr A, unsigned long k) /*{{{*/
{
    mpq_mat_realloc(B, A->m, A->n); /* no-op if A == B */
    for (unsigned int i = 0; i < B->m; i++) {
        for (unsigned int j = 0; j < B->n; j++) {
            mpq_srcptr aij = mpq_mat_entry_const(A, i, j);
            mpq_ptr bij = mpq_mat_entry(B, i, j);
            mpz_mul_ui(mpq_denref(bij), mpq_denref(aij), k);
            mpz_set(mpq_numref(bij), mpq_numref(aij));
            mpq_canonicalize(bij);
        }
    }
}
/*}}}*/
/* }}} */
/*{{{ powering */
// Returns A^n for n >= 0, A being a square matrix
void mpz_mat_pow_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n) /*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_realloc(B, A->m, A->n);
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
    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for (; k >>= 1;) {
        mpz_mat_mul(Q, B, B);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul(Q, B, A);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
} /*}}}*/

void mpz_mat_pow_ui_mod_ui(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n,
                           unsigned long p) /*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_realloc(B, A->m, A->n);
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
    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for (; k >>= 1;) {
        mpz_mat_mul_mod_ui(Q, B, B, p);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul_mod_ui(Q, B, A, p);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
} /*}}}*/

void mpz_mat_pow_ui_mod_mpz(mpz_mat_ptr B, mpz_mat_srcptr A, unsigned long n,
                            mpz_srcptr p) /*{{{*/
{
    ASSERT_ALWAYS(A->n == A->m);
    if (n == 0) {
        mpz_mat_realloc(B, A->m, A->n);
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
    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    mpz_mat_set(B, A);
    mpz_mat Q;
    mpz_mat_init(Q, A->m, A->n);
    for (; k >>= 1;) {
        mpz_mat_mul_mod_mpz(Q, B, B, p);
        mpz_mat_swap(Q, B);
        if (n & k) {
            mpz_mat_mul_mod_mpz(Q, B, A, p);
            mpz_mat_swap(Q, B);
        }
    }
    mpz_mat_clear(Q);
} /*}}}*/
/*}}}*/
/* {{{ polynomial evaluation */

void mpz_poly_eval_mpz_mat(mpz_mat_ptr D, mpz_poly_srcptr f,
        mpz_mat_srcptr M) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int const n = M->n;
    int const d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat(X, f, M);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    mpz_mat_set_mpz(D, mpz_poly_coeff_const(f, f->deg));
    for (int i = f->deg - 1; i >= 0; i--) {
        mpz_mat_mul(D, M, D);
        mpz_mat_add_mpz(D, mpz_poly_coeff_const(f, i));
    }
}
/*}}}*/
void mpz_poly_eval_mpz_mat_mod_ui(mpz_mat_ptr D, mpz_poly_srcptr f,
        mpz_mat_srcptr M, unsigned long p) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int const n = M->n;
    int const d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat_mod_ui(X, f, M, p);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    mpz_mat_set_ui(D, mpz_fdiv_ui(mpz_poly_coeff_const(f, f->deg), p));
    for (int i = f->deg - 1; i >= 0; i--) {
        mpz_mat_mul_mod_ui(D, M, D, p);
        mpz_mat_add_ui(D, mpz_fdiv_ui(mpz_poly_coeff_const(f, i), p));
        mpz_mat_mod_ui(D, D, p);
    }
}
/*}}}*/
void mpz_poly_eval_mpz_mat_mod_mpz(mpz_mat_ptr D, mpz_poly_srcptr f,
        mpz_mat_srcptr M, mpz_srcptr p) /*{{{*/
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned int const n = M->n;
    int const d = f->deg;
    if (d < 0) {
        mpz_mat_realloc(D, n, n);
        mpz_mat_set_ui(D, 0);
        return;
    }
    if (D == M) {
        mpz_mat X;
        mpz_mat_init(X, 0, 0);
        mpz_poly_eval_mpz_mat_mod_mpz(X, f, M, p);
        mpz_mat_swap(D, X);
        mpz_mat_clear(X);
        return;
    }
    mpz_mat_realloc(D, n, n);
    cxx_mpz tmp;
    mpz_fdiv_r(tmp, mpz_poly_coeff_const(f, f->deg), p);
    mpz_mat_set_mpz(D, tmp);
    for (int i = f->deg - 1; i >= 0; i--) {
        mpz_mat_mul_mod_mpz(D, M, D, p);
        mpz_fdiv_r(tmp, mpz_poly_coeff_const(f, i), p);
        mpz_mat_add_mpz(D, tmp);
        mpz_mat_mod_mpz(D, D, p);
    }
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
    unsigned int const m = M->m;
    unsigned int const n = M->n;
    ASSERT_ALWAYS(M != T);
    mpq_mat_realloc(T, m, m);
    mpq_mat_set_ui(T, 1);
    unsigned int rank = 0;
    mpq_t tmp1, tmp2;
    mpq_init(tmp1);
    mpq_init(tmp2);
    for (unsigned int j = 0; j < n; j++) {
        unsigned int i1;
        for (i1 = rank; i1 < M->m; i1++) {
            if (mpq_cmp_ui(mpq_mat_entry(M, i1, j), 0, 1) != 0)
                break;
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        if (i1 > rank) {
            mpq_mat_swaprows(M, rank, i1);
            mpq_mat_swaprows(T, rank, i1);
        }
        /* canonicalize this row */
        mpq_inv(tmp1, mpq_mat_entry(M, rank, j));
        mpq_mat_mulrow(M, rank, tmp1);
        mpq_mat_mulrow(T, rank, tmp1);
        for (unsigned int i0 = 0; i0 < m; i0++) {
            if (i0 == rank)
                continue;
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
    unsigned int const m = M->m;
    unsigned int const n = M->n;
    ASSERT_ALWAYS(M != T);
    if (T)
        mpz_mat_realloc(T, m, m);
    if (T)
        mpz_mat_set_ui(T, 1);
    unsigned int rank = 0;
    cxx_mpz gcd, tmp2, tmp3;
    for (unsigned int j = 0; j < n; j++) {
        unsigned int i1 = M->m;
        mpz_set(gcd, p);
        for (unsigned int i = rank; i < M->m; i++) {
            mpz_gcd(tmp2, mpz_mat_entry(M, i, j), p);
            if (mpz_cmp(tmp2, gcd) < 0) {
                i1 = i;
                mpz_set(gcd, tmp2);
                if (mpz_cmp_ui(gcd, 1) == 0)
                    break;
            }
        }
        if (i1 == M->m) /* this column is zero */
            continue;
        if (i1 > rank) {
            mpz_mat_swaprows(M, rank, i1);
            if (T)
                mpz_mat_swaprows(T, rank, i1);
        }
        // i1 = rank;
        /* canonicalize this row */
        /* gcd is gcd(M[rank, j], p) */
        mpz_divexact(tmp2, mpz_mat_entry(M, rank, j), gcd);
        mpz_divexact(tmp3, p, gcd);
        mpz_invert(tmp2, tmp2, tmp3);
        mpz_mat_mulrow_mod(M, rank, tmp2, p);
        if (T)
            mpz_mat_mulrow_mod(T, rank, tmp2, p);
        ASSERT_ALWAYS(mpz_cmp(mpz_mat_entry(M, rank, j), gcd) == 0);
        for (unsigned int i0 = 0; i0 < m; i0++) {
            if (i0 == rank)
                continue;
            /* we've normalized row rank, so all it takes is to multiply by
             * M[i0, j]
             */
            mpz_fdiv_q(tmp2, mpz_mat_entry(M, i0, j), gcd);
            /* for i0 > rank, we should have exact quotient */
            mpz_mat_submulrow_mod(M, i0, rank, tmp2, p);
            if (T)
                mpz_mat_submulrow_mod(T, i0, rank, tmp2, p);
            ASSERT_ALWAYS(mpz_cmp_ui(mpz_mat_entry(M, i0, j), 0) == 0);
        }
        rank++;
    }
}
/* }}} */

void mpz_mat_gauss_backend_mod_ui(mpz_mat_ptr M, mpz_mat_ptr T, unsigned long p)
{
    mpz_mat_gauss_backend_mod_mpz(M, T, cxx_mpz(p));
}

/* {{{ Hermite normal form over the integers
 * The algorithm here is very naive and straightforward, but at least
 * it's not embarassingly slow on tiny matrices.  Our intent is not to
 * run this code on anything larger than 30*30, in which case it's
 * sufficient.
 */
/*{{{ hnf helper algorithm (for column matrices) */

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
static int mpz_mat_hnf_helper(mpz_mat_ptr dT,
                              cxx_mpz_mat & a, unsigned int n0)
{
    unsigned int const n = a->m;
    int signdet = 1;
    mpz_mat_realloc(dT, n, n);
    mpz_mat_set_ui(dT, 1);

    ASSERT_ALWAYS(n0 <= n);

    if (n == n0)
        return signdet;

    /* row n0 is always considered active. If it is zero (which might
     * happen on entry), then this will be fixed on the first round if
     * there is another non-zero coefficient. Otherwise the column is
     * reduced and there is nothing to do. */
    std::vector<unsigned int> active;
    for(unsigned int i = n0 + 1 ; i < n ; i++)
        active.push_back(i);

    do {
        {
            unsigned int i0 = n0;
            for(const unsigned int i : active) {
                if (mpz_sgn(a(i, 0)) == 0)
                    continue;
                if (mpz_sgn(a(i0, 0)) == 0 || mpz_cmpabs(a(i, 0), a(i0, 0)) < 0)
                    i0 = i;
            }
            if (i0 != n0) {
                mpz_mat_swaprows(a, n0, i0);
                mpz_mat_swaprows(dT, n0, i0);
                signdet *= -1;
            }
        }

        active.clear();
        cxx_mpz q, r;

        if (mpz_sgn(a(n0, 0)) == 0)
            return signdet;

        /* It's not mandatory to do it right here if we do ndiv, but for
         * fdiv it is important.
         */
        if (mpz_sgn(a(n0, 0)) < 0) {
            signdet *= -1;
            mpz_neg(a(n0, 0), a(n0, 0));
            for(unsigned int j = 0 ; j < dT->n ; j++)
                mpz_neg(mpz_mat_entry(dT, n0, j), mpz_mat_entry(dT, n0, j));
        }

        /* reduce A[0..n[ wrt r0 == A[head] */
        for (unsigned int i = 0; i < n; i++) {
            if (i == n0 || mpz_cmp_ui(a(i, 0), 0) == 0)
                continue;

            /* use ndiv_qr to get centered coefficients, otherwise use
             * fdiv_qr */
            mpz_fdiv_qr(q, r, a(i, 0), a(n0, 0));
            if (i > n0 && r != 0)
                active.push_back(i);
            mpz_swap(r, a(i, 0));
            mpz_mat_submulrow(dT, i, n0, q);
        }

        /* break here so that we always get at least one run of the loop */
    } while (!active.empty());

    return signdet;
}

/*}}}*/

/*
 * M is put into HNF form.
 * return +1 or -1, which is the determinant of the transformation matrix
 * T */
int mpz_mat_hermite_form(mpz_mat_ptr M) /* , mpz_mat_ptr T) */
{
    int signdet = 1;
    unsigned int const m = M->m;
    unsigned int const n = M->n;
    unsigned int rank = 0;
    cxx_mpz_mat dT, Mx, My;
    cxx_mpz_mat colm(m, 1);

    for (unsigned int j = 0; j < n && rank < m; j++) {

        mpz_mat_submat_swap(M, 0, j, colm, 0, 0, m, 1);
        signdet *= mpz_mat_hnf_helper(/* T, */ dT, colm, rank);
        rank += mpz_cmp_ui(colm(rank, 0), 0) != 0;
        mpz_mat_submat_swap(M, 0, j, colm, 0, 0, m, 1);

        /* apply dT to the right part. By construction, our
         * transformation matrix dT ha already been applied to T itself,
         * as well as to column j. We also know that it has trivial
         * action on the left part (column of indices below j, which
         * contain data only for i < rank <= j). So the only thing that
         * is left is the right submatrix of size m * (n - 1 - j) of M */
        mpz_mat_realloc(Mx, m, n - 1 - j);
        mpz_mat_submat_swap(Mx, 0, 0, M, 0, j + 1, m, n - 1 - j);
        mpz_mat_mul(My, dT, Mx);
        mpz_mat_submat_swap(M, 0, j + 1, My, 0, 0, m, n - 1 - j);
    }

    return signdet;
}
/* }}} */

/* This is almost like mpz_mat_hermite_form, except that we do it in a
 * different order, which is more suitable for displaying number
 * field elements in a way which ends up being similar to magma's
 *
 * M is put into HNF form.
 */
int mpz_mat_hermite_form_rev(mpz_mat_ptr M) /*, mpz_mat_ptr T */ // {{{
{
    mpz_mat_reverse_rows(M, M);
    mpz_mat_reverse_columns(M, M);
    int s = mpz_mat_hermite_form(M);
    mpz_mat_reverse_rows(M, M);
    mpz_mat_reverse_columns(M, M);
    if (M->m > M->n) {
        /* we need some swaps... */
        mpz_mat sM;
        mpz_mat_init(sM, M->m, M->n);
        mpz_mat_submat_swap(sM, 0, 0, M, M->m - M->n, 0, M->n, M->n);
        mpz_mat_submat_swap(sM, M->n, 0, M, 0, 0, M->m - M->n, M->n);
        mpz_mat_swap(sM, M);
        mpz_mat_clear(sM);
        /* While the transformations above had no effect on s (because
         * they compensate), this one has.
         * we have n circular shifts on length m, plus a reversal on m-n.
         * a circular shift on length k is exactly k-1 inversions, so
         * that sums up to n*(m-1) inversions. Then we add
         * (m-n)*(m-n-1)/2 inversions. This is, in total,
         * (m(m-1)+n(n-1))/2 inversions.
         * m*(m-1) is congruent to 2 mod 4 when m is 2 or 3 mod 4
         */
        int const ninvs = ((M->m & 2) + (M->n & 2)) / 2;
        if (ninvs)
            s = -s;
    }
    return s;
} //}}}

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

    mpz_mat_init(T, M->m, M->m);
    mpz_mat_init(H, M->m, M->n);
    mpz_mat_set(H, M);

    // Storing a reduced matrix of M in H with gauss
    mpz_mat_gauss_backend_mod_mpz(H, T, p);

    // Finding the rank of M and H
    r = H->m;
    while ((r > 0) && mpz_mat_isnull_row(H, r - 1)) {
        r--;
    }
    // Kernel is of dimension n-r, and a basis of the kernel is given in the n-r
    // last rows of T We shall keep the convention of magma

    // Reallocating K with n-r columns
    mpz_mat_realloc(K, H->m - r, H->m);
    mpz_mat_submat_swap(K, 0, 0, T, r, 0, H->m - r, H->m);

    mpz_mat_clear(T);
    mpz_mat_clear(H);
}
void mpz_mat_kernel_mod_ui(mpz_mat_ptr K, mpz_mat_srcptr M, unsigned long p)
{
    mpz_mat_kernel_mod_mpz(K, M, cxx_mpz(p));
}
/* }}} */
/*{{{ inversion*/
void mpq_mat_inv(mpq_mat_ptr dst, mpq_mat_srcptr src)
{
    ASSERT_ALWAYS(src->m == src->n);

    mpq_mat aux;
    mpq_mat_init(aux, src->m, src->n);
    mpq_mat_set(aux, src);

    /* This is self-assignment compatible (done reading src before
     * dst is touched). */
    mpq_mat_gauss_backend(aux, dst);

    mpq_mat_clear(aux);
}

void mpz_mat_inv_mod_mpz(mpz_mat_ptr dst, mpz_mat_srcptr src, mpz_srcptr p)
{
    ASSERT_ALWAYS(src->m == src->n);

    mpz_mat aux;
    mpz_mat_init(aux,src->m,src->n);
    mpz_mat_set(aux,src);

    mpz_mat_gauss_backend_mod_mpz(aux, dst, p);

    mpz_mat_clear(aux);
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
            mpz_set(M_tmp.coeff[row][col],
                    mpz_mat_entry_const(M, row - 1, col - 1));
        }
    }

    if (U) {
        mat_Z U_tmp;
        LLL_init(&U_tmp, U->m, U->m);

        LLL(det, M_tmp, &U_tmp, a, b);

        for (unsigned int row = 1; row < U->m + 1; row++) {
            for (unsigned int col = 1; col < U->n + 1; col++) {
                mpz_set(mpz_mat_entry(U, row - 1, col - 1),
                        U_tmp.coeff[row][col]);
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



namespace fmt {

template<typename matrix_type, fmt_helper_sagemath_types custom_format>
inline constexpr const char * cado_matrix_ring_name = "";

template<> inline constexpr const char * cado_matrix_ring_name<cxx_mpz_mat, fmt_helper_sagemath_types::SAGEMATH> = "ZZ";
template<> inline constexpr const char * cado_matrix_ring_name<cxx_mpq_mat, fmt_helper_sagemath_types::SAGEMATH> = "QQ";
template<> inline constexpr const char * cado_matrix_ring_name<cxx_mpz_mat, fmt_helper_sagemath_types::MAGMA> = "Integers()";
template<> inline constexpr const char * cado_matrix_ring_name<cxx_mpq_mat, fmt_helper_sagemath_types::MAGMA> = "Rationals()";


template<typename matrix_type>
static auto format_cado_matrix(matrix_type const & M,
        fmt_helper_sagemath_types custom_format,
        format_context& ctx) -> format_context::iterator
{
    /* using enum is one of the few bits of c++20 that aren't supported
     * by g++-10, and we do have a few of these around.
     */
    using fmt_helper_sagemath_types::SAGEMATH;
    using fmt_helper_sagemath_types::TEXT;
    using fmt_helper_sagemath_types::MACHINE;
    using fmt_helper_sagemath_types::MAGMA;

    if (custom_format == SAGEMATH) {
        format_to(ctx.out(), "matrix({}, {}, {}, [\n",
                cado_matrix_ring_name<matrix_type, SAGEMATH>,
                M->m, M->n);
    } else if (custom_format == MAGMA) {
        format_to(ctx.out(), "Matrix({}, {}, {}, [\n",
                cado_matrix_ring_name<matrix_type, MAGMA>,
                M->m, M->n);
    } else if (custom_format == TEXT) {
        format_to(ctx.out(), "[");
    }
    for(unsigned int i = 0 ; i < M->m ; i++) {
        for(unsigned int j = 0 ; j < M->n ; j++) {
            std::string after;
            if (custom_format != MACHINE)
                if (j < M->n - 1 || i < M->m - 1)
                    after = ",";
            if (custom_format == SAGEMATH || custom_format == MAGMA) {
                if (j < M->n - 1)
                    after +=  " ";
                else if (i < M->m - 1)
                    after += "\n";
            } else if (custom_format == TEXT || custom_format == MACHINE) {
                if (j < M->n - 1 || i < M->m - 1)
                    after +=  " ";
            }
            typename matrix_type::cxx_coeff_type Mij = M(i, j);
            format_to(ctx.out(), "{}", Mij);
            format_to(ctx.out(), "{}", after);
        }
    }
    if (custom_format == SAGEMATH || custom_format == MAGMA) {
        format_to(ctx.out(), "])");
    } else if (custom_format == TEXT) {
        format_to(ctx.out(), "]");
    }
    return ctx.out();
}

auto formatter<cxx_mpz_mat>::format(cxx_mpz_mat const & M, format_context& ctx) const -> format_context::iterator
{
    return format_cado_matrix(M, custom_format, ctx);
}
auto formatter<cxx_mpq_mat>::format(cxx_mpq_mat const & M, format_context& ctx) const -> format_context::iterator
{
    return format_cado_matrix(M, custom_format, ctx);
}

} /* namespace fmt */

std::ostream & operator<<(std::ostream & os, cxx_mpz_mat const & M) /*{{{*/
{
    fmt::print(os, "{:S}", M);
    return os;
}
/*}}}*/
std::ostream & operator<<(std::ostream & os, cxx_mpq_mat const & M) /*{{{*/
{
    fmt::print(os, "{:S}", M);
    return os;
}
/*}}}*/

/*{{{ I/O*/
void mpz_mat_fprint(FILE * stream, mpz_mat_srcptr M)
{
    cxx_mpz_mat tM;
    mpz_mat_set(tM, M);
    fmt::print(stream, "{:S}\n", tM);
}

void mpq_mat_fprint(FILE * stream, mpq_mat_srcptr M)
{
    cxx_mpq_mat tM;
    mpq_mat_set(tM, M);
    fmt::print(stream, "{:S}\n", tM);
}
/*}}}*/
