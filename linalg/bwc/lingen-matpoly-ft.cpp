#include "cado.h"
#include <stdlib.h>
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.hpp"
#include "lingen-matpoly-ft.hpp"
#include "logline.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* timings made on cochon, rev 6877b97 (buggy; fixed in 76dde6c) */
#define MP_FTI_DEPTH_ADJ_24_36_36 { { 1, 6 }, { 2, 4 }, { 3, 3 }, { 4, 2 }, { 10, 1 }, { 11, 2 }, { 13, 1 }, { 22, 2 }, { 28, 1 }, { 32, 2 }, { 33, 1 }, { 38, 0 }, { 39, 1 }, { 54, 0 }, { 55, 1 }, { 64, 0 }, { 65, 1 }, { 66, 0 }, { 103, 1 }, { 104, 0 }, { 107, 1 }, { 114, 0 }, { 115, 1 }, { 129, 0 }, }

#define MUL_FTI_DEPTH_ADJ_36_36_36 { { 1, 6 }, { 2, 3 }, { 3, 2 }, { 6, 1 }, { 7, 2 }, { 14, 1 }, { 23, 0 }, { 26, 1 }, { 44, 0 }, { 46, 1 }, { 54, 0 }, { 61, 1 }, { 62, 0 }, }

void matpoly_ft_init(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, unsigned int m, unsigned int n, const struct fft_transform_info * fti)/*{{{*/
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    t->m = m;
    t->n = n;
    t->data = malloc(m * n * fft_alloc_sizes[0]);
    memset(t->data, 0, m * n * fft_alloc_sizes[0]);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_prepare(tij, fti);
        }
    }
}/*}}}*/

void matpoly_ft_clear(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, const struct fft_transform_info * fti MAYBE_UNUSED)/*{{{*/
{
    free(t->data);
    memset(t, 0, sizeof(*t));
}/*}}}*/

void matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr a, const struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(fft_alloc_sizes[1]);
        unsigned int m = t->m;     /* for icc ... */
        unsigned int n = t->n;
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++) {
                size_t offset = (i*t->n + j) * fft_alloc_sizes[0];
                void * tij = pointer_arith(t->data, offset);
                absrc_vec aij = matpoly_part_const(ab, a, i, j, 0);
                /* ok, casting like this is a crude hack ! */
                fft_do_dft_fppol(tij, (const mp_limb_t *) aij, a->size, tt, fti, ab->p);
            }
        free(tt);
    }
}

void matpoly_ft_zero(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, const struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_zero(tij, fti);
        }
    }
}

void matpoly_ft_export(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, const struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_export(tij, fti);
        }
    }
}

void matpoly_ft_import(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, const struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_import(tij, fti);
        }
    }
}

void matpoly_ft_add(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t0->m == t1->m);
    ASSERT_ALWAYS(t0->n == t1->n);
    ASSERT_ALWAYS(t0->m == u->m);
    ASSERT_ALWAYS(t0->n == u->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < t0->m ; i++) {
        for(unsigned int j = 0 ; j < t0->n ; j++) {
            size_t offset = (i*t0->n+j) * fft_alloc_sizes[0];
            void * t0ij = pointer_arith(t0->data, offset);
            void * t1ij = pointer_arith(t1->data, offset);
            void * uij  = pointer_arith(u->data, offset);
            fft_add(uij, t0ij, t1ij, fti);
        }
    }
}

void matpoly_ft_addmul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti)
{
    ASSERT_ALWAYS(t0->n == t1->m);
    ASSERT_ALWAYS(t0->m == u->m);
    ASSERT_ALWAYS(t1->n == u->n);
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * qt = malloc(fft_alloc_sizes[1]);
        void * tt = malloc(fft_alloc_sizes[2]);
        memset(qt, 0, fft_alloc_sizes[1]);

        unsigned int m = t0->m;     /* for icc ... */
        unsigned int n = t1->n;
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++) {
                size_t tsize = fft_alloc_sizes[0];
                memset(tt, 0, fft_alloc_sizes[2]);
                for(unsigned int k = 0 ; k < t0->n ; k++) {
                    void * uij  = pointer_arith(u->data, (i*u->n+j) * tsize);
                    void * t0ik = pointer_arith(t0->data, (i*t0->n+k) * tsize);
                    void * t1kj = pointer_arith(t1->data, (k*t1->n+j) * tsize);
                    fft_addmul(uij, t0ik, t1kj, tt, qt, fti);
                }
            }
        free(tt);
        free(qt);
    }
}

void matpoly_ft_mul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti)
{
    matpoly_ft_zero(ab, u, fti);
    return matpoly_ft_addmul(ab, u, t0, t1, fti);
}

/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
*/

void matpoly_ft_ift(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, const struct fft_transform_info * fti)
{
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(fft_alloc_sizes[1]);
        unsigned int m = t->m;     /* for icc ... */
        unsigned int n = t->n;
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_ift_fppol((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p);
            }
        free(tt);
    }
}

void matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, unsigned int shift, const struct fft_transform_info * fti)
{
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(fft_alloc_sizes[1]);

        unsigned int m = t->m;     /* for icc ... */
        unsigned int n = t->n;
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < m ; i++)
            for(unsigned int j = 0 ; j < n ; j++) {
                size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
                void * tij = pointer_arith(t->data, offset);
                abvec aij = matpoly_part(ab, a, i, j, 0);
                /* ok, casting like this is a crude hack ! */
                fft_do_ift_fppol_mp((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p, shift);
            }
        free(tt);
    }

}

void matpoly_mul_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, const struct lingen_substep_schedule * S MAYBE_UNUSED)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);

    matpoly_ft tc, ta, tb;
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol(fti, p, a->size, b->size, a->n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }

    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, csize);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    matpoly_ft_dft(ab, ta, a, fti);
    matpoly_ft_dft(ab, tb, b, fti);
    matpoly_ft_mul(ab, tc, ta, tb, fti);
    c->size = csize;
    ASSERT_ALWAYS(c->size <= c->alloc);
    matpoly_ft_ift(ab, c, tc, fti);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);
}/*}}}*/

void matpoly_mp_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, const struct lingen_substep_schedule * S MAYBE_UNUSED)/*{{{*/
{
    matpoly_ft tc, ta, tb;
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol_mp(fti, p, MIN(a->size, b->size), MAX(a->size, b->size), a->n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }

    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, MAX(a->size, b->size) - MIN(a->size, b->size) + 1);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    matpoly_ft_dft(ab, ta, a, fti);
    matpoly_ft_dft(ab, tb, b, fti);
    matpoly_ft_mul(ab, tc, ta, tb, fti);
    c->size = MAX(a->size, b->size) - MIN(a->size, b->size) + 1;
    ASSERT_ALWAYS(c->size <= c->alloc);
    matpoly_ft_ift_mp(ab, c, tc, MIN(a->size, b->size) - 1, fti);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);
}/*}}}*/
