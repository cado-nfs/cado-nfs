#include "cado.h"
#include <stdlib.h>
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.h"
#include "lingen-matpoly-ft.h"
#include "logline.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* timings made on cochon, rev 6877b97 (buggy; fixed in 76dde6c) */
#define MP_FTI_DEPTH_ADJ_24_36_36 { { 1, 6 }, { 2, 4 }, { 3, 3 }, { 4, 2 }, { 10, 1 }, { 11, 2 }, { 13, 1 }, { 22, 2 }, { 28, 1 }, { 32, 2 }, { 33, 1 }, { 38, 0 }, { 39, 1 }, { 54, 0 }, { 55, 1 }, { 64, 0 }, { 65, 1 }, { 66, 0 }, { 103, 1 }, { 104, 0 }, { 107, 1 }, { 114, 0 }, { 115, 1 }, { 129, 0 }, }

#define MUL_FTI_DEPTH_ADJ_36_36_36 { { 1, 6 }, { 2, 3 }, { 3, 2 }, { 6, 1 }, { 7, 2 }, { 14, 1 }, { 23, 0 }, { 26, 1 }, { 44, 0 }, { 46, 1 }, { 54, 0 }, { 61, 1 }, { 62, 0 }, }


void matpoly_ft_init(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, unsigned int m, unsigned int n, struct fft_transform_info * fti)
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
}

void matpoly_ft_clear(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti MAYBE_UNUSED)
{
    free(t->data);
    memset(t, 0, sizeof(*t));
}

/* define some machinery to have loops on matrices that can do both the
 * draft and the non-draft case. The annoying thing is that the inherent
 * break-early logic of the draft mode is completely orthogonal to the
 * openmp philosophy, and we'd rather not blend the two at runtime
 * because that would probably be a performance hit (plus that would mean
 * duplicating some non-trivial code).
 */

template<bool draft_mode, typename T>/*{{{*/
struct ft_wrap
{
    double tt0=0;
    int draft_cutoff = 0;

    ft_wrap(ft_wrap const &) = delete;
    ft_wrap(int draft) : tt0(wct_seconds()), draft_cutoff(draft) {}

    bool isdone() const {
        return wct_seconds() - tt0 >= draft_cutoff;
    }

    double operator()(T & t) {
        double yield = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel reduction(+:yield)
#endif
        {
            /* instantiate the temporaries */
            typename T::template instance<T> x(t);
            unsigned int ndone = 0;
#ifdef HAVE_OPENMP
            unsigned int k0 = omp_get_thread_num();
            unsigned int dk = omp_get_num_threads();
#else
            unsigned int k0 = 0;
            unsigned int dk = 1;
#endif
            for(unsigned int k = k0; k < t.m() * t.n() ; k += dk) {
                unsigned int i = k / t.n();
                unsigned int j = k % t.n();
                x.doij(i,j);
                ++ndone;
                if (isdone()) break;
            }
            yield += ndone / (wct_seconds() - tt0);
        }
        return (t.m() * t.n() / yield) - (wct_seconds() - tt0);
    }
};/*}}}*/

template<typename T>/*{{{*/
struct ft_wrap<false, T>
{
    ft_wrap(ft_wrap const &) = delete;
    ft_wrap(int draft) {
        ASSERT_ALWAYS(draft == 0);
    }

    double operator()(T & t) {
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            /* instantiate the temporaries */
            typename T::template instance<T> x(t);
            unsigned int m = t.m();     /* for icc ... */
            unsigned int n = t.n();
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int i = 0 ; i < m ; i++)
                for(unsigned int j = 0 ; j < n ; j++)
                    x.doij(i, j);
        }
        return 0;
    }
};/*}}}*/

struct ft_base {/*{{{*/
    abdst_field ab;
    struct fft_transform_info * fti;
    size_t fft_alloc_sizes[3];
    int draft;

    ft_base(abdst_field ab, struct fft_transform_info * fti, int draft) : ab(ab), fti(fti), draft(draft) 
    {
        fft_get_transform_allocs(fft_alloc_sizes, fti);
    }

    template<typename T>
    double choose(T & t) {
        if (draft) {
            return ft_wrap<true, T>(draft)(t);
        } else {
            return ft_wrap<false, T>(draft)(t);
        }
    }
};/*}}}*/

struct matpoly_ft_dft_code : public ft_base {
    matpoly_ft_ptr t;
    matpoly_srcptr a;
    unsigned int m() const { return a->m; }
    unsigned int n() const { return a->n; }
    template<typename T>
    struct instance : public T {
        using T::fft_alloc_sizes;
        void * tt;
        instance(T & P) : T(P) {
            tt = malloc(fft_alloc_sizes[1]);
        }
        instance(instance const&) = delete;
        ~instance() { free(tt); }
        using T::t;
        using T::a;
        using T::ab;
        using T::fti;
        void doij(unsigned int i, unsigned int j) {
            size_t offset = (i*t->n + j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            absrc_vec aij = matpoly_part_const(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_dft_fppol(tij, (const mp_limb_t *) aij, a->size, tt, fti, ab->p);
        }
    };
    
    matpoly_ft_dft_code(abdst_field ab, matpoly_ft_ptr t, matpoly_srcptr a, struct fft_transform_info * fti, int draft) : ft_base(ab, fti, draft), t(t), a(a) {
        ASSERT_ALWAYS(t->m == a->m);
        ASSERT_ALWAYS(t->n == a->n);
    }
    double operator()() { return choose(*this); }
};

double matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr a, struct fft_transform_info * fti, int draft)
{
    return matpoly_ft_dft_code(ab, t, a, fti, draft)();
}

void matpoly_ft_zero(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
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

void matpoly_ft_export(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
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

void matpoly_ft_import(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
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

void matpoly_ft_add(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti)
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

struct matpoly_ft_addmul_code : public ft_base {
    matpoly_ft_ptr u;
    matpoly_ft_srcptr t0;
    matpoly_ft_srcptr t1;
    unsigned int m() const { return t0->m; }
    unsigned int n() const { return t1->n; }
    template<typename T>
    struct instance : public T {
        using T::fft_alloc_sizes;
        void * qt;
        void * tt;
        instance(T & P) : T(P) {
            qt = malloc(fft_alloc_sizes[1]);
            tt = malloc(fft_alloc_sizes[2]);
            memset(qt, 0, fft_alloc_sizes[1]);
        }
        instance(instance const&) = delete;
        ~instance() { free(tt); free(qt); }
        using T::u;
        using T::t0;
        using T::t1;
        using T::ab;
        using T::fti;
        void doij(unsigned int i, unsigned int j) {
            size_t tsize = fft_alloc_sizes[0];
            memset(tt, 0, fft_alloc_sizes[2]);
            for(unsigned int k = 0 ; k < t0->n ; k++) {
                void * uij  = pointer_arith(u->data, (i*u->n+j) * tsize);
                void * t0ik = pointer_arith(t0->data, (i*T::t0->n+k) * tsize);
                void * t1kj = pointer_arith(t1->data, (k*T::t1->n+j) * tsize);
                fft_addmul(uij, t0ik, t1kj, tt, qt, fti);
            }
        }
    };
    
    matpoly_ft_addmul_code(abdst_field ab, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti, int draft) : ft_base(ab, fti, draft), u(u), t0(t0), t1(t1) {
        ASSERT_ALWAYS(t0->n == t1->m);
        ASSERT_ALWAYS(t0->m == u->m);
        ASSERT_ALWAYS(t1->n == u->n);
    }
    double operator()() { return choose(*this); }
};

double matpoly_ft_addmul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti, int draft)
{
    return matpoly_ft_addmul_code(ab, u, t0, t1, fti, draft)();
}

double matpoly_ft_mul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti, int draft)
{
    matpoly_ft_zero(ab, u, fti);
    return matpoly_ft_addmul(ab, u, t0, t1, fti, draft);
}

/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
*/

struct matpoly_ft_ift_code : public ft_base {
    matpoly_ptr a;
    matpoly_ft_ptr t;
    unsigned int m() const { return t->m; }
    unsigned int n() const { return t->n; }
    template<typename T>
    struct instance : public T {
        using T::fft_alloc_sizes;
        void * tt;
        instance(T & P) : T(P) {
            tt = malloc(fft_alloc_sizes[1]);
        }
        instance(instance const&) = delete;
        ~instance() { free(tt); }
        using T::t;
        using T::a;
        using T::ab;
        using T::fti;
        void doij(unsigned int i, unsigned int j) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_ift_fppol((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p);
        }
    };
    
    matpoly_ft_ift_code(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, struct fft_transform_info * fti, int draft) : ft_base(ab, fti, draft), a(a), t(t) {
        ASSERT_ALWAYS(t->m == a->m);
        ASSERT_ALWAYS(t->n == a->n);
    }
    double operator()() { return choose(*this); }
};

double matpoly_ft_ift(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, struct fft_transform_info * fti, int draft)
{
    return matpoly_ft_ift_code(ab, a, t, fti, draft)();
}

struct matpoly_ft_ift_mp_code : public ft_base {
    matpoly_ptr a;
    matpoly_ft_ptr t;
    unsigned int shift;
    unsigned int m() const { return t->m; }
    unsigned int n() const { return t->n; }
    template<typename T>
    struct instance : public T {
        using T::fft_alloc_sizes;
        void * tt;
        instance(T & P) : T(P) {
            tt = malloc(fft_alloc_sizes[1]);
        }
        instance(instance const&) = delete;
        ~instance() { free(tt); }
        using T::t;
        using T::a;
        using T::shift;
        using T::ab;
        using T::fti;
        void doij(unsigned int i, unsigned int j) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_ift_fppol_mp((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p, shift);
        }
    };
    
    matpoly_ft_ift_mp_code(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti, int draft) : ft_base(ab, fti, draft), a(a), t(t), shift(shift) {
        ASSERT_ALWAYS(t->m == a->m);
        ASSERT_ALWAYS(t->n == a->n);
    }
    double operator()() { return choose(*this); }
};

double matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti, int draft)
{
    return matpoly_ft_ift_mp_code(ab, a, t, shift, fti, draft)();
}

double matpoly_mul_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, int draft)/*{{{*/
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

    double x0 = 0;

    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, csize);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    x0 += matpoly_ft_dft(ab, ta, a, fti, draft);
    x0 += matpoly_ft_dft(ab, tb, b, fti, draft);
    x0 += matpoly_ft_mul(ab, tc, ta, tb, fti, draft);
    c->size = csize;
    ASSERT_ALWAYS(c->size <= c->alloc);
    x0 += matpoly_ft_ift(ab, c, tc, fti, draft);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);

    return x0;
}/*}}}*/
double matpoly_mp_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj, int draft)/*{{{*/
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

    double x0 = 0;

    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, MAX(a->size, b->size) - MIN(a->size, b->size) + 1);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    x0 += matpoly_ft_dft(ab, ta, a, fti, draft);
    x0 += matpoly_ft_dft(ab, tb, b, fti, draft);
    x0 += matpoly_ft_mul(ab, tc, ta, tb, fti, draft);
    c->size = MAX(a->size, b->size) - MIN(a->size, b->size) + 1;
    ASSERT_ALWAYS(c->size <= c->alloc);
    x0 += matpoly_ft_ift_mp(ab, c, tc, MIN(a->size, b->size) - 1, fti, draft);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);

    return x0;
}/*}}}*/
