#ifndef CADO_LINGEN_MATPOLY_FT_HPP
#define CADO_LINGEN_MATPOLY_FT_HPP

#include "cado_config.h"

#include <climits>
#include <cstring>

#include <algorithm>
#include <array>

#include <gmp.h>

#include "gmp_aux.h"
#include "arith-hard.hpp"
#include "lingen_call_companion.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_matpoly_select.hpp"
#include "lingen_memory_pool.hpp"
#include "macros.h"
#include "misc.h"
#include "omp_proxy.h"
#include "submatrix_range.hpp"

class tree_stats;


template<typename fft_type>
class matpoly_ft {
    static constexpr bool is_binary = is_binary_fft<fft_type>::value;
    typedef typename fft_type::ptr ptr;
    typedef typename fft_type::srcptr srcptr;
    typedef memory_pool_wrapper<ptr, true> memory_pool_type;
    static memory_pool_type memory;
public:
    struct memory_guard : private memory_pool_type::guard_base {
        memory_guard(size_t s) : memory_pool_type::guard_base(memory, s) {}
        ~memory_guard() { memory_pool_type::guard_base::pre_dtor(memory); }
    };

    fft_type fti;
    unsigned int m = 0;
    unsigned int n = 0;
    std::array<size_t, 3> fft_alloc_sizes;
    ptr data = NULL;
    unsigned int nrows() const { return m; }
    unsigned int ncols() const { return n; }

    bool check_pre_init() const { return data == NULL; }
    /* {{{ ctor / dtor */
    matpoly_ft(fft_type const & fti)
        : fti(fti)
        , fft_alloc_sizes(fti.get_alloc_sizes())
    { }
    matpoly_ft(unsigned int m, unsigned int n, fft_type const & fti)
        : matpoly_ft(fti)
    {
        this->m = m;
        this->n = n;
        data = memory.alloc(m * n * fft_alloc_sizes[0]);
        memset(data, 0, m * n * fft_alloc_sizes[0]);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                fti.prepare(part(i,j));
            }
        }
    }
    matpoly_ft() = default;
    matpoly_ft(matpoly_ft const&) = delete;
    matpoly_ft& operator=(matpoly_ft const&) = delete;
    matpoly_ft(matpoly_ft && a)
        : fti(a.fti)
        , m(a.m)
        , n(a.n)
        , fft_alloc_sizes(a.fft_alloc_sizes)
    {
        data=a.data;
        a.m=a.n=0;
        a.data=NULL;
    }
    matpoly_ft& operator=(matpoly_ft && a) {
        if (data) memory.free(data, m * n * fft_alloc_sizes[0]);
        fft_alloc_sizes = a.fft_alloc_sizes;
        fti = a.fti;
        m=a.m;
        n=a.n;
        data=a.data;
        a.m=a.n=0;
        a.data=NULL;
        return *this;
    }
    ~matpoly_ft() {
        if (data) memory.free(data, m * n * fft_alloc_sizes[0]);
    }
    /* }}} */
    /* {{{ direct access interface */
    ptr part(unsigned int i, unsigned int j) {
        return (ptr) pointer_arith(data, (i*n+j) * fft_alloc_sizes[0]);
    }
    srcptr part(unsigned int i, unsigned int j) const {
        return (ptr) pointer_arith(data, (i*n+j) * fft_alloc_sizes[0]);
    }
    /* }}} */
    /* {{{ views -- some of the modifiers are implemented here only */
    struct view_t : public submatrix_range {/*{{{*/
        matpoly_ft & M;
        view_t(matpoly_ft & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly_ft & M) : submatrix_range(M), M(M) {}
        ptr part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
        void zero() { /*{{{*/
            unsigned int nrows = this->nrows();
            unsigned int ncols = this->ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    M.fti.zero(part(i,j));
                }
            }
        }/*}}}*/
#if 0
        void fill_random(cxx_gmp_randstate & rstate) { /*{{{*/
            unsigned int nrows = this->nrows();
            unsigned int ncols = this->ncols();
            // the rstate is shared: it is not safe to openmp-it.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    M.fti.fill_random(part(i,j), rstate);
                }
            }
        }/*}}}*/
#endif

        void to_export() { /*{{{*/
            unsigned int nrows = this->nrows();
            unsigned int ncols = this->ncols();
            ASSERT(check());
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    M.fti.to_export(part(i,j));
                }
            }
        } /*}}}*/
        void to_import() { /*{{{*/
            unsigned int nrows = this->nrows();
            unsigned int ncols = this->ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    M.fti.to_import(part(i,j));
                }
            }
            ASSERT(check());
        }/*}}}*/
        int check() { return const_view_t(*this).check(); }
    };/*}}}*/
    struct const_view_t : public submatrix_range {/*{{{*/
        matpoly_ft const & M;
        const_view_t(matpoly_ft const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly_ft const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
        bool check() {/*{{{*/
            for(unsigned int i = 0 ; i < nrows() ; i++) {
                for(unsigned int j = 0 ; j < ncols() ; j++) {
                    if (!M.fti.check(part(i,j), 1))
                        return false;
                }
            }
            return true;
        }/*}}}*/
    };/*}}}*/
    /* }}} */
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }

    void zero(submatrix_range R) { view(R).zero(); }
    void zero() { view().zero(); }

#if 0
    void fill_random(submatrix_range const & R, cxx_gmp_randstate & rstate) {
        view(R).fill_random(rstate);
    }
    void fill_random(cxx_gmp_randstate & rstate) {
        view().fill_random(rstate);
    }
#endif

    void to_import(submatrix_range const & R) { view(R).to_import(); }
    void to_import() { view().to_import(); }

    void to_export(submatrix_range const & R) { view(R).to_import(); }
    void to_export() { view().to_export(); }

    int check(submatrix_range const & R) const { return view(R).check(); }
    int check() const { return view().check(); }

    static void dft(view_t t, typename matpoly<is_binary>::const_view_t a)/*{{{*/
    {
        unsigned int nrows = a.nrows();
        unsigned int ncols = a.ncols();
        ASSERT_ALWAYS(t.nrows() == nrows);
        ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
        unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
        {
            ptr tt = memory.alloc(t.M.fft_alloc_sizes[1]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    ptr tij = t.part(i, j);
                    typename matpoly<is_binary>::srcptr aij = a.part(i, j);
                    /* ok, casting like this is a crude hack ! */
                    t.M.fti.dft(tij, (const mp_limb_t *) aij, a.M.get_size(), tt);
                }
            }
            memory.free(tt, t.M.fft_alloc_sizes[1]);
        }
    }/*}}}*/
    static void ift(typename matpoly<is_binary>::view_t a, view_t t)/*{{{*/
    {
        unsigned int nrows = a.nrows();
        unsigned int ncols = a.ncols();
        ASSERT_ALWAYS(t.nrows() == nrows);
        ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
        unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
        {
            ptr tt = memory.alloc(t.M.fft_alloc_sizes[1]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    ptr tij = t.part(i,j);
                    typename matpoly<is_binary>::ptr aij = a.part(i, j);
                    /* ok, casting like this is a crude hack ! */
                    t.M.fti.ift((mp_limb_t *) aij, a.M.get_size(), tij, tt);
                }
            }
            memory.free(tt, t.M.fft_alloc_sizes[1]);
        }
    }/*}}}*/
    static void addcompose(view_t t, const_view_t t0, const_view_t t1)/*{{{*/
    {
        unsigned int nrows = t.nrows();
        unsigned int ncols = t.ncols();
        unsigned int nadd = t0.ncols();
        ASSERT_ALWAYS(t0.nrows() == nrows);
        ASSERT_ALWAYS(t1.ncols() == ncols);
        ASSERT_ALWAYS(t0.ncols() == nadd);
        ASSERT_ALWAYS(t1.nrows() == nadd);
        ASSERT(t0.check());
        ASSERT(t1.check());
        ASSERT(t.check());
#ifdef HAVE_OPENMP
        unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
        {
            ptr qt = memory.alloc(t.M.fft_alloc_sizes[1]);
            ptr tt = memory.alloc(t.M.fft_alloc_sizes[2]);
                
            memset(qt, 0, t.M.fft_alloc_sizes[1]);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int i = 0 ; i < nrows ; i++) {
                for(unsigned int j = 0 ; j < ncols ; j++) {
                    memset(tt, 0, t.M.fft_alloc_sizes[2]);
                    for(unsigned int k = 0 ; k < nadd ; k++) {
                        t.M.fti.addcompose(t.part(i,j), t0.part(i,k), t1.part(k,j), tt, qt);
                    }
                }
            }
            memory.free(tt, t.M.fft_alloc_sizes[2]);
            memory.free(qt, t.M.fft_alloc_sizes[1]);
        }
    }/*}}}*/

    static void mul(view_t t, const_view_t t0, const_view_t t1)
    {
        t.zero();
        addcompose(t, t0, t1);
    }

    /*
void add(matpoly_ft::view_t t, matpoly_ft::const_view_t t0, matpoly_ft::const_view_t t1)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
    ASSERT_ALWAYS(t0.nrows() == nrows);
    ASSERT_ALWAYS(t0.ncols() == ncols);
    ASSERT_ALWAYS(t1.nrows() == nrows);
    ASSERT_ALWAYS(t1.ncols() == ncols);

#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            fft_add(t.part(i,j), t0.part(i,j), t1.part(i,j), t.M.fti);
        }
    }
}
*/

    /* In a way, this is the only real API exported by this module */

    static auto mp_caching_adj(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M) -> matpoly<is_binary>;
    static auto mul_caching_adj(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M) -> matpoly<is_binary>;
    static matpoly<is_binary> mp_caching(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M) {
        return mp_caching_adj(stats, a, b, UINT_MAX, M);
    }
    static matpoly<is_binary> mul_caching(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M) {
        return mul_caching_adj(stats, a, b, UINT_MAX, M);
    }
};

#ifdef LINGEN_BINARY
extern template class matpoly_ft<gf2x_fake_fft_info>;
extern template class matpoly_ft<gf2x_cantor_fft_info>;
extern template class matpoly_ft<gf2x_ternary_fft_info>;
#else
extern template class matpoly_ft<fft_transform_info>;
#endif


#endif	/* LINGEN_MATPOLY_FT_HPP_ */
