#include "cado.h"
#include <cstdlib>
#include <sstream>
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.hpp"
#include "lingen-matpoly-ft.hpp"
#include "logline.h"
#include "sha1.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* timings made on cochon, rev 6877b97 (buggy; fixed in 76dde6c) */
#define MP_FTI_DEPTH_ADJ_24_36_36 { { 1, 6 }, { 2, 4 }, { 3, 3 }, { 4, 2 }, { 10, 1 }, { 11, 2 }, { 13, 1 }, { 22, 2 }, { 28, 1 }, { 32, 2 }, { 33, 1 }, { 38, 0 }, { 39, 1 }, { 54, 0 }, { 55, 1 }, { 64, 0 }, { 65, 1 }, { 66, 0 }, { 103, 1 }, { 104, 0 }, { 107, 1 }, { 114, 0 }, { 115, 1 }, { 129, 0 }, }

#define MUL_FTI_DEPTH_ADJ_36_36_36 { { 1, 6 }, { 2, 3 }, { 3, 2 }, { 6, 1 }, { 7, 2 }, { 14, 1 }, { 23, 0 }, { 26, 1 }, { 44, 0 }, { 46, 1 }, { 54, 0 }, { 61, 1 }, { 62, 0 }, }

matpoly_ft::memory_pool matpoly_ft::memory;

matpoly_ft::memory_pool_guard::memory_pool_guard(size_t s)
{
    ASSERT_ALWAYS(memory.allocated == 0);
    memory.peak = 0;
    memory.allowed = s;
}
matpoly_ft::memory_pool_guard::~memory_pool_guard() {
    ASSERT_ALWAYS(matpoly_ft::memory.allocated == 0);
}

void * matpoly_ft::memory_pool::alloc(size_t s)
{
    std::lock_guard<std::mutex> dummy(mm);
    ASSERT_ALWAYS(allocated + s <= allowed);
    allocated += s;
    if (allocated > peak) peak = allocated;
    return malloc(s);
}
void matpoly_ft::memory_pool::free(void * p, size_t s)
{
    std::lock_guard<std::mutex> dummy(mm);
    ASSERT_ALWAYS(allocated >= s);
    allocated -= s;
    ::free(p);
}


matpoly_ft::matpoly_ft(abdst_field ab, unsigned int m, unsigned int n, const struct fft_transform_info * fti)/*{{{*/
    : ab(ab)
    , m(m)
    , n(n)
    , fti(fti)
{
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    data = memory.alloc(m * n * fft_alloc_sizes[0]);
    memset(data, 0, m * n * fft_alloc_sizes[0]);
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            fft_transform_prepare(part(i,j), fti);
        }
    }
}/*}}}*/

matpoly_ft::~matpoly_ft() /* {{{ */
{
    if (data)
        memory.free(data, m * n * fft_alloc_sizes[0]);
}/*}}}*/

matpoly_ft::matpoly_ft(matpoly_ft && a)
{
    fti = a.fti;
    ab = a.ab;
    memcpy(fft_alloc_sizes, a.fft_alloc_sizes, sizeof(fft_alloc_sizes));
    m=a.m;
    n=a.n;
    data=a.data;
    a.m=a.n=0;
    a.data=NULL;
}
matpoly_ft& matpoly_ft::operator=(matpoly_ft&& a)
{
    if (data)
        memory.free(data, m * n * fft_alloc_sizes[0]);
    fti = a.fti;
    ab = a.ab;
    memcpy(fft_alloc_sizes, a.fft_alloc_sizes, sizeof(fft_alloc_sizes));
    m=a.m;
    n=a.n;
    data=a.data;
    a.m=a.n=0;
    a.data=NULL;
    return *this;
}

void dft(matpoly_ft::view_t t, matpoly::const_view_t a)
{
    ASSERT_ALWAYS(t.M.fti);
    unsigned int nrows = a.nrows();
    unsigned int ncols = a.ncols();
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(t.M.fft_alloc_sizes[1]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                void * tij = t.part(i, j);
                absrc_vec aij = a.part(i, j);
                /* ok, casting like this is a crude hack ! */
                fft_do_dft_fppol(tij, (const mp_limb_t *) aij, a.M.size, tt, t.M.fti, t.M.ab->p);
            }
        }
        free(tt);
    }
}

void zero(matpoly_ft::view_t t)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            fft_zero(t.part(i,j), t.M.fti);
        }
    }
}

void fill_random(matpoly_ft::view_t t, gmp_randstate_t rstate)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            fft_fill_random(t.part(i,j), t.M.fti, rstate);
        }
    }
}

void to_export(matpoly_ft::view_t t)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            fft_transform_export(t.part(i,j), t.M.fti);
        }
    }
}

void to_import(matpoly_ft::view_t t)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            fft_transform_import(t.part(i,j), t.M.fti);
        }
    }
}

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

void addmul(matpoly_ft::view_t t, matpoly_ft::const_view_t t0, matpoly_ft::const_view_t t1)
{
    unsigned int nrows = t.nrows();
    unsigned int ncols = t.ncols();
    unsigned int nadd = t0.ncols();
    ASSERT_ALWAYS(t0.nrows() == nrows);
    ASSERT_ALWAYS(t1.ncols() == ncols);
    ASSERT_ALWAYS(t0.ncols() == nadd);
    ASSERT_ALWAYS(t1.nrows() == nadd);
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * qt = malloc(t.M.fft_alloc_sizes[1]);
        void * tt = malloc(t.M.fft_alloc_sizes[2]);
        memset(qt, 0, t.M.fft_alloc_sizes[1]);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                memset(tt, 0, t.M.fft_alloc_sizes[2]);
                for(unsigned int k = 0 ; k < nadd ; k++) {
                    fft_addmul(t.part(i,j), t0.part(i,k), t1.part(k,j), tt, qt, t.M.fti);
                }
            }
        }
        free(tt);
        free(qt);
    }
}

void mul(matpoly_ft::view_t t, matpoly_ft::const_view_t t0, matpoly_ft::const_view_t t1)
{
    zero(t);
    addmul(t, t0, t1);
}

/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, const struct fft_transform_info * fti);
*/

void ift(matpoly::view_t a, matpoly_ft::view_t t)
{
    ASSERT_ALWAYS(t.M.fti);
    unsigned int nrows = a.nrows();
    unsigned int ncols = a.ncols();
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(t.M.fft_alloc_sizes[1]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                void * tij = t.part(i,j);
                abdst_vec aij = a.part(i, j);
                /* ok, casting like this is a crude hack ! */
                fft_do_ift_fppol((mp_limb_t *) aij, a.M.size, tij, tt, t.M.fti, t.M.ab->p);
            }
        }
        free(tt);
    }
}

void ift_mp(matpoly::view_t a, matpoly_ft::view_t t, unsigned int shift)
{
    ASSERT_ALWAYS(t.M.fti);
    unsigned int nrows = a.nrows();
    unsigned int ncols = a.ncols();
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
        void * tt = malloc(t.M.fft_alloc_sizes[1]);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                void * tij = t.part(i, j);
                abvec aij = a.part(i, j);
                /* ok, casting like this is a crude hack ! */
                fft_do_ift_fppol_mp((mp_limb_t *) aij, a.M.size, tij, tt, t.M.fti, t.M.ab->p, shift);
            }
        }
        free(tt);
    }

}

#if 0
std::string sha1sum_block(matpoly const & X, size_t s)
{
    std::stringstream ss;
    for(unsigned int i = 0 ; i < X.m ; i++) {
        ss << "[";
        for(unsigned int j = 0 ; j < X.n ; j++) {
            sha1_checksumming_stream S;
            S.write((const char *) X.part(i,j), X.size*sizeof(mp_limb_t));
            char checksum[41];
            S.checksum(checksum);
            if (s < 41) checksum[s]='\0';
            ss << " " << checksum;
        }
        ss << " ]\n";
    }
    return ss.str();
}

static inline mp_size_t fti_rsize0(const struct fft_transform_info * fti)
{
    mp_size_t w = fti->w;
    mp_size_t n = 1 << fti->depth;
    mp_bitcnt_t nw = (mp_bitcnt_t) n * w;
    mp_size_t rsize0 = nw/FLINT_BITS;  /* need rsize0+1 words for x\in R */
    return rsize0;
}

std::string sha1sum_block(matpoly_ft const & X, size_t s)
{
    std::ostringstream ss;
    for(unsigned int i = 0 ; i < X.m ; i++) {
        ss << "[";
        for(unsigned int j = 0 ; j < X.n ; j++) {
            const void * T = X.part(i,j);
            char checksum[41];
            if (!fft_transform_check(T, X.fti, 0)) {
                for(int i = 0 ; i < 40 ; checksum[i++]='!');
                checksum[40]='\0';
            } else {
                mp_size_t rsize0 = fti_rsize0(X.fti);
                mp_limb_t ** Tptrs = (mp_limb_t**) T;
                sha1_checksumming_stream S;
                unsigned int depth = X.fti->depth;
                unsigned int n = 1 << depth;
                for(unsigned int k = 0 ; k < 4*n+2 ; k++) {
                    S.write((const char *) Tptrs[k], (rsize0+1)*sizeof(mp_limb_t));
                }
                S.checksum(checksum);
            }
            if (s < 41) checksum[s]='\0';
            ss << " " << checksum;
        }
        ss << " ]\n";
    }
    return ss.str();
}
#endif

/* middle product and multiplication are really the same thing, so better
 * avoid code duplication */

struct op_mul {/*{{{*/
    size_t csize;
    op_mul(matpoly const & a, matpoly const & b, unsigned int adj, fft_transform_info * fti)
    {
        csize = a.size + b.size; csize -= (csize > 0);
        fft_get_transform_info_fppol(fti, a.ab->p, a.size, b.size, a.n);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t)
    {
        ::ift(a, t);
    }
};/*}}}*/
struct op_mp {/*{{{*/
    size_t csize;
    unsigned int shift;
    op_mp(matpoly const & a, matpoly const & b, unsigned int adj, fft_transform_info * fti)
    {
        csize = MAX(a.size, b.size) - MIN(a.size, b.size) + 1;
        shift = MIN(a.size, b.size) - 1;
        fft_get_transform_info_fppol_mp(fti, a.ab->p, MIN(a.size, b.size), MAX(a.size, b.size), a.n);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t)
    {
        ::ift_mp(a, t, shift);
    }
};/*}}}*/

template<typename T>
static void mp_or_mul(T& OP, matpoly & c, matpoly const & a, matpoly const & b, const struct fft_transform_info * fti, const struct lingen_substep_schedule * S)/*{{{*/
{
    if (c.m != a.m || c.n != a.n || c.alloc != OP.csize)
        c = matpoly(a.ab, a.m, b.n, OP.csize);

    const unsigned int r = 1; // only for notational consistency w/ bigmatpoly

    subdivision mpi_split(a.n, r);
    subdivision shrink0_split(a.m, S->shrink0);
    subdivision shrink2_split(b.n, S->shrink2);
    /* The order in which we do the transforms is not really our main
     * concern at this point. If sharing makes sense, then probably
     * shrink0 and shrink2 do not. So they're serving opposite purposes.
     */
    /* Declare ta, tb, tc early on so that we don't malloc/free n times.
     */
    matpoly_ft ta,tb,tc;
    const unsigned int nr1 = mpi_split.block_size_upper_bound();
    const unsigned int rank = 0;
    ASSERT_ALWAYS(rank < r);
    bool inner_is_row_major;
    {
        /* first, upper bounds on nrs0 and nrs2 */
        unsigned int nrs0 = iceildiv(a.m, S->shrink0);
        unsigned int nrs2 = iceildiv(b.n, S->shrink2);
        tc = matpoly_ft (c.ab, nrs0, nrs2, fti);
        /* We must decide on an ordering beforehand. We cannot do this
         * dynamically because of rounding issues: e.g. for 13=7+6, we
         * will do both 7*6 and 6*7 in the inner loops.
         */
        inner_is_row_major = nrs0 < nrs2;
        if (inner_is_row_major) {
            ta = matpoly_ft(a.ab, nrs0, r * S->batch, fti);
            tb = matpoly_ft(a.ab, r * S->batch, 1, fti);
        } else {
            ta = matpoly_ft(a.ab, 1, r * S->batch, fti);
            tb = matpoly_ft(a.ab, r * S->batch, nrs2, fti);
        }
    }
    unsigned int k0mpi,k1mpi;
    std::tie(k0mpi, k1mpi) = mpi_split.nth_block(rank);
    for(unsigned int round0 = 0 ; round0 < S->shrink0 ; round0++) {
        unsigned int i0,i1;
        std::tie(i0, i1) = shrink0_split.nth_block(round0);
        unsigned int nrs0 = i1-i0;
        for(unsigned int round2 = 0 ; round2 < S->shrink2 ; round2++) {
            unsigned int j0,j1;
            std::tie(j0, j1) = shrink2_split.nth_block(round2);
            unsigned int nrs2 = j1-j0;
            submatrix_range Rc(i0,j0,i1-i0,j1-j0);
            submatrix_range Rct(0,0,i1-i0,j1-j0);

            /* Now do a subblock */
            tc.zero();
            for(unsigned int k = 0 ; k < nr1 ; k += S->batch) {
                unsigned int k0 = k0mpi + k;
                unsigned int k1 = MIN(k1mpi, k0 + S->batch);
                if (inner_is_row_major) {
                    submatrix_range Ra(i0,k0,nrs0,k1-k0);
                    submatrix_range Rat(0,rank*S->batch, nrs0,k1-k0);
                    ta.zero();  // for safety because of rounding.
                    dft(ta.view(Rat), a.view(Ra));
                    // allgather ta among r nodes.
                    for(unsigned int j = 0 ; j < nrs2 ; j++) {
                        tb.zero();
                        dft(tb.view(submatrix_range(rank*S->batch,0,k1-k0,1)),
                                b.view(submatrix_range(k0,j0+j,k1-k0,1)));
                        // allgather tb among r nodes
                        // rounding might surprise us.
                        addmul(tc.view(submatrix_range(0,j,nrs0,1)),
                                ta.view(submatrix_range(0,0,nrs0,ta.ncols())),
                                tb.view(submatrix_range(0,0,tb.nrows(),1)));
                    }
                } else {
                    submatrix_range Rb(k0,j0,k1-k0,nrs2);
                    submatrix_range Rbt(rank*S->batch,0,k1-k0,nrs2);
                    tb.zero();
                    dft(tb.view(Rbt), b.view(Rb));
                    // allgather tb among r nodes
                    for(unsigned int i = 0 ; i < nrs0 ; i++) {
                        ta.zero();
                        dft(ta.view(submatrix_range(0,rank*S->batch,1,k1-k0)),
                                a.view(submatrix_range(i0+i,k0,1,k1-k0)));
                        // allgather ta among r nodes
                        addmul(tc.view(submatrix_range(i,0,1,nrs2)),
                                ta.view(submatrix_range(0,0,1,ta.ncols())),
                                tb.view(submatrix_range(0,0,tb.nrows(),nrs2)));
                    }
                }
            }
            c.size = OP.csize;
            ASSERT_ALWAYS(c.size <= c.alloc);
            OP.ift(c.view(Rc), tc.view(Rct));
        }
    }
}/*}}}*/

void matpoly_mp_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S)/*{{{*/
{
    struct fft_transform_info fti[1];
    op_mp OP(a, b, adj, fti);
    mp_or_mul(OP, c, a, b, fti, S);
} /* }}} */
void matpoly_mul_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S)/*{{{*/
{
    struct fft_transform_info fti[1];
    op_mul OP(a, b, adj, fti);
    mp_or_mul(OP, c, a, b, fti, S);
} /* }}} */

