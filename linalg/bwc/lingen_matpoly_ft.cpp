#include "cado.h"
#include <cstdlib>
#include <sstream>
#include "macros.h"
#include "utils.h"
#include "lingen_matpoly.hpp"
#include "lingen_matpoly_ft.hpp"
#include "logline.h"
#include "sha1.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

/* timings made on cochon, rev 6877b97 (buggy; fixed in 76dde6c) */
#define MP_FTI_DEPTH_ADJ_24_36_36 { { 1, 6 }, { 2, 4 }, { 3, 3 }, { 4, 2 }, { 10, 1 }, { 11, 2 }, { 13, 1 }, { 22, 2 }, { 28, 1 }, { 32, 2 }, { 33, 1 }, { 38, 0 }, { 39, 1 }, { 54, 0 }, { 55, 1 }, { 64, 0 }, { 65, 1 }, { 66, 0 }, { 103, 1 }, { 104, 0 }, { 107, 1 }, { 114, 0 }, { 115, 1 }, { 129, 0 }, }

#define MUL_FTI_DEPTH_ADJ_36_36_36 { { 1, 6 }, { 2, 3 }, { 3, 2 }, { 6, 1 }, { 7, 2 }, { 14, 1 }, { 23, 0 }, { 26, 1 }, { 44, 0 }, { 46, 1 }, { 54, 0 }, { 61, 1 }, { 62, 0 }, }

matpoly_ft::memory_pool matpoly_ft::memory;

matpoly_ft::memory_pool_guard::memory_pool_guard(size_t s) : mysize(s)
{
    oldsize = memory.allowed;
    if (oldsize == SIZE_MAX || s == SIZE_MAX)
        memory.allowed = SIZE_MAX;
    else
        memory.allowed += s;
    if (oldsize == 0)
        ASSERT_ALWAYS(memory.allocated == 0);
        memory.peak = 0;
}
matpoly_ft::memory_pool_guard::~memory_pool_guard() {
    if (oldsize == 0)
        ASSERT_ALWAYS(memory.allocated == 0);
    memory.allowed = oldsize;
    ASSERT_ALWAYS(memory.allocated <= memory.allowed);
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

int matpoly_ft::check() const { return ::check(view()); }

int check(matpoly_ft::const_view_t t)
{
    int ok = 1;
    for(unsigned int i = 0 ; i < t.nrows() ; i++) {
        for(unsigned int j = 0 ; j < t.ncols() ; j++) {
            if (!fft_transform_check(t.part(i,j), t.M.fti, 1))
                ok = 0;
        }
    }
    return ok;
}
int check(matpoly_ft::view_t t)
{
    return check((matpoly_ft::const_view_t) t);
}

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
// the rstate is shared: it is not safe to openmp-it.
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
    ASSERT(check(t));
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
    ASSERT(check(t));
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
    ASSERT(check(t0));
    ASSERT(check(t1));
    ASSERT(check(t));
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

#include "lingen_matpoly_bigmatpoly_ft_common.hpp"

template<> struct OP_CTX<matpoly> : public OP_CTX_base<matpoly> {
    typedef matpoly T;
    template<typename... Args>
    OP_CTX(Args&&... args) : OP_CTX_base<T>(args...) {}
    inline int a_irank() const { return 0; }
    inline int b_irank() const { return 0; }
    inline int a_jrank() const { return 0; }
    inline int b_jrank() const { return 0; }
    inline int mesh_size() const { return 1; }
    static const bool uses_mpi = false;
    inline void mesh_checks() const { }
    void alloc_c_if_needed(size_t size) {
        if (c.m != a.m || c.n != a.n || c.alloc != size)
            c = T(a.ab, a.m, b.n, size);
    }
    inline matpoly const & a_local()const  { return a; }
    inline matpoly const & b_local() const { return b; }
    inline matpoly & c_local() { return c; }
    inline void do_allgather(void *, int) const {}
    inline void begin_smallstep(std::string const &, unsigned int) const { }
    inline void end_smallstep() const {}
    inline void skip_smallstep(std::string const &, unsigned int) const { }
    inline bool local_smallsteps_done() const { return true; }
    template<typename OP> void doit(OP & op, lingen_call_companion::mul_or_mp_times * M) {
        if (M && op.get_transform_ram() > M->per_transform_ram) {
            fprintf(stderr, "Transform size for %s with input operand sizes (%zu, %zu) is %zu, which exceeds expected %zu (anticipated for operand sizes (%zu, %zu). Updating\n",
                    OP::name,
                    a.size,
                    b.size,
                    op.get_transform_ram(),
                    M->per_transform_ram,
                    M->asize,
                    M->bsize
                   );
            size_t ntransforms = M->ram / M->per_transform_ram;
            ASSERT_ALWAYS(M->ram % M->per_transform_ram == 0);
            M->per_transform_ram = op.get_transform_ram();
            M->ram = ntransforms * M->per_transform_ram;
        }
        matpoly_ft::memory_pool_guard dummy(M ? M->ram : SIZE_MAX);
        mp_or_mul(*this, op, op.fti, M ? & M->S : NULL);
    }
};


void matpoly_mp_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M)/*{{{*/
{
    op_mp op(a, b, adj);
    OP_CTX<matpoly>(c, a, b).doit(op, M);
} /* }}} */
void matpoly_mul_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M)/*{{{*/
{
    op_mul op(a, b, adj);
    OP_CTX<matpoly>(c, a, b).doit(op, M);
} /* }}} */

