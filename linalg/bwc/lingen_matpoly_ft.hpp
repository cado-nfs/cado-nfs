#ifndef LINGEN_MATPOLY_FT_H_
#define LINGEN_MATPOLY_FT_H_

#include <mutex>
#include "lingen_matpoly.hpp"
#include "flint-fft/fft.h"
#include "lingen_substep_schedule.h"
#include "tree_stats.hpp"
#include "misc.h"

struct matpoly_ft {
private:
    class memory_pool {
        std::mutex mm;
        public:
        size_t allowed=0;
        size_t allocated=0;
        size_t peak=0;
        void * alloc(size_t);
        void free(void *, size_t);
    };
    static memory_pool memory;
public:
    struct memory_pool_guard {
        size_t oldsize;
        size_t mysize;
        memory_pool_guard(size_t s);
        ~memory_pool_guard();
    };
    abdst_field ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
    const struct fft_transform_info * fti = NULL;
    size_t fft_alloc_sizes[3];
    void * data = NULL;
    inline unsigned int nrows() const { return m; }
    inline unsigned int ncols() const { return n; }

    bool check_pre_init() const { return data == NULL; }
    matpoly_ft(abdst_field ab, unsigned int m, unsigned int n, const struct fft_transform_info * fti);
    matpoly_ft(abdst_field ab, const struct fft_transform_info * fti) : ab(ab), fti(fti) {
        fft_get_transform_allocs(fft_alloc_sizes, fti);
    }
    matpoly_ft() = default;
    matpoly_ft(matpoly_ft const&) = delete;
    matpoly_ft& operator=(matpoly_ft const&) = delete;
    matpoly_ft(matpoly_ft &&);
    matpoly_ft& operator=(matpoly_ft &&);
    ~matpoly_ft();
    int check() const;
    inline void * part(unsigned int i, unsigned int j) {
        return pointer_arith(data, (i*n+j) * fft_alloc_sizes[0]);
    }
    inline const void * part(unsigned int i, unsigned int j) const {
        return pointer_arith(data, (i*n+j) * fft_alloc_sizes[0]);
    }
    void zero(submatrix_range const & R);
    void fill_random(submatrix_range const & R, gmp_randstate_t rstate);
    inline void zero() { zero(view()); }
    inline void fill_random(gmp_randstate_t rstate) { fill_random(view(), rstate); }
    void to_import(submatrix_range const & R);
    void to_export(submatrix_range const & R);
    void to_import() { to_import(view()); }
    void to_export() { to_export(view()); }
#if 0
    void dft(matpoly const & p, submatrix_range);
    void ift(matpoly & p, submatrix_range);
    void ift_mp(matpoly & p, unsigned int shift, submatrix_range);
    void add(matpoly_ft const & t0, matpoly_ft const & t1, submatrix_range);
    // void sub(matpoly_ft const & t0, matpoly_ft const & t1, submatrix_range);
    void mul(matpoly_ft const & t0, matpoly_ft const & t1, submatrix_range R0, submatrix_range R1);
    void addmul(matpoly_ft const & t0, matpoly_ft const & t1, submatrix_range R0, submatrix_range R1);

    inline void dft(matpoly const & p) { dft(p, submatrix_range(*this)); }
    inline void ift(matpoly & p) { ift(p, submatrix_range(*this)); }
    inline void ift_mp(matpoly & p, unsigned int shift) { ift_mp(p, shift, submatrix_range(*this)); }
    inline void add(matpoly_ft const & t0, matpoly_ft const & t1) {
        add(t0, t1, submatrix_range(*this));
    }
    /*
    inline void sub(matpoly_ft const & t0, matpoly_ft const & t1) {
        sub(t0, t1, submatrix_range(*this));
    }
    */
    inline void mul(matpoly_ft const & t0, matpoly_ft const & t1) {
        mul(t0, t1, submatrix_range(t0), submatrix_range(t1));
    }
    inline void addmul(matpoly_ft const & t0, matpoly_ft const & t1) {
        addmul(t0, t1, submatrix_range(t0), submatrix_range(t1));
    }
    inline void to_import() { to_import(submatrix_range(*this)); }
    inline void to_export() { to_export(submatrix_range(*this)); }
#endif

    struct view_t;
    struct const_view_t;

    struct view_t : public submatrix_range {
        matpoly_ft & M;
        view_t(matpoly_ft & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly_ft & M) : submatrix_range(M), M(M) {}
        inline void * part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        inline const void * part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    struct const_view_t : public submatrix_range {
        matpoly_ft const & M;
        const_view_t(matpoly_ft const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly_ft const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline const void * part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }
};

/* In a way, this is the only real API exported by this module */
void matpoly_mul_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S = NULL);

static inline void matpoly_mul_caching(matpoly & c, matpoly const & a, matpoly const & b, const struct lingen_substep_schedule * S = NULL) { return matpoly_mul_caching_adj(c, a, b, UINT_MAX, S); }

void matpoly_mp_caching_adj(matpoly & c, matpoly const & a, matpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S = NULL);

static inline void matpoly_mp_caching(matpoly & c, matpoly const & a, matpoly const & b, const struct lingen_substep_schedule * S = NULL) { return matpoly_mp_caching_adj(c, a, b, UINT_MAX, S); }

void zero(matpoly_ft::view_t);
void fill_random(matpoly_ft::view_t, gmp_randstate_t);
void dft(matpoly_ft::view_t, matpoly::const_view_t);
void ift(matpoly::view_t, matpoly_ft::view_t);
void ift_mp(matpoly::view_t, matpoly_ft::view_t, unsigned int shift);
void add(matpoly_ft::view_t, matpoly_ft::const_view_t, matpoly_ft::const_view_t);
// void sub(matpoly_ft::view_t, matpoly_ft::const_view_t, matpoly_ft::const_view_t);
void mul(matpoly_ft::view_t, matpoly_ft::const_view_t, matpoly_ft::const_view_t);
void addmul(matpoly_ft::view_t, matpoly_ft::const_view_t, matpoly_ft::const_view_t);
void to_import(matpoly_ft::view_t);
void to_export(matpoly_ft::view_t);

int check(matpoly_ft::view_t);
int check(matpoly_ft::const_view_t);
inline void matpoly_ft::zero(submatrix_range const & R) { ::zero(view(R)); }
inline void matpoly_ft::fill_random(submatrix_range const & R, gmp_randstate_t rstate) { ::fill_random(view(R), rstate); }
inline void matpoly_ft::to_import(submatrix_range const & R) { ::to_import(view(R)); }
inline void matpoly_ft::to_export(submatrix_range const & R) { ::to_export(view(R)); }
#endif	/* LINGEN_MATPOLY_FT_H_ */
