#ifndef CADO_LINGEN_POLYMAT_HPP
#define CADO_LINGEN_POLYMAT_HPP

#ifdef LINGEN_BINARY
#error                                                                         \
    "lingen_polymat does not work with binary (maybe use bblas instead -- could end up being the same interface, IDK)"
#endif

#include <cstddef>
#include <climits>
#include <cstdint>

#include <utility>
#include <vector>

#include "arith-hard.hpp"
#include "gmp_aux.h"
#include "macros.h"

template<bool is_binary>
class matpoly;

/* This is used only for lingen. */
struct polymat {
    arith_hard * ab = nullptr;
    unsigned int m = 0;
    unsigned int n = 0;

  private:
    size_t size = 0;
    size_t alloc = 0;

  public:
    size_t capacity() const { return alloc; }
    size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }
    typedef arith_hard::elt elt;
    typedef elt * ptr;
    typedef elt const * srcptr;
    ptr x = nullptr;

    polymat()
        : polymat(nullptr, 0, 0, 0)
    {
    }
    polymat(polymat const &) = delete;
    polymat & operator=(polymat const &) = delete;
    polymat(polymat &&) noexcept;
    polymat & operator=(polymat &&) noexcept;
    ~polymat();
    polymat(arith_hard * ab, unsigned int m, unsigned int n, size_t len);
    int check_pre_init() const;
    void realloc(size_t newalloc);
    void zero();
#if 0
    void swap(polymat & b);
#endif
    void fill_random(unsigned int nsize, cxx_gmp_randstate & rstate);
    int cmp(polymat const & b) const;

    /* {{{ access interface for polymat */
    ptr part(unsigned int i, unsigned int j, size_t k)
    {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        ASSERT_ALWAYS(size);
        return ab->vec_subvec(x, (k * m + i) * n + j);
    }
    elt & coeff(unsigned int i, unsigned int j, size_t k)
    {
        return ab->vec_item(part(i, j, k), 0);
    }
    srcptr part(unsigned int i, unsigned int j, size_t k) const
    {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        ASSERT_ALWAYS(size);
        return ab->vec_subvec(x, (k * m + i) * n + j);
    }
    elt const & coeff(unsigned int i, unsigned int j, size_t k) const
    {
        return ab->vec_item(part(i, j, k), 0);
    }
    /* }}} */

    void addmat(size_t kc, polymat const & a, size_t ka,
                polymat const & b, size_t kb);
    void submat(size_t kc, polymat const & a, size_t ka,
                polymat const & b, size_t kb);
    void mulmat(size_t kc, polymat const & a, size_t ka,
                polymat const & b, size_t kb);
    void addmulmat(size_t kc, polymat const & a, size_t ka,
                   polymat const & b, size_t kb);
    void multiply_column_by_x(unsigned int j, size_t size);
    void truncate(polymat const & src, size_t nsize);
    void extract_column(unsigned int jdst, size_t kdst,
                        polymat const & src, unsigned int jsrc,
                        size_t ksrc);
    void extract_row_fragment(unsigned int i1, unsigned int j1,
                              polymat const & src, unsigned int i0,
                              unsigned int j0, size_t n);
    void mul(polymat const & a, polymat const & b);
    void addmul(polymat const & a, polymat const & b);
    void mp(polymat const & a, polymat const & b);
    void addmp(polymat const & a, polymat const & b);

    void set_matpoly(matpoly<false> const &);
    void rshift(polymat const & src, size_t k);
};

/* {{{ cut-off structures */
/* This structure is used to decide which algorithm to use for a given
 * input length. This essentially is a function from N*N to a finite set of
 * choices (so far, {0,1} only). The value returned for a balanced input
 * length x*x is:
 * if x >= cut: 1
 * if x < cut:
 *      if table == NULL:
 *              return 0
 *      find last (s,a) in table such that s <= x
 *      return a
 * For unbalanced input length x*y, we compare MIN(x,y) with subdivide.
 * At the moment there is no additional table designed to improve the
 * choice.
 *
 * Cutoffs are used only by lingen-polymat.c and lingen-bigpolymat.c -- however
 * the tuning program also uses it, so we expose it here too.
 */
/* Such a cutoff structure is used for finding which algorithm to used
 * for:
 *  - a balanced x*x product. -> basecase(0) or karatsuba(1)
 *  - an unbalanced x*y product. -> basecase(0) or split into balanced
 *  sizes(1)
 */

struct polymat_cutoff_info {
    size_t cut = SIZE_MAX;
    size_t subdivide = 0;
    std::vector<std::pair<size_t, int>> table;
    polymat_cutoff_info() = default;
    polymat_cutoff_info(size_t cut, size_t subdivide)
        : cut(cut)
        , subdivide(subdivide)
    {}
};
void polymat_cutoff_add_step(struct polymat_cutoff_info * c, size_t size,
                             int alg);
void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff,
                                 struct polymat_cutoff_info * old_cutoff);
void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff,
                                struct polymat_cutoff_info * old_cutoff);
/* }}} */

#endif /* LINGEN_POLYMAT_HPP_ */
