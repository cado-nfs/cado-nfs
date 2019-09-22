#ifndef LINGEN_MATPOLY_HPP_
#define LINGEN_MATPOLY_HPP_

#include <mutex>
#include <list>
#include <tuple>
#include <functional>   /* reference_wrapper */

#include "mpfq_layer.h"

class matpoly;
struct polymat;

#include "lingen_polymat.hpp"
#include "lingen_submatrix.hpp"
#include "lingen_memory_pool.hpp"

/* This is used only for lingen. */

/* We use abvec because this offers the possibility of having flat data
 *
 * Note that this ends up being exactly the same data type as polymat.
 * The difference here is that the stride is not the same.
 */

class matpoly {
    /* It's only exposed when we compile the mpi-enabled code, of course.
     * But on the other hand it's harmless to keep the friend declaration
     * in all cases.
     */
    friend class bigmatpoly;

    typedef abdst_vec ptr;
    typedef memory_pool_wrapper<ptr, true> memory_pool_type;
    static memory_pool_type memory;
public:
    struct memory_guard : private memory_pool_type::guard_base {
        memory_guard(size_t s) : memory_pool_type::guard_base(memory, s) {}
        ~memory_guard() { memory_pool_type::guard_base::pre_dtor(memory); }
    };

    static constexpr bool over_gf2 = false;
    /* if we ever want the check binary to make sure that the
     * specification works correctly also wrt. pre-init state. Not sure
     * it's terribly useful.
     *
    struct must_be_pre_init : public std::runtime_error {
        must_be_pre_init() : std::runtime_error("this data should be in pre-init state") {}
    };
    void make_sure_pre_init() const { if (!check_pre_init()) throw must_be_pre_init(); }
    void make_sure_not_pre_init() const { if (!check_pre_init()) throw must_be_pre_init(); }
    */
    // static void add_to_main_memory_pool(size_t s);
    abdst_field ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
private:
    size_t size = 0;
    size_t alloc = 0;
    abvec x = NULL;
public:
    inline size_t capacity() const { return alloc; }
    inline size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }

    inline unsigned int nrows() const { return m; }
    inline unsigned int ncols() const { return n; }
    size_t alloc_size() const;
    const void * data_area() const { return x; }
    size_t data_size() const { return m * n * abvec_elt_stride(ab, size); }
    bool is_tight() const { return alloc == size; }

    matpoly() { m=n=0; size=alloc=0; ab=NULL; x=NULL; }
    matpoly(abdst_field ab, unsigned int m, unsigned int n, int len);
    matpoly(matpoly const&) = delete;
    matpoly& operator=(matpoly const&) = delete;
    matpoly& set(matpoly const&);
    matpoly(matpoly &&);
    matpoly& operator=(matpoly &&);
    ~matpoly();
    bool check_pre_init() const ATTRIBUTE_WARN_UNUSED_RESULT { return x == NULL; }
    void realloc(size_t);
    inline void shrink_to_fit() { realloc(size); }
    void zero();
    void clear() { *this = matpoly(); }

    /* {{{ access interface for matpoly */
    inline abdst_vec part(unsigned int i, unsigned int j, unsigned int k=0) {
        return abvec_subvec(ab, x, (i*n+j)*alloc+k);
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k=0) {
        return abvec_coeff_ptr(ab, part(i,j,k), 0);
    }
    struct coeff_accessor_proxy {
        abdst_field ab;
        abdst_elt p;
        coeff_accessor_proxy(matpoly& F, unsigned int i,
                unsigned int j, unsigned int k)
            : ab(F.ab), p(F.coeff(i, j, k))
        {
        }
        coeff_accessor_proxy& operator+=(absrc_elt x) {
            abadd(ab, p, p, x);
            return *this;
        }
    };
    inline coeff_accessor_proxy coeff_accessor(unsigned int i, unsigned int j, unsigned int k = 0) {
        return coeff_accessor_proxy(*this, i, j, k);
    }
    inline absrc_vec part(unsigned int i, unsigned int j, unsigned int k=0) const {
        return abvec_subvec_const(ab, x, (i*n+j)*alloc+k);
    }
    inline absrc_elt coeff(unsigned int i, unsigned int j, unsigned int k=0) const {
        return abvec_coeff_ptr_const(ab, part(i,j,k), 0);
    }
    /* }}} */
    void set_constant_ui(unsigned long e);
    void set_constant(absrc_elt e);
    void fill_random(unsigned int size, gmp_randstate_t rstate);
    int cmp(matpoly const & b) const;
    void multiply_column_by_x(unsigned int j, unsigned int size);
    void divide_column_by_x(unsigned int j, unsigned int size);
    void truncate(matpoly const & src, unsigned int size);
    void truncate(unsigned int size) { truncate(*this, size); }
    int tail_is_zero(unsigned int size);
    void zero_pad(unsigned int nsize); /* changes size to nsize */
    void extract_column(
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc);
    void zero_column(unsigned int jdst, unsigned int kdst);
#if 0
    /* These two are implemented, but unused and untested anyway */
    void transpose_dumb(matpoly const & src);
    void extract_row_fragment(unsigned int i1, unsigned int j1,
        matpoly const & src, unsigned int i0, unsigned int j0,
        unsigned int n);
#endif
    void rshift(matpoly const &, unsigned int k);

    /* It is probably wise to avoid the mul and mp functions below. The
     * first-class citizens are the caching alternatives.
     */
    void add(matpoly const & a, matpoly const & b);
    void add(matpoly const & a) { add(*this, a); }
    void sub(matpoly const & a, matpoly const & b);
    void sub(matpoly const & a) { sub(*this, a); }
    void addmul(matpoly const & a, matpoly const & b);
    static matpoly mul(matpoly const & a, matpoly const & b);
    void addmp(matpoly const & a, matpoly const & c);
    static matpoly mp(matpoly const & a, matpoly const & c);

    void set_polymat(polymat const & src);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero(unsigned int k);
    struct view_t;
    struct const_view_t;

    struct view_t : public submatrix_range {
        matpoly & M;
        view_t(matpoly & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly & M) : submatrix_range(M), M(M) {}
        inline abdst_vec part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        inline absrc_vec part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };

    struct const_view_t : public submatrix_range {
        matpoly const & M;
        const_view_t(matpoly const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline absrc_vec part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }
    matpoly truncate_and_rshift(unsigned int truncated_size, unsigned int rshift);
};

#endif	/* LINGEN_MATPOLY_HPP_ */
