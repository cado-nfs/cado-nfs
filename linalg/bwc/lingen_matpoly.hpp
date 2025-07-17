#ifndef CADO_LINGEN_MATPOLY_HPP
#define CADO_LINGEN_MATPOLY_HPP

#include <cstddef>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_call_companion.hpp"
#include "lingen_memory_pool.hpp"
#include "macros.h"
#include "arith-hard.hpp"
#include "submatrix_range.hpp"
#include "tree_stats.hpp"
#include "runtime_numeric_cast.hpp"

struct polymat;

/* This is used only for lingen. */

/* We use abvec because this offers the possibility of having flat data
 *
 * Note that this ends up being exactly the same data type as polymat.
 * The difference here is that the stride is not the same.
 */

template<bool is_binary> class matpoly;
template<bool is_binary> class bigmatpoly;

template<>
class matpoly<false> {
    /* It's only exposed when we compile the mpi-enabled code, of course.
     * But on the other hand it's harmless to keep the friend declaration
     * in all cases.
     */
    friend class bigmatpoly<false>;

public:
    static constexpr bool is_binary = false;
    typedef ::arith_hard arith_hard;
    typedef arith_hard::elt elt;
    typedef elt * ptr;
    typedef elt const * srcptr;

private:
    typedef memory_pool_wrapper<ptr, true> memory_pool_type;
    static memory_pool_type memory;

public:
    /* The matpoly memory guard, as long as it remains in scope,
     * registers an allowed allocation size for the singleton object
     * matpoly::memory.
     *
     * Copying the guard is a no-op, since doing so does not add up to
     * the allowed area.
     *
     * On the other hand, moving a guard does transfer the allowed
     * allocation size to the new scope
     */
    struct memory_guard : private memory_pool_type::guard_base {
        explicit memory_guard(size_t s) : memory_pool_type::guard_base(memory, s) {}
        memory_guard(memory_guard const &) : memory_pool_type::guard_base(memory, 0) {}
        memory_guard& operator=(memory_guard const &) = delete;

        /* See memory_pool_type::guard_base for the reason why we can't
         * move the guard objects (at least yet)
         */
        memory_guard(memory_guard &&) = delete;
        memory_guard& operator=(memory_guard &&) = delete;
        // we're in a dtor, hence nothrow, yet we have
        // ASSERT_ALWAYS...
        // coverity[exn_spec_violation]
        ~memory_guard() { memory_pool_type::guard_base::pre_dtor(memory); }
    };

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
    arith_hard * ab = nullptr;
    unsigned int m = 0;
    unsigned int n = 0;
private:
    size_t size = 0;
    size_t alloc = 0;
    ptr x = nullptr;
public:
    size_t capacity() const { return alloc; }
    size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }
    size_t get_true_nonzero_size() const;

    unsigned int nrows() const { return m; }
    unsigned int ncols() const { return n; }
    const void * data_area() const { return x; }
    size_t data_entry_size_in_bytes() const { return ab->vec_elt_stride(size); }
    size_t data_size_in_bytes() const { return m * n * data_entry_size_in_bytes(); }
private:
    size_t data_entry_alloc_size_in_bytes(size_t a) const { return ab->vec_elt_stride(a); }
    size_t data_alloc_size_in_bytes(size_t a) const { return m * n * data_entry_alloc_size_in_bytes(a); }
public:
    size_t data_entry_alloc_size_in_bytes() const { return data_entry_alloc_size_in_bytes(alloc); }
    size_t data_alloc_size_in_bytes() const { return data_alloc_size_in_bytes(alloc); }
    bool is_tight() const { return alloc == size; }

    matpoly() = default;
    matpoly(arith_hard * ab, unsigned int m, unsigned int n, size_t len);
    matpoly(matpoly const&) = delete;
    matpoly& operator=(matpoly const&) = delete;
    matpoly& set(matpoly const&);
    matpoly(matpoly &&) noexcept;
    matpoly& operator=(matpoly &&) noexcept;
    ~matpoly();
    matpoly similar_shell() const { return { ab, m, n, 0 }; }
    bool check_pre_init() const ATTRIBUTE_WARN_UNUSED_RESULT {
        return x == nullptr;
    }
    void realloc(size_t);
    void shrink_to_fit() { realloc(size); }
    void zero();
    void clear() { *this = matpoly(); }

    /* {{{ access interface for matpoly */
    /* part_head does not _promise_ to point to coefficient k exactly. In
     * the simd case (lingen_matpoly_binary.hpp) it points to the word
     * where coefficient k can be found */
    ptr part(unsigned int i, unsigned int j) {
        return ab->vec_subvec(x, (i*n+j)*alloc);
    }
    ptr part_head(unsigned int i, unsigned int j, size_t k) {
        return ab->vec_subvec(part(i, j), k);
    }
    elt & coeff(unsigned int i, unsigned int j, size_t k=0) {
        return ab->vec_item(part(i, j), k);
    }
    struct coeff_accessor_proxy {
        arith_hard * ab;
        arith_hard::elt * p;
        coeff_accessor_proxy(matpoly& F, unsigned int i,
                unsigned int j, size_t k)
            : ab(F.ab), p(F.part_head(i, j, k))
        {
        }
        coeff_accessor_proxy& operator+=(elt const & x) {
            ab->add_and_reduce(*p, x);
            return *this;
        }
    };
    coeff_accessor_proxy coeff_accessor(unsigned int i, unsigned int j, size_t k = 0) {
        return { *this, i, j, k };
    }
    srcptr part(unsigned int i, unsigned int j) const {
        return ab->vec_subvec(x, (i*n+j)*alloc);
    }
    srcptr part_head(unsigned int i, unsigned int j, size_t k) const {
        return ab->vec_subvec(part(i, j), k);
    }
    elt const & coeff(unsigned int i, unsigned int j, size_t k=0) const {
        return ab->vec_item(part(i, j), k);
    }
    /* }}} */
    void set_constant_ui(unsigned long e);
    void set_constant(elt const & e);
    /* Note that this method does not change the size field */
    void fill_random(size_t k0, size_t k1, cxx_gmp_randstate & rstate);
    void clear_and_set_random(size_t len, cxx_gmp_randstate & rstate)
    {
        if (len > capacity())
            zero_pad(len);
        zero_pad(capacity());
        fill_random(0, len, rstate);
        set_size(len);
    }

    int cmp(matpoly const & b) const;
    void multiply_column_by_x(unsigned int j, size_t size);
    void divide_column_by_x(unsigned int j, size_t size);
    void truncate(matpoly const & src, size_t size);
    void truncate(size_t size) { truncate(*this, size); }
    int tail_is_zero(size_t size);
    /* not to be confused with the former. the following two are in fact
     * relevant only to the binary interface. They're just noops here.
     */
    static bool high_word_is_clear() { return true; }
    void clear_high_word() {}

    /* This changes size to nsize, and fills [size..nsize[ with zeroes */
    void zero_pad(size_t nsize); /* changes size to nsize */

    void zero_with_size(size_t size) { set_size(0); zero_pad(size); }

    void extract_column(
        unsigned int jdst, size_t kdst,
        matpoly const & src, unsigned int jsrc, size_t ksrc);
    void zero_column(unsigned int jdst, size_t kdst);
#if 0
    /* These two are implemented, but unused and untested anyway */
    void transpose_dumb(matpoly const & src);
    void extract_row_fragment(unsigned int i1, unsigned int j1,
        matpoly const & src, unsigned int i0, unsigned int j0,
        size_t n);
#endif
    void rshift(matpoly const &, size_t k);
    void rshift(size_t k);

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

    /* We have nothing terrific to identify as substeps here, and there's
     * no schedule information that we intend to exploit */
    static matpoly mp(tree_stats & stats, matpoly const & a, matpoly const & b, lingen_call_companion::mul_or_mp_times * M)
    {
        if (!M) return mp(a, b);
        const tree_stats::smallstep_sentinel dummy(stats, M->step_name());
        return mp(a, b);
    }

    static matpoly mul(tree_stats & stats, matpoly const & a, matpoly const & b, lingen_call_companion::mul_or_mp_times * M) {
        if (!M) return mul(a, b);
        const tree_stats::smallstep_sentinel dummy(stats, M->step_name());
        return mul(a, b);
    }
 
    void set_polymat(polymat const & src);
    int coeff_is_zero(size_t k) const;
    void coeff_set_zero(size_t k);

    struct view_t : public submatrix_range {
        matpoly & M;
        view_t(matpoly & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly & M) : submatrix_range(M), M(M) {}
        ptr part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
        void zero();
    };

    struct const_view_t : public submatrix_range {
        matpoly const & M;
        const_view_t(matpoly const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return { *this, S }; }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return { *this, S }; }
    view_t view() { return { *this }; }
    const_view_t view() const { return { *this }; }

    static void copy(view_t t, const_view_t a);
    static void addmul(view_t t, const_view_t t0, const_view_t t1);
    static void addmp(view_t t, const_view_t t0, const_view_t t1);

    matpoly truncate_and_rshift(size_t truncated_size, size_t rshift);
};

extern template class matpoly<false>;

#endif	/* LINGEN_MATPOLY_HPP_ */
