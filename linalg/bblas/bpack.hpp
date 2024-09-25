#ifndef BPACK_HPP_
#define BPACK_HPP_

#include <cstdint>                       // for uint64_t, uint8_t
#include <vector>
#include <type_traits>
#include <algorithm>
#include <gmp.h>                          // for gmp_randstate_t
#include "bblas_bitmat.hpp"  // for bitmat
// wee need these to resolve the member constexpr
#include "bblas_mat64.hpp" // IWYU pragma: keep
#include "bblas_mat8.hpp" // IWYU pragma: keep
#include "macros.h"                       // for ASSERT_ALWAYS

/* a bpack_view is *non-owning* view on a bit bitmat<T>, made of bit
 * matrices stored in row major order, the base type being the template
 * parameter "bitmat<T>".
 */

template<typename T> struct bpack; // IWYU pragma: keep
template<typename T> struct bpack_const_view; // IWYU pragma: keep
template<typename T> struct bpack_view; // IWYU pragma: keep

template<typename T>
struct bpack_ops {
    static void mul(bpack_view<T> C, bpack_const_view<T> A, bpack_const_view<T> B);
    static void mul(bpack<T> & C, bpack<T> const & A, bpack<T> const & B) {
        mul(C.view(), A.view(), B.view());
    }
    /* Do X <- A * X
     * This works in place on the bitmat<T> X. A is considered "implicitly
     * lower triangular". */
    static void mul_lt_ge(bpack_const_view<T> A, bpack_view<T> X);
    static void mul_lt_ge(bpack<T> const & A, bpack<T> & X) {
        mul_lt_ge(A.view(), X.view());
    }

    /*
    static void add(bpack<T> & C, bpack<T> const & A, bpack<T> const & B);
    static void transpose(bpack<T> & C, bpack<T> const & A);
    static void mul_lt_ge(bpack<T> & C, bpack<T> const & A, bpack<T> const & B) {
        mul(C, A, B);
    }
    static void addmul(bpack<T> & C, bpack<T> const & A, bpack<T> const & B);
    static void addmul(bpack<T> & C,
            bpack<T> const & A,
            bpack<T> const & B,
            unsigned int i0,
            unsigned int i1,
            unsigned int yi0,
            unsigned int yi1);
    static void trsm(bpack<T> const & L,
            bpack<T> & U,
            unsigned int yi0,
            unsigned int yi1);
    static void trsm(bpack<T> const & L, bpack<T> & U);
    */
    /* Keeps only the upper triangular part in U, and copy the lower
     * triangular, together with a unit block, to L */
    static void extract_LU(bpack_view<T> L, bpack_view<T> U);
    static void extract_uppertriangular(bpack_view<T> a, bpack_const_view<T> const b);
    static void extract_lowertriangular(bpack_view<T> a, bpack_const_view<T> const b);

    static void extract_LU(bpack<T> & L, bpack<T> & U)
    {
        extract_LU(L.view(), U.view());
    }
    static void extract_uppertriangular(bpack<T> & a, bpack<T> const & b)
    {
        extract_uppertriangular(a.view(), b.view());
    }
    static void extract_lowertriangular(bpack<T> & a, bpack<T> const & b)
    {
        extract_lowertriangular(a.view(), b.view());
    }
};

template<typename matrix_pointer> struct bpack_view_base {
    typedef typename std::remove_pointer<matrix_pointer>::type matrix;
    static constexpr const unsigned int B = matrix::width;
    typedef typename matrix::datatype U;
    protected:
    matrix_pointer X;
    public:
    matrix & cell(unsigned int bi, unsigned int bj) {
        return X[bi * nblocks + bj];
    }
    matrix const & cell(unsigned int bi, unsigned int bj) const {
        return X[bi * nblocks + bj];
    }
    unsigned int mblocks;
    unsigned int nblocks;
    inline unsigned int nrows() const { return mblocks * B; }
    inline unsigned int ncols() const { return nblocks * B; }
    inline unsigned int nrowblocks() const { return mblocks; }
    inline unsigned int ncolblocks() const { return nblocks; }
    bpack_view_base(matrix_pointer X, unsigned int mblocks, unsigned int nblocks) : X(X), mblocks(mblocks), nblocks(nblocks) { }
};

template<typename T>
struct bpack_const_view : public bpack_view_base<bitmat<T> const *>
{
    typedef bpack_view_base<bitmat<T> const *> super;
    using super::B;
    using super::mblocks;
    using super::nblocks;
    private:
    using super::X;
    public:
    using super::cell;
    typedef bpack_const_view<T> const_view_t;
    typedef bpack_view<T> view_t;
    bpack_const_view(bitmat<T> const * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    bool is_lowertriangular() const;
    bool is_uppertriangular() const;
    bool triangular_is_unit() const;
    bool operator==(int a) const;
    inline bool overlaps(bpack_const_view<T> v) {
        if (v.X <= X && (v.X + v.mblocks + v.nblocks) > X) return true;
        if (X <= v.X && (X + mblocks + nblocks) > v.X) return true;
        return false;
    }
    bool operator==(const_view_t) const;
    inline bool operator==(view_t v) const { return *this == v.const_view(); }
};

template<typename T>
struct bpack_view : bpack_view_base<bitmat<T> *> {
    typedef bpack_view_base<bitmat<T> *> super;
    using super::B;
    using super::mblocks;
    using super::nblocks;
    private:
    using super::X;
    public:
    using super::cell;
    bpack_view(bitmat<T> * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    typedef bpack_const_view<T> const_view_t;
    typedef bpack_view<T> view_t;
    const_view_t const_view() const { return const_view_t(super::X, mblocks, nblocks); }
    const_view_t view() const { return const_view_t(super::X, mblocks, nblocks); }
    view_t view() { return *this; }
    /* Goal: obtain a PLE decomposition of the bitmat<T> X, together with a list
     * of the pivot rows.
     *
     * we assume that X has size 64*m * 64*n, and that blocks are stored
     * row-major (i.e. we have m lists of n consecutive mat64's).
     *
     * The L part is stored inside the bitmat<T> X.
     *
     * The permutations are stored implicitly. We know that permutations are
     * formed as (current index i, other index >= i). Hence it is sufficient
     * to store the list of other indices, up to the rank. This information
     * is sufficient to recover the list of pivot rows.
     */
    std::vector<unsigned int> ple();
    std::vector<unsigned int> ple(std::vector<unsigned int> const &);
    void propagate_row_permutations(std::vector<unsigned int> const &);
    void invert_lower_triangular();

    inline bool is_lowertriangular() const { return const_view().is_lowertriangular(); }
    inline bool is_uppertriangular() const { return const_view().is_uppertriangular(); }
    inline bool triangular_is_unit() const { return const_view().triangular_is_unit(); }

    void fill_random(gmp_randstate_t rstate);
    void make_uppertriangular();
    void make_lowertriangular();
    void make_unit_uppertriangular();
    void make_unit_lowertriangular();
    void triangular_make_unit();
    bpack_view<T>& set(int a);
    bpack_view<T>& set(bpack_const_view<T> v);
    inline bpack_view<T>& set(bpack_view<T> v) { return set(v.const_view()); }
    inline bpack_view<T>& set(bpack<T> const & v) { return set(v.view()); }
    inline bool operator==(int a) const { return const_view() == a; }
    inline bool overlaps(bpack_view<T> v) { return const_view().overlaps(v.const_view()); }
    inline bool overlaps(bpack_const_view<T> v) { return const_view().overlaps(v); }
    inline bool operator==(const_view_t v) const { return const_view() == v; }
    inline bool operator==(view_t v) const { return const_view() == v.const_view(); }
};

template<typename T>
struct bpack : public bpack_ops<T> {
    typedef bpack_ops<T> ops;
    typedef bitmat<T> matrix;
    static constexpr const unsigned int B = matrix::width;
    typename matrix::vector_type X;
    unsigned int mblocks;
    unsigned int nblocks;
    inline unsigned int nrows() const { return mblocks * B; }
    inline unsigned int ncols() const { return nblocks * B; }
    inline unsigned int nrowblocks() const { return mblocks; }
    inline unsigned int ncolblocks() const { return nblocks; }
    bpack(unsigned int m, unsigned int n) : X((m/B)*(n/B)), mblocks(m/B), nblocks(n/B) {
        ASSERT_ALWAYS(m % B == 0);
        ASSERT_ALWAYS(n % B == 0);
        *this = 0;
    }
    typedef bpack_view<T> view_t;
    view_t view() { return view_t(X.data(), mblocks, nblocks); }
    typedef bpack_const_view<T> const_view_t;
    const_view_t view() const { return const_view_t(X.data(), mblocks, nblocks); }
    const_view_t const_view() const { return const_view_t(X.data(), mblocks, nblocks); }

    matrix & cell(unsigned int bi, unsigned int bj) { return view().cell(bi, bj); }
    matrix const & cell(unsigned int bi, unsigned int bj) const { return view().cell(bi, bj); }
    inline bpack<T>& operator=(int a) { view().set(a); return *this; }
    inline bool operator==(int a) const { return view() == a; }
    inline bpack(const_view_t a)
        : bpack(a.nrows(), a.ncols())
    {
        std::copy_n(&a.cell(0,0), a.nrowblocks() * a.ncolblocks(), &cell(0,0));
    }

    /* Most member functions are done at the view() or const_view() level
     */
    std::vector<unsigned int> ple() { return view().ple(); }
    std::vector<unsigned int> ple(std::vector<unsigned int> const & d) { return view().ple(d); }
    void propagate_row_permutations(std::vector<unsigned int> const & p){
        return view().propagate_row_permutations(p);
    }
    void invert_lower_triangular() { return view().invert_lower_triangular(); }

    inline bool is_lowertriangular() const { return view().is_lowertriangular(); }
    inline bool is_uppertriangular() const { return view().is_uppertriangular(); }
    inline bool triangular_is_unit() const { return view().triangular_is_unit(); }

    inline void fill_random(gmp_randstate_t rstate) { view().fill_random(rstate); }
    inline void make_uppertriangular() { view().make_uppertriangular(); }
    inline void make_lowertriangular() { view().make_lowertriangular(); }
    inline void make_unit_uppertriangular() { view().make_unit_uppertriangular(); }
    inline void make_unit_lowertriangular() { view().make_unit_lowertriangular(); }
    inline void triangular_make_unit() { view().triangular_make_unit(); }
    inline bool operator==(const_view_t v) const { return const_view() == v; }
    inline bool operator==(view_t v) const { return const_view() == v.const_view(); }
    inline bool operator==(bpack<T> const & v) const { return const_view() == v.const_view(); }
};

/* The code is in bpack.cpp ; presently there are no specializations, but
 * if the need arises, we may define a few of them.
 *
 * The ple() code uses another internal structure, of which ple() is only
 * a front-end.
 */
extern template struct bpack_ops<uint64_t>;
extern template struct bpack_ops<uint8_t>;
extern template struct bpack_const_view<uint64_t>;
extern template struct bpack_const_view<uint8_t>;
extern template struct bpack_view<uint64_t>;
extern template struct bpack_view<uint8_t>;
extern template struct bpack<uint64_t>;
extern template struct bpack<uint8_t>;

#endif	/* BPACK_HPP_ */
