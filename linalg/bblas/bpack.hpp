#ifndef BPACK_HPP_
#define BPACK_HPP_

#include "bblas.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include <vector>
#include <type_traits>
#include <algorithm>

/* a bpack_view is *non-owning* view on a bit matrix, made of bit
 * matrices stored in row major order, the base type being the template
 * parameter "matrix".
 */

template<typename matrix> struct bpack;
template<typename matrix> struct bpack_const_view;
template<typename matrix> struct bpack_view;

template<typename matrix>
struct bpack_ops {
    static void mul(bpack_view<matrix> C, bpack_const_view<matrix> A, bpack_const_view<matrix> B);
    static void mul(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B) {
        mul(C.view(), A.view(), B.view());
    }
    /* This works in place on the matrix X. A is considered "implicitly
     * lower triangular". */
    static void mul_lt_ge(bpack_const_view<matrix> A, bpack_view<matrix> X);
    static void mul_lt_ge(bpack<matrix> const & A, bpack<matrix> & X) {
        mul_lt_ge(A.view(), X.view());
    }

    /*
    static void add(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B);
    static void transpose(bpack<matrix> & C, bpack<matrix> const & A);
    static void mul_lt_ge(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B) {
        mul(C, A, B);
    }
    static void addmul(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B);
    static void addmul(bpack<matrix> & C,
            bpack<matrix> const & A,
            bpack<matrix> const & B,
            unsigned int i0,
            unsigned int i1,
            unsigned int yi0,
            unsigned int yi1);
    static void trsm(bpack<matrix> const & L,
            bpack<matrix> & U,
            unsigned int yi0,
            unsigned int yi1);
    static void trsm(bpack<matrix> const & L, bpack<matrix> & U);
    static void extract_uppertriangular(bpack<matrix> & a, bpack<matrix> const & b);
    static void extract_lowertriangular(bpack<matrix> & a, bpack<matrix> const & b);
    */
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

template<typename matrix>
struct bpack_const_view : public bpack_view_base<matrix const *>
{
    typedef bpack_view_base<matrix const *> super;
    using super::B;
    using super::mblocks;
    using super::nblocks;
    private:
    using super::X;
    public:
    using super::cell;
    bpack_const_view(matrix const * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    bool is_lowertriangular() const;
    bool is_uppertriangular() const;
    bool triangular_is_unit() const;
    bool operator==(int a) const;
    inline bool overlaps(bpack_const_view<matrix> v) {
        if (v.X <= X && (v.X + v.mblocks + v.nblocks) > X) return true;
        if (X <= v.X && (X + mblocks + nblocks) > v.X) return true;
        return false;
    }
};

template<typename matrix>
struct bpack_view : bpack_view_base<matrix *> {
    typedef bpack_view_base<matrix *> super;
    using super::B;
    using super::mblocks;
    using super::nblocks;
    private:
    using super::X;
    public:
    using super::cell;
    bpack_view(matrix * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    typedef bpack_const_view<matrix> const_view_t;
    typedef bpack_view<matrix> view_t;
    const_view_t const_view() const { return const_view_t(super::X, mblocks, nblocks); }
    const_view_t view() const { return const_view_t(super::X, mblocks, nblocks); }
    view_t view() { return *this; }
    /* Goal: obtain a PLE decomposition of the matrix X, together with a list
     * of the pivot rows.
     *
     * we assume that X has size 64*m * 64*n, and that blocks are stored
     * row-major (i.e. we have m lists of n consecutive mat64's).
     *
     * The L part is stored inside the matrix X.
     *
     * The permutations are stored implicitly. We know that permutations are
     * formed as (current index i, other index >= i). Hence it is sufficient
     * to store the list of other indices, up to the rank. This information
     * is sufficient to recover the list of pivot rows.
     */
    std::vector<unsigned int> ple();
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
    bpack_view<matrix>& set(int a);
    bpack_view<matrix>& set(bpack_const_view<matrix> v);
    inline bpack_view<matrix>& set(bpack_view<matrix> v) { return set(v.const_view()); }
    inline bpack_view<matrix>& set(bpack<matrix> const & v) { return set(v.view()); }
    inline bool operator==(int a) const { return const_view() == a; }
    inline bool overlaps(bpack_view<matrix> v) { return const_view().overlaps(v.const_view()); }
    inline bool overlaps(bpack_const_view<matrix> v) { return const_view().overlaps(v); }
};

template<typename matrix>
struct bpack : public bpack_ops<matrix> {
    typedef bpack_ops<matrix> ops;
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
    typedef bpack_view<matrix> view_t;
    view_t view() { return view_t(&X[0], mblocks, nblocks); }
    typedef bpack_const_view<matrix> const_view_t;
    const_view_t view() const { return const_view_t(&X[0], mblocks, nblocks); }

    matrix & cell(unsigned int bi, unsigned int bj) { return view().cell(bi, bj); }
    matrix const & cell(unsigned int bi, unsigned int bj) const { return view().cell(bi, bj); }
    inline bpack<matrix>& operator=(int a) { view().set(a); return *this; }
    inline bool operator==(int a) const { return view() == a; }
    inline bpack<matrix>(const_view_t a)
        : bpack(a.nrows(), a.ncols())
    {
        std::copy_n(&a.cell(0,0), a.nrowblocks() * a.ncolblocks(), &cell(0,0));
    }

    /* Most member functions are done at the view() or const_view() level
     */
    std::vector<unsigned int> ple() { return view().ple(); }
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
};

/* The code is in bpack.cpp ; presently there are no specializations, but
 * if the need arises, we may define a few of them.
 *
 * The ple() code uses another internal structure, of which ple() is only
 * a front-end.
 */
extern template struct bpack_ops<mat64>;
extern template struct bpack_ops<mat8>;
extern template struct bpack_const_view<mat64>;
extern template struct bpack_const_view<mat8>;
extern template struct bpack_view<mat64>;
extern template struct bpack_view<mat8>;
extern template struct bpack<mat64>;
extern template struct bpack<mat8>;

#endif	/* BPACK_HPP_ */
