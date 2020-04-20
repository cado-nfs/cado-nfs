#ifndef BPACK_HPP_
#define BPACK_HPP_

#include "bblas.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include <vector>
#include <type_traits>

/* a bpack_view is *non-owning* view on a bit matrix, made of bit
 * matrices stored in row major order, the base type being the template
 * parameter "matrix".
 */

template<typename matrix> struct bpack;
template<typename matrix> struct bpack_const_view;
template<typename matrix> struct bpack_view;

template<typename matrix>
struct bpack_ops {
    /*
    static void add(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B);
    static void transpose(bpack<matrix> & C, bpack<matrix> const & A);
    static void mul(bpack<matrix> & C, bpack<matrix> const & A, bpack<matrix> const & B);
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
    matrix_pointer X;
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
    using super::X;
    bpack_const_view(matrix const * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    bool is_lowertriangular() const;
    bool is_uppertriangular() const;
    bool triangular_is_unit() const;
};

template<typename matrix>
struct bpack_view : bpack_view_base<matrix *> {
    typedef bpack_view_base<matrix *> super;
    using super::B;
    using super::mblocks;
    using super::nblocks;
    using super::X;
    bpack_view(matrix * X, unsigned int mblocks, unsigned int nblocks) : super(X, mblocks, nblocks) {}
    typedef bpack_const_view<matrix> const_view_t;
    const_view_t const_view() const { return const_view_t(X, mblocks, nblocks); }
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
};

template<typename matrix>
struct bpack : public bpack_ops<matrix> {
    typedef bpack_ops<matrix> ops;
    static constexpr const unsigned int B = matrix::width;
    std::vector<matrix> X;
    unsigned int mblocks;
    unsigned int nblocks;
    inline unsigned int nrows() const { return mblocks * B; }
    inline unsigned int ncols() const { return nblocks * B; }
    inline unsigned int nrowblocks() const { return mblocks; }
    inline unsigned int ncolblocks() const { return nblocks; }
    bpack(unsigned int m, unsigned int n) : X((m/B)*(n/B)), mblocks(m/B), nblocks(n/B) {
        ASSERT_ALWAYS(m % B == 0);
        ASSERT_ALWAYS(n % B == 0);
    }
    void set_zero();
    typedef bpack_view<matrix> view_t;
    view_t view() { return view_t(&X[0], mblocks, nblocks); }
    typedef bpack_const_view<matrix> const_view_t;
    const_view_t view() const { return const_view_t(&X[0], mblocks, nblocks); }

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
