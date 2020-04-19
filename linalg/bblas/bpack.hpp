#ifndef BPACK_HPP_
#define BPACK_HPP_

#include "bblas.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include <vector>

/* a bpack is *non-owning* view on a bit matrix, made of bit matrices
 * stored in row major order, the base type being the template parameter
 * "matrix".
 */
template<typename matrix>
struct bpack {
    static constexpr const unsigned int B = matrix::width;
    typedef typename matrix::datatype U;
    matrix * X;
    unsigned int m;
    unsigned int n;
    inline unsigned int nrows() const { return m * B; }
    inline unsigned int ncols() const { return n * B; }
    inline unsigned int nrowblocks() const { return m; }
    inline unsigned int ncolblocks() const { return n; }
    bpack(matrix * X, unsigned int m, unsigned int n) : X(X), m(m), n(n) { }

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
};

/* The code is in bpack.cpp ; presently there are no specializations, but
 * if the need arises, we may define a few of them.
 *
 * The ple() code uses another internal structure, of which ple() is only
 * a front-end.
 */
extern template struct bpack<mat64>;
extern template struct bpack<mat8>;

#endif	/* BPACK_HPP_ */
