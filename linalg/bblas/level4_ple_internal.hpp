#ifndef LEVEL4_PLE_INTERNAL_HPP_
#define LEVEL4_PLE_INTERNAL_HPP_

#include "level4.hpp"
#include <vector>

struct PLE {/*{{{*/
    mat64_ptr X;
    unsigned int m;
    unsigned int n;

    int find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const;

    void propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const;

    void propagate_permutations(unsigned int ii0, unsigned int bj0, unsigned int const * q0, unsigned int const * q1) const;

    void move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q) const;

    void trsm(unsigned int bi,
            unsigned int bj,
            unsigned int yi0,
            unsigned int yi1) const;

    void sub(unsigned int bi,
            unsigned int bj,
            unsigned int yi0,
            unsigned int yi1,
            unsigned int ii) const;

    PLE(mat64_ptr X, unsigned int m, unsigned int n) : X(X), m(m), n(n) {}

    int operator()(unsigned int * p0);
};/*}}}*/

#endif	/* LEVEL4_PLE_INTERNAL_HPP_ */
