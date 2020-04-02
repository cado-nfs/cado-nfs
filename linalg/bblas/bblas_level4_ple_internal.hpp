#ifndef LEVEL4_PLE_INTERNAL_HPP_
#define LEVEL4_PLE_INTERNAL_HPP_

#include "bblas_level4.hpp"
#include <vector>

template<typename matrix>
struct PLE {/*{{{*/
    static constexpr const unsigned int B = matrix::width;
    typedef typename matrix::datatype U;
    matrix * X;
    unsigned int m;
    unsigned int n;

    int find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const;

    void propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const;

    void propagate_permutations(unsigned int ii1, unsigned int bj0, std::vector<unsigned int>::const_iterator q0, std::vector<unsigned int>::const_iterator q1) const;

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

    struct debug_stuff {
    unsigned int m;
    unsigned int n;
        typename matrix::vector_type X_orig;
        typename matrix::vector_type X, X_target;
        debug_stuff(PLE const & ple)
            : m(ple.m)
            , n(ple.n)
            , X_orig(ple.X, ple.X + m * n)
        {}
        void start_check(matrix const * X0)
        {
            X = typename matrix::vector_type(X0, X0 + m * n);
            X_target = X_orig;
        }
        void start_check(typename matrix::vector_type const & X0)
        {
            ASSERT_ALWAYS(X0.size() == m * n);
            start_check(&X0[0]);
        }
        void apply_permutations(std::vector<unsigned int>::const_iterator, std::vector<unsigned int>::const_iterator);
        void apply_permutations(std::vector<unsigned int> const & V)
        {
            apply_permutations(V.begin(), V.end());
        }
        typename matrix::vector_type get_LL(unsigned int rr);
        typename matrix::vector_type get_UU(unsigned int rr);
        bool complete_check(typename matrix::vector_type const & LL, typename matrix::vector_type const & UU) ;
        bool check(matrix const * X0, std::vector<unsigned int>::const_iterator p0, unsigned int ii);
    };


    PLE(matrix * X, unsigned int m, unsigned int n) : X(X), m(m), n(n) {}

    std::vector<unsigned int> operator()(debug_stuff * D = NULL);
};/*}}}*/

#include "bblas_level4_ple_internal_inl.hpp"

#endif	/* LEVEL4_PLE_INTERNAL_HPP_ */
