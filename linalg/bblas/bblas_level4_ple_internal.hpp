#ifndef LEVEL4_PLE_INTERNAL_HPP_
#define LEVEL4_PLE_INTERNAL_HPP_

#include "bblas_level4.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include <vector>

/* Activate this to get extra timings for the lower-level steps of PLE.
 * The summary isn't very surprising.
 *  - for PLE on square matrices, sub dominates
 *  - for PLE on tall, skinny matrices, propagate_pivot dominates
 *  - for PLE on long flat matrices, trsm dominates.
 */
#define xxxTIME_PLE

template<typename matrix>
struct PLE {/*{{{*/
    static constexpr const unsigned int B = matrix::width;
    typedef typename matrix::datatype U;
    matrix * X;
    unsigned int m;
    unsigned int n;
#ifdef TIME_PLE
    static unsigned long ncalls;
    static double t_find_pivot;
    static double t_propagate_pivot;
    static double t_propagate_permutation;
    static double t_move_l_fragments;
    static double t_trsm;
    static double t_sub;
    static double t_total;
    static void print_and_flush_stats();
#endif

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


    PLE(matrix * X, unsigned int m, unsigned int n) : X(X), m(m), n(n) { }

    std::vector<unsigned int> operator()(debug_stuff * D = NULL);
};/*}}}*/

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
template<typename matrix> constexpr const unsigned int PLE<matrix>::B;
#ifdef TIME_PLE
template<typename matrix> unsigned long PLE<matrix>::ncalls = 0;
template<typename matrix> double PLE<matrix>::t_find_pivot = 0;
template<typename matrix> double PLE<matrix>::t_propagate_pivot = 0;
template<typename matrix> double PLE<matrix>::t_propagate_permutation = 0;
template<typename matrix> double PLE<matrix>::t_move_l_fragments = 0;
template<typename matrix> double PLE<matrix>::t_trsm = 0;
template<typename matrix> double PLE<matrix>::t_sub = 0;
template<typename matrix> double PLE<matrix>::t_total = 0;
#endif

#include "bblas_level4_ple_internal_inl.hpp"

// extern template struct PLE<mat64>;
// extern template struct PLE<mat8>;

#endif	/* LEVEL4_PLE_INTERNAL_HPP_ */
