#ifndef LEVEL4_PLE_INTERNAL_HPP_
#define LEVEL4_PLE_INTERNAL_HPP_

#include "bblas_level4.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include <algorithm>
#include <queue>
#include <vector>


/* Activate this to get extra timings for the lower-level steps of PLE.
 * The summary isn't very surprising.
 *  - for PLE on square matrices, sub dominates
 *  - for PLE on tall, skinny matrices, propagate_pivot dominates
 *  - for PLE on long flat matrices, trsm dominates.
 */
#define xxxTIME_PLE


template<typename matrix>
struct PLE : public bpack_view<matrix> {/*{{{*/
    using typename bpack_view<matrix>::U;
    using bpack_view<matrix>::B;
    using bpack_view<matrix>::mblocks;
    using bpack_view<matrix>::nblocks;
    // using bpack_view<matrix>::X;
    using bpack_view<matrix>::cell;
    using bpack_view<matrix>::view;
    using bpack_view<matrix>::const_view;
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
    std::vector<unsigned int> weights;
#if 0
    struct prio_cmp {
        typedef std::vector<unsigned int>::iterator T;
        /* priority queue selects the "largest" in the sense of the
         * comparison operator */
        bool operator()(T a, T b) const {
            if (*a != *b) return *a > *b;
            return a > b;
        }
    };
    /* We're managing the make_heap etc by ourselves */
    std::vector<std::vector<unsigned int>::iterator> prio;
#elif 0
    indexed_priority_queue<unsigned int, unsigned int, std::greater<unsigned int>> prio;
#else
    std::vector<unsigned int> prio_to_data;
    std::vector<unsigned int> data_to_prio;
#endif

    int find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);

    void propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);

    void propagate_row_permutations(unsigned int ii1, unsigned int bj0, std::vector<unsigned int>::const_iterator q0, std::vector<unsigned int>::const_iterator q1);

    void move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q);

    void trsm(unsigned int bi,
            unsigned int bj,
            unsigned int yi0,
            unsigned int yi1);

    void sub(unsigned int bi,
            unsigned int bj,
            unsigned int yi0,
            unsigned int yi1,
            unsigned int ii);

    struct debug_stuff : public bpack<matrix> {
        /* orig == target ? */
        using bpack<matrix>::mblocks;
        using bpack<matrix>::nblocks;
        using bpack<matrix>::nrows;
        using bpack<matrix>::ncols;
        using bpack<matrix>::cell;
        bpack<matrix> orig, target;
        debug_stuff(PLE const & ple)
            : bpack<matrix>(ple.nrows(), ple.ncols())
            , orig(ple.view())
            , target(orig)
        {}
        void apply_permutations(std::vector<unsigned int>::const_iterator, std::vector<unsigned int>::const_iterator);
        bpack<matrix> get_LL(unsigned int rr);
        bpack<matrix> get_UU(unsigned int rr);
        bool complete_check(bpack<matrix> const & LL, bpack<matrix> const & UU) ;
        bool check(bpack_const_view<matrix> X0, std::vector<unsigned int>::const_iterator p0, unsigned int ii) {
            (bpack<matrix>&) *this = X0;
            target = orig;
            apply_permutations(p0, p0 + ii);
            auto LL = get_LL(ii);
            auto UU = get_UU(ii);
            return complete_check(LL, UU);
        }
    };


    PLE(bpack_view<matrix> b, std::vector<unsigned int> d) : bpack_view<matrix>(b), weights(d) {
        ASSERT_ALWAYS(d.size() == b.nrows());
#if 0
        prio.reserve(d.size());
        for(auto it = weights.begin() ; it != weights.end() ; ++it)
            prio.push_back(it);
        std::sort(prio.begin(), prio.end(), prio_cmp());
#else
        std::is_sorted(weights.begin(), weights.end());
        for(unsigned int ii = 0 ; ii < d.size() ; ii++) {
            prio_to_data.push_back(ii);
            data_to_prio.push_back(ii);
        }
#endif
    }
    PLE(bpack_view<matrix> b) : PLE(b, std::vector<unsigned int>(b.nrows(), 0)) {}

    std::vector<unsigned int> operator()(debug_stuff * D = NULL);
};/*}}}*/

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
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

// extern template struct PLE<mat64>;
// extern template struct PLE<mat8>;

// must forward-declare the few specializations that we have.
template<> void PLE<mat64>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);
template<> void PLE<mat8>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);
template<> void PLE<mat64>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q);

extern template struct PLE<mat64>;
extern template struct PLE<mat8>;

#endif	/* LEVEL4_PLE_INTERNAL_HPP_ */
