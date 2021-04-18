#ifndef LEVEL4_PLE_INTERNAL_HPP_
#define LEVEL4_PLE_INTERNAL_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "bblas_level4.hpp"
#include <cstdint>
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


template<typename T>
struct PLE : public bpack_view<T> {/*{{{*/
    using typename bpack_view<T>::U;
    using bpack_view<T>::B;
    using bpack_view<T>::mblocks;
    using bpack_view<T>::nblocks;
    // using bpack_view<T>::X;
    using bpack_view<T>::cell;
    using bpack_view<T>::view;
    using bpack_view<T>::const_view;
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
    std::vector<unsigned int> prio_to_data;
    std::vector<unsigned int> data_to_prio;

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

    struct debug_stuff : public bpack<T> {
        /* orig == target ? */
        using bpack<T>::mblocks;
        using bpack<T>::nblocks;
        using bpack<T>::nrows;
        using bpack<T>::ncols;
        using bpack<T>::cell;
        bpack<T> orig, target;
        debug_stuff(PLE const & ple)
            : bpack<T>(ple.nrows(), ple.ncols())
            , orig(ple.view())
            , target(orig)
        {}
        void apply_permutations(std::vector<unsigned int>::const_iterator, std::vector<unsigned int>::const_iterator);
        bpack<T> get_LL(unsigned int rr);
        bpack<T> get_UU(unsigned int rr);
        bool complete_check(bpack<T> const & LL, bpack<T> const & UU) ;
        bool check(bpack_const_view<T> X0, std::vector<unsigned int>::const_iterator p0, unsigned int ii) {
            (bpack<T>&) *this = X0;
            target = orig;
            apply_permutations(p0, p0 + ii);
            auto LL = get_LL(ii);
            auto UU = get_UU(ii);
            return complete_check(LL, UU);
        }
    };


    PLE(bpack_view<T> b, std::vector<unsigned int> d);
    PLE(bpack_view<T> b);

    std::vector<unsigned int> operator()(debug_stuff * D = NULL);
};/*}}}*/

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
#ifdef TIME_PLE
template<typename T> unsigned long PLE<T>::ncalls = 0;
template<typename T> double PLE<T>::t_find_pivot = 0;
template<typename T> double PLE<T>::t_propagate_pivot = 0;
template<typename T> double PLE<T>::t_propagate_permutation = 0;
template<typename T> double PLE<T>::t_move_l_fragments = 0;
template<typename T> double PLE<T>::t_trsm = 0;
template<typename T> double PLE<T>::t_sub = 0;
template<typename T> double PLE<T>::t_total = 0;
#endif

// extern template struct PLE<uint64_t>;
// extern template struct PLE<uint8_t>;

// must forward-declare the few specializations that we have.
#if defined(HAVE_AVX2) || defined(HAVE_SSE41)
template<> void PLE<uint64_t>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);
#endif
#ifdef HAVE_MMX
template<> void PLE<uint8_t>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j);
#endif
#if defined(HAVE_AVX2) || defined(HAVE_SSE41)
template<> void PLE<uint64_t>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q);
#endif

extern template struct PLE<uint64_t>;
extern template struct PLE<uint8_t>;

#endif	/* LEVEL4_PLE_INTERNAL_HPP_ */
