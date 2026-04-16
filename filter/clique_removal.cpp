#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <functional>
#include <set>
#include <thread>
#include <utility>
#include <vector>

#include "fmt/base.h"

#include "macros.h"
#include "purge_matrix.hpp"
#include "timing.h"
#include "typedefs.h"

/* Contribution of each column of weight w to the weight of the connected
   component.
   Reference [1]: "The filtering step of discrete logarithm and integer
   factorization algorithms", Cyril Bouvier, https://hal.inria.fr/hal-00734654,
   2013, 27 pages.
   Each weight is identified by LAMBDA (0 to 6) and NU (0 to 3),
   the default one is \Omega_{31} (LAMBDA=3, NU=1).
   Cavallar's weight function is \Omega_{23} (LAMBDA=2, NU=3).
*/

template<int lambda> float w_function_helper(float w);
template<> float w_function_helper<0>(float) { return 1; }
template<> float w_function_helper<1>(float w) { return powf(2.0 / 3.0, w-2); }
template<> float w_function_helper<2>(float w) { return powf(0.5, w-2); }
template<> float w_function_helper<3>(float w) { return powf(0.8, w-2); }
template<> float w_function_helper<4>(float w) { return 1.0F / log2f (w); }
template<> float w_function_helper<5>(float w) { return 2.0F / w; }
template<> float w_function_helper<6>(float w) { return 4.0F / (w * w); }
template<int lambda, int nu>
struct w_function {
    float operator()(weight_t w) const {
        if (w >= 3)
            return w_function_helper<lambda>(w);
        else if (w == 2) {
            static_assert(nu >= 0 && nu <= 3);
            if constexpr (nu == 0) {
                return 0;
            } else {
                return 1.0 / (1 << (4 - nu));
            }
        } else {
            return 0;
        }
    }
};

#ifndef USE_WEIGHT_LAMBDA
#define USE_WEIGHT_LAMBDA 3
#endif
#ifndef USE_WEIGHT_NU
#define USE_WEIGHT_NU 1
#endif

/* print info on the weight function that is used */
void purge_matrix::print_clique_removal_weight_function()
{
    const w_function<USE_WEIGHT_LAMBDA,USE_WEIGHT_NU> f;
    fmt::print("Weight function used during clique removal:\n"
            "  0     1     2     3     4     5     6     7\n"
            "0.000 0.000 ");
    for (weight_t k = 2; k < 8; k++)
        fmt::print("{:0.3f} ", f(k));
    fmt::print("\n");
}

/*********** Functions to compute and delete connected component **************/

/* Compute connected component beginning at row i0.
 *
 * Return the number of rows in the component and a pair (i0, w)
 * where w is the weight computed by the function F on all column
 * weights.
 *
 * Multithread safe (if each thread has its own row_buffer).
 *
 * If a row is found in the component with index less than i0, return
 * a 0-row component. It indicates that the same component was found
 * earlier.
 *
 * XXX
 * there's a design decision regarding the underlying type for the
 * row_buffer. Namely:
 *  - an std::set allows clean code and avoids quadratic behavior in
 *  the connected component size. Note that it isn't clear we're
 *  going to have large connected components anyway, so maybe it's
 *  not that much of a benefit.
 *  - an std::vector has reserved storage and avoid reallocations (at
 *  least when amortized), which is not the case with std::set (it
 *  would be cumbersome to do this and keep iterator validity
 *  guarantees).
 *
 * We have both, at least for now.
 */

template <weight_function F>
auto purge_matrix::compute_connected_component_with_set(size_t i0,
        std::vector<size_t> const & sum2,
        F const & f) const
    -> std::pair<size_t, connected_component>
{
    std::set<size_t> buf;
    buf.insert(i0);

    float w = 0.; /* initial weight */

    /* Loop on all connected rows */
    for (auto const i : buf) {
        auto const * p = rows[i];

        /* Loop on all columns of the current row */
        for (auto const * q = p; *q != END_OF_ROW; q++) {
            index_t const h = *q;
            weight_t const wh = column_weights[h];
            w += f(wh);
            /* When we just _compute_ the connected components, we don't
             * have to check that sum2[h] != 0, since it's pretty much
             * guaranteed. Such isn't the case when we delete the
             * connected components, since wh==2 can appear for columns
             * that used to be heavier before. We don't have sum2 data
             * for them.
             */
            if (UNLIKELY(wh == 2 && sum2[h])) {
                /* TODO: this assert is actually not useful. By
                 * construction we always have sums of two things
                 * in this table. We could even make sum2 be a
                 * table of XORs, if we want.
                 */
                ASSERT_ALWAYS(i <= sum2[h]);
                size_t const i1 = sum2[h] - i;
                /* See if the connected component was found earlier.  */
                if (i1 < i0)
                    return {0, {}};

                if (std::ranges::find(buf, i1) == buf.end())
                    buf.insert(i1);
            }
        }
    }

    return {buf.size(), {.i=i0, .w=w}};
}

template <weight_function F>
auto purge_matrix::compute_connected_component(size_t i0,
        std::vector<size_t> const & sum2,
        std::vector<size_t> & buf, F const & f) const
    -> std::pair<size_t, connected_component>
{
    buf.clear();
    buf.push_back(i0);

    float w = 0.; /* initial weight */

    /* Loop on all connected rows */
    for (size_t k = 0; k < buf.size(); k++) {
        auto const i = buf[k];
        auto const * p = rows[i];

        /* Loop on all columns of the current row */
        for (auto const * q = p; *q != END_OF_ROW; q++) {
            index_t const h = *q;
            weight_t const wh = column_weights[h];
            w += f(wh);
            /* When we just _compute_ the connected components, we don't
             * have to check that sum2[h] != 0, since it's pretty much
             * guaranteed. Such isn't the case when we delete the
             * connected components, since wh==2 can appear for columns
             * that used to be heavier before. We don't have sum2 data
             * for them.
             */
            if (UNLIKELY(wh == 2 && sum2[h])) {
                /* TODO: this assert is actually not useful. By
                 * construction we always have sums of two things
                 * in this table. We could even make sum2 be a
                 * table of XORs, if we want.
                 */
                ASSERT_ALWAYS(i <= sum2[h]);
                size_t const i1 = sum2[h] - i;
                /* See if the connected component was found earlier.  */
                if (i1 < i0)
                    return {0, {}};

                /* As of c++17, std::set<T>::insert does *NOT*
                 * invalidate iterators to the set. So we could do
                 * the thing below with a set in log time. Alas, this has
                 * its downsides of frequents allocations and
                 * deallocations.
                 *
                 * The vector version is linear here, hence
                 * quadratic overall. For this reason we use fall back to
                 * the set version, which has better complexity, if we
                 * ever detect that some long components exist.
                 */
                if (buf.size() >= ccc_quadratic_abort)
                    return compute_connected_component_with_set(i0, sum2, f);

                if (std::ranges::find(buf, i1) == buf.end())
                    buf.push_back(i1);
            }
        }
    }

    return {buf.size(), {.i=i0, .w=w}};
}


/* Delete the connected component containing the row "current_row" */
void
purge_matrix::delete_connected_component(size_t i0,
                                std::vector<size_t> const & sum2,
                                std::vector<size_t> & buf)
{
    /* reuse the code of compute_connected_component, but with a
     * trivial weight function.  */
    const auto f = [](auto) { return float(); };
    compute_connected_component(i0, sum2, buf, f);

    /* This updates the different data fields atomically *BUT* it doesn't
     * update sum2! New columns with weight 2 can appear, and we won't
     * have their sum2 data right away. */
    for (auto const i : buf)
        delete_row(i);
}


/************* Code for clique (= connected component) removal ***************/

/* Print stats on the length of the cliques (connected components). The stats
 * are expensive to compute, so these functions should not be called by default.
 */
void purge_matrix::print_stats_on_cliques (FILE *out,
        std::vector<size_t> const & sum2,
        int verbose) const
{
    std::vector<size_t> component_sizes;
    component_sizes.reserve(rows.size());

    /* use a single buffer, avoid constant reallocations */
    std::vector<size_t> buf;

    for (size_t i = 0; i < rows.size(); i++) {
        if (is_active(i)) {
            const auto f = [](auto) { return float(); };
            auto [ size, C ] = compute_connected_component(i, sum2, buf, f);
            component_sizes.push_back(size);
        }
    }

    print_stats_generic (out, component_sizes, "cliques", "length", verbose);
}


/* Code for clique_removal.
 *
 * A connected component "belongs" to an index range I if its
 * least-indexed row is within I.
 *
 * The single-thread version of this code considers only one index range
 * I = [0, rows.size()).
 *
 * An n-thread instance behaves differently. The range [0, rows.size())
 * is split in chunks of size C, and thread k goes
 * through the ranges [(k+ell*n)*C, (k+1+ell*n)*C), one by one. This is
 * done so that the different threads have roughly similar work load (as
 * we start the search with rows of incresasing index, the probability of
 * discovering a lower-indexed row in the component is of course larger).
 */
size_t purge_matrix::clique_removal_core (
        sum2_type const & sum2,
        int64_t target_excess, size_t max_nb_comp,
        size_t nthreads,
        int verbose)
{
    /* We want the max_nb_comp components with the largest computed
     * score. The first step towards this is a min heap.
     *
     * Each thread will compute its own min heap.
     */
    std::vector<std::vector<connected_component>> Q(nthreads);

    static constexpr size_t batch_size_for_connected_components = 1024;
    auto fragment = [&](std::vector<connected_component> & Q,
            size_t i0, size_t i1, size_t stride)
    {
        std::vector<size_t> buf;

        for( ; i0 < rows.size() ; i0 += stride, i1 += stride) {
            const w_function<USE_WEIGHT_LAMBDA,USE_WEIGHT_NU> f;
            const std::greater<connected_component> G;
            i1 = std::min(i1, rows.size());
            for (size_t i = i0 ; i < i1 ; i++) {
                if (is_active(i)) {
                    auto [ size, C ] = compute_connected_component(i, sum2, buf, f);
                    if (!size)
                        continue;

                    if (Q.size() < max_nb_comp) {
                        /* insert, and make sure we still have a min-heap */
                        Q.push_back(C);
                        std::ranges::push_heap(Q, G);
                    } else {
                        /* If this connected component is even lighter
                         * than our min, there's no point in adding it.
                         */
                        if (C.w < Q.front().w)
                            continue;

                        /* the lightest connected component can be
                         * dropped, leaving only the max_nb_comp-1
                         * heaviest ones. Then we add our new component.
                         */
                        std::ranges::pop_heap(Q, G);
                        Q.back() = C;
                        std::ranges::push_heap(Q, G);
                    }
                }
            }
        }
    };

    if (nthreads == 1) {
        fragment(Q[0], 0, rows.size(), rows.size());
    } else {
        std::vector<std::thread> T;
        T.reserve(nthreads);
        for (size_t k = 0; k < nthreads; k++) {
            const size_t i0 = batch_size_for_connected_components * k;
            const size_t i1 = batch_size_for_connected_components * (k + 1);
            const size_t stride = batch_size_for_connected_components * nthreads;
            T.emplace_back(fragment, std::ref(Q[k]), i0, i1, stride);
        }
        for (auto & t: T)
            t.join();
    }

    /* At this point, Q[k] contains the max_nb_comp heaviest
     * connected component computed by thread k.
     *
     * We'll then proceed through them from heaviest to lightest, but we
     * won't do any further inserts.
     *
     * Two possible approaches here. It is fine if each thread does its
     * preferred removals based on its local order. But we prefer to do
     * it globally, which has the advantage of yielding more
     * deterministic results.
     *
     * To do so, the easiest approach is to sort the concatenation of all
     * Q and proceed from the top. Normal comparison works, here, since
     * the weights are stored with the components.
     */
    std::vector<connected_component> all_Q;
    for(auto const & q : Q)
        all_Q.insert(all_Q.end(), q.begin(), q.end());
    std::ranges::sort(all_Q);

    fmt::print("Cliq. rem.: computed heaviest connected components at "
            "{:2.2f}\n", seconds());
    fflush (stdout);

    if (verbose > 0)
        print_stats_on_cliques (stdout, sum2, verbose);

    size_t nb_clique_deleted = 0;
    float last_weight = 0.0;

    /* We're not doing the removals in a multithreaded way, on purpose
     * because we want determinism. But could it be a bit of an
     * opinionated, expensive, and ultimately unnecessary decision?
     */
    std::vector<size_t> buf;
    for( ; excess() > target_excess ; ) {
        if (all_Q.empty()) {
            fmt::print(stderr, "Cliq. rem.:"
                    " Warning, the list of connected components is empty\n");
            break;
        }

        auto [ i, w ] = all_Q.back();
        last_weight = w;
        delete_connected_component (i, sum2, buf);
        all_Q.pop_back();
        nb_clique_deleted++;
    }

    if (nb_clique_deleted)
        fmt::print("Cliq. rem.: weight of last removed connected component: {}\n",
                last_weight);

    return nb_clique_deleted;
}

/***************************** Clique removal ********************************/

void purge_matrix::clique_removal (int64_t target_excess,
                 size_t nthreads, int verbose)
{
    if (excess() <= target_excess)
        return;

    /* max_nb_comp is the maximum number of connected components that
     * each threads is going to store.
     *
     * Each connected component reduces the excess by at most one, and
     * perhaps zero (albeit rarely).
     *
     * Each thread will compute max_nb_comp, for determinism (in case one
     * thread actually owns all the heavy components!). Therefore for
     * nthreads>=2, we'll compute at least 2*max_nb_comp components
     * overall.
     *
     * If nthreads==1, we increase max_nb_comp by 25% to take into account the
     * fact that some components (we expect not many) won't reduce the
     * excess.
     */
    size_t max_nb_comp = excess() - target_excess;
    if (nthreads == 1)
        max_nb_comp+= max_nb_comp/4;

    /* Then, call the core function that computes and deletes connected
     * components until excess is equal to target_excess.
     */
    auto sum2 = compute_sum2(nthreads);
    fmt::print("Cliq. rem.: computed sum2 at {:2.2f}\n", seconds());

    size_t nb_cliques_del = clique_removal_core (sum2, target_excess,
            max_nb_comp, nthreads, verbose);

    fmt::print("Cliq. rem.: deleted {} heaviest connected "
            "components at {:2.2f}\n", nb_cliques_del, seconds());
    if (verbose > 0)
        fmt::print("# INFO: max_nb_comp_per_thread={} target_excess="
                "{}\n", max_nb_comp, target_excess);
    fflush (stdout);
}
