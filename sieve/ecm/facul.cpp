#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include <gmp.h>

#include "facul_method.hpp"
#include "facul_strategies.hpp"
#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "modset.hpp" // for FaculModulusBase
#include "macros.h"

/*
 * Apply a bunch of ECM curves to find a factor with high probability.
 *
 * There are cases where we apparently are interested in seeing the
 * cofactors, so we might as well try to return them.
 */
static facul_status facul_aux(
        std::vector<cxx_mpz> & factors,
        std::vector<std::unique_ptr<FaculModulusBase>> & todo,
        facul_strategies_base const & strategies,
        std::vector<facul_method_side> const & methods,
        size_t start_method,
        int side)
{
    for (size_t i = start_method ; i < methods.size() && !todo.empty() ; i++) {
        auto const & meth = methods[i];
        if (meth.side != side)
            continue; /* this method is not for this side */

#ifdef ENABLE_UNSAFE_FACUL_STATS
        if (i < STATS_LEN)
            stats_called_aux[i]++;
#endif /* ENABLE_UNSAFE_FACUL_STATS */

        std::vector<std::unique_ptr<FaculModulusBase>> next_todo;
        for(auto & n : todo) {
            std::vector<std::unique_ptr<FaculModulusBase>> more_composites;

            facul_status const res = n->facul_doit_onefm(
                    factors, *meth.method, more_composites,
                    strategies.lpb[side],
                    strategies.BB[side],
                    strategies.BBB[side]);

            if (res == FACUL_NOT_SMOOTH)
                /* The cofactor m is not smooth. Abort */
                return FACUL_NOT_SMOOTH;

            if (res == FACUL_MAYBE) {
                if (more_composites.empty()) {
                    /* Zero factors found. If it was the last method for this
                       side, then one stops the cofactorization. Otherwise, one
                       tries with an other method */
                    next_todo.emplace_back(std::move(n));
                } else {
                    for(auto & c : more_composites)
                        next_todo.emplace_back(std::move(c));
                }
            }
        }
        std::swap(todo, next_todo);
    }
    return todo.empty() ? FACUL_SMOOTH : FACUL_MAYBE;
}

/* This function tries to factor all the cofactors in the vector N with
 * strategies. The returned facul_result objects (one for each input number)
 * will have the prime factors.
 *
 * If we find composites, but not to the point of obtaining a complete
 * factorization, then they're not returned.
 */

std::vector<facul_result>
facul_all(std::vector<cxx_mpz> const & N, facul_strategies const & strategies)
{
    const size_t nsides = N.size();

    std::vector<unsigned int> sizes(N.size());
    std::transform(N.cbegin(), N.cend(), sizes.begin(),
                   [] (cxx_mpz const & n) { return n.bits(); });

    auto const & methods = strategies(sizes);

    std::vector<facul_result> res(nsides, FACUL_NOT_SMOOTH);

    std::vector<std::vector<std::unique_ptr<FaculModulusBase>>> composites(nsides);

#ifdef PARI
    fmt::print(stderr, join(N, " "));
#endif

    ASSERT_ALWAYS(nsides == strategies.B.size());
    ASSERT_ALWAYS(nsides == strategies.BB.size());

    for(size_t side = 0 ; side < nsides ; side++) {
        ASSERT_ALWAYS(mpz_sgn(N[side]) >= 0);
        if (N[side] == 1) {
            res[side].status = FACUL_SMOOTH;
            continue;
        }

        if (mpz_get_d(N[side]) < strategies.BB[side]) {
            res[side].status = FACUL_SMOOTH;
            res[side].primes.emplace_back(N[side]);
            continue;
        } else if (strategies.BB[side] == 0) {
            /* We encounter this in the sublat case: since we don't claim
             * to have perfectly sieved the factor base, strat.BB is
             * actually {0, 0} in that case. It means that we don't have
             * this same quick primality check. Yet we used to flatly
             * compare against B*B. It's not very consistent with the
             * fact that (of course) facul_doit_onefm obeys strat.BB.
             *
             * A slight improvement is to indeed check against B*B, but
             * check for primality on top of that. (we might want to
             * defer the primality test to after the point where the
             * other side gets sieved, though).
             */
            auto const Bd = double(strategies.B[side]);
            if (mpz_get_d(N[side]) < Bd * Bd && mpz_probab_prime_p(N[side], 1)) {
                res[side].status = FACUL_SMOOTH;
                res[side].primes.emplace_back(N[side]);
                continue;
            }
        }
        /* If none of the previous tests passed, then we hope to find a
         * factor with our many rounds of facul_doit_onefm.
         */

        std::unique_ptr<FaculModulusBase> n(FaculModulusBase::init_mpz(N[side]));

        if (!n) {
            /* This one is out of bounds, too bad. Abort. */
            res[side].status = FACUL_NOT_SMOOTH;
            return res;
        }
        composites[side].emplace_back(std::move(n));
    }

    /* Note that it's important that we do the above tests _before_ we
     * possibly bail out because there are no methods to be tried. We do
     * have test cases (F9_cofactest) which for some weird reason try
     * with an empty method list (why?), and yet we'd like to see factors
     * in the cases that can be detected easily.
     */
    if (methods.empty())
        return res;

#ifdef ENABLE_UNSAFE_FACUL_STATS
    int stats_nb_side = 0, stats_index_transition = 0;
#endif             /* ENABLE_UNSAFE_FACUL_STATS */

    /* last_i[s] is the index of last method tried on side s */
    std::vector<size_t> last_i(nsides, 0);

    for (size_t i = 0 ; i < methods.size() ; i++) {
        auto const & meth = methods[i];

#ifdef ENABLE_UNSAFE_FACUL_STATS
        stats_current_index = i - stats_nb_side * stats_index_transition;
        if (meth.is_last) {
            stats_nb_side = 1;
            stats_index_transition = i + 1;
        }
#endif /* ENABLE_UNSAFE_FACUL_STATS */
        int const side = meth.side;

        if (!res[side].primes.empty() || composites[side].size() != 1) {
            /* We just achieved a split on this side. We'll proceed with
             * facul_aux unless the other cofactor is identified as
             * non-smooth before that. */
            continue;
        }

        /* If all sides are smooth, we can exit the loop,
           otherwise we must continue with the next methods,
           since methods might be interleaved between side 0 and 1,
           thus we don't have an easy way to skip all methods for this side.
           We could do this with another representation, say methods[0][i]
           for side 0, 0 <= i < m, methods[1][j] for side 1, 0 <= j < n,
           and which_method[k] = {0, 1} for 0 <= k < m+n. */
        bool all_smooth = std::all_of(res.cbegin(), res.cend(),
                                      [](facul_result const & r) {
                                          return r.status == FACUL_SMOOTH;
                                      });
        if (all_smooth)
            return res;

        if (res[side].status == FACUL_SMOOTH)
            continue;

#ifdef ENABLE_UNSAFE_FACUL_STATS
        if (stats_current_index < STATS_LEN)
            stats_called[stats_current_index]++;
#endif /* ENABLE_UNSAFE_FACUL_STATS */

        last_i[side] = i;

        std::vector<std::unique_ptr<FaculModulusBase>> more_composites;

        res[side].status = composites[side].front()->facul_doit_onefm(
            res[side].primes, *meth.method, more_composites,
            strategies.lpb[side], strategies.BB[side], strategies.BBB[side]);

        switch(res[side].status) {
            case FACUL_NOT_SMOOTH:
                /* The cofactor is not smooth. Abort. */
                return res;
            case FACUL_MAYBE:
                if (more_composites.empty()) {
                    /* No factor found. If it was the last method for this
                     * side, then one stops the cofactorization. Otherwise, one
                     * tries with an other method.
                     */
                    if (meth.is_last) {
                        res[side].status = FACUL_NOT_SMOOTH;
                        // res[side].composites.emplace_back(N[side]);
                        return res;
                    }
                } else {
                    /* a composite split! */
                    std::swap(composites[side], more_composites);
                }
                continue;
            case FACUL_SMOOTH:
                continue;
        }
    }

    for(size_t side = 0 ; side < nsides ; side++) {
        if (res[side].status != FACUL_MAYBE)
            continue;

        res[side].status = facul_aux(res[side].primes, composites[side],
                strategies, methods, last_i[side] + 1, side);

        if (res[side].status == FACUL_NOT_SMOOTH)
            return res;
    }

    for (size_t side = 0; side < nsides; side++) {
        std::ranges::sort(res[side].primes);
    }

    return res;
}

facul_result facul(cxx_mpz const & N, facul_strategy_oneside const & strategy)
{
    facul_result res;

#ifdef PARI
    gmp_fprintf(stderr, "%Zd", N);
#endif

    if (N <= 0)
        return FACUL_NOT_SMOOTH;
    if (N == 1)
        return FACUL_SMOOTH;

    /* Use the fastest modular arithmetic that's large enough for this input */
    std::unique_ptr<FaculModulusBase> m(FaculModulusBase::init_mpz(N));

    /* If the composite does not fit into our modular arithmetic, return
       no factor */
    if (!m)
        return FACUL_NOT_SMOOTH;

    std::vector<std::unique_ptr<FaculModulusBase>> todo;
    todo.emplace_back(std::move(m));

    /* We used to call facul_doit. In fact, facul_aux does exactly the
     * same thing, so let's use that instead.
     *
     * The catch is that facul_aux takes a vector of facul_method_side
     * objects. Before we refactor this with facul_method, just do a
     * copy.
     */
    std::vector<facul_method_side> methods;
    methods.reserve(strategy.methods.size());
    for(auto const & m : strategy.methods) {
        methods.emplace_back(&m, 0);
    }
    const facul_strategies_base strat_base { strategy };

    res.status = facul_aux(res.primes, todo, strat_base, methods, 0, 0);

    /* Sort the factors we found */
    std::ranges::sort(res.primes);

    return res;
}

