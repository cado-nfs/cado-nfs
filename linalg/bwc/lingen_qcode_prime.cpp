#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <utility>
#include <sstream>
#include <tuple>
#include <vector>
#include <algorithm>

#include "gmp_aux.h"
#include "cxx_mpz.hpp"
#include "arith-hard.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_bw_dimensions.hpp"
#include "lingen_call_companion.hpp"
#include "lingen_expected_pi_length.hpp"
#include "lingen_qcode_prime.hpp"
#include "macros.h"
#include "tree_stats.hpp"

/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

/*{{{ basecase */

/* {{{ col sorting */
/* Note on col sorting: we sort only with respect to the global delta[]
 * parameter. As it turns out, we also access the column index in the
 * same aray and sort with respect to it, but this is only cosmetic.
 *
 * Note that unlike what is asserted in other copies of the code, sorting
 * w.r.t. the local delta[] value is completely useless. Code which
 * relies on this should be fixed.
 */

/* }}} */

/*}}}*/

struct bw_lingen_basecase_raw_object {
    bmstatus<false> & bm;
    matpoly<false> const & E;
    bool generator_found = false;
    bw_dimensions<false> & d;
    unsigned int const m;
    unsigned int const n;
    unsigned int const b;
    matpoly<false>::arith_hard * ab;
    unsigned int const pi_room_base;
    /* This keeps tracks of the weights of the columns of pi, both
     * globally and locally */
    std::vector<unsigned int> pi_length_global;
    std::vector<unsigned int> pi_length_local;
    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    std::vector<bool> is_pivot;

    bw_lingen_basecase_raw_object(bmstatus<false> & bm, matpoly<false> const & E)
        : bm(bm)
        , E(E)
        , d(bm.d)
        , m(d.m)
        , n(d.n)
        , b(m+n)
        , ab(&d.ab)
        , pi_room_base(expected_pi_length(d, bm.delta, E.get_size()))
        , pi_length_global(b, 1)
        , pi_length_local(b, 1)
        , is_pivot(b, false)
    {   /* {{{ */
        ASSERT(E.m == m);
        ASSERT(E.n == b);
        unsigned int mi, ma;
        std::tie(mi, ma) = get_minmax_delta(bm.delta);

        for(unsigned int i = 0 ; i < b ; i++) {
            /* see bug 21744-bis. Our column priority will allow adding
             * to a row if its delta is above the average. But then this
             * might surprise us, because in truth we dont start with
             * something of large degree here. So we're bumping pi_length a
             * little bit to accomodate for this situation.
             */
            pi_length_global[i] = 1 + bm.delta[i] - mi;
        }

    }   /* }}} */


    std::vector<unsigned int> columns_needing_recomputation()   /* {{{ */
    {
        std::vector<unsigned int> todo;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            /* We should never have to recompute from pi using discarded
             * columns. Discarded columns should always correspond to
             * pivots */
            ASSERT_ALWAYS(bm.lucky[j] >= 0);
            todo.push_back(j);
        }
        return todo;
    }   /* }}} */

    /* This recomputes the columns of degree T of e*pi, because the
     * previous computation ended with a complete cancellation at this
     * column.
     */
    void recompute_columns(     /* {{{ */
            unsigned int t,
            matpoly<false> & e,
            matpoly<false> const & pi,
            std::vector<unsigned int> const & todo)
    {
        /* icc openmp doesn't grok todo.size() as being a constant
         * loop bound */
        unsigned int const nj = todo.size();

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            auto * e_ur = ab->alloc<arith_hard::elt_ur_for_addmul>(1);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int jl = 0 ; jl < nj ; ++jl) {
                for(unsigned int i = 0 ; i < m ; ++i) {
                    unsigned int const j = todo[jl];
                    unsigned int const lj = std::min(pi_length_local[j], t + 1);
                    ab->set_zero(*e_ur);
                    for(unsigned int k = 0 ; k < b ; ++k) {
                        for(unsigned int s = 0 ; s < lj ; s++) {
                            ab->addmul_ur(*e_ur,
                                    E.coeff(i, k, t - s),
                                    pi.coeff(k, j, s));
                        }
                    }
                    ab->reduce(e.coeff(i, j, 0), *e_ur);
                }
            }

            ab->free(e_ur);
        }
    } /* }}} */


    unsigned int check_cancellations( /* {{{ */
            matpoly<false> const & e,
            std::vector<unsigned int> const & todo)
    {

        unsigned int newluck = 0;
        for(auto const j : todo) {
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                nz += ab->is_zero(e.coeff(i, j, 0));
            }
            if (nz == m) {
                newluck++, bm.lucky[j]++;
            } else if (bm.lucky[j] > 0) {
                bm.lucky[j] = 0;
            }
        }


        if (newluck) {
            /* If newluck == n, then we probably have a generator. We add an
             * extra guarantee. newluck==n, for a total of k iterations in a
             * row, means m*n*k coefficients cancelling magically. We would
             * like this to be impossible by mere chance. Thus we want n*k >
             * luck_mini, which can easily be checked */

            int const luck_mini = expected_pi_length<false>(d);
            unsigned int luck_sure = 0;

            printf("t=%d, canceled columns:", bm.t);
            for(unsigned int j = 0 ; j < b ; j++) {
                if (bm.lucky[j] > 0) {
                    printf(" %u", j);
                    luck_sure += bm.lucky[j] >= luck_mini;
                }
            }

            if (newluck == n && luck_sure == n) {
                if (!generator_found) {
                    printf(", complete generator found, for sure");
                }
                generator_found = 1;
            }
            printf(".\n");
        }

        return generator_found;
    } /* }}} */

    std::vector<unsigned int> get_column_order() /* {{{ */
    {
        /* Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        std::vector<std::pair<int, int>> ctable_pre;
        ctable_pre.reserve(b);
        for(unsigned int j = 0; j < b; j++)
            ctable_pre.emplace_back(bm.delta[j], j);
        std::ranges::sort(ctable_pre);
        std::vector<unsigned int> ctable;
        ctable.reserve(ctable_pre.size());
        for(auto const & j : ctable_pre)
            ctable.emplace_back(j.second);
        return ctable;
    } /* }}} */

    std::pair<matpoly<false>, std::vector<unsigned int>> compute_transvections(
            matpoly<false> & e,
            std::vector<unsigned int> const & ctable)
    /* {{{ */
    {
        // return the list of tranvections as well as the list of pivot
        // columns

        /*
         * The matrix T is *not* used for actually storing the product of
         * the transvections, just the *list* of transvections. Then,
         * instead of applying them row-major, we apply them column-major
         * (abiding by the ordering of pivots), so that we get a better
         * opportunity to do lazy reductions.
         */
        matpoly<false> T(ab, b, b, 1);

        T.set_constant_ui(1);

        std::vector<unsigned int> pivot_columns;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int const j = ctable[jl];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                if (!ab->is_zero(e.coeff(u, j, 0)))
                    break;
            }
            if (u == m) continue;
            ASSERT_ALWAYS(pivot_columns.size() < m);
            /* }}} */
            pivot_columns.push_back(j);
            /* {{{ Cancel this coeff in all other columns. */
            auto * inv = ab->alloc();
            int const rc = ab->inverse(*inv, e.coeff(u, j, 0));
            if (!rc) {
                std::ostringstream os;
                ab->cxx_out(os, *inv);
                fprintf(stderr, "Error, found a factor of the modulus: %s\n",
                        os.str().c_str());
                exit(EXIT_FAILURE);
            }
            ab->neg(*inv, *inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int const k = ctable[kl];
                if (ab->is_zero(e.coeff(u, k, 0)))
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                auto * lambda = ab->alloc();
                ab->mul(*lambda, *inv, e.coeff(u, k, 0));
                ASSERT(bm.delta[j] <= bm.delta[k]);
                /* {{{ Apply on both e and pi */
                for(unsigned int i = 0 ; i < m ; i++) {
                    ab->addmul_and_reduce(e.coeff(i, k, 0), *lambda, e.coeff(i, j, 0));
                }
                if (bm.lucky[k] < 0) {
                    /* This column is already discarded, don't bother */
                    continue;
                }
                if (bm.lucky[j] < 0) {
                    /* This column is discarded. This is going to
                     * invalidate another column of pi. Not a problem,
                     * unless it's been marked as lucky previously ! */
                    ASSERT_ALWAYS(bm.lucky[k] <= 0);
                    printf("Column %u discarded from now on (through addition from column %u)\n", k, j);
                    bm.lucky[k] = -1;
                    continue;
                }
                /* We do *NOT* really update T. T is only used as
                 * storage!
                 */
                ab->set(T.coeff(j, k, 0), *lambda);
                /* }}} */
                ab->free(lambda);
            }
            ab->free(inv);
            /* }}} */
        }

        return { std::move(T), std::move(pivot_columns) };
    }    /* }}} */


    void apply_transvections(matpoly<false> & pi, matpoly<false> const & T,
            std::vector<unsigned int> const & columns)
    {
        /* {{{ apply the transformations, using the transvection
         * reordering trick */


#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            auto * tmp = ab->alloc<arith_hard::elt_ur_for_addmul>();
            auto * tmp_pi = ab->alloc<arith_hard::elt_ur_for_addmul>();

            for(unsigned int jl = 0 ; jl < b ; ++jl) {
                unsigned int const j = columns[jl];
                /* compute column j completely. We may put this interface in
                 * matpoly, but it's really special-purposed, to the point
                 * that it really makes little sense IMO
                 *
                 * Beware: operations on the different columns are *not*
                 * independent, here ! Operations on the different degrees,
                 * on the other hand, are. As well of course as the
                 * operations on the different entries in each column.
                 */

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
                {
                    for(unsigned int kl = m ; kl < b ; kl++) {
                        unsigned int const k = columns[kl];
                        ASSERT_ALWAYS(ab->cmp(T.coeff(k, j, 0), k==j) == 0);
                    }
                    for(unsigned int kl = 0 ; kl < std::min(m,jl) ; kl++) {
                        unsigned int const k = columns[kl];
                        if (ab->is_zero(T.coeff(k, j, 0))) continue;
                        // note that pi_length_global is for the _global_
                        // deltas!
                        ASSERT_ALWAYS(pi_length_global[k] <= pi_length_global[j]);
                        // This line is not a debug check! It's super
                        // important! See #30105
                        pi_length_local[j] = std::max(pi_length_local[k], pi_length_local[j]);
                    }
                }

                /* Icc 2019 synthetizes a pragma omp single around
                 * accesses to pi_length_local. I don't think it makes
                 * sense in that particular case, it's fine enough here
                 * to read the data now after the critical section above.
                 */
                unsigned int const dummy = pi_length_local[j];

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
                for(unsigned int i = 0 ; i < b ; i++) {
                    for(unsigned int s = 0 ; s < dummy ; s++) {
                        arith_hard::elt & piijs = pi.coeff(i, j, s);

                        ab->set(*tmp_pi, piijs);

                        for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                            unsigned int const k = columns[kl];
                            /* TODO: if column k was already a pivot on previous
                             * turn (which could happen, depending on m and n),
                             * then the corresponding entry is probably zero
                             * (exact condition needs to be written more
                             * accurately).
                             */

                            arith_hard::elt const & Tkj = T.coeff(k, j, 0);
                            if (ab->is_zero(Tkj)) continue;
                            /* pi[i,k] has length pi_length_global[k]. Multiply
                             * that by T[k,j], which is a constant. Add
                             * to the unreduced thing. We don't have an
                             * mpfq api call for that operation.
                             */
                            arith_hard::elt const & piiks = pi.coeff(i, k, s);
                            ab->mul_ur(*tmp, piiks, Tkj);
                            ab->add(*tmp_pi, *tmp);
                        }
                        ab->reduce(piijs, *tmp_pi);
                    }
                }
            }
            ab->free(tmp);
            ab->free(tmp_pi);
        }
    } /* }}} */

    std::vector<unsigned int> expand_to_all_columns_and_update_flags(
            std::vector<unsigned int> const & pivot_columns,
            std::vector<unsigned int> const & ctable)
    {
        auto all = pivot_columns;

        is_pivot.assign(b, false);
        for(auto const j : pivot_columns)
            is_pivot[j] = true;

        /* non-pivot columns are only added to and never read, so it does
         * not really matter where we put their computation, provided
         * that the columns that we do read are done at this point.
         */
        for(auto const j : ctable)
            if (!is_pivot[j])
                all.push_back(j);
        return all;
    }

    void shift_pivot_columns_in_pi(matpoly<false> & pi)
    {
        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_length_local[j] >= pi.capacity()) {
                if (!generator_found) {
                    pi.realloc(pi.capacity() + std::max(pi.capacity() / (m+n),
                                size_t(1)));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm.t,
                            100 * pi.capacity() / pi_room_base);
                    for(unsigned int j = 0; j < b; j++)
                        printf(" %u", pi_length_local[j]);
                    printf("\n");
                } else {
                    ASSERT_ALWAYS(bm.lucky[j] <= 0);
                    if (bm.lucky[j] == 0)
                        printf("t=%u, column %u discarded from now on\n",
                                bm.t, j);
                    bm.lucky[j] = -1;
                    pi_length_global[j]++;
                    pi_length_local[j]++;
                    bm.delta[j]++;
                    continue;
                }
            }
            pi.multiply_column_by_x(j, pi_length_local[j]);
            pi_length_local[j]++;
            pi_length_global[j]++;
            bm.delta[j]++;
        }
    }
    /* }}} */

    matpoly<false> get_pi()// {{{
    {
        matpoly<false> pi(ab, b, b, int(pi_room_base));
        pi.set_size(pi_room_base);

        /* NOTE that this sets pi.size=1 */
        pi.set_constant_ui(1);

        matpoly<false> e(ab, m, b, 1);
        e.set_size(1);

        for (unsigned int t = 0; t < E.get_size() ; t++, bm.t++) {

            /*  Update the columns of e for degree t. Save computation
             * time by not recomputing those which can easily be derived from
             * previous iteration. Notice that the columns of e are exactly
             * at the physical positions of the corresponding columns of pi.
             */

            auto todo = columns_needing_recomputation();
            recompute_columns(t, e, pi, todo);

            if (check_cancellations(e, todo)) break;

            auto ctable = get_column_order();

            matpoly<false> T;
            std::vector<unsigned int> pivot_columns;

            std::tie(T, pivot_columns) = compute_transvections(e, ctable);

            unsigned int r = pivot_columns.size();
            ASSERT_ALWAYS(r == m);

            auto all = expand_to_all_columns_and_update_flags(
                    pivot_columns, ctable);
            apply_transvections(pi, T, all);

            shift_pivot_columns_in_pi(pi);

            unsigned int pisize = 0;
            for(auto ell : pi_length_local)
                pisize = std::max(pisize, ell);
            /* Given the structure of the computation, there's no reason for the
             * initial estimate to go wrong.
             */
            ASSERT_ALWAYS(pisize <= pi.capacity());
            pi.set_size(pisize);
        }


        for(unsigned int j = 0; j < b; j++) {
            for(unsigned int k = pi_length_local[j] ; k < pi.get_size() ; k++) {
                for(unsigned int i = 0 ; i < b ; i++) {
                    ASSERT_ALWAYS(ab->is_zero(pi.coeff(i, j, k)));
                }
            }
        }

        bm.done = generator_found;
        return pi;
    }// }}}
};

matpoly<false> bw_lingen_basecase_raw(bmstatus<false> & bm, matpoly<false> const & E)
{
    return bw_lingen_basecase_raw_object(bm, E).get_pi();
}
/* }}} */

/* wrap this up */
matpoly<false>
bw_lingen_basecase(bmstatus<false> & bm, matpoly<false> & E)
{
    lingen_call_companion const & C = bm.companion(bm.depth, E.get_size());
    tree_stats::sentinel const dummy(bm.stats, "basecase", E.get_size(), C.total_ncalls, true);
    bmstatus<false>::depth_sentinel ddummy(bm);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    matpoly<false> pi = bw_lingen_basecase_raw(bm, E);
    bm.stats.end_smallstep();
    E = matpoly<false>();
    return pi;
}

void test_basecase(matpoly<false>::arith_hard * ab, unsigned int m, unsigned int n, size_t L, cxx_gmp_randstate & rstate)/*{{{*/
{
    /* used by testing code */
    bmstatus<false> bm(m,n,ab->characteristic());
    unsigned int const t0 = iceildiv(m,n);
    bm.set_t0(t0);
    matpoly<false> E(ab, m, m+n, L);
    E.zero_pad(L);
    E.fill_random(0, L, rstate);
    bw_lingen_basecase_raw(bm, E);
}/*}}}*/
