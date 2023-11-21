#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <sys/param.h>
#include <cstdio>                        // for printf, fprintf, stderr
#include <cstdlib>                       // for exit, qsort, EXIT_FAILURE
#include <sstream>
#include <tuple>                          // for tie, tuple
#include <vector>                         // for vector
#include "cxx_mpz.hpp"
#include "arith-hard.hpp" // IWYU pragma: keep
#include "lingen_bmstatus.hpp"            // for bmstatus
#include "lingen_bw_dimensions.hpp"
#include "lingen_call_companion.hpp"      // for lingen_call_companion
#include "lingen_expected_pi_length.hpp"
#include "lingen_qcode_prime.hpp"
#include "macros.h"                       // for ASSERT_ALWAYS, ASSERT, icei...
#include "tree_stats.hpp"                 // for tree_stats, tree_stats::sen...

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
/* We sort only with respect to the global delta[] parameter. As it turns
 * out, we also access the column index in the same aray and sort with
 * respect to it, but this is only cosmetic.
 *
 * Note that unlike what is asserted in other coiped of the code, sorting
 * w.r.t. the local delta[] value is completely useless. Code which
 * relies on this should be fixed.
 */

typedef int (*sortfunc_t) (const void*, const void*);

static int lexcmp2(const int x[2], const int y[2])
{
    for(int i = 0 ; i < 2 ; i++) {
        int d = x[i] - y[i];
        if (d) return d;
    }
    return 0;
}

/* }}} */

/*}}}*/

matpoly bw_lingen_basecase_raw(bmstatus & bm, matpoly const & E) /*{{{*/
{
    int generator_found = 0;

    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    matpoly::arith_hard * ab = & d.ab;
    ASSERT(E.m == m);
    ASSERT(E.n == b);

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int mi, ma;
    std::tie(mi, ma) = get_minmax_delta(bm.delta);

    unsigned int pi_room_base = expected_pi_length(d, bm.delta, E.get_size());

    matpoly pi(ab, b, b, pi_room_base);
    pi.set_size(pi_room_base);

    /* Also keep track of the
     * number of coefficients for the columns of pi. Set pi to Id */

    std::vector<unsigned int> pi_lengths(b, 1);
    std::vector<unsigned int> pi_real_lengths(b, 1);
    pi.set_constant_ui(1);

    for(unsigned int i = 0 ; i < b ; i++) {
        pi_lengths[i] = 1;
        /* Fix for check 21744-bis. Our column priority will allow adding
         * to a row if its delta is above the average. But then this
         * might surprise us, because in truth we dont start with
         * something of large degree here. So we're bumping pi_length a
         * little bit to accomodate for this situation.
         */
        pi_lengths[i] += bm.delta[i] - mi;
    }

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    std::vector<unsigned int> pivots(m, 0);
    std::vector<int> is_pivot(b, 0);

    matpoly e(ab, m, b, 1);
    e.set_size(1);

    matpoly T(ab, b, b, 1);
    int (*ctable)[2] = new int[b][2];

    for (unsigned int t = 0; t < E.get_size() ; t++, bm.t++) {

        /* {{{ Update the columns of e for degree t. Save computation
         * time by not recomputing those which can easily be derived from
         * previous iteration. Notice that the columns of e are exactly
         * at the physical positions of the corresponding columns of pi.
         */

        std::vector<unsigned int> todo;
        for(unsigned int j = 0 ; j < b ; j++) {
            if (is_pivot[j]) continue;
            /* We should never have to recompute from pi using discarded
             * columns. Discarded columns should always correspond to
             * pivots */
            ASSERT_ALWAYS(bm.lucky[j] >= 0);
            todo.push_back(j);
        }
        /* icc openmp doesn't grok todo.size() as being a constant
         * loop bound */
        unsigned int nj = todo.size();

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            auto e_ur = ab->alloc<arith_hard::elt_ur_for_addmul>(1);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int jl = 0 ; jl < nj ; ++jl) {
                for(unsigned int i = 0 ; i < m ; ++i) {
                    unsigned int j = todo[jl];
                    unsigned int lj = MIN(pi_real_lengths[j], t + 1);
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
        /* }}} */

        /* {{{ check for cancellations */

        unsigned int newluck = 0;
        for(unsigned int jl = 0 ; jl < todo.size() ; ++jl) {
            unsigned int j = todo[jl];
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

            int luck_mini = expected_pi_length(d);
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
        /* }}} */

        if (generator_found) break;

        /* {{{ Now see in which order I may look at the columns of pi, so
         * as to keep the nominal degrees correct. In contrast with what
         * we used to do before, we no longer apply the permutation to
         * delta. So the delta[] array keeps referring to physical
         * indices, and we'll tune this in the end. */
        for(unsigned int j = 0; j < b; j++) {
            ctable[j][0] = bm.delta[j];
            ctable[j][1] = j;
        }
        qsort(ctable, b, 2 * sizeof(int), (sortfunc_t) & lexcmp2);
        /* }}} */

        /* {{{ Now do Gaussian elimination */

        /*
         * The matrix T is *not* used for actually storing the product of
         * the transvections, just the *list* of transvections. Then,
         * instead of applying them row-major, we apply them column-major
         * (abiding by the ordering of pivots), so that we get a better
         * opportunity to do lazy reductions.
         */

        T.set_constant_ui(1);

        is_pivot.assign(b, 0);
        unsigned int r = 0;

        std::vector<unsigned int> pivot_columns;
        /* Loop through logical indices */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            unsigned int u = 0;
            /* {{{ Find the pivot */
            for( ; u < m ; u++) {
                if (!ab->is_zero(e.coeff(u, j, 0)))
                    break;
            }
            if (u == m) continue;
            ASSERT(r < m);
            /* }}} */
            pivots[r++] = j;
            is_pivot[j] = 1;
            pivot_columns.push_back(j);
            /* {{{ Cancel this coeff in all other columns. */
            auto inv = ab->alloc();
            int rc = ab->inverse(*inv, e.coeff(u, j, 0));
            if (!rc) {
                std::ostringstream os;
                ab->cxx_out(os, *inv);
                fprintf(stderr, "Error, found a factor of the modulus: %s\n",
                        os.str().c_str());
                exit(EXIT_FAILURE);
            }
            ab->neg(*inv, *inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                if (ab->is_zero(e.coeff(u, k, 0)))
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                auto lambda = ab->alloc();
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
        /* }}} */

        /* {{{ apply the transformations, using the transvection
         * reordering trick */

        /* non-pivot columns are only added to and never read, so it does
         * not really matter where we put their computation, provided
         * that the columns that we do read are done at this point.
         */
        for(unsigned int jl = 0; jl < b; jl++) {
            unsigned int j = ctable[jl][1];
            if (!is_pivot[j])
                pivot_columns.push_back(j);
        }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            auto tmp = ab->alloc<arith_hard::elt_ur_for_addmul>();
            auto tmp_pi = ab->alloc<arith_hard::elt_ur_for_addmul>();

            for(unsigned int jl = 0 ; jl < b ; ++jl) {
                unsigned int j = pivot_columns[jl];
                /* compute column j completely. We may put this interface in
                 * matpoly, but it's really special-purposed, to the point
                 * that it really makes little sense IMO
                 *
                 * Beware: operations on the different columns are *not*
                 * independent, here ! Operations on the different degrees,
                 * on the other hand, are. As well of course as the
                 * operations on the different entries in each column.
                 */

#ifndef NDEBUG
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
                {
                    for(unsigned int kl = m ; kl < b ; kl++) {
                        unsigned int k = pivot_columns[kl];
                        ASSERT_ALWAYS(ab->cmp(T.coeff(k, j, 0), k==j) == 0);
                    }
                    for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                        unsigned int k = pivot_columns[kl];
                        if (ab->is_zero(T.coeff(k, j, 0))) continue;
                        ASSERT_ALWAYS(pi_lengths[k] <= pi_lengths[j]);
                        pi_real_lengths[j] = std::max(pi_real_lengths[k], pi_real_lengths[j]);
                    }
                }
#endif

                /* Icc 2019 synthetizes a pragma omp single around
                 * accesses to pi_real_lengths. I don't think it makes
                 * sense in that particular case, it's fine enough here
                 * to read the data now after the critical section above.
                 */
                unsigned int dummy = pi_real_lengths[j];

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
                for(unsigned int i = 0 ; i < b ; i++) {
                    for(unsigned int s = 0 ; s < dummy ; s++) {
                        arith_hard::elt & piijs = pi.coeff(i, j, s);

                        ab->set(*tmp_pi, piijs);

                        for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                            unsigned int k = pivot_columns[kl];
                            /* TODO: if column k was already a pivot on previous
                             * turn (which could happen, depending on m and n),
                             * then the corresponding entry is probably zero
                             * (exact condition needs to be written more
                             * accurately).
                             */

                            arith_hard::elt const & Tkj = T.coeff(k, j, 0);
                            if (ab->is_zero(Tkj)) continue;
                            /* pi[i,k] has length pi_lengths[k]. Multiply
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
        /* }}} */

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_real_lengths[j] >= pi.capacity()) {
                if (!generator_found) {
                    pi.realloc(pi.capacity() + MAX(pi.capacity() / (m+n), 1));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm.t,
                            100 * pi.capacity() / pi_room_base);
                    for(unsigned int j = 0; j < b; j++)
                        printf(" %u", pi_real_lengths[j]);
                    printf("\n");
                } else {
                    ASSERT_ALWAYS(bm.lucky[j] <= 0);
                    if (bm.lucky[j] == 0)
                        printf("t=%u, column %u discarded from now on\n",
                                bm.t, j);
                    bm.lucky[j] = -1;
                    pi_lengths[j]++;
                    pi_real_lengths[j]++;
                    bm.delta[j]++;
                    continue;
                }
            }
            pi.multiply_column_by_x(j, pi_real_lengths[j]);
            pi_real_lengths[j]++;
            pi_lengths[j]++;
            bm.delta[j]++;
        }
        /* }}} */
    }

    delete[] ctable;

    unsigned int pisize = 0;
    for(unsigned int j = 0; j < b; j++) {
        if (pi_real_lengths[j] > pisize)
            pisize = pi_real_lengths[j];
    }
    /* Given the structure of the computation, there's no reason for the
     * initial estimate to go wrong.
     */
    ASSERT_ALWAYS(pisize <= pi.capacity());
    pi.set_size(pisize);

    for(unsigned int j = 0; j < b; j++) {
        for(unsigned int k = pi_real_lengths[j] ; k < pi.get_size() ; k++) {
            for(unsigned int i = 0 ; i < b ; i++) {
                ASSERT_ALWAYS(ab->is_zero(pi.coeff(i, j, k)));
            }
        }
    }

    bm.done = generator_found;
    return pi;
}
/* }}} */

/* wrap this up */
matpoly
bw_lingen_basecase(bmstatus & bm, matpoly & E)
{
    lingen_call_companion const & C = bm.companion(bm.depth(), E.get_size());
    tree_stats::sentinel dummy(bm.stats, "basecase", E.get_size(), C.total_ncalls, true);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    matpoly pi = bw_lingen_basecase_raw(bm, E);
    bm.stats.end_smallstep();
    E = matpoly();
    return pi;
}

void test_basecase(matpoly::arith_hard * ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate)/*{{{*/
{
    /* used by testing code */
    bmstatus bm(m,n,ab->characteristic());
    unsigned int t0 = iceildiv(m,n);
    bm.set_t0(t0);
    matpoly E(ab, m, m+n, L);
    E.zero_pad(L);
    E.fill_random(0, L, rstate);
    bw_lingen_basecase_raw(bm, E);
}/*}}}*/
