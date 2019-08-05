/* Copyright (C) 1999--2007 Emmanuel Thom'e --- see LICENSE file */
#include "cado.h"

#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <stdexcept>
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif
#include <cassert>
#include <fstream>

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "mpfq_layer.h"
#include "memusage.h"

/* lingen-matpoly is the default code. */
#include "lingen-matpoly.hpp"

#define ENABLE_MPI_LINGEN

#ifdef ENABLE_MPI_LINGEN
#include "lingen-bigmatpoly.hpp"
#endif

#include "lingen-bigmatpoly-ft.hpp"

#include "bw-common.h"		/* Handy. Allows Using global functions
                                 * for recovering parameters */
#include "plingen.hpp"
#include "plingen-tuning.hpp"
#include "logline.h"
#include "tree_stats.hpp"
#include "sha1.h"

/* Call tree for methods within this program:
 *
 * Two separate entry points.
 *
 * bw_biglingen_collective        when need to go collective
 *     |<--> bw_biglingen_recursive , loops to bw_biglingen_collective
 *     \->bw_lingen_single               when it makes sense to do this locally again
 * bw_lingen_single                      when computation can be done locally
 *    |<->bw_lingen_recursive
 *    |->bw_lingen_basecase
 *
 */

#define MPI_MY_SIZE_T   MPI_UNSIGNED_LONG

static unsigned int display_threshold = 10;
static int with_timings = 0;

/* This is an indication of the number of bytes we read at a time for A
 * (input) and F (output) */
static unsigned int io_block_size = 1 << 20;

/* If non-zero, then reading from A is actually replaced by reading from
 * a random generator */
static unsigned int random_input_length = 0;

static int split_input_file = 0;  /* unsupported ; do acollect by ourselves */
static int split_output_file = 0; /* do split by ourselves */

gmp_randstate_t rstate;

static const char * checkpoint_directory;
static unsigned int checkpoint_threshold = 100;
static int save_gathered_checkpoints = 0;

static int allow_zero_on_rhs = 0;

int rank0_exit_code = EXIT_SUCCESS;


int global_flag_ascii = 0;
int global_flag_tune = 0;

struct bmstatus {/*{{{*/
    bw_dimensions d;
    unsigned int t;
    std::vector<int> lucky;

    double t_basecase;
    double t_mp;
    double t_mul;
    double t_cp_io;

    unsigned int lingen_threshold;
    unsigned int lingen_mpi_threshold;
    int mpi_dims[2]; /* mpi_dims[0] = mpi[0] * thr[0] */
    MPI_Comm com[3]; /* [0]: MPI_COMM_WORLD, reordered.
                        [1]: row-wise
                        [2]: column-wise */

    lingen_hints_t hints;

    tree_stats stats;

    int depth() const { return stats.non_transition_depth(); }

    bmstatus(unsigned int m, unsigned int n)/*{{{*/
    {
        memset(&d, 0, sizeof(bw_dimensions));
        d.m = m;
        d.n = n;
        lucky.assign(m+n, 0);
    }/*}}}*/
    lingen_call_companion const& companion(int depth, size_t L)/*{{{*/
    {
        lingen_hints_t::key_type K { depth, L };

        if (hints.find(K) != hints.end())
            return hints[K];

        fprintf(stderr, "# No tuned configuration for"
                " depth=%d input_length=%zu\n",
                depth, L);

        struct same_depth {
            int d;
            bool operator()(lingen_hints_t::value_type const & kv) const {
                return kv.first.depth == d;
            }
        };

        auto it = std::find_if(hints.begin(), hints.end(), same_depth { depth } );
        if (it == hints.end()) {
            fprintf(stderr, "# No tuned configuration for"
                    " depth=%d !!!\n", depth);
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "# Using nearby configuration for"
                " depth=%d input_length=%zu\n",
                it->first.depth, it->first.L);

        hints[K]=it->second;

        return it->second;
    }/*}}}*/
    bool recurse(int depth, size_t L) {/*{{{*/
        return companion(depth, L).recurse;
    }/*}}}*/
    bool recurse(size_t L) {/*{{{*/
        return companion(depth(), L).recurse;
    }/*}}}*/
};/*}}}*/

void plingen_decl_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "ascii",
            "read and write data in ascii");
    param_list_decl_usage(pl, "timings",
            "provide timings on all output lines");
    param_list_decl_usage(pl, "tune",
            "activate tuning mode");
    param_list_decl_usage(pl, "allow_zero_on_rhs",
            "do not cry if the generator corresponds to a zero contribution on the RHS vectors");

    /* we must be square ! And thr is not supported. */
    param_list_decl_usage(pl, "mpi", "number of MPI nodes across which the execution will span, with mesh dimensions");
    param_list_decl_usage(pl, "thr", "number of threads (on each node) for the program, with mesh dimensions");

    param_list_decl_usage(pl, "nrhs",
            "number of columns to treat differently, as corresponding to rhs vectors");
    param_list_decl_usage(pl, "rhs",
            "file with rhs vectors (only the header is read)");

    param_list_decl_usage(pl, "afile",
            "input sequence file");
    param_list_decl_usage(pl, "random-input-with-length",
            "use surrogate for input");
    param_list_decl_usage(pl, "split-input-file",
            "work with split files on input");
    param_list_decl_usage(pl, "split-output-file",
            "work with split files on output");
    param_list_decl_usage(pl, "random_seed",
            "seed the random generator");
    param_list_decl_usage(pl, "ffile",
            "output generator file");

    param_list_decl_usage(pl, "checkpoint-directory",
            "where to save checkpoints");
    param_list_decl_usage(pl, "checkpoint-threshold",
            "threshold for saving checkpoints");
    param_list_decl_usage(pl, "display-threshold",
            "threshold for outputting progress lines");
    param_list_decl_usage(pl, "io-block-size",
            "chunk size for reading the input or writing the output");

    param_list_decl_usage(pl, "lingen_mpi_threshold",
            "use MPI matrix operations above this size");
    param_list_decl_usage(pl, "lingen_threshold",
            "use recursive algorithm above this size");
    param_list_decl_usage(pl, "save_gathered_checkpoints",
            "save global checkpoints files, instead of per-job files");

    param_list_configure_switch(pl, "--tune", &global_flag_tune);
    param_list_configure_switch(pl, "--ascii", &global_flag_ascii);
    param_list_configure_switch(pl, "--timings", &with_timings);
    param_list_configure_alias(pl, "seed", "random_seed");

    plingen_tuning_decl_usage(pl);
    tree_stats::declare_usage(pl);
}/*}}}*/

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

std::tuple<unsigned int, unsigned int> get_minmax_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(auto x : delta) {
        if (x > maxdelta) maxdelta = x;
        if (x < mindelta) mindelta = x;
    }
    return std::make_tuple(mindelta, maxdelta);
}/*}}}*/
unsigned int get_min_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta(delta);
    return mindelta;
}/*}}}*/
unsigned int get_max_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta(delta);
    return maxdelta;
}/*}}}*/
std::tuple<unsigned int, unsigned int> get_minmax_delta_on_solutions(bmstatus & bm, std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(unsigned int j = 0 ; j < bm.d.m + bm.d.n ; j++) {
        if (bm.lucky[j] <= 0) continue;
        if (delta[j] > maxdelta) maxdelta = delta[j];
        if (delta[j] < mindelta) mindelta = delta[j];
    }
    return std::make_tuple(mindelta, maxdelta);
}/*}}}*/
unsigned int get_max_delta_on_solutions(bmstatus & bm, std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta_on_solutions(bm, delta);
    return maxdelta;
}/*}}}*/



static inline unsigned int expected_pi_length(bw_dimensions & d, unsigned int len = 0)/*{{{*/
{
    /* The idea is that we want something which may account for something
     * exceptional, bounded by probability 2^-64. This corresponds to a
     * column in e (matrix of size m*b) to be spontaneously equal to
     * zero. This happens with probability (#K)^-m.
     * The solution to
     * (#K)^(-m*x) > 2^-64
     * is m*x*log_2(#K) < 64
     *
     * We thus need to get an idea of the value of log_2(#K).
     *
     * (Note that we know that #K^abgroupsize(ab) < 2^64, but that bound
     * might be very gross).
     *
     * The same esitmate can be used to appreciate what we mean by
     * ``luck'' in the end. If a column happens to be zero more than
     * expected_pi_length(d,0) times in a row, then the cause must be
     * more than sheer luck, and we use it to detect generating rows.
     */

    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    abdst_field ab MAYBE_UNUSED = d.ab;
    unsigned int res = 1 + iceildiv(len * m, b);
    cxx_mpz p;
    abfield_characteristic(ab, (mpz_ptr) p);
    unsigned int l;
    if (mpz_cmp_ui(p, 1024) >= 0) {
        l = mpz_sizeinbase(p, 2);
        l *= abfield_degree(ab);    /* roughly log_2(#K) */
    } else {
        mpz_pow_ui(p, p, abfield_degree(ab));
        l = mpz_sizeinbase(p, 2);
    }
    // unsigned int safety = iceildiv(abgroupsize(ab), m * sizeof(abelt));
    unsigned int safety = iceildiv(64, m * l);
    return res + safety;
}/*}}}*/
static inline unsigned int expected_pi_length(bw_dimensions & d, std::vector<unsigned int> const & delta, unsigned int len)/*{{{*/
{
    // see comment above.

    unsigned int mi, ma;
    std::tie(mi, ma) = get_minmax_delta(delta);

    return expected_pi_length(d, len) + ma - mi;
}/*}}}*/

static inline unsigned int expected_pi_length_lowerbound(bw_dimensions & d, unsigned int len)/*{{{*/
{
    /* generically we expect that len*m % (m+n) columns have length
     * 1+\lfloor(len*m/(m+n))\rfloor, and the others have length one more.
     * For one column to have a length less than \lfloor(len*m/(m+n))\rfloor,
     * it takes probability 2^-(m*l) using the notations above. Therefore
     * we can simply count 2^(64-m*l) accidental zero cancellations at
     * most below the bound.
     * In particular, it is sufficient to derive from the code above!
     */
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    abdst_field ab MAYBE_UNUSED = d.ab;
    unsigned int res = 1 + (len * m) / b;
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    unsigned int l;
    if (mpz_cmp_ui(p, 1024) >= 0) {
        l = mpz_sizeinbase(p, 2);
        l *= abfield_degree(ab);    /* roughly log_2(#K) */
    } else {
        mpz_pow_ui(p, p, abfield_degree(ab));
        l = mpz_sizeinbase(p, 2);
    }
    mpz_clear(p);
    unsigned int safety = iceildiv(64, m * l);
    return safety < res ? res - safety : 0;
}/*}}}*/


/* This destructively cancels the first len coefficients of E, and
 * computes the appropriate matrix pi which achieves this. The
 * elimination is done in accordance with the nominal degrees found in
 * delta.
 *
 * The result is expected to have degree ceil(len*m/b) coefficients, so
 * that E*pi is divisible by X^len.
 */

/* TODO: adapt for GF(2) */
int
bw_lingen_basecase_raw(bmstatus & bm, matpoly & pi, matpoly const & E, std::vector<unsigned int> & delta) /*{{{*/
{
    int generator_found = 0;

    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    abdst_field ab = d.ab;
    ASSERT(E.m == m);
    ASSERT(E.n == b);

    ASSERT(pi.m == 0);
    ASSERT(pi.n == 0);
    ASSERT(pi.alloc == 0);

    /* Allocate something large enough for the result. This will be
     * soon freed anyway. Set it to identity. */
    unsigned int mi, ma;
    std::tie(mi, ma) = get_minmax_delta(delta);

    unsigned int pi_room_base = expected_pi_length(d, delta, E.size);

    pi = matpoly(ab, b, b, pi_room_base);
    pi.size = pi_room_base;

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
        pi_lengths[i] += delta[i] - mi;
    }

    /* Keep a list of columns which have been used as pivots at the
     * previous iteration */
    std::vector<unsigned int> pivots(m, 0);
    std::vector<int> is_pivot(b, 0);

    matpoly e(ab, m, b, 1);
    e.size = 1;

    matpoly T(ab, b, b, 1);
    int (*ctable)[2] = new int[b][2];

    for (unsigned int t = 0; t < E.size ; t++, bm.t++) {

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
            abelt_ur tmp_ur;
            abelt_ur_init(ab, &tmp_ur);

            abelt_ur e_ur;
            abelt_ur_init(ab, &e_ur);
            
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for(unsigned int jl = 0 ; jl < nj ; ++jl) {
                for(unsigned int i = 0 ; i < m ; ++i) {
                    unsigned int j = todo[jl];
                    unsigned int lj = MIN(pi_real_lengths[j], t + 1);
                    abelt_ur_set_zero(ab, e_ur);
                    for(unsigned int k = 0 ; k < b ; ++k) {
                        for(unsigned int s = 0 ; s < lj ; s++) {
                            abmul_ur(ab, tmp_ur,
                                    E.coeff(i, k, t - s),
                                    pi.coeff(k, j, s));
                            abelt_ur_add(ab, e_ur, e_ur, tmp_ur);
                        }
                    }
                    abreduce(ab, e.coeff(i, j, 0), e_ur);
                }
            }

            abelt_ur_clear(ab, &tmp_ur);
            abelt_ur_clear(ab, &e_ur);
        }
        /* }}} */

        /* {{{ check for cancellations */

        unsigned int newluck = 0;
        for(unsigned int jl = 0 ; jl < todo.size() ; ++jl) {
            unsigned int j = todo[jl];
            unsigned int nz = 0;
            for(unsigned int i = 0 ; i < m ; i++) {
                nz += abcmp_ui(ab, e.coeff(i, j, 0), 0) == 0;
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
            ctable[j][0] = delta[j];
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
                if (abcmp_ui(ab, e.coeff(u, j, 0), 0) != 0)
                    break;
            }
            if (u == m) continue;
            assert(r < m);
            /* }}} */
            pivots[r++] = j;
            is_pivot[j] = 1;
            pivot_columns.push_back(j);
            /* {{{ Cancel this coeff in all other columns. */
            abelt inv;
            abinit(ab, &inv);
            int rc = abinv(ab, inv, e.coeff(u, j, 0));
            if (!rc) {
                fprintf(stderr, "Error, found a factor of the modulus: ");
                abfprint(ab, stderr, inv);
                fprintf(stderr, "\n");
                exit(EXIT_FAILURE);
            }
            abneg(ab, inv, inv);
            for (unsigned int kl = jl + 1; kl < b ; kl++) {
                unsigned int k = ctable[kl][1];
                if (abcmp_ui(ab, e.coeff(u, k, 0), 0) == 0)
                    continue;
                // add lambda = e[u,k]*-e[u,j]^-1 times col j to col k.
                abelt lambda;
                abinit(ab, &lambda);
                abmul(ab, lambda, inv, e.coeff(u, k, 0));
                assert(delta[j] <= delta[k]);
                /* {{{ Apply on both e and pi */
                abelt tmp;
                abinit(ab, &tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    /* TODO: Would be better if mpfq had an addmul */
                    abmul(ab, tmp, lambda, e.coeff(i, j, 0));
                    abadd(ab,
                            e.coeff(i, k, 0),
                            e.coeff(i, k, 0),
                            tmp);
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
                abset(ab, T.coeff(j, k, 0), lambda);
                abclear(ab, &tmp);
                /* }}} */
                abclear(ab, &lambda);
            }
            abclear(ab, &inv); /* }}} */
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
            abelt_ur tmp_pi;
            abelt_ur tmp;
            abelt_ur_init(ab, &tmp);
            abelt_ur_init(ab, &tmp_pi);

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
                        absrc_elt Tkj = T.coeff(k, j, 0);
                        ASSERT_ALWAYS(abcmp_ui(ab, Tkj, k==j) == 0);
                    }
                    for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                        unsigned int k = pivot_columns[kl];
                        absrc_elt Tkj = T.coeff(k, j, 0);
                        if (abcmp_ui(ab, Tkj, 0) == 0) continue;
                        ASSERT_ALWAYS(pi_lengths[k] <= pi_lengths[j]);
                        pi_real_lengths[j] = std::max(pi_real_lengths[k], pi_real_lengths[j]);
                    }
                }
#endif

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
                for(unsigned int i = 0 ; i < b ; i++) {
                    for(unsigned int s = 0 ; s < pi_real_lengths[j] ; s++) {
                        abdst_elt piijs = pi.coeff(i, j, s);

                        abelt_ur_set_elt(ab, tmp_pi, piijs);

                        for(unsigned int kl = 0 ; kl < MIN(m,jl) ; kl++) {
                            unsigned int k = pivot_columns[kl];
                            /* TODO: if column k was already a pivot on previous
                             * turn (which could happen, depending on m and n),
                             * then the corresponding entry is probably zero
                             * (exact condition needs to be written more
                             * accurately).
                             */

                            absrc_elt Tkj = T.coeff(k, j, 0);
                            if (abcmp_ui(ab, Tkj, 0) == 0) continue;
                            /* pi[i,k] has length pi_lengths[k]. Multiply
                             * that by T[k,j], which is a constant. Add
                             * to the unreduced thing. We don't have an
                             * mpfq api call for that operation.
                             */
                            absrc_elt piiks = pi.coeff(i, k, s);
                            abmul_ur(ab, tmp, piiks, Tkj);
                            abelt_ur_add(ab, tmp_pi, tmp_pi, tmp);
                        }
                        abreduce(ab, piijs, tmp_pi);
                    }
                }
            }
            abelt_ur_clear(ab, &tmp);
            abelt_ur_clear(ab, &tmp_pi);
        }
        /* }}} */

        ASSERT_ALWAYS(r == m);

        /* {{{ Now for all pivots, multiply column in pi by x */
        for (unsigned int j = 0; j < b ; j++) {
            if (!is_pivot[j]) continue;
            if (pi_real_lengths[j] >= pi.alloc) {
                if (!generator_found) {
                    pi.realloc(pi.alloc + MAX(pi.alloc / (m+n), 1));
                    printf("t=%u, expanding allocation for pi (now %zu%%) ; lengths: ",
                            bm.t,
                            100 * pi.alloc / pi_room_base);
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
                    delta[j]++;
                    continue;
                }
            }
            pi.multiply_column_by_x(j, pi_real_lengths[j]);
            pi_real_lengths[j]++;
            pi_lengths[j]++;
            delta[j]++;
        }
        /* }}} */
    }

    pi.size = 0;
    for(unsigned int j = 0; j < b; j++) {
        if (pi_real_lengths[j] > pi.size)
            pi.size = pi_real_lengths[j];
    }
    /* Given the structure of the computation, there's no reason for the
     * initial estimate to go wrong.
     */
    ASSERT_ALWAYS(pi.size <= pi.alloc);
    for(unsigned int j = 0; j < b; j++) {
        for(unsigned int k = pi_real_lengths[j] ; k < pi.size ; k++) {
            for(unsigned int i = 0 ; i < b ; i++) {
                ASSERT_ALWAYS(abis_zero(ab, pi.coeff(i, j, k)));
            }
        }
    }
    pi.size = MIN(pi.size, pi.alloc);

    return generator_found;
}
int
bw_lingen_basecase(bmstatus & bm, matpoly & pi, matpoly & E, std::vector<unsigned int> & delta)
{
    lingen_call_companion const & C = bm.companion(bm.depth(), E.size);
    tree_stats::sentinel dummy(bm.stats, __func__, E.size, C.total_ncalls, true);
    bm.stats.plan_smallstep("basecase", C.ttb);
    bm.stats.begin_smallstep("basecase");
    int done = bw_lingen_basecase_raw(bm, pi, E, delta);
    bm.stats.end_smallstep();
    E = matpoly();
    return done;
}
/*}}}*/

/*}}}*/

/* TODO: adapt for GF(2) */
double avg_matsize(abdst_field ab, unsigned int m, unsigned int n, int ascii)/*{{{*/
{
    if (!ascii) {
        /* Easy case first. If we have binary input, then we know a priori
         * that the input data must have size a multiple of the element size.
         */
        size_t elemsize = abvec_elt_stride(ab, 1);
        size_t matsize = elemsize * m * n;
        return matsize;
    }

    /* Ascii is more complicated. We're necessarily fragile here.
     * However, assuming that each coefficient comes with only one space,
     * and each matrix with an extra space (this is how the GPU program
     * prints data -- not that this ends up having a considerable impact
     * anyway...), we can guess the number of bytes per matrix. */

    /* Formula for the average number of digits of an integer mod p,
     * written in base b:
     *
     * (k-((b^k-1)/(b-1)-1)/p)  with b = Ceiling(Log(p)/Log(b)).
     */
    double avg;
    mpz_t p, a;
    mpz_init(p);
    mpz_init(a);
    abfield_characteristic(ab, p);
    unsigned long k = ceil(log(mpz_get_d(p))/log(10));
    unsigned long b = 10;
    mpz_ui_pow_ui(a, b, k);
    mpz_sub_ui(a, a, 1);
    mpz_fdiv_q_ui(a, a, b-1);
    avg = k - mpz_get_d(a) / mpz_get_d(p);
    mpz_clear(p);
    mpz_clear(a);
    // printf("Expect roughly %.2f decimal digits for integers mod p.\n", avg);
    double matsize = (avg + 1) * m * n + 1;
    // printf("Expect roughly %.2f bytes for each sequence matrix.\n", matsize);
    return matsize;
}/*}}}*/

/* {{{ I/O helpers */
/* TODO: adapt for GF(2) */
/* {{{ matpoly_write
 * writes some of the matpoly data to f, either in ascii or binary
 * format. This can be used to write only part of the data (degrees
 * [k0..k1[). Returns the number of coefficients (i.e., matrices, so at
 * most k1-k0) successfully written, or
 * -1 on error (e.g. when some matrix was only partially written).
 */
int matpoly_write(abdst_field ab, std::ostream& os, matpoly const & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    unsigned int m = transpose ? M.n : M.m;
    unsigned int n = transpose ? M.m : M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.size && k1 <= M.size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                absrc_elt x;
                x = transpose ? M.coeff(j, i, k)
                              : M.coeff(i, j, k);
                if (ascii) {
                    if (j) err = !(os << " ");
                    if (!err) err = !(abcxx_out(ab, os, x));
                } else {
                    err = !(os.write((const char *) x, (size_t) abvec_elt_stride(ab, 1)));
                }
                if (!err) matnb++;
            }
            if (!err && ascii) err = !(os << "\n");
        }
        if (ascii) err = err || !(os << "\n");
        if (err) {
            abclear(ab, &tmp);
            return (matnb == 0) ? (int) (k - k0) : -1;
        }
    }
    abclear(ab, &tmp);
    return k1 - k0;
}

/* fw must be an array of FILE* pointers of exactly the same size as the
 * matrix to be written.
 */
template<typename Ostream>
int matpoly_write_split(abdst_field ab, std::vector<Ostream> & fw, matpoly const & M, unsigned int k0, unsigned int k1, int ascii)
{
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.size && k1 <= M.size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < M.m ; i++) {
            for(unsigned int j = 0 ; !err && j < M.n ; j++) {
                std::ostream& os = fw[i*M.n+j];
                absrc_elt x = M.coeff(i, j, k);
                if (ascii) {
                    err = !(abcxx_out(ab, os, x));
                    if (!err) err = !(os << "\n");
                } else {
                    err = !(os.write((const char *) x, (size_t) abvec_elt_stride(ab, 1)));
                }
                if (!err) matnb++;
            }
        }
        if (err) {
            abclear(ab, &tmp);
            return (matnb == 0) ? (int) (k - k0) : -1;
        }
    }
    abclear(ab, &tmp);
    return k1 - k0;
}
/* }}} */

/* {{{ matpoly_read
 * reads some of the matpoly data from f, either in ascii or binary
 * format. This can be used to parse only part of the data (degrees
 * [k0..k1[, k1 being an upper bound). Returns the number of coefficients
 * (i.e., matrices, so at most k1-k0) successfully read, or
 * -1 on error (e.g. when some matrix was only partially read).
 *
 * Note that the matrix must *not* be in pre-init state. It must have
 * been already allocated.
 */

int matpoly_read(abdst_field ab, FILE * f, matpoly & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    ASSERT_ALWAYS(!M.check_pre_init());
    unsigned int m = transpose ? M.n : M.m;
    unsigned int n = transpose ? M.m : M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.size && k1 <= M.size));
    abelt tmp;
    abinit(ab, &tmp);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                abdst_elt x;
                x = transpose ? M.coeff(j, i, k)
                              : M.coeff(i, j, k);
                if (ascii) {
                    err = abfscan(ab, f, x) == 0;
                } else {
                    err = fread(x, abvec_elt_stride(ab, 1), 1, f) < 1;
                }
                if (!err) matnb++;
            }
        }
        if (err) return (matnb == 0) ? (int) (k - k0) : -1;
    }
    abclear(ab, &tmp);
    return k1 - k0;
}/* }}} */

/* }}} */

/*{{{ Checkpoints */

/* There's much copy-paste here */

struct cp_info {
    bmstatus & bm;
    int level;
    unsigned int t0;
    unsigned int t1;
    int mpi;
    int rank;
    char * auxfile;
    char * sdatafile;
    char * gdatafile;
    const char * datafile;
    FILE * aux;
    FILE * data;
    cp_info(bmstatus & bm, unsigned int t0, unsigned int t1, int mpi);
    ~cp_info();
    int save_aux_file(size_t pi_size, std::vector<unsigned int> const & delta, int done);
    int load_aux_file(size_t * p_pi_size, std::vector<unsigned int> & delta, int * p_done);
    int load_data_file(matpoly & pi, size_t pi_size);
    int save_data_file(matpoly const & pi, size_t pi_size);
};

cp_info::cp_info(bmstatus & bm, unsigned int t0, unsigned int t1, int mpi)
    : bm(bm), t0(t0), t1(t1), mpi(mpi)
{
    if (mpi)
        MPI_Comm_rank(bm.com[0], &(rank));
    else
        rank = 0;
    int rc;
    level = bm.depth();
    rc = asprintf(&auxfile, "%s/pi.%d.%u.%u.aux",
            checkpoint_directory, level, t0, t1);
    ASSERT_ALWAYS(rc >= 0);
    rc = asprintf(&gdatafile, "%s/pi.%d.%u.%u.single.data",
            checkpoint_directory, level, t0, t1);
    ASSERT_ALWAYS(rc >= 0);
    rc = asprintf(&sdatafile, "%s/pi.%d.%u.%u.%d.data",
            checkpoint_directory, level, t0, t1, rank);
    ASSERT_ALWAYS(rc >= 0);
    datafile = mpi ? sdatafile : gdatafile;
}

cp_info::~cp_info()
{
    free(sdatafile);
    free(gdatafile);
    free(auxfile);
}

int cp_info::save_aux_file(size_t pi_size, std::vector<unsigned int> const & delta, int done)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    if (rank) return 1;
    FILE * aux = fopen(auxfile, "w");
    int rc;
    if (aux == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", auxfile);
        return 0;
    }
    rc = fprintf(aux, "%zu\n", pi_size);
    if (rc <= 0) goto cp_info_save_aux_file_bailout;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fprintf(aux, "%s%u", i?" ":"", delta[i]);
        if (rc <= 0) goto cp_info_save_aux_file_bailout;
    }
    rc = fprintf(aux, "\n");
    if (rc <= 0) goto cp_info_save_aux_file_bailout;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fprintf(aux, "%s%d", i?" ":"", bm.lucky[i]);
        if (rc <= 0) goto cp_info_save_aux_file_bailout;
    }
    rc = fprintf(aux, "\n");
    if (rc <= 0) goto cp_info_save_aux_file_bailout;
    rc = fprintf(aux, "%d\n", done);
    if (rc <= 0) goto cp_info_save_aux_file_bailout;
    rc = fclose(aux);
    if (rc == 0) return 1;
cp_info_save_aux_file_bailout:
    fclose(aux);
    unlink(auxfile);
    return 0;
}/*}}}*/

int cp_info::load_aux_file(size_t * p_pi_size, std::vector<unsigned int> & delta, int * p_done)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    if (rank) return 1;
    FILE * aux = fopen(auxfile, "r");
    int rc;
    if (aux == NULL) {
        // fprintf(stderr, "Warning: cannot open %s\n", auxfile);
        return 0;
    }
    rc = fscanf(aux, "%zu", p_pi_size);
    if (rc != 1) { fclose(aux); return 0; }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fscanf(aux, "%u", &(delta[i]));
        if (rc != 1) { fclose(aux); return 0; }
    }
    for(unsigned int i = 0 ; i < m + n ; i++) {
        rc = fscanf(aux, "%d", &(bm.lucky[i]));
        if (rc != 1) { fclose(aux); return 0; }
    }
    rc = fscanf(aux, "%d", p_done);
    if (rc != 1) { fclose(aux); return 0; }
    rc = fclose(aux);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
int cp_info::load_data_file(matpoly & pi, size_t pi_size)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    FILE * data = fopen(datafile, "rb");
    int rc;
    if (data == NULL) {
        fprintf(stderr, "Warning: cannot open %s\n", datafile);
        return 0;
    }
    pi = matpoly(ab, m+n, m+n, pi_size);
    pi.size = pi_size;
    rc = matpoly_read(ab, data, pi, 0, pi.size, 0, 0);
    if (rc != (int) pi.size) { fclose(data); return 0; }
    rc = fclose(data);
    return rc == 0;
}/*}}}*/

/* TODO: adapt for GF(2) */
/* I think we always have pi_size == pi.size, the only questionable
 * situation is when we're saving part of a big matrix */
int cp_info::save_data_file(matpoly const & pi, size_t pi_size)/*{{{*/
{
    abdst_field ab = bm.d.ab;
    std::ofstream data(datafile, std::ios_base::out | std::ios_base::binary);
    int rc;
    if (!data) {
        fprintf(stderr, "Warning: cannot open %s\n", datafile);
        unlink(auxfile);
        return 0;
    }
    rc = matpoly_write(ab, data, pi, 0, pi_size, 0, 0);
    if (rc != (int) pi.size) goto cp_info_save_data_file_bailout;
    data.close();
    if (!data)
        return 1;
cp_info_save_data_file_bailout:
    unlink(datafile);
    unlink(auxfile);
    return 0;
}/*}}}*/

int load_checkpoint_file(bmstatus & bm, matpoly & pi, unsigned int t0, unsigned int t1, std::vector<unsigned int> & delta, int * p_done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;

    cp_info cp(bm, t0, t1, 0);

    ASSERT_ALWAYS(pi.check_pre_init());
    size_t pi_size;
    /* Don't output a message just now, since after all it's not
     * noteworthy if the checkpoint file does not exist. */
    int ok = cp.load_aux_file(&pi_size, delta, p_done);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s", cp.datafile);
        ok = cp.load_data_file(pi, pi_size);
        logline_end(&bm.t_cp_io,"");
        if (!ok)
            fprintf(stderr, "Warning: I/O error while reading %s\n", cp.datafile);
    }
    if (ok) bm.t = t1;
    return ok;
}/*}}}*/

int save_checkpoint_file(bmstatus & bm, matpoly & pi, unsigned int t0, unsigned int t1, std::vector<unsigned int> const & delta, int done)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    cp_info cp(bm, t0, t1, 0);
    logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s%s",
            cp.datafile,
            cp.mpi ? " (MPI, scattered)" : "");
    int ok = cp.save_aux_file(pi.size, delta, done);
    if (ok) ok = cp.save_data_file(pi, pi.size);
    logline_end(&bm.t_cp_io,"");
    if (!ok && !cp.rank)
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
    return ok;
}/*}}}*/

#ifdef ENABLE_MPI_LINGEN
int load_mpi_checkpoint_file_scattered(bmstatus & bm, bigmatpoly & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> & delta, int * p_done)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, t0, t1, 1);
    ASSERT_ALWAYS(xpi.check_pre_init());
    size_t pi_size;
    int ok = cp.load_aux_file(&pi_size, delta, p_done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&pi_size, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(p_done, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s (MPI, scattered)",
                cp.datafile);
        do {
            FILE * data = fopen(cp.datafile, "rb");
            int rc;
            ok = data != NULL;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                if (data) free(data);
                break;
            }
            xpi.finish_init(ab, m+n, m+n, pi_size);
            xpi.size = pi_size;
            rc = matpoly_read(ab, data, xpi.my_cell(), 0, xpi.size, 0, 0);
            ok = ok && rc == (int) xpi.size;
            rc = fclose(data);
            ok = ok && (rc == 0);
        } while (0);
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        logline_end(&bm.t_cp_io,"");
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp.datafile);
        }
    }
    if (ok) bm.t = t1;
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_scattered(bmstatus & bm, bigmatpoly const & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> const & delta, int done)/*{{{*/
{
    /* corresponding t is bm.t - E.size ! */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    cp_info cp(bm, t0, t1, 1);
    int ok = cp.save_aux_file(xpi.size, delta, done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (!ok && !rank) unlink(cp.auxfile);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s (MPI, scattered)",
                cp.datafile);
        ok = cp.save_data_file(xpi.my_cell(), xpi.size);
        logline_end(&bm.t_cp_io,"");
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
        if (!ok) {
            if (cp.datafile) unlink(cp.datafile);
            if (!rank) unlink(cp.auxfile);
        }
    }
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
    }
    return ok;
}/*}}}*/

int load_mpi_checkpoint_file_gathered(bmstatus & bm, bigmatpoly & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> & delta, int * p_done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    size_t pi_size;
    int ok = cp.load_aux_file(&pi_size, delta, p_done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(&pi_size, 1, MPI_MY_SIZE_T, 0, bm.com[0]);
    MPI_Bcast(&delta[0], m + n, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&bm.lucky[0], m + n, MPI_INT, 0, bm.com[0]);
    MPI_Bcast(p_done, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        logline_begin(stdout, SIZE_MAX, "Reading checkpoint %s (MPI, gathered)",
                cp.datafile);
        do {
            FILE * data = NULL;
            if (!rank) ok = (data = fopen(cp.datafile, "rb")) != NULL;
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                if (data) free(data);
                break;
            }

            xpi.finish_init(ab, m+n, m+n, pi_size);
            xpi.set_size(pi_size);

            double avg = avg_matsize(ab, m + n, m + n, 0);
            unsigned int B = iceildiv(io_block_size, avg);

            /* This is only temp storage ! */
            matpoly pi(ab, m + n, m + n, B);
            pi.size = B;

            for(unsigned int k = 0 ; ok && k < xpi.size ; k += B) {
                unsigned int nc = MIN(B, xpi.size - k);
                if (!rank)
                    ok = matpoly_read(ab, data, pi, 0, nc, 0, 0) == (int) nc;
                MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
                xpi.scatter_mat_partial(pi, k, nc);
            }

            if (!rank) {
                int rc = fclose(data);
                ok = ok && (rc == 0);
            }
        } while (0);
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
        logline_end(&bm.t_cp_io,"");
        if (!ok && !rank) {
            fprintf(stderr, "Warning: I/O error while reading %s\n",
                    cp.datafile);
        }
    }
    if (ok) bm.t = t1;
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file_gathered(bmstatus & bm, bigmatpoly const & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> const & delta, int done)/*{{{*/
{
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    cp_info cp(bm, t0, t1, 1);
    cp.datafile = cp.gdatafile;
    logline_begin(stdout, SIZE_MAX, "Saving checkpoint %s (MPI, gathered)",
            cp.datafile);
    int ok = cp.save_aux_file(xpi.size, delta, done);
    MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
    if (ok) {
        do {
            std::ofstream data;
            if (!rank) {
                data.open(cp.datafile, std::ios_base::out | std::ios_base::binary);
                ok = (bool) data;
            }
            MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            if (!ok) {
                if (!rank)
                    fprintf(stderr, "Warning: cannot open %s\n", cp.datafile);
                break;
            }

            double avg = avg_matsize(ab, m + n, m + n, 0);
            unsigned int B = iceildiv(io_block_size, avg);

            /* This is only temp storage ! */
            matpoly pi(ab, m + n, m + n, B);
            pi.size = B;

            for(unsigned int k = 0 ; ok && k < xpi.size ; k += B) {
                unsigned int nc = MIN(B, xpi.size - k);
                xpi.gather_mat_partial(pi, k, nc);
                if (!rank)
                    ok = matpoly_write(ab, data, pi, 0, nc, 0, 0) == (int) nc;
                MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
            }

            if (!rank) {
                data.close();
                ok = ok && (bool) data;
            }
        } while (0);
        MPI_Bcast(&ok, 1, MPI_INT, 0, bm.com[0]);
        if (!ok && !rank) {
            if (cp.datafile) unlink(cp.datafile);
            unlink(cp.auxfile);
        }
    }
    logline_end(&bm.t_cp_io,"");
    if (!ok && !rank) {
        fprintf(stderr, "Warning: I/O error while saving %s\n", cp.datafile);
    }
    return ok;
}/*}}}*/

int load_mpi_checkpoint_file(bmstatus & bm, bigmatpoly & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> & delta, int * p_done)/*{{{*/
{
    /* read scattered checkpoint with higher priority if available,
     * because we like distributed I/O. Otherwise, read gathered
     * checkpoint if we could find one.
     */
    if (!checkpoint_directory) return 0;
    if ((t1 - t0) < checkpoint_threshold) return 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    cp_info cp(bm, t0, t1, 1);
    int ok = 0;
    int aux_ok = rank || access(cp.auxfile, R_OK) == 0;
    int sdata_ok = access(cp.sdatafile, R_OK) == 0;
    int scattered_ok = aux_ok && sdata_ok;
    MPI_Allreduce(MPI_IN_PLACE, &scattered_ok, 1, MPI_INT, MPI_MIN, bm.com[0]);
    if (scattered_ok) {
        ok = load_mpi_checkpoint_file_scattered(bm, xpi, t0, t1, delta, p_done);
        if (ok) return ok;
    }
    int gdata_ok = rank || access(cp.gdatafile, R_OK) == 0;
    int gathered_ok = aux_ok && gdata_ok;
    MPI_Bcast(&gathered_ok, 1, MPI_INT, 0, bm.com[0]);
    if (gathered_ok) {
        ok = load_mpi_checkpoint_file_gathered(bm, xpi, t0, t1, delta, p_done);
    }
    return ok;
}/*}}}*/

int save_mpi_checkpoint_file(bmstatus & bm, bigmatpoly const & xpi, unsigned int t0, unsigned int t1, std::vector<unsigned int> const & delta, int done)/*{{{*/
{
    if (save_gathered_checkpoints) {
        return save_mpi_checkpoint_file_gathered(bm, xpi, t0, t1, delta, done);
    } else {
        return save_mpi_checkpoint_file_scattered(bm, xpi, t0, t1, delta, done);
    }
}/*}}}*/
#endif  /* ENABLE_MPI_LINGEN */

/*}}}*/

/**********************************************************************/

/*{{{ Main entry points and recursive algorithm (with and without MPI) */

/* Forward declaration, it's used by the recursive version */
int bw_lingen_single(bmstatus & bm, matpoly & pi, matpoly & E, std::vector<unsigned int> & delta);

#ifdef ENABLE_MPI_LINGEN
int bw_biglingen_collective(bmstatus & bm, bigmatpoly & pi, bigmatpoly & E, std::vector<unsigned int> & delta);
#endif

std::string sha1sum(matpoly const & X)
{
    sha1_checksumming_stream S;
    S.write((const char *) X.x, X.m*(X.n)*X.size*sizeof(mp_limb_t));
    char checksum[41];
    S.checksum(checksum);
    return std::string(checksum);
}

int
bw_lingen_recursive(bmstatus & bm, matpoly & pi, matpoly & E, std::vector<unsigned int> & delta) /*{{{*/
{
    int depth = bm.depth();
    size_t z = E.size;

    lingen_call_companion const & C = bm.companion(depth, z);

    tree_stats::sentinel dummy(bm.stats, __func__, E.size, C.total_ncalls);

    {
        bm.stats.begin_plan_smallstep("MP", C.mp.tt);
        /*
        bm.stats.plan_smallstep("dft_A", C.mp.t_dft_A);
        bm.stats.plan_smallstep("dft_A_comm", C.mp.t_dft_A_comm);
        bm.stats.plan_smallstep("dft_B", C.mp.t_dft_B);
        bm.stats.plan_smallstep("dft_B_comm", C.mp.t_dft_B_comm);
        bm.stats.plan_smallstep("addmul", C.mp.t_conv);
        bm.stats.plan_smallstep("ift_C", C.mp.t_ift_C);
        */
        bm.stats.end_plan_smallstep();

        bm.stats.begin_plan_smallstep("MUL", C.mul.tt);
        /*
        bm.stats.plan_smallstep("dft_A", C.mul.t_dft_A);
        bm.stats.plan_smallstep("dft_A_comm", C.mul.t_dft_A_comm);
        bm.stats.plan_smallstep("dft_B", C.mul.t_dft_B);
        bm.stats.plan_smallstep("dft_B_comm", C.mul.t_dft_B_comm);
        bm.stats.plan_smallstep("addmul", C.mul.t_conv);
        bm.stats.plan_smallstep("ift_C", C.mul.t_ift_C);
        */
        bm.stats.end_plan_smallstep();
    }

    bw_dimensions & d = bm.d;
    int done;

    /* we have to start with something large enough to get all
     * coefficients of E_right correct */
    size_t half = E.size - (E.size / 2);
    unsigned int pi_expect = expected_pi_length(d, delta, E.size);
    unsigned int pi_expect_lowerbound = expected_pi_length_lowerbound(d, E.size);
    unsigned int pi_left_expect = expected_pi_length(d, delta, half);
    unsigned int pi_left_expect_lowerbound = expected_pi_length_lowerbound(d, half);
    unsigned int pi_left_expect_used_for_shift = MIN(pi_left_expect, half + 1);

    /* declare an lazy-alloc all matrices */
    matpoly E_left;
    matpoly pi_left;
    matpoly pi_right;
    matpoly E_right;

    E_left = E.truncate_and_rshift(half, half + 1 - pi_left_expect_used_for_shift);

    // this (now) consumes E_left entirely.
    done = bw_lingen_single(bm, pi_left, E_left, delta);

    ASSERT_ALWAYS(pi_left.size);

    if (done) {
        pi = std::move(pi_left);
        return done;
    }

    bm.stats.begin_smallstep("MP");
    ASSERT_ALWAYS(pi_left.size <= pi_left_expect);
    ASSERT_ALWAYS(done || pi_left.size >= pi_left_expect_lowerbound);

    /* XXX I don't understand why I need to do this. It seems to me that
     * MP(XA, B) and MP(A, B) should be identical whenever deg A > deg B.
     */
    ASSERT_ALWAYS(pi_left_expect_used_for_shift >= pi_left.size);
    if (pi_left_expect_used_for_shift != pi_left.size) {
        E.rshift(E, pi_left_expect_used_for_shift - pi_left.size);
        /* Don't shrink_to_fit at this point, because we've only made a
         * minor adjustment. */
    }

    logline_begin(stdout, z, "t=%u %*sMP(%zu, %zu) -> %zu",
            bm.t, depth,"",
            E.size, pi_left.size, E.size - pi_left.size + 1);

    {
        E_right = matpoly(d.ab, d.m, d.m+d.n, E.size - pi_left.size + 1);
        matpoly_ft::memory_pool_guard dummy(C.mp.ram);
        matpoly_mp_caching(E_right, E, pi_left, &C.mp.S);
        E = matpoly();
    }

    logline_end(&(bm.t_mp), "");
    bm.stats.end_smallstep();

    unsigned int pi_right_expect = expected_pi_length(d, delta, E_right.size);
    unsigned int pi_right_expect_lowerbound = expected_pi_length_lowerbound(d, E_right.size);

    done = bw_lingen_single(bm, pi_right, E_right, delta);
    ASSERT_ALWAYS(pi_right.size <= pi_right_expect);
    ASSERT_ALWAYS(done || pi_right.size >= pi_right_expect_lowerbound);

    /* stack is now pi_left, pi_right */

    bm.stats.begin_smallstep("MUL");
    logline_begin(stdout, z, "t=%u %*sMUL(%zu, %zu) -> %zu",
            bm.t, depth, "",
            pi_left.size, pi_right.size, pi_left.size + pi_right.size - 1);

    {
        pi = matpoly(d.ab, d.m+d.n, d.m+d.n, pi_left.size + pi_right.size - 1);
        matpoly_ft::memory_pool_guard dummy(C.mul.ram);
        matpoly_mul_caching(pi, pi_left, pi_right, &C.mul.S);
    }

    /* Note that the leading coefficients of pi_left and pi_right are not
     * necessarily full-rank, so that we have to fix potential zeros. If
     * we don't, the degree of pi artificially grows with the recursive
     * level.
     */
#if 1
    /* In fact, it's not entirely impossible that pi grows more than
     * what we had expected on entry, e.g. if we have one early
     * generator. So we can't just do this. Most of the time it will
     * work, but we can't claim that it will always work.
     *
     * One possible sign is when the entry deltas are somewhat even, and
     * the result deltas are unbalanced.
     */
    for(; pi.size > pi_expect ; pi.size--) {
        /* These coefficients really must be zero */
        ASSERT_ALWAYS(pi.coeff_is_zero(pi.size - 1));
    }
    ASSERT_ALWAYS(pi.size <= pi_expect);
#endif
    /* Now below pi_expect, it's not impossible to have a few
     * cancellations as well.
     */
    for(; pi.size ; pi.size--) {
        if (!pi.coeff_is_zero(pi.size - 1)) break;
    }
    ASSERT_ALWAYS(done || pi.size >= pi_expect_lowerbound);


    logline_end(&bm.t_mul, "");
    bm.stats.end_smallstep();

    return done;
}/*}}}*/

int bw_lingen_single(bmstatus & bm, matpoly & pi, matpoly & E, std::vector<unsigned int> & delta) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(!rank);
    unsigned int t0 = bm.t;
    unsigned int t1 = bm.t + E.size;

    int done;

    lingen_call_companion const & C = bm.companion(bm.depth(), E.size);

    if (load_checkpoint_file(bm, pi, t0, t1, delta, &done))
        return done;

    // ASSERT_ALWAYS(E.size < bm.lingen_mpi_threshold);

    // fprintf(stderr, "Enter %s\n", __func__);
    if (!bm.recurse(E.size)) {
        tree_stats::transition_sentinel dummy(bm.stats, "recursive threshold", E.size, C.total_ncalls);
        bm.t_basecase -= seconds();
        done = bw_lingen_basecase(bm, pi, E, delta);
        bm.t_basecase += seconds();
    } else {
        done = bw_lingen_recursive(bm, pi, E, delta);
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    save_checkpoint_file(bm, pi, t0, t1, delta, done);

    return done;
}/*}}}*/

#ifdef ENABLE_MPI_LINGEN
int bw_biglingen_recursive(bmstatus & bm, bigmatpoly & pi, bigmatpoly & E, std::vector<unsigned int> & delta) /*{{{*/
{
    int depth = bm.depth();
    size_t z = E.size;

    lingen_call_companion const & C = bm.companion(depth, z);

    tree_stats::sentinel dummy(bm.stats, __func__, E.size, C.total_ncalls);

    {
        bm.stats.begin_plan_smallstep("MP", C.mp.tt);
        bm.stats.plan_smallstep("dft_A", C.mp.t_dft_A);
        bm.stats.begin_plan_smallstep("dft_A_comm", C.mp.t_dft_A_comm);
        bm.stats.plan_smallstep("export", C.mp.t_dft_A_comm);
        bm.stats.plan_smallstep("comm", C.mp.t_dft_A_comm);
        bm.stats.plan_smallstep("import", C.mp.t_dft_A_comm);
        bm.stats.end_plan_smallstep();
        bm.stats.plan_smallstep("dft_B", C.mp.t_dft_B);
        bm.stats.begin_plan_smallstep("dft_B_comm", C.mp.t_dft_B_comm);
        bm.stats.plan_smallstep("export", C.mp.t_dft_B_comm);
        bm.stats.plan_smallstep("comm", C.mp.t_dft_B_comm);
        bm.stats.plan_smallstep("import", C.mp.t_dft_B_comm);
        bm.stats.end_plan_smallstep();
        bm.stats.plan_smallstep("addmul", C.mp.t_conv);
        bm.stats.plan_smallstep("ift_C", C.mp.t_ift_C);
        bm.stats.end_plan_smallstep();

        bm.stats.begin_plan_smallstep("MUL", C.mul.tt);
        bm.stats.plan_smallstep("dft_A", C.mul.t_dft_A);
        bm.stats.begin_plan_smallstep("dft_A_comm", C.mul.t_dft_A_comm);
        bm.stats.plan_smallstep("export", C.mul.t_dft_A_comm);
        bm.stats.plan_smallstep("comm", C.mul.t_dft_A_comm);
        bm.stats.plan_smallstep("import", C.mul.t_dft_A_comm);
        bm.stats.end_plan_smallstep();
        bm.stats.plan_smallstep("dft_B", C.mul.t_dft_B);
        bm.stats.begin_plan_smallstep("dft_B_comm", C.mul.t_dft_B_comm);
        bm.stats.plan_smallstep("export", C.mul.t_dft_B_comm);
        bm.stats.plan_smallstep("comm", C.mul.t_dft_B_comm);
        bm.stats.plan_smallstep("import", C.mul.t_dft_B_comm);
        bm.stats.end_plan_smallstep();
        bm.stats.plan_smallstep("addmul", C.mul.t_conv);
        bm.stats.plan_smallstep("ift_C", C.mul.t_ift_C);
        bm.stats.end_plan_smallstep();
    }


    bw_dimensions & d = bm.d;
    int done;

    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    /* we have to start with something large enough to get
     * all coefficients of E_right correct */
    size_t half = E.size - (E.size / 2);
    unsigned int pi_expect = expected_pi_length(d, delta, E.size);
    unsigned int pi_expect_lowerbound = expected_pi_length_lowerbound(d, E.size);
    unsigned int pi_left_expect = expected_pi_length(d, delta, half);
    unsigned int pi_left_expect_lowerbound = expected_pi_length_lowerbound(d, half);
    unsigned int pi_left_expect_used_for_shift = MIN(pi_left_expect, half + 1);

    /* declare an lazy-alloc all matrices */
    bigmatpoly_model const& model(E);
    bigmatpoly E_left(model);
    bigmatpoly E_right(model);
    bigmatpoly pi_left(model);
    bigmatpoly pi_right(model);

    E_left = E.truncate_and_rshift(half, half + 1 - pi_left_expect_used_for_shift);

    done = bw_biglingen_collective(bm, pi_left, E_left, delta);

    ASSERT_ALWAYS(pi_left.size);
    E_left = bigmatpoly(model);

    if (done) {
        pi = std::move(pi_left);
        return done;
    }

    bm.stats.begin_smallstep("MP");
    ASSERT_ALWAYS(pi_left.size <= pi_left_expect);
    ASSERT_ALWAYS(done || pi_left.size >= pi_left_expect_lowerbound);

    /* XXX I don't understand why I need to do this. It seems to me that
     * MP(XA, B) and MP(A, B) should be identical whenever deg A > deg B.
     */
    ASSERT_ALWAYS(pi_left_expect_used_for_shift >= pi_left.size);
    if (pi_left_expect_used_for_shift != pi_left.size) {
        E.rshift(E, pi_left_expect_used_for_shift - pi_left.size);
        /* Don't shrink_to_fit at this point, because we've only made a
         * minor adjustment. */
    }

    logline_begin(stdout, z, "t=%u %*sMPI-MP(%zu, %zu) -> %zu",
            bm.t, depth, "",
            E.size, pi_left.size, E.size - pi_left.size + 1);

    {
        ASSERT_ALWAYS(pi_left.ab);
        ASSERT_ALWAYS(E.ab);
        matpoly_ft::memory_pool_guard dummy(C.mp.ram);
        /* XXX should we pre-alloc ? We do that in the non-mpi case, but
         * that seems to be useless verbosity */
        bigmatpoly_mp_caching(bm.stats, E_right, E, pi_left, &C.mp.S);
        E = bigmatpoly(model);
        ASSERT_ALWAYS(E_right.ab);
        MPI_Barrier(bm.com[0]);
    }

    logline_end(&bm.t_mp, "");
    bm.stats.end_smallstep();

    unsigned int pi_right_expect = expected_pi_length(d, delta, E_right.size);
    unsigned int pi_right_expect_lowerbound = expected_pi_length_lowerbound(d, E_right.size);

    done = bw_biglingen_collective(bm, pi_right, E_right, delta);
    ASSERT_ALWAYS(pi_right.size <= pi_right_expect);
    ASSERT_ALWAYS(done || pi_right.size >= pi_right_expect_lowerbound);
    
    E_right = bigmatpoly(model);

    bm.stats.begin_smallstep("MUL");
    logline_begin(stdout, z, "t=%u %*sMPI-MUL(%zu, %zu) -> %zu",
            bm.t, depth, "",
            pi_left.size, pi_right.size, pi_left.size + pi_right.size - 1);

    {
        ASSERT_ALWAYS(pi_left.ab);
        ASSERT_ALWAYS(pi_right.ab);
        matpoly_ft::memory_pool_guard dummy(C.mul.ram);
        /* XXX should we pre-alloc ? We do that in the non-mpi case, but
         * that seems to be useless verbosity */
        bigmatpoly_mul_caching(bm.stats, pi, pi_left, pi_right, &C.mul.S);
        ASSERT_ALWAYS(pi.ab);
        MPI_Barrier(bm.com[0]);
    }

    /* Note that the leading coefficients of pi_left and pi_right are not
     * necessarily full-rank, so that we have to fix potential zeros. If
     * we don't, the degree of pi artificially grows with the recursive
     * level.
     */
#if 1
    /* In fact, it's not entirely impossible that pi grows more than
     * what we had expected on entry, e.g. if we have one early
     * generator. So we can't just do this. Most of the time it will
     * work, but we can't claim that it will always work.
     *
     * One possible sign is when the entry deltas are somewhat even, and
     * the result deltas are unbalanced.
     */
    for(; pi.size > pi_expect ; pi.size--) {
        /* These coefficients really must be zero */
        ASSERT_ALWAYS(pi.coeff_is_zero(pi.size - 1));
    }
    ASSERT_ALWAYS(pi.size <= pi_expect);
#endif
    /* Now below pi_expect, it's not impossible to have a few
     * cancellations as well.
     */
    for(; pi.size ; pi.size--) {
        if (!pi.coeff_is_zero(pi.size - 1)) break;
    }
    ASSERT_ALWAYS(done || pi.size >= pi_expect_lowerbound);

    logline_end(&bm.t_mul, "");
    bm.stats.end_smallstep();

    return done;
}/*}}}*/

int bw_biglingen_collective(bmstatus & bm, bigmatpoly & pi, bigmatpoly & E, std::vector<unsigned int> & delta)/*{{{*/
{
    /* as for bw_lingen_single, we're tempted to say that we're just a
     * trampoline. In fact, it's not really satisfactory: we're really
     * doing stuff here. In a sense though, it's not *that much* of a
     * trouble, because the mpi threshold will be low enough that doing
     * our full job here is not too much of a problem.
     */
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    int done;
    int rank;
    int size;
    MPI_Comm_rank(bm.com[0], &rank);
    MPI_Comm_size(bm.com[0], &size);
    unsigned int t0 = bm.t;
    unsigned int t1 = bm.t + E.size;

    lingen_call_companion const & C = bm.companion(bm.depth(), E.size);
    bool go_mpi = C.go_mpi;
    // bool go_mpi = E.size >= bm.lingen_mpi_threshold;

    if (load_mpi_checkpoint_file(bm, pi, t0, t1, delta, &done))
        return done;

    // fprintf(stderr, "Enter %s\n", __func__);
    if (go_mpi) {
        done = bw_biglingen_recursive(bm, pi, E, delta);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        tree_stats::transition_sentinel dummy(bm.stats, "mpi threshold", E.size, C.total_ncalls);

        matpoly sE(ab, m, b, E.size);
        matpoly spi;

        double expect0 = bm.hints.tt_gather_per_unit * E.size;
        bm.stats.plan_smallstep("gather(L+R)", expect0);
        bm.stats.begin_smallstep("gather(L+R)");
        E.gather_mat(sE);
        bm.stats.end_smallstep();

        /* Only the master node does the local computation */
        if (!rank)
            done = bw_lingen_single(bm, spi, sE, delta);

        double expect1 = bm.hints.tt_scatter_per_unit * E.size;
        bm.stats.plan_smallstep("scatter(L+R)", expect1);
        bm.stats.begin_smallstep("scatter(L+R)");
        pi = bigmatpoly(ab, E.get_model(), b, b, 0);
        pi.scatter_mat(spi);
        MPI_Bcast(&done, 1, MPI_INT, 0, bm.com[0]);
        MPI_Bcast(&delta[0], b, MPI_UNSIGNED, 0, bm.com[0]);
        MPI_Bcast(&bm.lucky[0], b, MPI_UNSIGNED, 0, bm.com[0]);
        MPI_Bcast(&(bm.t), 1, MPI_UNSIGNED, 0, bm.com[0]);
        /* Don't forget to broadcast delta from root node to others ! */
        bm.stats.end_smallstep();
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    save_mpi_checkpoint_file(bm, pi, t0, t1, delta, done);

    MPI_Barrier(bm.com[0]);

    return done;
}/*}}}*/
#endif  /* ENABLE_MPI_LINGEN */

/*}}}*/

/**********************************************************************/

/**********************************************************************/
/* {{{ reading A and writing F ... */
struct bm_io {/*{{{*/
    bmstatus & bm;
    unsigned int t0 = 0;
    FILE ** fr = NULL; /* array of n files when split_input_file is supported
                   (which is not the case as of now),
                   or otherwise just a 1-element array */
    char * iobuf = NULL;
    const char * input_file = NULL;
    const char * output_file = NULL;
    int ascii = 0;
    /* This is only a rolling window ! */
    matpoly A, F;
    unsigned int (*fdesc)[2] = NULL;
    /* This k is the coefficient in A(X) div X of the next coefficient to
     * be read. This is thus the total number of coefficients of A(X) div
     * X which have been read so far.
     * In writing mode, k is the number of coefficients of F which have
     * been written so far.
     */
    unsigned int k = 0;
    unsigned int guessed_length = 0;

    bool leader() const {
        int rank;
        MPI_Comm_rank(bm.com[0], &rank);
        return rank == 0;
    }
    unsigned int set_write_behind_size(std::vector<unsigned int> const & delta);
    void zero1(unsigned int deg);
    int read1(unsigned int io_window);
    bm_io(bm_io const&)=delete;
    bm_io(bmstatus & bm, const char * input_file, const char * output_file, int ascii);
    ~bm_io();
    void begin_read();
    void end_read();
    void guess_length();
    void compute_initial_F() ;

    template<class Consumer, class Sink>
        void compute_final_F(Sink & S, Consumer& pi, std::vector<unsigned int> & delta);
    template<class Producer>
        void compute_E(Producer& E, unsigned int expected, unsigned int allocated);
    template<typename T, typename Sink>
        void output_flow(T & pi, std::vector<unsigned int> & delta);
};
/*}}}*/

/* The reading mode of bm_io is streaming, but with a look-back
 * functionality: we want to be able to access coefficients a few places
 * earlier, so we keep them in memory.
 *
 * The writing mode has a write-ahead feature. Coefficients of the result
 * are written at various times, some earlier than others. The time span
 * is the same as for the reading mode.
 */


/* We write the coefficients of the reversed polynomial \hat{F*\pi}, in
 * increasing degree order. Thus the first coefficients written
 * corresponds to high degree coefficients in \pi.
 * This is mostly for historical reasons, since in fact, we would prefer
 * having the coefficients in the reverse order.
 */

/* let mindelta and maxdelta be (as their name suggests) the minimum and
 * maximum over all deltas corresponding to solution columns.
 *
 * For 0<=i<n, coeff i,j,k of pi becomes coeff i,j,k' of the final f,
 * with k'=k-(maxdelta-delta[j]).
 *
 * For 0<=i<m, coeff n+i,j,k of pi becomes coeff c[i],j,k' of the final f,
 * with k'=k-(maxdelta-delta[j])-(t0-e[j]).
 *
 * Therefore the maximum write-behind distance is (maxdelta-mindelta)+t0.
 * We need one coeff more (because offset goes from 0 to maxoffset).
 */
unsigned int bm_io::set_write_behind_size(std::vector<unsigned int> const & delta)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta_on_solutions(bm, delta);
    unsigned int window = maxdelta - mindelta + t0 + 1;
    if (d.nrhs) {
        /* in sm-outside-matrix mode for DLP, we form the matrix F
         * slightly differently, as some *rows* are shifted out before
         * writing */
        window++;
    }
    if (leader()) {
        F.realloc(window);
        F.zero();
        F.size = window;
    }
    return window;
}/*}}}*/

/* {{{ bm_output_* classes: the link from the in-memory structures to the
 * filesystem.
 * Note that coefficients are written transposed */
struct bm_output_singlefile {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    std::ofstream f;
    char * iobuf;
    char * filename;
    bm_output_singlefile(bm_io &aa, matpoly const & P, const char * suffix = "")
        : aa(aa), P(P)
    {
        if (!aa.leader()) return;
        int rc = asprintf(&filename, "%s%s", aa.output_file, suffix);
        ASSERT_ALWAYS(rc >= 0);
        std::ios_base::openmode mode = std::ios_base::out;
        if (!aa.ascii) mode |= std::ios_base::binary;  
        f.open(filename, mode);
        DIE_ERRNO_DIAG(!f, "fopen", filename);
        iobuf = (char*) malloc(2 * io_block_size);
        f.rdbuf()->pubsetbuf(iobuf, 2 * io_block_size);
    }
    ~bm_output_singlefile()
    {
        if (!aa.leader()) return;
        printf("Saved %s\n", filename);
        free(filename);
        f.close();
        free(iobuf);
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.alloc;
        matpoly_write(aa.bm.d.ab, f, P, deg, deg + 1, aa.ascii, 1);
    }
};/*}}}*/
struct bm_output_splitfile {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    std::vector<std::ofstream> fw;
    bm_output_splitfile(bm_io & aa, matpoly const & P, const char * suffix = "")
        : aa(aa), P(P)
    {
        if (!aa.leader()) return;
        std::ios_base::openmode mode = std::ios_base::out;
        if (!aa.ascii) mode |= std::ios_base::binary;  
        for(unsigned int i = 0 ; i < P.m ; i++) {
            for(unsigned int j = 0 ; j < P.n ; j++) {
                char * str;
                int rc = asprintf(&str, "%s.sols%d-%d.%d-%d%s",
                        aa.output_file,
                        j, j + 1,
                        i, i + 1, suffix);
                ASSERT_ALWAYS(rc >= 0);
                fw.emplace_back(std::ofstream { str, mode } );
                DIE_ERRNO_DIAG(!fw.back(), "fopen", str);
                free(str);
            }
        }
        /* Do we want specific caching bufs here ? I doubt it */
        /*
           iobuf = (char*) malloc(2 * io_block_size);
           setbuffer(f, iobuf, 2 * io_block_size);
           */
    }
    ~bm_output_splitfile() {
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.alloc;
        matpoly_write_split(aa.bm.d.ab, fw, P, deg, deg + 1, aa.ascii);
    }

};/*}}}*/
struct bm_output_checksum {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    sha1_checksumming_stream f;
    const char * suffix;
    bm_output_checksum(bm_io & aa, matpoly const & P, const char * suffix = NULL)
        : aa(aa), P(P), suffix(suffix) { }
    ~bm_output_checksum() {
        char out[41];
        f.checksum(out);
        if (suffix)
            printf("checksum(%s): %s\n", suffix, out);
        else
            printf("checksum: %s\n", out);
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.alloc;
        matpoly_write(aa.bm.d.ab, f, P, deg, deg + 1, aa.ascii, 1);
    }
};/*}}}*/
/* }}} */

void bm_io::zero1(unsigned int deg)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int n = d.n;
    deg = deg % F.alloc;
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            abset_zero(ab, F.coeff(i, j, deg));
        }
    }
}/*}}}*/

/* Transfers of matrix entries, from memory to memory ; this is possibly
 * almost a transparent layer, and the matpoly_* instances really _are_
 * transparent layers. However we want a unified interface for the case
 * where we read an MPI-scattered matrix by chunks.
 *
 * Given that we are talking memory-to-memory, the classes below are used
 * as follows. The *_write_task classes are at the beginning of the
 * computation, where the matrix A is read from disk, and the matrix E is
 * written to (possibly distributed) memory. The one which matters is E
 * here. Symmetrically, the *_read_task are for reading the computed
 * matrix pi from (possibly distributed) memory, and writing the final
 * data F to disk.
 *
 * Pay attention to the fact that the calls to these structures are
 * collective calls.
 */
#ifdef ENABLE_MPI_LINGEN
class bigmatpoly_consumer_task { /* {{{ */
    /* This reads a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable. */
    bigmatpoly const & xpi;
    matpoly pi; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    int rank;

    public:
    bigmatpoly_consumer_task(bm_io & aa, bigmatpoly const & xpi) : xpi(xpi) {
        bmstatus & bm = aa.bm;
        bw_dimensions & d = bm.d;
        unsigned int m = d.m;
        unsigned int n = d.n;
        unsigned int b = m + n;
        MPI_Comm_rank(bm.com[0], &rank);

        /* Decide on the temp storage size */
        double avg = avg_matsize(d.ab, n, n, aa.ascii);
        B = iceildiv(io_block_size, avg);
        if (!rank && !random_input_length) {
            printf("Writing F to %s\n", aa.output_file);
            printf("Writing F by blocks of %u coefficients"
                    " (%.1f bytes each)\n", B, avg);
        }
        pi = matpoly(d.ab, b, b, B);

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    inline unsigned int size() { return xpi.size; }

    inline absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(!rank);
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        return pi.coeff(i, j, k - k0);
    }

    absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX || k - k0 >= B) {
            k0 = k - (k % B);
            pi.zero();
            pi.size = B;
            xpi.gather_mat_partial(pi, k0, MIN(xpi.size - k0, B));
        }
        if (rank) return NULL;
        return coeff_const_locked(i, j, k);
    }
};      /* }}} */
class bigmatpoly_producer_task { /* {{{ */
    /* This writes a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable.
     * Note that in any case, the coefficient indices must be progressive
     * in the write.
     */
    bigmatpoly & xE;
    matpoly E; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    int rank;

    /* forbid copies */
    bigmatpoly_producer_task(bigmatpoly_producer_task const &) = delete;

    public:
    bigmatpoly_producer_task(bm_io & aa, bigmatpoly & xE) : xE(xE) {
        bmstatus & bm = aa.bm;
        bw_dimensions & d = bm.d;
        unsigned int m = d.m;
        unsigned int n = d.n;
        unsigned int b = m + n;
        MPI_Comm_rank(aa.bm.com[0], &rank);


        /* Decide on the MPI chunk size */
        double avg = avg_matsize(d.ab, m, n, 0);
        B = iceildiv(io_block_size, avg);

        if (!rank) {
            /* TODO: move out of here */
            if (aa.input_file)
                printf("Reading A from %s\n", aa.input_file);
            printf("Computing E by chunks of %u coefficients (%.1f bytes each)\n",
                    B, avg);
        }

        E = matpoly(d.ab, m, b, B);

        /* Setting E.size is rather artificial, since we use E essentially
         * as an area where we may write coefficients freely. The only aim is
         * to escape some safety checks involving ->size in matpoly_part */
        E.zero();

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    // inline unsigned int size() { return xE.size; }
    inline void set_size(unsigned int s) {
        xE.set_size(s);
    }

    inline abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        ASSERT(!rank);

        E.size += (E.size == k - k0);
        ASSERT(E.size == k - k0 + 1);

        return E.coeff(i, j, k - k0);
    }

    abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX) {
            E.zero();
            E.size = 0;
            ASSERT(k == 0);
            k0 = 0;
        } else {
            if (k >= (k0 + B)) {
                /* We require progressive reads */
                ASSERT(k == k0 + B);
                xE.scatter_mat_partial(E, k0, B);
                E.zero();
                k0 += B;
            }
        }
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        if (rank) return NULL;
        return coeff_locked(i, j, k);
    }

    void finalize(unsigned int length) {
        if (length > k0) {
            /* Probably the last chunk hasn't been dispatched yet */
            ASSERT(length < k0 + B);
            xE.scatter_mat_partial(E, k0, length - k0);
        }
        set_size(length);
    }
    friend void matpoly_extract_column(
        bigmatpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc);
};

void matpoly_extract_column(
        bigmatpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc)
{
    dst.coeff_locked(0, jdst, kdst);
    dst.E.extract_column(jdst, kdst - dst.k0, src, jsrc, ksrc);
}

/* }}} */
#endif  /* ENABLE_MPI_LINGEN */
class matpoly_consumer_task {/*{{{*/
    /* This does the same, but for a simple matpoly. Of course this is
     * much simpler ! */
    matpoly const & pi;

    public:
    matpoly_consumer_task(bm_io & aa, matpoly const & pi) :
            pi(pi)
    {
        if (!random_input_length) {
            printf("Writing F to %s\n", aa.output_file);
        }
    }

    inline unsigned int chunk_size() const { return 1; }

    inline unsigned int size() { return pi.size; }

    absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        return pi.coeff(i, j, k);
    }
    inline absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_const_locked(i, j, k);
    }

    ~matpoly_consumer_task() { }
};/*}}}*/
class matpoly_producer_task { /* {{{ */
    matpoly & E;

    /* forbid copies */
    matpoly_producer_task(matpoly_producer_task const &) = delete;

    public:
    matpoly_producer_task(bm_io & aa, matpoly & E) :
        E(E)
    {
            /* TODO: move out of here */
            if (aa.input_file)
                printf("Reading A from %s\n", aa.input_file);
    }

    inline unsigned int chunk_size() const { return 1; }
    inline void set_size(unsigned int s) { E.size = s; }
    // inline unsigned int size() { return E.size; }

    abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        E.size += (E.size == k);
        ASSERT(E.size == k + 1);
        return E.coeff(i, j, k);
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_locked(i, j, k);
    }
    void finalize(unsigned int length) { set_size(length); }
    ~matpoly_producer_task() { }
    friend void matpoly_extract_column(
        matpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc);
};
void matpoly_extract_column(
        matpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc)
{
    dst.E.extract_column(jdst, kdst, src, jsrc, ksrc);
}

/* }}} */


template<class Consumer, class Sink>
void bm_io::compute_final_F(Sink & S, Consumer& pi, std::vector<unsigned int> & delta)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    int leader = rank == 0;


    /* We are not interested by pi.size, but really by the number of
     * coefficients for the columns which give solutions. */
    unsigned int maxdelta = get_max_delta_on_solutions(bm, delta);

    if (leader) printf("Final f(X)=f0(X)pi(X) has degree %u\n", maxdelta);

    unsigned int window = F.alloc;

    /* Which columns of F*pi will make the final generator ? */
    std::vector<unsigned int> sols(n, 0);
    for(unsigned int j = 0, jj=0 ; j < m + n ; j++) {
        if (bm.lucky[j] <= 0)
            continue;
        sols[jj++]=j;
    }

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_s = pi.size() / 100;

    /*
     * first compute the rhscontribs. Use that to decide on a renumbering
     * of the columns, because we'd like to get the "primary" solutions
     * first, in a sense. Those are the ones with fewer zeroes in the
     * rhscontrib part. So we would want that matrix to have its columns
     * sorted in decreasing weight order
     *
     * An alternative, possibly easier, is to have a function which
     * decides the solution ordering precisely based on the inspection of
     * this rhscoeffs matrix (?). But how should we spell that info when
     * we give it to mksol ??
     */

    /* This **modifies** the "sols" array */
    if (d.nrhs) {
        /* The data is only gathered at rank 0. So the stuff we compute
         * can only be meaningful there. On the other hand, it is
         * necessary that we call *collectively* the pi.coeff() routine,
         * because we need synchronisation of all ranks.
         */

        matpoly rhs(ab, d.nrhs, n, 1);
        rhs.size = 1;

        {
            Sink Srhs(*this, rhs, ".rhs");

            /* Now redo the exact same loop as above, this time
             * adding the contributions to the rhs matrix. */
            for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int iF, offset;
                    if (ipi < n) {
                        iF = ipi;
                        offset = 0;
                    } else {
                        iF = fdesc[ipi-n][1];
                        offset = t0 - fdesc[ipi-n][0];
                    }
                    if (iF >= d.nrhs) continue;
                    ASSERT_ALWAYS(delta[jpi] >= offset);
                    unsigned kpi = delta[jpi] - offset;

                    ASSERT_ALWAYS(d.nrhs);
                    ASSERT_ALWAYS(iF < d.nrhs);
                    ASSERT_ALWAYS(jF < n);
                    absrc_elt src = pi.coeff_const(ipi, jpi, kpi);

                    if (leader) {
                        abdst_elt dst = rhs.coeff(iF, jF, 0);
                        abadd(ab, dst, dst, src);
                    }
                }
            }

            if (leader) {
                /* Now comes the time to prioritize the different solutions. Our
                 * goal is to get the unessential solutions last ! */
                int (*sol_score)[2] = new int[n][2];
                memset(sol_score, 0, n * 2 * sizeof(int));
                /* score per solution is the number of non-zero coefficients,
                 * that's it. Since we have access to lexcmp2, we want to use it.
                 * Therefore, desiring the highest scoring solutions first, we
                 * negate the hamming weight.
                 */
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    sol_score[jF][1] = jF;
                    for(unsigned int iF = 0 ; iF < d.nrhs ; iF++) {
                        int z = !abis_zero(ab, rhs.coeff(iF, jF, 0));
                        sol_score[jF][0] -= z;
                    }
                }
                qsort(sol_score, n, 2 * sizeof(int), (sortfunc_t) & lexcmp2);

                if (leader) {
                    printf("Reordered solutions:\n");
                    for(unsigned int i = 0 ; i < n ; i++) {
                        printf(" %d (col %d in pi, weight %d on rhs vectors)\n", sol_score[i][1], sols[sol_score[i][1]], -sol_score[i][0]);
                    }
                }

                /* We'll now modify the sols[] array, so that we get a reordered
                 * F, too (and mksol/gather don't have to care about our little
                 * tricks */
                {
                    matpoly rhs2(ab, d.nrhs, n, 1);
                    rhs2.size = 1;
                    for(unsigned int i = 0 ; i < n ; i++) {
                        rhs2.extract_column(i, 0, rhs, sol_score[i][1], 0);
                    }
                    rhs = std::move(rhs2);
                    if (sol_score[0][0] == 0) {
                        if (allow_zero_on_rhs) {
                            printf("Note: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (ok because of allow_zero_on_rhs=1)\n");
                        } else {
                            fprintf(stderr, "ERROR: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (maybe use allow_zero_on_rhs ?)\n");
                            rank0_exit_code = EXIT_FAILURE;
                        }
                    }
                    /* ugly: use sol_score[i][0] now to provide the future
                     * "sols" array. We'll get rid of sol_score right afterwards
                     * anyway.
                     */
                    for(unsigned int i = 0 ; i < n ; i++) {
                        sol_score[i][0] = sols[sol_score[i][1]];
                    }
                    for(unsigned int i = 0 ; i < n ; i++) {
                        sols[i] = sol_score[i][0];
                    }
                }
                delete[] sol_score;

                Srhs.write1(0);
            }
        }
    }

    /* we need to read pi backwards. The number of coefficients in pi is
     * pilen = maxdelta + 1 - t0. Hence the first interesting index is
     * maxdelta - t0. However, for notational ease, we'll access
     * coefficients from index pi.size downwards. The latter is always
     * large enough.
     */

    ASSERT_ALWAYS(pi.size() >= maxdelta + 1 - t0);

    for(unsigned int s = 0 ; s < pi.size() ; s++) {
        unsigned int kpi = pi.size() - 1 - s;

        /* as above, we can't just have ranks > 0 do nothing. We need
         * synchronization of the calls to pi.coeff()
         */

        /* This call is here only to trigger the gather call. This one
         * must therefore be a collective call. Afterwards we'll use a
         * _locked call which is okay to call only at rank 0.
         */
        pi.coeff_const(0, 0, kpi);

        if (rank) continue;

        /* Coefficient kpi + window of F has been totally computed,
         * because of previous runs of this loop (which reads the
         * coefficients of pi).
         */
        if (kpi + window == maxdelta) {
            /* Avoid writing zero coefficients. This can occur !
             * Example:
             * tt=(2*4*1200) mod 1009, a = (tt cat tt cat * tt[0..10])
             */
            for(unsigned int j = 0 ; j < n ; j++) {
                int z = 1;
                for(unsigned int i = 0 ; z && i < n ; i++) {
                    absrc_elt src = F.coeff(i, j, 0);
                    z = abis_zero(ab, src);
                }

                if (z) {
                    /* This is a bit ugly. Given that we're going to
                     * shift one column of F, we'll have a
                     * potentially deeper write-back buffer. Columns
                     * which seemed to be ready still are, but they
                     * will now be said so only at the next step.
                     */
                    printf("Reduced solution column #%u from"
                            " delta=%u to delta=%u\n",
                            sols[j], delta[sols[j]], delta[sols[j]]-1);
                    window++;
                    F.realloc(window);
                    F.size = window;
                    delta[sols[j]]--;
                    /* shift this column */
                    for(unsigned int k = 1 ; k < window ; k++) {
                        F.extract_column(j, k-1, F, j, k);
                    }
                    F.zero_column(j, window - 1);
                    break;
                }
            }
        }
        /* coefficient of degree maxdelta-kpi-window is now complete */
        if (kpi + window <= maxdelta) {
            S.write1((maxdelta-kpi) - window);
            zero1((maxdelta-kpi) - window);
        }
        /* Now use our coefficient. This might tinker with
         * coefficients up to degree kpi-(window-1) in the file F */

        if (kpi > maxdelta + t0 ) {
            /* this implies that we always have kF > delta[jpi],
             * whence we expect a zero contribution */
            for(unsigned int i = 0 ; i < m + n ; i++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    absrc_elt src = pi.coeff_const_locked(i, jpi, kpi);
                    ASSERT_ALWAYS(abis_zero(ab, src));
                }
            }
            continue;
        }

        for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
            for(unsigned int jF = 0 ; jF < n ; jF++) {
                unsigned int jpi = sols[jF];
                unsigned int iF, offset;
                if (ipi < n) {
                    /* Left part of the initial F is x^0 times
                     * identity. Therefore, the first n rows of pi
                     * get multiplied by this identity matrix, this
                     * is pretty simple.
                     */
                    iF = ipi;
                    offset = 0;
                } else {
                    /* next m rows of the initial F are of the form
                     * x^(some value) times some canonical basis
                     * vector. Therefore, the corresponding row in pi
                     * ends up contributing to some precise row in F,
                     * and with an offset which is dictated by the
                     * exponent of x.
                     */
                    iF = fdesc[ipi-n][1];
                    offset = t0 - fdesc[ipi-n][0];
                }
                unsigned int subtract = maxdelta - delta[jpi] + offset;
                ASSERT(subtract < window);
                if (maxdelta < kpi + subtract) continue;
                unsigned int kF = (maxdelta - kpi) - subtract;
                unsigned int kF1 = kF - (iF < d.nrhs);
                abdst_elt dst;
                if (kF1 == UINT_MAX) {
                    /* this has been addressed in the first pass,
                     * earlier.
                     */
                    continue;
                } else {
                    dst = F.coeff(iF, jF, kF1 % window);
                }
                absrc_elt src = pi.coeff_const_locked(ipi, jpi, kpi);
                ASSERT_ALWAYS(kF <= delta[jpi] || abis_zero(ab, src));
                abadd(ab, dst, dst, src);
            }
        }

        if (leader && s > next_report_s) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Written %u coefficients (%.1f%%) in %.1f s\n",
                        s, 100.0 * s / pi.size(), tt-tt0);
                next_report_t = tt + 10;
                next_report_s = s + pi.size() / 100;
            }
        }
    }
    /* flush the pipe */
    if (leader && window <= maxdelta) {
        for(unsigned int s = window ; s-- > 0 ; )
            S.write1(maxdelta - s);
    }
}/*}}}*/

/* read 1 coefficient into the sliding window of input coefficients of
 * the input series A. The io_window parameter controls the size of the
 * sliding window. There are in fact two behaviours:
 *  - io_window == 0: there is no sliding window, really, and the new
 *    coefficient is appended as the last coefficient of A.
 *  - io_window > 0: there, we really have a sliding window. Coeffs
 *    occupy places in a circular fashion within the buffer.
 */
int bm_io::read1(unsigned int io_window)/*{{{*/
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    static unsigned int generated_random_coefficients = 0;

    unsigned int pos = k;
    if (io_window) {
        pos = pos % io_window;
        ASSERT_ALWAYS(A.size == io_window);
    } else {
        ASSERT_ALWAYS(A.size == k);
        if (k >= A.alloc)
            A.realloc(A.alloc + 1);
        ASSERT_ALWAYS(k < A.alloc);
        A.size++;
    }
    if (random_input_length) {
        if (generated_random_coefficients >= random_input_length)
            return 0;

        for (unsigned int i = 0; i < m ; i++) {
            for (unsigned int j = 0; j < n ; j++) {
                abdst_elt x = A.coeff(i, j, pos);
                abrandom(ab, x, rstate);
            }
        }
        generated_random_coefficients++;
        k++;
        return 1;
    }

    for (unsigned int i = 0; i < m ; i++) {
        for (unsigned int j = 0; j < n ; j++) {
            abdst_elt x = A.coeff(i, j, pos);
            int rc;
            if (ascii) {
                rc = abfscan(ab, fr[0], x);
                /* rc is the number of bytes read -- non-zero on success */
            } else {
                size_t elemsize = abvec_elt_stride(ab, 1);
                rc = fread(x, elemsize, 1, fr[0]);
                rc = rc == 1;
                abnormalize(ab, x);
            }
            if (!rc) {
                if (i == 0 && j == 0) {
                    return 0;
                }
                fprintf(stderr,
                        "Parse error while reading coefficient (%d,%d,%d)%s\n",
                        i, j, 1 + k,
                        ascii ? "" : " (forgot --ascii?)");
                exit(EXIT_FAILURE);
            }
        }
    }
    k++;
    return 1;
}/*}}}*/

bm_io::bm_io(bmstatus & bm, const char * input_file, const char * output_file, int ascii)/*{{{*/
    : bm(bm)
    , input_file(input_file)
    , output_file(output_file)
    , ascii(ascii)
    , A(bm.d.ab, bm.d.m, bm.d.n, 1)
{
}/*}}}*/

bm_io::~bm_io()/*{{{*/
{
    if (fdesc) free(fdesc);
}/*}}}*/

void bm_io::begin_read()/*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    if (rank) return;

    if (random_input_length) {
        /* see below. I think it would be a bug to not do that */
        read1(0);
        return;
    }

    if (split_input_file) {
        fprintf(stderr, "--split-input-file not supported yet\n");
        exit(EXIT_FAILURE);
    }
    fr = (FILE**) malloc(sizeof(FILE*));
    fr[0] = fopen(input_file, ascii ? "r" : "rb");

    DIE_ERRNO_DIAG(fr[0] == NULL, "fopen", input_file);
    iobuf = (char*) malloc(2 * io_block_size);
    setbuffer(fr[0], iobuf, 2 * io_block_size);

    /* read the first coefficient ahead of time. This is because in most
     * cases, we'll discard it. Only in the DL case, we will consider the
     * first coefficient as being part of the series. Which means that
     * the coefficient reads in the I/O loop will sometimes correspond to
     * the coefficient needed at that point in time, while we will also
     * (in the DL case) need data from the previous read.
     */
    if (!read1(0)) {
        fprintf(stderr, "Read error from %s\n", input_file);
        exit(EXIT_FAILURE);
    }
}/*}}}*/

void bm_io::end_read()/*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    if (rank) return;
    if (random_input_length) return;
    fclose(fr[0]);
    free(fr);
    fr = NULL;
    free(iobuf);
    iobuf = 0;
}/*}}}*/

void bm_io::guess_length()/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    if (random_input_length) {
        guessed_length = random_input_length;
        return;
    }

    if (!rank) {
        struct stat sbuf[1];
        int rc = fstat(fileno(fr[0]), sbuf);
        if (rc < 0) {
            perror(input_file);
            exit(EXIT_FAILURE);
        }

        size_t filesize = sbuf->st_size;

        if (!ascii) {
            size_t avg = avg_matsize(ab, m, n, ascii);
            if (filesize % avg) {
                fprintf(stderr, "File %s has %zu bytes, while its size should be amultiple of %zu bytes (assuming binary input; perhaps --ascii is missing ?).\n", input_file, filesize, avg);
                exit(EXIT_FAILURE);
            }
            guessed_length = filesize / avg;
        } else {
            double avg = avg_matsize(ab, m, n, ascii);
            double expected_length = filesize / avg;
            if (!rank)
                printf("# Expect roughly %.2f items in the sequence.\n", expected_length);

            /* First coefficient is always lighter, so we add a +1. */
            guessed_length = 1 + expected_length;
        }
    }
    MPI_Bcast(&(guessed_length), 1, MPI_UNSIGNED, 0, bm.com[0]);
}/*}}}*/

void bm_io::compute_initial_F() /*{{{ */
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    int rank;
    fdesc = (unsigned int(*)[2])malloc(2 * m * sizeof(unsigned int));
    MPI_Comm_rank(bm.com[0], &rank);
    if (!rank) {
        /* read the first few coefficients. Expand A accordingly as we are
         * doing the read */

        ASSERT(A.m == m);
        ASSERT(A.n == n);

        abelt tmp;
        abinit(ab, &tmp);

        /* First try to create the initial F matrix */
        printf("Computing t0\n");

        /* We want to create a full rank m*m matrix M, by extracting columns
         * from the first coefficients of A */

        matpoly M(ab, m, m, 1);
        M.size = 1;

        /* For each integer i between 0 and m-1, we have a column, picked
         * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
         * the other ones, has coefficient at row pivots[i] unequal to zero.
         */
        std::vector<unsigned int> pivots(m, 0);
        std::vector<unsigned int> exponents(m, 0);
        std::vector<unsigned int> cnum(m, 0);
        unsigned int r = 0;

        for (unsigned int k = 0; r < m ; k++) {
            /* read a new coefficient */
            read1(0);

            for (unsigned int j = 0; r < m && j < n; j++) {
                /* Extract a full column into M (column j, degree k in A) */
                /* adjust the coefficient degree to take into account the
                 * fact that for SM columns, we might in fact be
                 * interested by the _previous_ coefficient */
                M.extract_column(r, 0, A, j, k + (j >= bm.d.nrhs));

                /* Now reduce it modulo all other columns */
                for (unsigned int v = 0; v < r; v++) {
                    unsigned int u = pivots[v];
                    /* the v-th column in the M is known to
                     * kill coefficient u (more exactly, to have a -1 as u-th
                     * coefficient, and zeroes for the other coefficients
                     * referenced in the pivots[0] to pivots[v-1] indices).
                     */
                    /* add M[u,r]*column v of M to column r of M */
                    for(unsigned int i = 0 ; i < m ; i++) {
                        if (i == u) continue;
                        abmul(ab, tmp,
                                  M.coeff(i, v, 0),
                                  M.coeff(u, r, 0));
                        abadd(ab, M.coeff(i, r, 0),
                                  M.coeff(i, r, 0),
                                  tmp);
                    }
                    abset_zero(ab,
                            M.coeff(u, r, 0));
                }
                unsigned int u = 0;
                for( ; u < m ; u++) {
                    if (abcmp_ui(ab, M.coeff(u, r, 0), 0) != 0)
                        break;
                }
                if (u == m) {
                    printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                           k + (j >= bm.d.nrhs), j, r);

                    /* we need at least m columns to get as starting matrix
                     * with full rank. Given that we have n columns per
                     * coefficient, this means at least m/n matrices.
                     */

                    if (k * n > m + 40) {
                        printf("The choice of starting vectors was bad. "
                               "Cannot find %u independent cols within A\n", m);
                        exit(EXIT_FAILURE);
                    }
                    continue;
                }

                /* Bingo, it's a new independent col. */
                pivots[r] = u;
                cnum[r] = j;
                exponents[r] = k;

                /* Multiply the column so that the pivot becomes -1 */
                int rc = abinv(ab, tmp, M.coeff(u, r, 0));
                if (!rc) {
                    fprintf(stderr, "Error, found a factor of the modulus: ");
                    abfprint(ab, stderr, tmp);
                    fprintf(stderr, "\n");
                    exit(EXIT_FAILURE);
                }
                abneg(ab, tmp, tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    abmul(ab, M.coeff(i, r, 0),
                              M.coeff(i, r, 0),
                              tmp);
                }

                r++;

                // if (r == m)
                    printf
                        ("[X^%d] A, col %d increases rank to %d (head row %d)%s\n",
                         k + (j >= bm.d.nrhs), j, r, u,
                         (j < bm.d.nrhs) ? " (column not shifted because of the RHS)":"");
            }
        }

        if (r != m) {
            printf("This amount of data is insufficient. "
                   "Cannot find %u independent cols within A\n", m);
            exit(EXIT_FAILURE);
        }

        t0 = exponents[r - 1] + 1;
        /* We always have one extra coefficient of back log */
        ASSERT_ALWAYS(t0 == k - 1);

        printf("Found satisfying init data for t0=%d\n", t0);

        /* We've also got some adjustments to make: room for one extra
         * coefficient is needed in A. Reading of further coefficients will
         * pickup where they stopped, and will always leave the last t0+2
         * coefficients readable. */
        A.realloc(t0 + 2);
        A.size++;


        for(unsigned int j = 0 ; j < m ; j++) {
            fdesc[j][0] = exponents[j];
            fdesc[j][1] = cnum[j];
            ASSERT_ALWAYS(exponents[j] < t0);
        }
        // free(pivots);
        // free(exponents);
        // free(cnum);
        // matpoly_clear(ab, M);
        abclear(ab, &tmp);
    }
    MPI_Bcast(fdesc, 2*m, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&(t0), 1, MPI_UNSIGNED, 0, bm.com[0]);
    bm.t = t0;
}				/*}}} */

template<class Writer>
void bm_io::compute_E(Writer& E, unsigned int expected, unsigned int allocated)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    unsigned int window = t0 + 2;

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_k = expected / 100;

    int over = 0;
    unsigned int final_length = UINT_MAX;

    for(unsigned int kE = 0 ; kE + t0 < allocated ; kE ++) {

        if (kE % E.chunk_size() || kE + t0 >= expected) {
            MPI_Bcast(&over, 1, MPI_INT, 0, bm.com[0]);
            if (over) break;
        }
        if (rank == 0 && !over) {
            over = !read1(window);
            if (over) {
                printf("EOF met after reading %u coefficients\n", k);
                final_length = kE;
            }
        }

        E.coeff(0, 0, kE);

        if (over || rank)
            continue;

        if (kE + t0 > allocated) {
            fprintf(stderr, "Going way past guessed length%s ???\n", ascii ? " (more than 5%%)" : "");
        }

        for(unsigned int j = 0 ; j < n ; j++) {
            /* If the first columns of F are the identity matrix, then
             * in E we get data from coefficient kE+t0 in A. More
             * generally, if it's x^q*identity, we read
             * coeficient of index kE + t0 - q.
             *
             * Note that we prefer to take q=0 anyway, since a
             * choice like q=t0 would create duplicate rows in E,
             * and that would be bad.
             */
            unsigned int kA = kE + t0 + (j >= bm.d.nrhs);
            matpoly_extract_column(E, j, kE, A, j, kA % window);
        }

        for(unsigned int jE = n ; jE < m + n ; jE++) {
            unsigned int jA = fdesc[jE-n][1];
            unsigned int offset = fdesc[jE-n][0];
            unsigned int kA = kE + offset + (jA >= bm.d.nrhs);
            matpoly_extract_column(E, jE, kE, A, jA, kA % window);
        }
        /* This is because k integrates some backlog because of
         * the SM / non-SM distinction (for DL) */
        ASSERT_ALWAYS(kE + 1 + t0 + 1 == k);
        if (k > next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Read %u coefficients (%.1f%%)"
                        " in %.1f s (%.1f MB/s)\n",
                        k, 100.0 * k / expected,
                        tt-tt0, k * avg_matsize(ab, m, n, ascii) / (tt-tt0)/1.0e6);
                next_report_t = tt + 10;
                next_report_k = k + expected / 100;
            }
        }
    }
    MPI_Bcast(&final_length, 1, MPI_UNSIGNED, 0, bm.com[0]);
    ASSERT_ALWAYS(final_length != UINT_MAX);
    E.finalize(final_length);
}/*}}}*/

template<typename T> struct matpoly_factory {};

#ifdef ENABLE_MPI_LINGEN
template<> struct matpoly_factory<bigmatpoly> {
    typedef bigmatpoly T;
    typedef bigmatpoly_producer_task producer_task;
    typedef bigmatpoly_consumer_task consumer_task;
    bigmatpoly_model model;
    matpoly_factory(MPI_Comm * comm, unsigned int m, unsigned int n) : model(comm, m, n) {}
    T init(abdst_field ab, unsigned int m, unsigned int n, int len) {
        return bigmatpoly(ab, model, m, n, len);
    }
    static int bw_lingen(bmstatus & bm, T & pi, T & E, std::vector<unsigned int> & delta) {
        return bw_biglingen_collective(bm, pi, E, delta);
    }
    static size_t alloc(T const & p) { return p.my_cell().alloc; }
};
#endif

template<> struct matpoly_factory<matpoly> {
    typedef matpoly T;
    typedef matpoly_producer_task producer_task;
    typedef matpoly_consumer_task consumer_task;
    matpoly_factory() {}
    T init(abdst_field ab, unsigned int m, unsigned int n, int len) {
        return matpoly(ab, m, n, len);
    }
    static int bw_lingen(bmstatus & bm, T & pi, T & E, std::vector<unsigned int> & delta) {
        return bw_lingen_single(bm, pi, E, delta);
    }
    static size_t alloc(T const & p) { return p.alloc; }

};

template<typename T, typename Sink>
void bm_io::output_flow(T & pi, std::vector<unsigned int> & delta)
{
    unsigned int n = bm.d.n;

    matpoly::memory_pool_guard dummy(SIZE_MAX);

    F = matpoly(bm.d.ab, n, n, t0 + 1);

    set_write_behind_size(delta);

    Sink S(*this, F);

    typename matpoly_factory<T>::consumer_task pi_consumer(*this, pi);

    compute_final_F(S, pi_consumer, delta);

    /* We need this because we want all our deallocation to happen before
     * the guard's dtor gets called */
    F = matpoly();
}


/*}}}*/

void usage()
{
    fprintf(stderr, "Usage: ./plingen [options, to be documented]\n");
    fprintf(stderr,
            "General information about bwc options is in the README file\n");
    exit(EXIT_FAILURE);
}

unsigned int count_lucky_columns(bmstatus & bm)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    int luck_mini = expected_pi_length(d);
    MPI_Bcast(&bm.lucky[0], b, MPI_UNSIGNED, 0, bm.com[0]);
    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; nlucky += bm.lucky[j++] >= luck_mini) ;
    return nlucky;
}/*}}}*/

int check_luck_condition(bmstatus & bm)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int nlucky = count_lucky_columns(bm);

    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    if (!rank) {
        printf("Number of lucky columns: %u (%u wanted)\n", nlucky, n);
    }

    if (nlucky == n)
        return 1;

    if (!rank) {
        fprintf(stderr, "Could not find the required set of solutions (nlucky=%u)\n", nlucky);
    }
    if (random_input_length) {
        static int once=0;
        if (once++) {
            if (!rank) {
                fprintf(stderr, "Solution-faking loop crashed\n");
            }
            MPI_Abort(bm.com[0], EXIT_FAILURE);
        }
        if (!rank) {
            printf("Random input: faking successful computation\n");
        }
        for(unsigned int j = 0 ; j < n ; j++) {
            bm.lucky[(j * 1009) % (m+n)] = expected_pi_length(d);
        }
        return check_luck_condition(bm);
    }

    return 0;
}/*}}}*/

void display_deltas(bmstatus & bm, std::vector<unsigned int> const & delta)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;

    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    if (!rank) {
        printf("Final, t=%u: delta =", bm.t);
        for(unsigned int j = 0; j < m + n; j++) {
            printf(" %u", delta[j]);
            if (bm.lucky[j] < 0) {
                printf("(*)");
            }
        }
        printf("\n");
    }
}/*}}}*/

void print_node_assignment(MPI_Comm comm)/*{{{*/
{
    int rank;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    struct utsname me[1];
    int rc = uname(me);
    if (rc < 0) { perror("uname"); MPI_Abort(comm, 1); }
    size_t sz = 1 + sizeof(me->nodename);
    char * global = (char*) malloc(size * sz);
    memset(global, 0, size * sz);
    memcpy(global + rank * sz, me->nodename, sizeof(me->nodename));

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            global, sz, MPI_BYTE, comm);
    if (rank == 0) {
        char name[80];
        int len=80;
        MPI_Comm_get_name(comm, name, &len);
        name[79]=0;
        for(int i = 0 ; i < size ; i++) {
            printf("# %s rank %d: %s\n", name, i, global + i * sz);
        }
    }
    free(global);
}/*}}}*/

void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate)/*{{{*/
{
    /* used by testing code */
    bmstatus bm(m,n);
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    abfield_specify(bm.d.ab, MPFQ_PRIME_MPZ, p);
    unsigned int t0 = bm.t = iceildiv(m,n);
    std::vector<unsigned int> delta(m+n, t0);
    matpoly E(ab, m, m+n, L);
    E.fill_random(L, rstate);
    matpoly pi;
    bw_lingen_basecase_raw(bm, pi, E, delta);
    mpz_clear(p);
}/*}}}*/

/* Counting memory usage in the recursive algorithm.
 *
 * The recursive algorithm is designed to allow the allocated memory for
 * the input to be reused for placing the output. Some memory might have
 * been saved by upper layers. We also have some local allocation.
 *
 * Notations: The algorithm starts at depth 0 with an
 * input length L, and the notation \ell_i denotes L/2^(i+1). We have
 * \ell_i=2\ell_{i+1}. The notation \alpha denotes m/(m+n). Note that the
 * input has size \alpha*(1-\alpha)*L times (m+n)^2*\log_2(p) (divided by
 * r^2 if relevant).
 *
 * We define five quantities. All are understood as multiples of
 * (m+n)^2*\log_2(p).
 *
 * MP(i) is the extra storage needed for the MP operation at depth i.
 *
 * MUL(i) is the extra storage needed for the MUL operation at depth i.
 *
 * IO(i) is the common size of the input and output data of the call at
 *       depth i. We have
 *              IO(i) = 2\alpha\ell_i
 *
 * ST(i) is the storage *at all levels above the current one* (i.e. with
 *    depth strictly less than i) for the data that is still live and
 *    need to exist until after we return. This count is maximized in the
 *    leftmost branch, where chopped E at all levels must be kept.
 *    chopped E at depth i (not counted in ST(i) !) is:
 *          \alpha(1+\alpha) \ell_i
 *    (counted as the degree it takes to make the necessary data that
 *    we want to use to compute E_right),
 *    so the cumulated cost above is twice the depth 0 value, minus the
 *    depth i value, i.e.
 *              ST(i) = \alpha(1+\alpha)(L-2\ell_i).
 * SP(i) is the "spike" at depth i: not counting allocation that is
 *    already reserved for IO or ST, this is the amount of extra memory
 *    that is required by the call at depth i. We have:
 *      SP(i) = max {
 *              \alpha\ell_i,
 *              \alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+SP(i+1),
 *             2\alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *             2\alpha\ell_i-2\alpha\ell_i+SP(i+1)
 *             4\alpha\ell_i-2\alpha\ell_i+MUL(i)
 *             }
 *            = max {
 *              \alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+SP(i+1),
 *             2\alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *             4\alpha\ell_i-2\alpha\ell_i+MUL(i)
 *                           }
 * 
 * Combining this together, and using
 * ST(i)+\alpha(1+\alpha)\ell_i=ST(i+1), we have:
 *
 * IO(i)+ST(i)+SP(i) = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              ST(i) + 2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *              ST(i) + 4\alpha\ell_i+MUL(i)
 *                      }
 *                   = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              \alpha(1+\alpha)(L-\ell_i)  + 2\alpha\ell_i + MP(i),
 *              \alpha(1+\alpha)(L-2\ell_i) + 4\alpha\ell_i + MUL(i)
 *                      }
 *                   = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              \alpha((1+\alpha)L+(1-\alpha)\ell_i) + MP(i),
 *              \alpha((1+\alpha)L+2(1-\alpha)\ell_i) + MUL(i),
 *                      }
 *
 * Let RMP(i) be the amount of memory that is reserved while we are doing
 * the MP operation, and define RMUL similarly. We have:
 *      RMP(i)  = \alpha(1+\alpha)(L-\ell_i)  + 2\alpha\ell_i
 *      RMUL(i) = \alpha(1+\alpha)(L-2\ell_i) + 4\alpha\ell_i
 * whence:
 *      RMP(i) = \alpha((1+\alpha)L+(1-\alpha)\ell_i)
 *      RMUL(i) = \alpha((1+\alpha)L+2(1-\alpha)\ell_i)
 *
 * We have RMP(i) <= RMUL(i) <= RMP(0) <= RMUL(0) = 2\alpha*L. We'll use
 * the un-simplified expression later.
 *
 * Furthermore IO(infinity)=SP(infinity)=0, and ST(infinity)=\alpha(1+\alpha)L
 *
 * So that eventually, the amount of reserved memory for the whole
 * algorithm is RMUL(0)=2\alpha*L (which is 2/(1-\alpha)=2*(1+m/n) times
 * the input size). On top of that we have the memory required
 * for the transforms.
 *
 *
 * When going MPI, matrices may be rounded with some inaccuracy.
 * Splitting in two a 3x3 matrix leads to a 2x2 chunk, which is 1.77
 * times more than the simplistic proportionality rule.
 *
 * Therefore it makes sense to distinguish between matrices of size
 * m*(m+n) and (m+n)*(m+n). If we recompute RMUL(i) by taking this into
 * account, we obtain:
 *      [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[4\alpha\ell_i]
 * where we only paid attention to the rounding issues with dimensions,
 * as those are more important than for degrees. Bottom line, the max is
 * expected to be for i=0, and that will be made only of pi matrices.
 */

/* Some of the early reading must be done before we even start, since
 * the code that we run depends on the input size.
 */
template<typename T>
void lingen_main_code(matpoly_factory<T> & F, abdst_field ab, bm_io & aa)
{
    int rank;
    bmstatus & bm(aa.bm);
    MPI_Comm_rank(bm.com[0], &rank);
    unsigned int m = bm.d.m;
    unsigned int n = bm.d.n;
    unsigned int guess = aa.guessed_length;
    size_t safe_guess = aa.ascii ? ceil(1.05 * guess) : guess;

    std::vector<unsigned int> delta(m+n, bm.t);

    
    /* c0 is (1+m/n) times the input size */
    size_t c0 = abvec_elt_stride(bm.d.ab,
                iceildiv(m+n, bm.mpi_dims[0]) *
                iceildiv(m+n, bm.mpi_dims[1]) *
                iceildiv(m*safe_guess, m+n));
    matpoly::memory_pool_guard main(2*c0);
    if (!rank) {
        char buf[20];
        printf("# Estimated memory for JUST transforms (per node): %s\n",
                size_disp(2*c0, buf));
        printf("# Estimated peak total memory (per node): max at depth %d: %s\n",
                bm.hints.ipeak,
                size_disp(bm.hints.peak, buf));
    }


    T E  = F.init(ab, m, m+n, safe_guess);
    T pi = F.init(ab, 0, 0, 0);   /* pre-init for now */

    typename matpoly_factory<T>::producer_task E_producer(aa, E);

    aa.compute_E(E_producer, guess, safe_guess);
    aa.end_read();

    matpoly_factory<T>::bw_lingen(bm, pi, E, delta);
    bm.stats.final_print();

    display_deltas(bm, delta);
    if (!rank) printf("(pi.alloc = %zu)\n", matpoly_factory<T>::alloc(pi));

    if (check_luck_condition(bm)) {
        if (random_input_length) {
            aa.output_flow<T, bm_output_checksum>(pi, delta);
        } else if (split_output_file) {
            aa.output_flow<T, bm_output_splitfile>(pi, delta);
        } else {
            aa.output_flow<T, bm_output_singlefile>(pi, delta);
        }
    }

    /* clear everything */
}

/* We don't have a header file for this one */
extern "C" void check_for_mpi_problems();

int main(int argc, char *argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

    bw_common_decl_usage(pl);
    plingen_decl_usage(pl);
    logline_decl_usage(pl);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    /* {{{ interpret our parameters */
    gmp_randinit_default(rstate);

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    logline_init_timer();

    param_list_parse_int(pl, "allow_zero_on_rhs", &allow_zero_on_rhs);
    param_list_parse_uint(pl, "random-input-with-length", &random_input_length);
    param_list_parse_int(pl, "split-output-file", &split_output_file);
    param_list_parse_int(pl, "split-input-file", &split_input_file);

    const char * afile = param_list_lookup_string(pl, "afile");

    if (bw->m == -1) {
	fprintf(stderr, "no m value set\n");
	exit(EXIT_FAILURE);
    }
    if (bw->n == -1) {
	fprintf(stderr, "no n value set\n");
	exit(EXIT_FAILURE);
    }
    if (!global_flag_tune && !(afile || random_input_length)) {
        fprintf(stderr, "No afile provided\n");
        exit(EXIT_FAILURE);
    }

    /* we allow ffile and ffile to be both NULL */
    const char * tmp = param_list_lookup_string(pl, "ffile");
    char * ffile = NULL;
    if (tmp) {
        ffile = strdup(tmp);
    } else if (afile) {
        int rc = asprintf(&ffile, "%s.gen", afile);
        ASSERT_ALWAYS(rc >= 0);
    }
    ASSERT_ALWAYS((afile==NULL) == (ffile == NULL));

    bmstatus bm(bw->m, bw->n);
    bw_dimensions & d = bm.d;

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    if (!global_flag_tune && !random_input_length) {
        if (!rhs_name) {
            fprintf(stderr, "When using plingen, you must either supply --random-input-with-length, or provide a rhs, or possibly provide rhs=none\n");
        } else if (strcmp(rhs_name, "none") == 0) {
            rhs_name = NULL;
        }
    }
    if ((rhs_name != NULL) && param_list_parse_uint(pl, "nrhs", &(bm.d.nrhs))) {
        fprintf(stderr, "the command line arguments rhs= and nrhs= are incompatible\n");
        exit(EXIT_FAILURE);
    }
    if (rhs_name && strcmp(rhs_name, "none") != 0) {
        if (!rank)
            get_rhs_file_header(rhs_name, NULL, &(bm.d.nrhs), NULL);
        MPI_Bcast(&bm.d.nrhs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    abfield_init(d.ab);
    abfield_specify(d.ab, MPFQ_PRIME_MPZ, bw->p);

    bm.lingen_threshold = 10;
    bm.lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen_threshold", &(bm.lingen_threshold));
    param_list_parse_uint(pl, "display-threshold", &(display_threshold));
    param_list_parse_uint(pl, "lingen_mpi_threshold", &(bm.lingen_mpi_threshold));
    param_list_parse_uint(pl, "io-block-size", &(io_block_size));
    gmp_randseed_ui(rstate, bw->seed);
    if (bm.lingen_mpi_threshold < bm.lingen_threshold) {
        bm.lingen_mpi_threshold = bm.lingen_threshold;
        fprintf(stderr, "Argument fixing: setting lingen_mpi_threshold=%u (because lingen_threshold=%u)\n",
                bm.lingen_mpi_threshold, bm.lingen_threshold);
    }
    checkpoint_directory = param_list_lookup_string(pl, "checkpoint-directory");
    param_list_parse_uint(pl, "checkpoint-threshold", &checkpoint_threshold);
    param_list_parse_int(pl, "save_gathered_checkpoints", &save_gathered_checkpoints);



#if defined(FAKEMPI_H_)
    bm.lingen_mpi_threshold = UINT_MAX;
#endif

    /* }}} */

    /* TODO: we should rather use lingen_platform.
     */
    /* {{{ Parse MPI args. Make bm.com[0] a better mpi communicator */
    bm.mpi_dims[0] = 1;
    bm.mpi_dims[1] = 1;
    param_list_parse_intxint(pl, "mpi", bm.mpi_dims);
    {
        /* Display node index wrt MPI_COMM_WORLD */
        print_node_assignment(MPI_COMM_WORLD);

        /* Reorder all mpi nodes so that each node gets the given number
         * of jobs, but close together.
         */
        int mpi[2] = { bm.mpi_dims[0], bm.mpi_dims[1], };
        int thr[2] = {1,1};
#ifdef  HAVE_OPENMP
        if (param_list_parse_intxint(pl, "thr", thr)) {
            if (!rank)
                printf("# Limiting number of openmp threads to %d\n",
                        thr[0] * thr[1]);
            omp_set_num_threads(thr[0] * thr[1]);
        }
#endif

#ifdef  FAKEMPI_H_
        if (mpi[0]*mpi[1] > 1) {
            fprintf(stderr, "non-trivial option mpi= can't be used with fakempi. Please do an MPI-enabled build (MPI=1)\n");
            exit(EXIT_FAILURE);
        }
#endif
        if (!rank)
            printf("# size=%d mpi=%dx%d thr=%dx%d\n", size, mpi[0], mpi[1], thr[0], thr[1]);
        ASSERT_ALWAYS(size == mpi[0] * mpi[1]);
        if (bm.mpi_dims[0] != bm.mpi_dims[1]) {
            if (!rank)
                fprintf(stderr, "The current plingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
                    bm.mpi_dims[0], bm.mpi_dims[1]);
            abort();
        }
        int irank = rank / mpi[1];
        int jrank = rank % mpi[1];
        bm.com[0] = MPI_COMM_WORLD;
        /* MPI Api has some very deprecated prototypes */
        MPI_Comm_set_name(bm.com[0], (char*) "world");

        char commname[32];
        snprintf(commname, sizeof(commname), "row%d\n", irank);
        MPI_Comm_split(MPI_COMM_WORLD, irank, jrank, &(bm.com[1]));
        MPI_Comm_set_name(bm.com[1], commname);

        snprintf(commname, sizeof(commname), "col%d\n", jrank);
        MPI_Comm_split(MPI_COMM_WORLD, jrank, irank, &(bm.com[2]));
        MPI_Comm_set_name(bm.com[2], commname);

        print_node_assignment(bm.com[0]);
    }
    /* }}} */


    /* plingen tuning accepts some arguments. We look them up so as to
     * avoid failures down the line */
    plingen_tuning_lookup_parameters(pl);
    
    tree_stats::interpret_parameters(pl);
    logline_interpret_parameters(pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* TODO: read the a files in scattered mode */

    /* This is our estimate of the *global amount of memory* that the
     * program will use.
     */
    size_t exp_lenA = 2 + iceildiv(bm.d.m, bm.d.n);
    matpoly::memory_pool_guard main_memory(
            /* bm_io */
            abvec_elt_stride(bm.d.ab,bm.d.m*bm.d.n*exp_lenA) +   /* bm_io A */
            abvec_elt_stride(bm.d.ab,bm.d.m*bm.d.m) +   /* bm_io M */
            0);

    /* We now have a protected structure for a_reading task which does
     * the right thing concerning parallelism among MPI nodes (meaning
     * that non-root nodes essentially do nothing while the master job
     * does the I/O stuff) */
    bm_io aa(bm, afile, ffile, global_flag_ascii);
    aa.begin_read();
    aa.guess_length();

    /* run the mpi problem detection only if we're certain that we're at
     * least close to the ballpark where this sort of checks make sense.
     */
    if ((size_t) aa.guessed_length * (size_t) (bm.d.m + bm.d.n) * (size_t) abvec_elt_stride(bm.d.ab, 1) >= (1 << 28)) {
        check_for_mpi_problems();
    }

    {
        matpoly::memory_pool_guard blanket(SIZE_MAX);
        matpoly_ft::memory_pool_guard blanket_ft(SIZE_MAX);
        try {
            bm.hints = plingen_tuning(bm.d, aa.guessed_length, bm.com[0], pl);
        } catch (std::overflow_error const & e) {
            fputs(e.what(), stderr);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    if (global_flag_tune) {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    aa.compute_initial_F();

    int go_mpi = aa.guessed_length >= bm.lingen_mpi_threshold;

    if (go_mpi && !rank) {
        if (size) {
            printf("Expected length %u exceeds MPI threshold %u, going MPI now.\n", aa.guessed_length, bm.lingen_mpi_threshold);
        } else {
            printf("Expected length %u exceeds MPI threshold %u, but the process is not running in an MPI context.\n", aa.guessed_length, bm.lingen_mpi_threshold);
        }
    }

    if (go_mpi && size > 1) {
        matpoly_factory<bigmatpoly> F(bm.com, bm.mpi_dims[0], bm.mpi_dims[1]);
#ifdef ENABLE_MPI_LINGEN
        lingen_main_code(F, d.ab, aa);
#else
        /* The ENABLE_MPI_LINGEN flag should be turned on for a proper
         * MPI run.
         */
        ASSERT_ALWAYS(0);
#endif
    } else if (!rank) {
        /* We don't want to bother with memory problems in the non-mpi
         * case when the tuning was done for MPI: this is because the
         * per-transform ram was computed in the perspective of an MPI
         * run, and not for a plain run.
         */
        matpoly_factory<matpoly> F;
        if (size > 1) {
            matpoly_ft::memory_pool_guard blanket_ft(SIZE_MAX);
            lingen_main_code(F, d.ab, aa);
        } else {
            /* on the other hand, plain non-mpi code should benefit from
             * that safety net, since the tuning is expected to have
             * computed the needed ram correctly.
             */
            lingen_main_code(F, d.ab, aa);
        }
    } else {
        /* we have go_mpi == 0 and rank > 0 : all we have to do is
         * wait...
         */
    }

    if (!rank && random_input_length) {
        printf("t_basecase = %.2f\n", bm.t_basecase);
        printf("t_mp = %.2f\n", bm.t_mp);
        printf("t_mul = %.2f\n", bm.t_mul);
        printf("t_cp_io = %.2f\n", bm.t_cp_io);
        long peakmem = PeakMemusage();
        if (peakmem > 0)
            printf("# PeakMemusage (MB) = %ld (VmPeak: can be misleading)\n", peakmem >> 10);
    }

    abfield_clear(d.ab);
    if (ffile) free(ffile);

    gmp_randclear(rstate);

    bw_common_clear(bw);

    return rank0_exit_code;
}

/* vim:set sw=4 sta et: */
