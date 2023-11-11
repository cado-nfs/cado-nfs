/* Copyright (C) 1999--2007 Emmanuel Thom'e --- see LICENSE file */
#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <sys/param.h>
#include <cstdint>                        // for SIZE_MAX
#include <cmath>                          // for ceil
#include <cstdio>                         // for printf, fprintf, size_t
#include <cstdlib>                        // for exit, free, EXIT_FAILURE
#include <cstring>                        // for strcmp, memcpy, memset, strdup

#include <algorithm>                      // for fill
#include <fstream>                        // for basic_ostream::write
#include <memory>                         // for unique_ptr
#include <stdexcept>                      // for runtime_error, overflow_error
#include <string>                         // for operator+, string
#include <vector>                         // for vector

#include <sys/utsname.h>                  // for uname, utsname
#include <gmp.h>                          // for gmp_randclear, gmp_randinit...

#include "bw-common.h"                    // for bw, bw_common_clear, bw_com...
#include "fmt/core.h"                     // for check_format_string
#include "fmt/format.h"                   // for basic_buffer::append, basic...
#include "fmt/printf.h" // IWYU pragma: keep
#ifdef LINGEN_BINARY
#include "gf2x-fft.h"                     // for gf2x_cantor_fft_info
#include "gf2x-ternary-fft.h"             // for gf2x_ternary_fft_info
#else
#include "flint-fft/transform_interface.h"          // fft_transform_info
#endif
#include "arith-hard.hpp"             // for abfield_clear, abfield_init
#include "lingen_bigmatpoly.hpp"          // for bigmatpoly, bigmatpoly_model
#include "lingen_bigmatpoly_ft.hpp"       // for bigmatpoly_ft
#include "lingen_bmstatus.hpp"            // for bmstatus
#include "lingen_bw_dimensions.hpp"       // for bw_dimensions
#include "lingen_call_companion.hpp"      // for lingen_call_companion, ling...
#include "lingen_checkpoints.hpp"         // for save_checkpoint_file, load_...
#include "lingen_expected_pi_length.hpp"  // for expected_pi_length, expecte...
#include "lingen_hints.hpp"               // for lingen_hints
#include "lingen_io_matpoly.hpp"          // for lingen_io_matpoly_decl_usage
#include "lingen_io_wrappers.hpp"         // for lingen_output_wrapper_base
#include "lingen_matpoly_ft.hpp"          // for matpoly_ft, matpoly_ft<>::m...
#include "lingen_matpoly_select.hpp"      // for matpoly, matpoly::memory_guard
#include "lingen_qcode_select.hpp"           // for bw_lingen_basecase
#include "lingen_substep_schedule.hpp"    // for lingen_substep_schedule
#include "lingen_tuning.hpp"              // for lingen_tuning, lingen_tunin...
#include "logline.h"                      // for logline_end, logline_init_t...
#include "macros.h"                       // for ASSERT_ALWAYS, iceildiv, MIN
#include "memusage.h"                     // for PeakMemusage
#include "misc.h"                         // for size_disp
#include "omp_proxy.h"                    // for omp_set_num_threads
#include "params.h"                       // for cxx_param_list, param_list_...
#include "select_mpi.h"                   // for MPI_Comm, MPI_Comm_rank
#include "sha1.h"                         // for sha1_checksumming_stream
#include "timing.h"                       // for seconds
#include "tree_stats.hpp"                 // for tree_stats, tree_stats::tra...
#include "portability.h" // strdup // IWYU pragma: keep


/* If non-zero, then reading from A is actually replaced by reading from
 * a random generator */
static unsigned int random_input_length = 0;
static unsigned int input_length = 0;

static int split_input_file = 0;  /* unsupported ; do acollect by ourselves */
static int split_output_file = 0; /* do split by ourselves */

gmp_randstate_t rstate;

static int allow_zero_on_rhs = 0;

int rank0_exit_code = EXIT_SUCCESS;

int global_flag_ascii = 0;
int global_flag_tune = 0;

void lingen_decl_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "ascii",
            "read and write data in ascii");
    param_list_decl_usage(pl, "timings",
            "provide timings on all output lines");
    param_list_decl_usage(pl, "tune",
            "activate tuning mode");
    param_list_decl_usage(pl, "allow_zero_on_rhs",
            "do not cry if the generator corresponds to a zero contribution on the RHS vectors");

    /* we must be square ! */
    param_list_decl_usage(pl, "mpi", "number of MPI nodes across which the execution will span, with mesh dimensions");
    param_list_decl_usage(pl, "thr", "number of threads (on each node) for the program, with mesh dimensions");

    param_list_decl_usage(pl, "nrhs",
            "number of columns that correspond to rhs vectors");
    param_list_decl_usage(pl, "rhs",
            "file with rhs vectors (only the header is read)");

    param_list_decl_usage(pl, "afile",
            "input sequence file");
    param_list_decl_usage(pl, "input_length",
            "input sequence length (defaults to auto-detect)");
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

#if 0
    param_list_decl_usage(pl, "lingen_mpi_threshold",
            "use MPI matrix operations above this size");
    param_list_decl_usage(pl, "lingen_threshold",
            "use recursive algorithm above this size");
#endif

    param_list_configure_switch(pl, "--tune", &global_flag_tune);
    param_list_configure_switch(pl, "--ascii", &global_flag_ascii);
    param_list_configure_alias(pl, "seed", "random_seed");

}/*}}}*/

/**********************************************************************/

/*{{{ Main entry points and recursive algorithm (with and without MPI) */

/* Forward declaration, it's used by the recursive version */
matpoly bw_lingen_single(bmstatus & bm, matpoly & E);

bigmatpoly bw_biglingen_collective(bmstatus & bm, bigmatpoly & E);

std::string sha1sum(matpoly const & X)
{
    sha1_checksumming_stream S;
    /*
    if (X.is_tight())
        S.write((const char *) X.data_area(), X.data_size_in_bytes());
        */
#if GMP_LIMB_BITS == 32
    size_t nbytes_buffer = iceildiv(X.data_entry_size_in_bytes(), 8);
    char buffer[nbytes_buffer];
    std::fill_n(buffer, nbytes_buffer, 0);
#endif
    for(unsigned int i = 0 ; i < X.nrows() ; i++)
        for(unsigned int j = 0 ; j < X.ncols() ; j++) {
#if GMP_LIMB_BITS == 32
            std::copy_n(reinterpret_cast<const char *>(X.part(i, j)), X.data_entry_size_in_bytes(), buffer);
            S.write(buffer, nbytes_buffer);
#else
            S.write((const char *) X.part(i, j), X.data_entry_size_in_bytes());
#endif
        }
    char checksum[41];
    S.checksum(checksum);
    return std::string(checksum);
}

template<typename matpoly_type, typename fft_type>
struct matching_ft_type {};

template<typename matpoly_type>
struct matpoly_diverter {};

template<typename fft_type>
struct matching_ft_type<matpoly, fft_type> {
    typedef matpoly_ft<fft_type> type;
};
template<> struct matpoly_diverter<matpoly> {
    static constexpr const char * prefix = "";
    static matpoly callback(bmstatus & bm, matpoly & E) {
        return bw_lingen_single(bm, E);
    }
};
constexpr const char * matpoly_diverter<matpoly>::prefix;

template<typename fft_type>
struct matching_ft_type<bigmatpoly, fft_type> {
    typedef bigmatpoly_ft<fft_type> type;
};
template<> struct matpoly_diverter<bigmatpoly> {
    static constexpr const char * prefix = "MPI-";
    static bigmatpoly callback(bmstatus & bm, bigmatpoly & E) {
        return bw_biglingen_collective(bm, E);
    }
};
constexpr const char * matpoly_diverter<bigmatpoly>::prefix;


template<typename matpoly_type>
matpoly_type generic_mp(matpoly_type & E, matpoly_type & pi_left, bmstatus & bm, lingen_call_companion & C)
{
    switch (C.mp.S.fft_type) {
        case lingen_substep_schedule::FFT_NONE:
            return matpoly_type::mp(bm.stats, E, pi_left, &C.mp);
        case lingen_substep_schedule::FFT_FLINT:
#ifndef LINGEN_BINARY
            return matching_ft_type<matpoly_type,
                    fft_transform_info>::type::mp_caching(
                            bm.stats, E, pi_left, & C.mp);
#else
            throw std::runtime_error("fft type \"flint\" does not make sense here");
#endif
#ifdef LINGEN_BINARY
        case lingen_substep_schedule::FFT_CANTOR:
            return matching_ft_type<matpoly_type,
                    gf2x_cantor_fft_info>::type::mp_caching(
                            bm.stats, E, pi_left, & C.mp);
        case lingen_substep_schedule::FFT_TERNARY:
            return matching_ft_type<matpoly_type,
                    gf2x_ternary_fft_info>::type::mp_caching(
                            bm.stats, E, pi_left, & C.mp);
#else
        case lingen_substep_schedule::FFT_TERNARY:
        case lingen_substep_schedule::FFT_CANTOR:
            throw std::runtime_error("fft types over GF(2)[x] do not make sense here");
#endif
    }
    throw std::runtime_error("invalid fft_type");
}

template<typename matpoly_type>
matpoly_type generic_mul(matpoly_type & pi_left, matpoly_type & pi_right, bmstatus & bm, lingen_call_companion & C)
{
    switch (C.mul.S.fft_type) {
        case lingen_substep_schedule::FFT_NONE:
            return matpoly_type::mul(bm.stats, pi_left, pi_right, & C.mul);
        case lingen_substep_schedule::FFT_FLINT:
#ifndef LINGEN_BINARY
            return matching_ft_type<matpoly_type,
                    fft_transform_info>::type::mul_caching(
                            bm.stats, pi_left, pi_right, & C.mul);
#else
            throw std::runtime_error("fft type \"flint\" does not make sense here");
#endif
#ifdef LINGEN_BINARY
        case lingen_substep_schedule::FFT_CANTOR:
            return matching_ft_type<matpoly_type,
                    gf2x_cantor_fft_info>::type::mul_caching(
                            bm.stats, pi_left, pi_right, & C.mul);
        case lingen_substep_schedule::FFT_TERNARY:
            return matching_ft_type<matpoly_type,
                    gf2x_ternary_fft_info>::type::mul_caching(
                            bm.stats, pi_left, pi_right, & C.mul);
#else
        case lingen_substep_schedule::FFT_TERNARY:
        case lingen_substep_schedule::FFT_CANTOR:
            throw std::runtime_error("fft types over GF(2)[x] do not make sense here");
#endif
    }
    throw std::runtime_error("invalid fft_type");
}

template<typename matpoly_type>
void truncate_overflow(bmstatus & bm, matpoly_type & pi, unsigned int pi_expect)
{
    if (pi.get_size() > pi_expect) {
        unsigned int nluck = 0;
        unsigned int maxdelta_lucky = 0;
        for(unsigned int j = 0 ; j < bm.d.m + bm.d.n ; j++) {
            if (!bm.lucky[j]) continue;
            nluck++;
            if (bm.delta[j] >= maxdelta_lucky)
                maxdelta_lucky = bm.delta[j];
        }
        printf("truncating excess cols\n"
                "   pi has length %zu, we expected %u at most.\n"
                "   max delta on the %u lucky columns: %u\n",
                pi.get_size(),
                pi_expect,
                nluck,
                maxdelta_lucky);
        bm.display_deltas();
        pi.set_size(pi_expect);
        pi.clear_high_word();
        /* in fact, there is really no reason to reach here. The expected
         * pi length should be large enough so that the "luck" condition
         * and eventually the "done"/"finished" flags are triggered
         * before we reach the point where we overflow. This means a
         * constant offset. Therefore, reaching here is a bug.
         *
         * Furthermore, while the idea of truncating is appealing, it
         * doesn't quote solve the problem when we wish to keep iterating
         * a bit more: because of the degree constraints, truncating will
         * make the next E matrix cancel, as we're almost knowingly
         * killing its rank...
         *
         * one of the better things to do could be to keep track of a
         * shift vector. but again, what for...
         */
        ASSERT_ALWAYS(0);
    }
}


template<typename matpoly_type>
matpoly_type bw_lingen_recursive(bmstatus & bm, matpoly_type & E) /*{{{*/
{
    int depth = bm.depth();
    size_t z = E.get_size();

    /* C0 is a copy. We won't use it for long anyway. We'll take a
     * reference _later_ */
    lingen_call_companion C0 = bm.companion(depth, z);

    tree_stats::sentinel dummy(bm.stats, fmt::sprintf("%srecursive", matpoly_diverter<matpoly_type>::prefix), z, C0.total_ncalls);

    bm.stats.plan_smallstep(C0.mp.step_name(), C0.mp.tt);
    bm.stats.plan_smallstep(C0.mul.step_name(), C0.mul.tt);

    bw_dimensions & d = bm.d;

    /* we have to start with something large enough to get all
     * coefficients of E_right correct */
    size_t half = E.get_size() - (E.get_size() / 2);
    // unsigned int pi_expect = expected_pi_length(d, bm.delta, E.get_size());
    unsigned int pi_expect_lowerbound = expected_pi_length_lowerbound(d, E.get_size());
    unsigned int pi_left_expect = expected_pi_length(d, bm.delta, half);
    unsigned int pi_left_expect_lowerbound = expected_pi_length_lowerbound(d, half);
    unsigned int pi_left_expect_used_for_shift = MIN(pi_left_expect, half + 1);


    matpoly_type E_left = E.truncate_and_rshift(half, half + 1 - pi_left_expect_used_for_shift);

    // this consumes E_left entirely.
    matpoly_type pi_left = matpoly_diverter<matpoly_type>::callback(bm, E_left);

    ASSERT_ALWAYS(pi_left.get_size());

    if (bm.done) {
        return pi_left;
    }

    truncate_overflow(bm, pi_left, pi_left_expect);

    ASSERT_ALWAYS(bm.done || pi_left.get_size() >= pi_left_expect_lowerbound);

    /* XXX I don't understand why I need to do this. It seems to me that
     * MP(XA, B) and MP(A, B) should be identical whenever deg A > deg B.
     */
    ASSERT_ALWAYS(pi_left_expect_used_for_shift >= pi_left.get_size());
    if (pi_left_expect_used_for_shift != pi_left.get_size()) {
        E.rshift(E, pi_left_expect_used_for_shift - pi_left.get_size());
        /* Don't shrink_to_fit at this point, because we've only made a
         * minor adjustment. */
    }

    matpoly_type E_right = E.similar_shell();

    if (!load_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E_right, bm.t, bm.t+E.get_size()-pi_left.get_size()+1)) {
        logline_begin(stdout, z, "t=%u %*s%sMP(%s, %zu, %zu) -> %zu",
                bm.t, depth,"", matpoly_diverter<matpoly_type>::prefix,
                C0.mp.fft_name(),
                E.get_size(), pi_left.get_size(), E.get_size() - pi_left.get_size() + 1);

        E_right = generic_mp(E, pi_left, bm, bm.companion(depth, z));
        logline_end(&bm.t_mp, "");
    }
    E.clear();


    unsigned int pi_right_expect = expected_pi_length(d, bm.delta, E_right.get_size());
    unsigned int pi_right_expect_lowerbound = expected_pi_length_lowerbound(d, E_right.get_size());

    matpoly_type pi_right = matpoly_diverter<matpoly_type>::callback(bm, E_right);

    truncate_overflow(bm, pi_right, pi_right_expect);

    ASSERT_ALWAYS(bm.done || pi_right.get_size() >= pi_right_expect_lowerbound);

    logline_begin(stdout, z, "t=%u %*s%sMUL(%s, %zu, %zu) -> %zu",
            bm.t, depth, "", matpoly_diverter<matpoly_type>::prefix,
            C0.mul.fft_name(),
            pi_left.get_size(), pi_right.get_size(), pi_left.get_size() + pi_right.get_size() - 1);

    matpoly_type pi = generic_mul(pi_left, pi_right, bm, bm.companion(depth, z));
    pi_left.clear();
    pi_right.clear();

    /* Note that the leading coefficients of pi_left and pi_right are not
     * necessarily full-rank, so that we have to fix potential zeros. If
     * we don't, the degree of pi artificially grows with the recursive
     * level.
     */
    unsigned int pisize = pi.get_size();
#if 0
    /* The asserts below are bogus. In fact, it's not entirely impossible
     * that pi grows more than what we had expected on entry, e.g. if we
     * have one early generator. So we can't just do this. Most of the
     * time it will work, but we can't claim that it will always work.
     *
     * One possible sign is when the entry deltas are somewhat even, and
     * the result deltas are unbalanced.
     *
     * It might also be worthwhile to do this check by column, and look
     * at both degree and valuation.
     */
    for(; pisize > pi_expect ; pisize--) {
        /* These coefficients really must be zero */
        ASSERT_ALWAYS(pi.coeff_is_zero(pisize - 1));
    }
    ASSERT_ALWAYS(pisize <= pi_expect);
#endif
    /* Now below pi_expect, it's not impossible to have a few
     * cancellations as well.
     */
    for(; pisize ; pisize--) {
        if (!pi.coeff_is_zero(pisize - 1)) break;
    }
    pi.set_size(pisize);
    ASSERT_ALWAYS(bm.done || pisize >= pi_expect_lowerbound);

    logline_end(&bm.t_mul, "");

    return pi;
}/*}}}*/


matpoly bw_lingen_single_nocp(bmstatus & bm, matpoly & E) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(!rank);
    matpoly pi;

    lingen_call_companion C = bm.companion(bm.depth(), E.get_size());

    // ASSERT_ALWAYS(E.size < bm.lingen_mpi_threshold);

    // fprintf(stderr, "Enter %s\n", __func__);
    if (!bm.recurse(E.get_size())) {
        tree_stats::transition_sentinel dummy(bm.stats, "recursive_threshold", E.get_size(), C.total_ncalls);
        bm.t_basecase -= seconds();
        E.clear_high_word();
        pi = bw_lingen_basecase(bm, E);
        bm.t_basecase += seconds();
    } else {
        pi = bw_lingen_recursive(bm, E);
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    return pi;
}/*}}}*/
matpoly bw_lingen_single(bmstatus & bm, matpoly & E) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(!rank);
    unsigned int t0 = bm.t;
    unsigned int t1 = bm.t + E.get_size();
    matpoly pi;

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E, t0, t1);

    if (load_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1))
        return pi;

    pi = bw_lingen_single_nocp(bm, E);

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1);

    return pi;
}/*}}}*/

bigmatpoly bw_biglingen_collective(bmstatus & bm, bigmatpoly & E)/*{{{*/
{
    /* as for bw_lingen_single, we're tempted to say that we're just a
     * trampoline. In fact, it's not really satisfactory: we're really
     * doing stuff here. In a sense though, it's not *that much* of a
     * trouble, because the mpi threshold will be low enough that doing
     * our full job here is not too much of a problem.
     */
    bw_dimensions & d = bm.d;
    matpoly::arith_hard * ab = & d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    int rank;
    int size;
    MPI_Comm_rank(bm.com[0], &rank);
    MPI_Comm_size(bm.com[0], &size);
    unsigned int t0 = bm.t;
    unsigned int t1 = bm.t + E.get_size();
    bigmatpoly_model const& model(E);
    int depth = bm.depth();
    size_t z = E.get_size();

    lingen_call_companion C = bm.companion(depth, z);
    bool go_mpi = C.go_mpi();
    // bool go_mpi = E.get_size() >= bm.lingen_mpi_threshold;

    bigmatpoly pi(model);

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E, t0, t1);

    if (load_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1))
        return pi;

    // fprintf(stderr, "Enter %s\n", __func__);
    if (go_mpi) {
        pi = bw_lingen_recursive(bm, E);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        tree_stats::transition_sentinel dummy(bm.stats, "mpi_threshold", E.get_size(), C.total_ncalls);

        matpoly sE(ab, m, b, E.get_size());
        matpoly spi;

        double expect0 = bm.hints.tt_gather_per_unit * E.get_size();
        bm.stats.plan_smallstep("gather(L+R)", expect0);
        bm.stats.begin_smallstep("gather(L+R)");
        E.gather_mat(sE);
        E = bigmatpoly(model);
        bm.stats.end_smallstep();

        /* Only the master node does the local computation */
        if (!rank)
            spi = bw_lingen_single_nocp(bm, sE);

        double expect1 = bm.hints.tt_scatter_per_unit * z;
        bm.stats.plan_smallstep("scatter(L+R)", expect1);
        bm.stats.begin_smallstep("scatter(L+R)");
        pi = bigmatpoly(ab, model, b, b, 0);
        pi.scatter_mat(spi);
        MPI_Bcast(&bm.done, 1, MPI_INT, 0, bm.com[0]);
        MPI_Bcast(bm.delta.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
        MPI_Bcast(bm.lucky.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
        MPI_Bcast(&(bm.t), 1, MPI_UNSIGNED, 0, bm.com[0]);
        /* Don't forget to broadcast delta from root node to others ! */
        bm.stats.end_smallstep();
    }
    // fprintf(stderr, "Leave %s\n", __func__);

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1);

    MPI_Barrier(bm.com[0]);

    return pi;
}/*}}}*/

/*}}}*/

/**********************************************************************/

/**********************************************************************/
unsigned int count_lucky_columns(bmstatus & bm)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    unsigned int b = m + n;
    int luck_mini = expected_pi_length(d);
    MPI_Bcast(bm.lucky.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
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
            unsigned int s = (j * 1009) % (m+n);
            bm.lucky[s]  = expected_pi_length(d);
            bm.delta[s] -= expected_pi_length(d);
        }
        return check_luck_condition(bm);
    }

    return 0;
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

/* We don't have a header file for this one */
extern "C" void check_for_mpi_problems();

int wrapped_main(int argc, char *argv[])
{

    cxx_param_list pl;

    bw_common_decl_usage(pl);
    lingen_decl_usage(pl);
    logline_decl_usage(pl);
    lingen_tuning_decl_usage(pl);
    lingen_checkpoint::decl_usage(pl);
    lingen_io_matpoly_decl_usage(pl);
    tree_stats::declare_usage(pl);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    /* {{{ interpret our parameters */
    gmp_randinit_default(rstate);

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    param_list_parse_int(pl, "allow_zero_on_rhs", &allow_zero_on_rhs);
    param_list_parse_uint(pl, "random-input-with-length", &random_input_length);
    param_list_parse_uint(pl, "input-length", &input_length);
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

    bmstatus bm(bw->m, bw->n, bw->p);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    if (!global_flag_tune && !random_input_length) {
        if (!rhs_name) {
            fprintf(stderr, "# When using lingen, you must either supply --random-input-with-length, or provide a rhs, or possibly provide rhs=none\n");
        } else if (strcmp(rhs_name, "none") == 0) {
            rhs_name = NULL;
        }
    }
    if (param_list_parse_uint(pl, "nrhs", &(bm.d.nrhs)) && rhs_name) {
        fprintf(stderr, "# the command line arguments rhs= and nrhs= are incompatible\n");
        exit(EXIT_FAILURE);
    }
    if (rhs_name && strcmp(rhs_name, "none") != 0) {
        if (!rank)
            get_rhs_file_header(rhs_name, NULL, &(bm.d.nrhs), NULL);
        MPI_Bcast(&bm.d.nrhs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    gmp_randseed_ui(rstate, bw->seed);

#if 0
    bm.lingen_threshold = 10;
    bm.lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen_threshold", &(bm.lingen_threshold));
    param_list_parse_uint(pl, "lingen_mpi_threshold", &(bm.lingen_mpi_threshold));
    if (bm.lingen_mpi_threshold < bm.lingen_threshold) {
        bm.lingen_mpi_threshold = bm.lingen_threshold;
        fprintf(stderr, "Argument fixing: setting lingen_mpi_threshold=%u (because lingen_threshold=%u)\n",
                bm.lingen_mpi_threshold, bm.lingen_threshold);
    }


#if defined(FAKEMPI_H_)
    bm.lingen_mpi_threshold = UINT_MAX;
#endif
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
            if (omp_get_max_threads() >= thr[0] * thr[1]) {
                if (!rank)
                    printf("# Limiting number of openmp threads to %d\n",
                            thr[0] * thr[1]);
                omp_set_num_threads(thr[0] * thr[1]);
            } else {
                if (!rank)
                    printf("# Number of openmp threads is capped at %d"
                            ", which is below thr=%dx%d. Keeping as it is\n",
                            omp_get_max_threads(),
                            thr[0], thr[1]);
            }
        }
#else
        if (param_list_parse_intxint(pl, "thr", thr)) {
            if (thr[0]*thr[1] != 1) {
                if (!rank) {
                    fprintf(stderr, "This program only wants openmp for multithreading. Ignoring thr argument.\n");
                }
                param_list_add_key(pl, "thr", "1x1", PARAMETER_FROM_CMDLINE);
            }
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
                fprintf(stderr, "The current lingen code is limited to square splits ; here, we received a %d x %d split, which will not work\n",
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

        constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
        if ((bm.d.m + bm.d.n) / simd < (unsigned int) mpi[0]) {
            printf("########################################################\n");
            printf("# Warning: this run will leave some resources idle:\n"
                   "# the matrices of size %u*%u and %u*%u can be split into\n"
                   "# chunks of minimal size %u, whence an mpi split over %d*%d is useless\n",
                   bm.d.m,
                   bm.d.m + bm.d.n,
                   bm.d.m + bm.d.n,
                   bm.d.m + bm.d.n,
                   simd,
                   mpi[0],
                   mpi[0]);
            printf("########################################################\n");
        }
    }
    /* }}} */

    /* lingen tuning accepts some arguments. We look them up so as to
     * avoid failures down the line */
    lingen_tuning_lookup_parameters(pl);
    
    tree_stats::interpret_parameters(pl);
    logline_interpret_parameters(pl);
    lingen_checkpoint::interpret_parameters(pl);
    lingen_io_matpoly_interpret_parameters(pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* TODO: read the a files in scattered mode */

    /* Don't police memory right now, we don't care */
    matpoly::memory_guard main_memory(SIZE_MAX);

    std::unique_ptr<lingen_input_wrapper_base> A_series;

    if (random_input_length) {
        A_series = std::unique_ptr<lingen_input_wrapper_base>(new lingen_random_input(&bm.d.ab, bm.d.m, bm.d.n, rstate, random_input_length));
    } else {
        A_series = std::unique_ptr<lingen_input_wrapper_base>(new lingen_file_input(&bm.d.ab, bm.d.m, bm.d.n, afile, global_flag_ascii, input_length));
    }

#ifdef LINGEN_BINARY
#define K_elts_to_bytes(x)      (iceildiv((x),ULONG_BITS) * sizeof(unsigned long))
#else
#define K_elts_to_bytes(x)      (bm.d.ab.vec_elt_stride((x)))
#endif

    /* run the mpi problem detection only if we're certain that we're at
     * least close to the ballpark where this sort of checks make sense.
     */
    if (K_elts_to_bytes((size_t) A_series->guessed_length() * (size_t) (bm.d.m + bm.d.n)) >= (1 << 28)) {
        check_for_mpi_problems();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* This will cause the initial read */
    std::unique_ptr<lingen_E_from_A> E_series = std::unique_ptr<lingen_E_from_A>(new lingen_E_from_A(bm.d, *A_series));

    bm.t = E_series->t0;

    size_t L = E_series->guessed_length();

    {
        matpoly::memory_guard blanket(SIZE_MAX);
#ifndef LINGEN_BINARY
        typename matpoly_ft<fft_transform_info>::memory_guard blanket_ft(SIZE_MAX);
#else
        typename matpoly_ft<gf2x_cantor_fft_info>::memory_guard blanket_ft(SIZE_MAX);
        typename matpoly_ft<gf2x_ternary_fft_info>::memory_guard blanket_ft2(SIZE_MAX);
#endif
        try {
            bm.hints = lingen_tuning(bm.d, L - bm.t, bm.com[0], pl);
        } catch (std::overflow_error const & e) {
            fputs(e.what(), stderr);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    if (global_flag_tune) {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    size_t safe_guess = global_flag_ascii ? ceil(1.05 * L) : L;

    /* c0 is (1+m/n) times the input size */
    size_t c0 = K_elts_to_bytes(
                iceildiv(bm.d.m + bm.d.n, bm.mpi_dims[0]) *
                iceildiv(bm.d.m + bm.d.n, bm.mpi_dims[1]) *
                iceildiv(bm.d.m*safe_guess, bm.d.m+bm.d.n));
    matpoly::memory_guard main(2*c0);

    if (!rank) {
        char buf[20];
        printf("# Estimated memory for JUST transforms (per node): %s\n",
                size_disp(2*c0, buf));
        printf("# Estimated peak total memory (per node): max at depth %d: %s\n",
                bm.hints.ipeak,
                size_disp(bm.hints.peak, buf));
    }

    int go_mpi = bm.companion(0, L).go_mpi();

    if (go_mpi) {
        if (!rank) {
            if (size > 1) {
                printf("Expected length %zu exceeds MPI threshold,"
                       " going MPI now.\n",
                       L);
            } else {
                printf("Expected length %zu exceeds MPI threshold, "
                       "but the process is not running in an MPI context.\n",
                       L);
            }
        }
        MPI_Barrier(bm.com[0]);
    }

    std::unique_ptr<lingen_output_wrapper_base> Fdst;
    std::unique_ptr<lingen_output_wrapper_base> Fdst_rhs;
    
    if (random_input_length) {
        Fdst = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_sha1sum(& bm.d.ab, bm.d.n, bm.d.n, "F"));
        Fdst_rhs = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_sha1sum(& bm.d.ab, bm.d.nrhs, bm.d.n, "Frhs"));
    } else if (split_output_file) {
        std::string pattern = ffile;
        Fdst = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_splitfile(& bm.d.ab, bm.d.n, bm.d.n, pattern + ".sols{2}-{3}.{0}-{1}", global_flag_ascii));
        Fdst_rhs = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_splitfile(& bm.d.ab, bm.d.nrhs, bm.d.n, pattern + ".sols{2}-{3}.{0}-{1}.rhs", global_flag_ascii));
    } else {
        Fdst = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_singlefile(& bm.d.ab, bm.d.n, bm.d.n, ffile, global_flag_ascii));
        Fdst_rhs = std::unique_ptr<lingen_output_wrapper_base>(new lingen_output_to_singlefile(& bm.d.ab, bm.d.nrhs, bm.d.n, std::string(ffile) + ".rhs", global_flag_ascii));
    }

    if (go_mpi && size > 1) {
        bigmatpoly_model model(bm.com, bm.mpi_dims[0], bm.mpi_dims[1]);
        bigmatpoly E(& bm.d.ab, model, bm.d.m, bm.d.m + bm.d.n, safe_guess);
        lingen_scatter<bigmatpoly> fill_E(E);
        lingen_F0 F0 = *E_series;
        pipe(*E_series, fill_E, "Read");
        bm.delta.assign(bm.d.m + bm.d.n, F0.t0);
        logline_init_timer();
        bigmatpoly pi = bw_biglingen_collective(bm, E);
        bm.stats.final_print();
        bm.display_deltas();
        if (!rank) printf("(pi.alloc = %zu)\n", pi.my_cell().capacity());
        constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
        pi.zero_pad(simd * iceildiv(pi.get_size(), simd));
        if (check_luck_condition(bm)) {
            lingen_gather_reverse<bigmatpoly> read_PI(pi);
            lingen_F_from_PI Fsrc(bm, read_PI, F0);
            pipe(Fsrc, *Fdst, "Written", true);
            Fsrc.write_rhs(*Fdst_rhs);
        }
    } else if (!rank) {
        /* We do this only in the rank==0 case, since we have really
         * nothing to do at the other ranks.
         */

        /* We don't want to bother with memory problems in the non-mpi
         * case when the tuning was done for MPI: this is because the
         * per-transform ram was computed in the perspective of an MPI
         * run, and not for a plain run.
         */
#ifndef LINGEN_BINARY
        typename matpoly_ft<fft_transform_info>::memory_guard blanket_ft(SIZE_MAX);
#else
        typename matpoly_ft<gf2x_cantor_fft_info>::memory_guard blanket_ft(SIZE_MAX);
        typename matpoly_ft<gf2x_ternary_fft_info>::memory_guard blanket_ft2(SIZE_MAX);
#endif
        matpoly E(& bm.d.ab, bm.d.m, bm.d.m + bm.d.n, safe_guess);
        lingen_scatter<matpoly> fill_E(E);
        lingen_F0 F0 = *E_series;
        pipe(*E_series, fill_E, "Read");
        bm.delta.assign(bm.d.m + bm.d.n, F0.t0);
        logline_init_timer();
        matpoly pi = bw_lingen_single(bm, E);
        bm.stats.final_print();
        bm.display_deltas();
        if (!rank) printf("(pi.alloc = %zu)\n", pi.capacity());
        constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
        pi.zero_pad(simd * iceildiv(pi.get_size(), simd));
        if (check_luck_condition(bm)) {
            lingen_gather_reverse<matpoly> read_PI(pi);
            lingen_F_from_PI Fsrc(bm, read_PI, F0);
            pipe(Fsrc, *Fdst, "Written", true);
            Fsrc.write_rhs(*Fdst_rhs);
        }
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

    if (ffile) free(ffile);

    gmp_randclear(rstate);

    return 0;   // ignored.
}

/* We do this so that the dtors of the data that gets allocated within
 * main are allowed to use MPI_Comm_rank.
 */
// coverity[root_function]
int main(int argc, char *argv[])
{
#ifdef  HAVE_OPENMP
    if (getenv("OMP_DYNAMIC") == NULL) {
        /* Change the default behavior with respect to dynamic thread
         * allocation, but do it with *lower* priority than the
         * environment variable (so that the possibility of changing the
         * behavior at runtime is retained).
         */
        omp_set_dynamic(true);
    }
#endif

    bw_common_init(bw, &argc, &argv);
    wrapped_main(argc, argv);
    bw_common_clear(bw);
    return rank0_exit_code;
}

/* vim:set sw=4 sta et: */
