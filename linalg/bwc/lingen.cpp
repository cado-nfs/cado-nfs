/* Copyright (C) 1999--2007 Emmanuel Thom'e --- see LICENSE file */
#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <sys/param.h>

#include <cstddef>

#include <stdexcept>                      // for runtime_error, overflow_error
#include <string>                         // for operator+, string
#include <vector>                         // for vector

#include <sys/utsname.h>                  // for uname, utsname
#include <gmp.h>

#include "fmt/format.h"                   // for basic_buffer::append, basic...

#ifdef LINGEN_BINARY
#include "gf2x-fft.h"                     // for gf2x_cantor_fft_info
#include "gf2x-ternary-fft.h"             // for gf2x_ternary_fft_info
#else
#include "flint-fft/transform_interface.h"          // fft_transform_info
#endif
#include "lingen_bigmatpoly.hpp"          // for bigmatpoly, bigmatpoly_model
#include "lingen_bigmatpoly_ft.hpp"       // for bigmatpoly_ft
#include "lingen_bmstatus.hpp"            // for bmstatus
#include "lingen_bw_dimensions.hpp"       // for bw_dimensions
#include "lingen_call_companion.hpp"      // for lingen_call_companion, ling...
#include "lingen_checkpoints.hpp"         // for save_checkpoint_file, load_...
#include "lingen_expected_pi_length.hpp"  // for expected_pi_length, expecte...
#include "lingen_hints.hpp"               // for lingen_hints
#include "lingen_io_wrappers.hpp"         // for lingen_output_wrapper_base
#include "lingen_matpoly_ft.hpp"          // for matpoly_ft, matpoly_ft<>::m...
#include "lingen_matpoly_select.hpp"      // for matpoly, matpoly::memory_guard
#include "lingen_qcode_select.hpp"           // for bw_lingen_basecase
#include "lingen_substep_schedule.hpp"    // for lingen_substep_schedule
#include "lingen_tuning.hpp"              // for lingen_tuning, lingen_tunin...
#include "logline.hpp"                      // for logline_end, logline_init_t...
#include "macros.h"                       // for ASSERT_ALWAYS, iceildiv, MIN
#include "select_mpi.h"                   // for MPI_Comm, MPI_Comm_rank
#include "sha1.h"                         // for sha1_checksumming_stream
#include "timing.h"                       // for seconds
#include "tree_stats.hpp"                 // for tree_stats, tree_stats::tra...

// some of our entry points are exported for tests
#include "lingen.hpp"


/**********************************************************************/

/* Main entry points and recursive algorithm (with and without MPI) */

// used for debugging
// NOLINTNEXTLINE(misc-use-internal-linkage)
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
static matpoly_type generic_mp(matpoly_type & E, matpoly_type & pi_left, bmstatus & bm, lingen_call_companion & C)
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
static matpoly_type generic_mul(matpoly_type & pi_left, matpoly_type & pi_right, bmstatus & bm, lingen_call_companion & C)
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
static void truncate_overflow(bmstatus & bm, matpoly_type & pi, unsigned int pi_expect)
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
        fmt::print("truncating excess cols\n"
                "   pi has length {}, we expected {} at most.\n"
                "   max delta on the {} lucky columns: {}\n",
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
         * doesn't quite solve the problem when we wish to keep iterating
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
static matpoly_type bw_lingen_recursive(bmstatus & bm, matpoly_type & E) /*{{{*/
{
    int const depth = bm.depth;
    size_t const z = E.get_size();

    /* C0 is a copy. We won't use it for long anyway. We'll take a
     * reference _later_ */
    lingen_call_companion const C0 = bm.companion(depth, z);

    tree_stats::sentinel const dummy(bm.stats, fmt::format("{}recursive", matpoly_diverter<matpoly_type>::prefix), z, C0.total_ncalls);
    bmstatus::depth_sentinel ddummy(bm);

    bm.stats.plan_smallstep(C0.mp.step_name(), C0.mp.tt);
    bm.stats.plan_smallstep(C0.mul.step_name(), C0.mul.tt);

    bw_dimensions & d = bm.d;

    /* we have to start with something large enough to get all
     * coefficients of E_right correct */
    size_t const half = E.get_size() - (E.get_size() / 2);
    // unsigned int pi_expect = expected_pi_length(d, bm.delta, E.get_size());
    unsigned int const pi_expect_lowerbound = expected_pi_length_lowerbound(d, E.get_size());
    unsigned int const pi_left_expect = expected_pi_length(d, bm.delta, half);
    unsigned int const pi_left_expect_lowerbound = expected_pi_length_lowerbound(d, half);
    unsigned int const pi_left_expect_used_for_shift = MIN(pi_left_expect, half + 1);


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


    unsigned int const pi_right_expect = expected_pi_length(d, bm.delta, E_right.get_size());
    unsigned int const pi_right_expect_lowerbound = expected_pi_length_lowerbound(d, E_right.get_size());

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


static matpoly bw_lingen_single_nocp(bmstatus & bm, matpoly & E) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(!rank);
    matpoly pi;

    lingen_call_companion const C = bm.companion(bm.depth, E.get_size());

    // ASSERT_ALWAYS(E.size < bm.lingen_mpi_threshold);

    // fmt::print(stderr, "Enter {}\n", __func__);
    if (!bm.recurse(E.get_size())) {
        tree_stats::transition_sentinel const dummy(bm.stats, "recursive_threshold", E.get_size(), C.total_ncalls);
        bm.t_basecase -= seconds();
        E.clear_high_word();
        pi = bw_lingen_basecase(bm, E);
        bm.t_basecase += seconds();
    } else {
        pi = bw_lingen_recursive(bm, E);
    }
    // fmt::print(stderr, "Leave {}\n", __func__);

    return pi;
}/*}}}*/
matpoly bw_lingen_single(bmstatus & bm, matpoly & E) /*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(!rank);
    unsigned int const t0 = bm.t;
    unsigned int const t1 = bm.t + E.get_size();
    matpoly pi;

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E, t0, t1);

    if (load_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1))
        return pi;

    // bm.display_deltas();
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
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const b = m + n;
    int rank;
    int size;
    MPI_Comm_rank(bm.com[0], &rank);
    MPI_Comm_size(bm.com[0], &size);
    unsigned int const t0 = bm.t;
    unsigned int const t1 = bm.t + E.get_size();
    bigmatpoly_model const& model(E);
    int const depth = bm.depth;
    size_t const z = E.get_size();

    lingen_call_companion const C = bm.companion(depth, z);
    bool const go_mpi = C.go_mpi();
    // bool go_mpi = E.get_size() >= bm.lingen_mpi_threshold;

    bigmatpoly pi(model);

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_E, E, t0, t1);

    if (load_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1))
        return pi;

    // fmt::print(stderr, "Enter {}\n", __func__);
    if (go_mpi) {
        pi = bw_lingen_recursive(bm, E);
    } else {
        /* Fall back to local code */
        /* This entails gathering E locally, computing pi locally, and
         * dispathing it back. */

        tree_stats::transition_sentinel const dummy(bm.stats, "mpi_threshold", E.get_size(), C.total_ncalls);

        matpoly sE(ab, m, b, E.get_size());
        matpoly spi;

        double const expect0 = bm.hints.tt_gather_per_unit * E.get_size();
        bm.stats.plan_smallstep("gather(L+R)", expect0);
        bm.stats.begin_smallstep("gather(L+R)");
        E.gather_mat(sE);
        E = bigmatpoly(model);
        bm.stats.end_smallstep();

        /* Only the master node does the local computation */
        if (!rank)
            spi = bw_lingen_single_nocp(bm, sE);

        double const expect1 = bm.hints.tt_scatter_per_unit * z;
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
    // fmt::print(stderr, "Leave {}\n", __func__);

    save_checkpoint_file(bm, LINGEN_CHECKPOINT_PI, pi, t0, t1);

    MPI_Barrier(bm.com[0]);

    return pi;
}/*}}}*/

/**/

/**********************************************************************/
/* vim:set sw=4 sta et: */
