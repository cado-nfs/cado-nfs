#include "cado.h"

#include <cstddef>      /* see https://gcc.gnu.org/gcc-4.9/porting_to.html */
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
#include <cassert>
#include <cfloat>
#include <stdexcept>
#include <algorithm>
#ifdef  HAVE_SIGHUP
#include <csignal>
#endif
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"
#include "macros.h"
#include "utils.h"
#ifndef SELECT_MPFQ_LAYER_u64k1
#include "lingen_qcode_prime.hpp"
// #include "lingen_polymat.hpp"
#else
#include "lingen_qcode_binary.hpp"
#endif
#include "lingen_matpoly_ft.hpp"
#include "lingen_mul_substeps.hpp"

#include "lingen.hpp"
#include "lingen_tuning.hpp"
#include "lingen_tuning_cache.hpp"
#include "lingen_platform.hpp"
#include "tree_stats.hpp"
#include "fmt/format.h"
#include "lingen_substep_schedule.hpp"
#include "lingen_substep_characteristics.hpp"

#include <vector>
#include <array>
#include <utility>
#include <map>
#include <string>
#include <ostream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

/* Given n>=1 return the list of all integers k (1<=k<=n) such
 * that it is possible to divide n into k blocks (not necessarily of
 * equal size) but such that all the splits that are obtained this way
 * have different maximal block sizes.
 *
 * (in the binary case we impose that block sizes are always multiples of
 * 64)
 *
 * The returned list is such that k->iceildiv(n, k) actually performs the
 * reversal of the list. And furthermore, for all k's such that k^2<=n,
 * the obtained splits are different. Hence it suffices to iterate over
 * all these k's
 */
std::vector<unsigned int> all_splits_of(unsigned int n)
{
#ifdef SELECT_MPFQ_LAYER_u64k1
    ASSERT_ALWAYS(n % 64 == 0);
    n /= 64;
#endif
    std::vector<unsigned int> res;
    for(unsigned int k = 1 ; k * k <= n ; k++) res.push_back(k);
    unsigned int j = res.size();
    if (res[j-1] * res[j-1] == n) j--;
    for( ; j-- ; ) res.push_back(iceildiv(n, res[j]));
    return res;
}

template<typename OP>
lingen_substep_schedule optimize(lingen_substep_characteristics<OP> const & U, lingen_platform const & P, lingen_tuning_cache & C, size_t reserved) { /* {{{ */
    unsigned int nr0 = U.mpi_split0(P).block_size_upper_bound();
    unsigned int nr1 = U.mpi_split1(P).block_size_upper_bound();
    unsigned int nr2 = U.mpi_split2(P).block_size_upper_bound();
    size_t min_my_ram = SIZE_MAX;
    lingen_substep_schedule S_lean;
    std::vector<lingen_substep_schedule> all_schedules;
    for(unsigned int shrink0 : all_splits_of(nr0)) {
        for(unsigned int shrink2 : all_splits_of(nr2)) {
            unsigned int nrs0 = U.shrink_split0(P, shrink0).block_size_upper_bound();
            unsigned int nrs2 = U.shrink_split2(P, shrink2).block_size_upper_bound();
#if 1
            /* first the splits with b0 == nrs0 */
            {
                unsigned int b0 = nrs0;
                for(unsigned int b1 : all_splits_of(nr1)) {
                    for(unsigned int b2 : all_splits_of(nrs2)) {
                        lingen_substep_schedule S;
                        S.shrink0 = shrink0;
                        S.shrink2 = shrink2;
                        S.batch = {{ b0, b1, b2 }};
                        size_t my_ram = U.get_peak_ram(P, S);
                        if (reserved + my_ram <= P.available_ram) {
                            all_schedules.push_back(S);
                        }
                        if (my_ram < min_my_ram) {
                            min_my_ram = my_ram;
                            S_lean = S;
                        }
                    }
                }
            }
            /* then the splits with b2 == nrs2 */
            {
                unsigned int b2 = nrs2;
                for(unsigned int b1 : all_splits_of(nr1)) {
                    for(unsigned int b0 : all_splits_of(nrs0)) {
                        lingen_substep_schedule S;
                        S.shrink0 = shrink0;
                        S.shrink2 = shrink2;
                        S.batch = {{ b0, b1, b2 }};
                        size_t my_ram = U.get_peak_ram(P, S);
                        if (reserved + my_ram <= P.available_ram) {
                            all_schedules.push_back(S);
                        }
                        if (my_ram < min_my_ram) {
                            min_my_ram = my_ram;
                            S_lean = S;
                        }
                    }
                }
            }
#else
            /* replicate the old choices. */
            for(unsigned int b1 : all_splits_of(nr1)) {
                lingen_substep_schedule S;
                S.shrink0 = shrink0;
                S.shrink2 = shrink2;
                if (nrs0 < nrs2) {
                    S.batch = { nrs0, b1, 1 };
                } else {
                    S.batch = { 1, b1, nrs2 };
                }
                size_t my_ram = U.get_peak_ram(P, S);
                if (reserved + my_ram <= P.available_ram) {
                    all_schedules.push_back(S);
                }
                if (my_ram < min_my_ram) {
                    min_my_ram = my_ram;
                    S_lean = S;
                }
            }
#endif
        }
    }
    std::sort(all_schedules.begin(), all_schedules.end());
    auto it = std::unique(all_schedules.begin(), all_schedules.end());
    all_schedules.erase(it, all_schedules.end());

    if (all_schedules.empty()) {
        char buf[20];
        std::ostringstream os;
        os << "Fatal error:"
            << " it is not possible to complete this calculation with only "
            << size_disp(P.available_ram, buf)
            << " of memory for intermediate transforms.\n";
        os << "Based on the cost for input length "
            << U.input_length
            << ", we need at the very least "
            << size_disp(U.get_peak_ram(P, S_lean), buf);
        os  << ", plus "
            << size_disp(reserved, buf)
            << " for reserved memory at upper levels";
        os << " [with shrink="
            << fmt::format("({},{})",
                    S_lean.shrink0, S_lean.shrink2)
            << ", batch="
            << fmt::format("({},{},{})",
                    S_lean.batch[0], S_lean.batch[1], S_lean.batch[2])
            << "]\n";
        fputs(os.str().c_str(), stderr);
        throw std::overflow_error(os.str());
    }

    for(auto & S : all_schedules) {
        /* This should ensure that all timings are obtained from cache */
        U.get_call_time(P, S, C);
    }

    U.sort_schedules(all_schedules, P, C);

    lingen_substep_schedule S = all_schedules.front();
    return S;
}
/* }}} */

template<typename fft_type>
struct lingen_tuner {
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;

    /* imported from the dims struct */
    abdst_field ab;
    cxx_mpz p;

    unsigned int m,n;

    size_t L;

    lingen_platform P;

    lingen_tuning_cache C;

    gmp_randstate_t rstate;

    const char * timing_cache_filename = NULL;

    /* stop measuring the time taken by the basecase when it is
     * more than this number times the time taken by the other
     * alternatives
     */
    double basecase_keep_until = 1.8;

    unsigned int forced_lingen_mpi_threshold = UINT_MAX;

    std::map<size_t, lingen_substep_schedule> schedules_mp, schedules_mul;

    static void declare_usage(cxx_param_list & pl) {/*{{{*/
        lingen_platform::declare_usage(pl);
        param_list_decl_usage(pl, "tuning_timing_cache_filename",
                "For --tune only: save (and re-load) timings for individual transforms in this file\n");
        param_list_decl_usage(pl, "basecase-keep-until",
                "For --tune only: stop measuring basecase timing when it exceeds the time of the recursive algorithm (counting its leaf calls) by this factor\n");
    }/*}}}*/

    static void lookup_parameters(cxx_param_list & pl) {/*{{{*/
        lingen_platform::lookup_parameters(pl);
        param_list_lookup_string(pl, "tuning_timing_cache_filename");
        param_list_lookup_string(pl, "basecase-keep-until");
    }/*}}}*/

    lingen_tuner(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl) :
        ab(d.ab), 
        m(d.m), n(d.n), L(L), P(comm, pl)
    {
#ifdef SELECT_MPFQ_LAYER_u64k1
        mpz_set_ui(p, 2);
#else
        mpz_set (p, abfield_characteristic_srcptr(ab));
#endif
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, 1);

        param_list_parse_double(pl, "basecase-keep-until", &basecase_keep_until);

        /* only the leader will do the tuning, so only the leader cares
         * about loading/saving it...
         */
        timing_cache_filename = param_list_lookup_string(pl, "tuning_timing_cache_filename");
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        if (rank == 0)
            C.load(timing_cache_filename);

        param_list_parse_uint(pl, "lingen_mpi_threshold", &forced_lingen_mpi_threshold);
    }

    ~lingen_tuner() {
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        if (rank == 0)
            C.save(timing_cache_filename);
        gmp_randclear(rstate);
    }

    std::tuple<size_t, double> mpi_threshold_comm_and_time() {/*{{{*/
        /* This is the time taken by gather() and scatter() right at the
         * threshold point. This total time is independent of the
         * threshold itself, in fact.
         */
        size_t data0 = m*(m+n)*(P.r*P.r-1);
        data0 = data0 * L * abvec_elt_stride(ab, 1);
        // We must **NOT** divide by r*r, because the problem is
        // precisely caused by the fact that gather() and scatter() all
        // imply one contention point which is the central node.
        // data0 = data0 / (r*r);
        std::tuple<size_t, double> vv { 2 * data0, 2 * data0 / P.mpi_xput};
        return vv;
    }/*}}}*/

    double compute_and_report_basecase(size_t length) { /*{{{*/
        double tt;

        lingen_tuning_cache::basecase_key K { mpz_sizeinbase(p, 2), m, n, length, P.openmp_threads };

        if (!C.has(K)) {
            tt = wct_seconds();
            test_basecase(ab, m, n, length, rstate);
            tt = wct_seconds() - tt;
            C[K] = { tt };
        }
        return C[K];
    }/*}}}*/

    typedef std::tuple<size_t, size_t, size_t, unsigned int> weighted_call_t;

    std::vector<weighted_call_t> calls_and_weights_at_depth(int i) {
        /*
         * Let L = (Q << (i+1)) + (u << i) + v, with u={0,1} and
         * 0<=v<(1<<i). We have L % (1 << i) = v. Note that u=(L>>i)%2.
         * Let U = u << i.
         *
         * input length at depth i, a.k.a. length(E), is either
         * L>>i=floor(L/2^i)=2Q+u or one above.  The number of calls is:
         *      
         * v times:
         *         length(E) = 2Q + u + 1
         *         length(E_left) = Q + 1
         *         length(E_right) = Q + u
         * otherwise:
         *         length(E) = 2Q + u
         *         length(E_left) = Q + u
         *         length(E_right) = Q
         * And in all cases:
         *         length(pi_left)  = ceiling(alpha * length(E_left))
         *         length(pi_right) = ceiling(alpha * length(E_right))
         *
         * The length of the chopped copy of E is:
         *         length(E_right) + length(pi_right) - 1
         *
         * The length of the product pi_left*pi_right is
         *      ceiling(alpha * length(E_left)) + ceiling(alpha * length(E_right)) - 1
         * This should also be less than or equal to
         *      ceiling(alpha * length(E))
         * So that in cases where the former bound is larger than the
         * latter, we're almost sure that the top coefficients multiply
         * to a zero matrix.
         *
         */                     

        size_t Q = L >> (i+1);
        size_t u = (L >> i) & 1;
        size_t v = L % (1 << i);

        if (v) {
            weighted_call_t w0 { 2*Q + u,     Q + u, Q,     (1 << i) - v };
            weighted_call_t w1 { 2*Q + u + 1, Q + 1, Q + u, v };
            std::vector<weighted_call_t> res {{ w0, w1 }};
            return res;
        } else {
            weighted_call_t w0 { 2*Q + u,     Q + u, Q,     (1 << i) - v };
            std::vector<weighted_call_t> res {{ w0 }};
            return res;
        }
    }

    lingen_substep_characteristics<op_mp<fft_type>> mp_substep(weighted_call_t const & cw) { /* {{{ */
        size_t length_E;
        size_t length_E_left;
        size_t length_E_right;
        unsigned int weight;

        std::tie(length_E, length_E_left, length_E_right, weight) = cw;

        ASSERT_ALWAYS(length_E >= 2);
        ASSERT_ALWAYS(weight);

        size_t csize = length_E_right;
        /* pi is the identity matrix for zero coefficients, but the
         * identity matrix already has length 1.
         */
        size_t bsize = 1 + iceildiv(m * length_E_left, m+n);
        size_t asize = csize + bsize - 1;

        ASSERT_ALWAYS(asize);
        ASSERT_ALWAYS(bsize);

        return lingen_substep_characteristics<op_mp<fft_type>>(
                ab, rstate, length_E,
                m, m+n, m+n,
                asize, bsize, csize);
    } /* }}} */
    void compute_schedules_for_mp(int i, bool print, size_t reserved=0) { /* {{{ */
        int printed_mem_once=0;
        for(auto const & cw : calls_and_weights_at_depth(i)) {
            size_t L = std::get<0>(cw);
            if (!recursion_makes_sense(L)) continue;
            auto step = mp_substep(cw);
            bool print_here = print && !printed_mem_once++;
            if (print_here)
                step.report_size_stats_human();

            lingen_substep_schedule S = optimize(step, P, C, reserved);

            if (print_here) {
                step.get_and_report_call_time(P, S, C);
            } else {
                step.get_call_time(P, S, C);
            }
            schedules_mp[L] = S;
        }
    } /* }}} */
    lingen_substep_characteristics<op_mul<fft_type>> mul_substep(weighted_call_t const & cw) { /* {{{ */
        size_t length_E;
        size_t length_E_left;
        size_t length_E_right;
        unsigned int weight;

        std::tie(length_E, length_E_left, length_E_right, weight) = cw;

        ASSERT_ALWAYS(length_E >= 2);
        ASSERT_ALWAYS(weight);

        size_t asize = 1 + iceildiv(m * length_E_left, m+n);
        size_t bsize = 1 + iceildiv(m * length_E_right, m+n);
        size_t csize = asize + bsize - 1;

        ASSERT_ALWAYS(asize);
        ASSERT_ALWAYS(bsize);

        return lingen_substep_characteristics<op_mul<fft_type>>(
                ab, rstate, length_E,
                m+n, m+n, m+n,
                asize, bsize, csize);

    } /* }}} */
    bool recursion_makes_sense(size_t L) const {
        return L >= 2;
    }

    void compute_schedules_for_mul(int i, bool print, size_t reserved) { /* {{{ */
        // we might want to consider the option of a single-node layer as
        // well. The effect of making this distinction between
        // lingen_threshold and lingen_mpi_threshold is not clear though,
        // since presently our formulas don't give rise to much
        // difference.
        // std::vector<lingen_platform> pps { P.single(), P };
        std::vector<lingen_platform> pps { P };
        for(auto const & pp : pps) {
            int printed_mem_once=0;
            for(auto const & cw : calls_and_weights_at_depth(i)) {
                size_t L = std::get<0>(cw);
                if (!recursion_makes_sense(L)) continue;
                auto step = mul_substep(cw);

                lingen_substep_schedule S = optimize(step, pp, C, reserved);
                if (print && !printed_mem_once++) {
                    step.report_size_stats_human();
                    step.get_and_report_call_time(pp, S, C);
                } else {
                    step.get_call_time(pp, S, C);
                }
                schedules_mul[L] = S;
            }
        }
    } /* }}} */

    lingen_hints tune_local() {
        size_t N = m*n*L/(m+n);
        char buf[20];
        printf("# Measuring lingen data for N ~ %zu m=%u n=%u for a %zu-bit prime p, using a %u*%u grid of %u-thread nodes [max target RAM = %s]\n",
                N, m, n, mpz_sizeinbase(p, 2),
                P.r, P.r, P.T,
                size_disp(P.available_ram, buf));
#ifdef HAVE_OPENMP
        printf("# Note: non-cached basecase measurements are done using openmp as it is configured for the running code, that is, with %d threads\n", P.openmp_threads);
#endif

        lingen_hints hints;

        int fl = log2(L) + 1;

        /* with basecase_keep_until == 0, then we never measure basecase */
        bool basecase_eliminated = basecase_keep_until == 0;
        std::map<size_t, std::tuple<bool, std::array<double, 3> >, lingen_tuning_cache::coarse_compare> best;
        size_t upper_threshold = SIZE_MAX;
        size_t peak = 0;
        int ipeak = -1;

        for(int i = fl ; i>=0 ; i--) {
            auto cws = calls_and_weights_at_depth(i);

            printf("####################### Measuring time at depth %d #######################\n", i);
            /* For input length L, the reserved
             * storage at depth i is
             *   RMP'(i)  = [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[2\alpha\ell_i]
             *   RMUL'(i) = [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[4\alpha\ell_i]
             * with the notations \alpha=m/(m+n), \ell_i=L/2^(i+1), and
             * [] denotes ceiling.
             * The details of the computation are in the comments in
             * lingen.cpp
             */
            size_t base_E  = iceildiv(m,P.r)*iceildiv(m+n,P.r)*mpz_size(p)*sizeof(mp_limb_t);
            size_t base_pi = iceildiv(m+n,P.r)*iceildiv(m+n,P.r)*mpz_size(p)*sizeof(mp_limb_t);
            size_t reserved_base = base_E * (L - (L >> i));
            size_t reserved_mp  = base_pi * iceildiv(m * iceildiv(L, 1<<i), m+n);
            size_t reserved_mul = base_pi * iceildiv(m * iceildiv(2*L, 1<<i), m+n);
            reserved_mp += reserved_base;
            reserved_mul += reserved_base;

            printf("# MP reserved storage = %s\n", size_disp(reserved_mp, buf));
            compute_schedules_for_mp(i, true, reserved_mp);

            printf("# MUL reserved storage = %s\n", size_disp(reserved_mul, buf));
            compute_schedules_for_mul(i, true, reserved_mul);

            double time_b = 0;
            double time_r = 0;
            double time_m = 0;
            double time_r_self = 0;
            double time_m_self = 0;
            size_t ram_mp = 0;
            size_t ram_mul = 0;

            bool basecase_was_eliminated = basecase_eliminated;

            ASSERT_ALWAYS(cws.size() <= 2);

            for(size_t idx = 0 ; idx < cws.size() ; idx++) {
                auto const & cw(cws[idx]);
                size_t L, Lleft, Lright;
                unsigned int weight;
                std::tie(L, Lleft, Lright, weight) = cw;

                /* the weight is the number of calls that must be made
                 * with this input length. The sum of weights at this
                 * depth must equal the total input length, whatever the
                 * level */
                ASSERT_ALWAYS(weight);

                if (!L) {
                    decltype(best)::mapped_type v { false, {{ 0, 0, 0 }}};
                    best[L] = v;
                    continue;
                }

                lingen_call_companion::key K { i, L };

                /* We **MUST** create hints[K], at this point */

                if (hints.find(K) == hints.end()) {
                    double ttb = DBL_MAX;
                    double ttr = DBL_MAX;
                    double ttrchildren = DBL_MAX;

                    lingen_call_companion U;
                    U.total_ncalls = 0;

                    bool forced = false;
                    bool rwin;

                    if (recursion_makes_sense(L) && forced_lingen_mpi_threshold != UINT_MAX) {
                        if (L >= forced_lingen_mpi_threshold) {
                            printf("# Forcing recursion at this level,"
                                    " since L=%zu>=lingen_mpi_threshold=%u\n",
                                    L, forced_lingen_mpi_threshold);
                            rwin = true;
                        } else {
                            printf("# Forcing basecase at this level,"
                                    " since L=%zu<lingen_mpi_threshold=%u\n",
                                    L, forced_lingen_mpi_threshold);
                            rwin = false;
                        }
                        forced = true;
                    }
                    if (!recursion_makes_sense(L) || (!(forced && rwin) && !basecase_eliminated))
                        ttb = compute_and_report_basecase(L);

                    if (recursion_makes_sense(L)) {
                        auto MP = mp_substep(cw);
                        U.mp = MP.get_companion(P, schedules_mp[L], C);
                        U.mp.reserved_ram = reserved_mp;

                        auto MUL = mul_substep(cw);
                        U.mul = MUL.get_companion(P, schedules_mul[L], C);
                        U.mul.reserved_ram = reserved_mp;

                        ttr = U.mp.tt.t + U.mul.tt.t;
                        ttrchildren = 0;
                        ttrchildren += std::get<1>(best[Lleft])[std::get<0>(best[Lleft])];
                        ttrchildren += std::get<1>(best[Lright])[std::get<0>(best[Lright])];

                        size_t m;
                        m = U.mp.ram + U.mp.reserved_ram;
                        if (m > ram_mp) ram_mp = m;
                        if (m > peak) { ipeak = i; peak = m; }

                        m = U.mul.ram + U.mul.reserved_ram;
                        if (m > ram_mul) ram_mul = m;
                        if (m > peak) { ipeak = i; peak = m; }
                    }

                    if (ttb >= basecase_keep_until * (ttr + ttrchildren))
                        basecase_eliminated = true;

                    if (!forced)
                        rwin = ttb >= std::min(1.0, basecase_keep_until) * (ttr + ttrchildren);

                    /* if basecase_keep_until < 1, then we probably want
                     * to prevent the basecase from being counted as
                     * winning at this point.
                     */
                    rwin = rwin || basecase_eliminated;
                    decltype(best)::mapped_type vv { rwin, {{ttb, ttr + ttrchildren, ttr}} };
                    best[L] = vv;

                    U.recurse = rwin;
                    /* See comment in compute_schedules_for_mul.
                     * Presently we don't identify cases where
                     * lingen_threshold makes sense at all */
                    U.go_mpi = rwin;
                    U.ttb = ttb;

                    ASSERT_ALWAYS(hints.find(K) == hints.end());
                    hints[K] = U;
                }
                ASSERT_ALWAYS(best.find(L) != best.end());
                hints[K].total_ncalls += weight;

                time_b += std::get<1>(best[L])[0] * weight;
                time_r += std::get<1>(best[L])[1] * weight;
                time_r_self += std::get<1>(best[L])[2] * weight;
                time_m += std::get<1>(best[L])[idx] * weight;
                time_m_self += std::get<1>(best[L])[2*idx] * weight;
            }

            size_t L0 = std::get<0>(cws.front());
            size_t L1 = std::get<0>(cws.back());
            /* calls_and_weights_at_depth must return a sorted list */
            ASSERT_ALWAYS(L0 <= L1);
            size_t L0r = lingen_round_operand_size(L0);
            size_t L1r = lingen_round_operand_size(L1);
            bool approx_same = L0r == L1r;
            bool rec0 = std::get<0>(best[L0]);
            bool rec1 = std::get<0>(best[L1]);

            std::ostringstream os;
            // os << "Depth " << i << ":";
            std::string oss = os.str();
            int pad = oss.size();
            const char * msg = oss.c_str();
            const char * msg2 = "";
            const char * strbest = " [BEST]";
            if (basecase_was_eliminated || !recursion_makes_sense(L1))
                strbest="";
            if (time_b < DBL_MAX) {
                const char * isbest = (!rec0 && !rec1) ? strbest : "";
                printf("#%*s basecase(threshold>%zu): %.2f [%.1fd]%s\n",
                        pad, msg, L1,
                        time_b, time_b / 86400, isbest);
                msg = msg2;
            }
            if (!approx_same && recursion_makes_sense(L1)) {
                const char * isbest = (rec1 && !rec0) ? strbest : "";
                std::ostringstream os2;
                os2 << " mixed(threshold=" << L1 << "): ";
                std::string ss2 = os2.str();
                printf("#%*s%s%.2f [%.1fd] (self: %.2f [%.1fd])%s\n",
                        pad, msg, ss2.c_str(),
                        time_m, time_m / 86400,
                        time_m_self, time_m_self / 86400, isbest);
                msg = msg2;
                int pad2 = pad + ss2.size();
                char buf[20];
                char buf2[20];
                if (ram_mp > ram_mul) {
                    printf("#%*s(memory(MP): %s, incl %s reserved)\n",
                            pad2, msg,
                            size_disp(ram_mp, buf),
                            size_disp(reserved_mp, buf2));
                } else {
                    printf("#%*s(memory(MUL): %s, incl %s reserved)\n",
                            pad2, msg,
                            size_disp(ram_mul, buf),
                            size_disp(reserved_mul, buf2));
                }

            }
            if (recursion_makes_sense(L0)) {
                const char * isbest = rec0 ? strbest : "";
                std::ostringstream os2;
                os2 << " recursive(threshold<=" << L0 << "): ";
                std::string ss2 = os2.str();
                printf("#%*s%s%.2f [%.1fd] (self: %.2f [%.1fd])%s\n", pad, msg, ss2.c_str(),
                        time_r, time_r / 86400, time_r_self, time_r_self / 86400, isbest);
                msg = msg2;
                int pad2 = pad + ss2.size();
                char buf[20];
                char buf2[20];
                if (ram_mp > ram_mul) {
                    printf("#%*s(memory(MP): %s, incl %s reserved)\n",
                            pad2, msg,
                            size_disp(ram_mp, buf),
                            size_disp(reserved_mp, buf2));
                } else {
                    printf("#%*s(memory(MUL): %s, incl %s reserved)\n",
                            pad2, msg,
                            size_disp(ram_mul, buf),
                            size_disp(reserved_mul, buf2));
                }
            }

            if (rec0) {
                // theshold is <= L0
                if (upper_threshold > L0) {
                    printf("# We expect lingen_mpi_threshold <= %zu\n", L0);
                    upper_threshold = L0;
                }
            } else if (rec1 && !rec0) {
                ASSERT_ALWAYS(cws.size() == 2);
                // threshold is =L1
                if (upper_threshold != L1) {
                    printf("# We expect lingen_mpi_threshold = %zu\n", L1);
                    upper_threshold = L1;
                }
            } else {
                // threshold is > L1
                if (upper_threshold <= L1) {
                    printf("# we expect lingen_mpi_threshold > %zu\n", L1);
                    upper_threshold = SIZE_MAX;
                }
            }
        }
        printf("################################# Total ##################################\n");
        if (forced_lingen_mpi_threshold != UINT_MAX) {
            printf("# Using explicit lingen_mpi_threshold=%zu (from command-line)\n", upper_threshold);
        } else {
            printf("# Automatically tuned lingen_mpi_threshold=%zu\n", upper_threshold);
        }
        size_t size_com0;
        double tt_com0;
        std::tie(size_com0, tt_com0) = mpi_threshold_comm_and_time();
        printf("# Communication time at lingen_mpi_threshold (%s): %.2f [%.1fd]\n", size_disp(size_com0, buf), tt_com0, tt_com0/86400);
        double time_best = std::get<1>(best[L])[std::get<0>(best[L])];
        time_best += tt_com0;
        printf("# Expected total time: %.2f [%.1fd], peak memory %s (at depth %d)\n", time_best, time_best / 86400, size_disp(peak, buf), ipeak);
        hints.ipeak=ipeak;
        hints.peak=peak;
        printf("(%u,%u,%u,%.1f,%1.f)\n",m,n,P.r,time_best,(double)peak/1024./1024./1024.);

        /* This one is strictly linear anyway */
        hints.tt_gather_per_unit = tt_com0 / 2 / L;
        hints.tt_scatter_per_unit = tt_com0 / 2 / L;

        return hints;
    }
    lingen_hints tune() {
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        lingen_hints hints;

        if (rank == 0)
            hints = tune_local();

        hints.share(0, P.comm);

        return hints;
    }
};

/* For the moment we're only hooking the fft_transform_info version */
lingen_hints lingen_tuning(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl)
{
#ifndef SELECT_MPFQ_LAYER_u64k1
    return lingen_tuner<fft_transform_info>(d, L, comm, pl).tune();
#else
    return lingen_tuner<gf2x_cantor_fft_info>(d, L, comm, pl).tune();
#endif
}

void lingen_tuning_decl_usage(cxx_param_list & pl)
{
#ifndef SELECT_MPFQ_LAYER_u64k1
    lingen_tuner<fft_transform_info>::declare_usage(pl);
#else
    lingen_tuner<gf2x_cantor_fft_info>::declare_usage(pl);
#endif
}

void lingen_tuning_lookup_parameters(cxx_param_list & pl)
{
#ifndef SELECT_MPFQ_LAYER_u64k1
    lingen_tuner<fft_transform_info>::lookup_parameters(pl);
#else
    lingen_tuner<gf2x_cantor_fft_info>::lookup_parameters(pl);
#endif
}

