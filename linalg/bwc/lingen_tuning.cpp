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
#include "lingen_qcode_select.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_mul_substeps.hpp"

#include "lingen.hpp"
#include "lingen_tuning.hpp"
#include "lingen_tuning_cache.hpp"
#include "lingen_platform.hpp"
#include "tree_stats.hpp"
#include "fmt/format.h"
#include "fmt/printf.h"
#include "lingen_substep_schedule.hpp"
#include "lingen_substep_characteristics.hpp"

#include <vector>
#include <array>
#include <utility>
#include <map>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

/* Given n>=1 return the list of all integers k (1<=k<=n) such
 * that it is possible to divide n into k blocks (not necessarily of
 * equal size) but such that all the splits that are obtained this way
 * have different maximal block sizes.
 *
 * (in the binary case we impose that block sizes are always multiples of
 * 8 -- this is done to limit the number of possible configurations)
 *
 * The returned list is such that k->iceildiv(n, k) actually performs the
 * reversal of the list. And furthermore, for all k's such that k^2<=n,
 * the obtained splits are different. Hence it suffices to iterate over
 * all these k's
 */
std::vector<unsigned int> all_splits_of(unsigned int n)
{
#ifdef SELECT_MPFQ_LAYER_u64k1
    n /= 8;
#endif
    std::vector<unsigned int> res;
    for(unsigned int k = 1 ; k * k <= n ; k++) res.push_back(k);
    unsigned int j = res.size();
    if (res[j-1] * res[j-1] == n) j--;
    for( ; j-- ; ) res.push_back(iceildiv(n, res[j]));
    return res;
}

std::vector<lingen_substep_schedule> 
optimize(std::ostream& os, lingen_substep_characteristics const & U, lingen_platform const & P, unsigned int mesh, lingen_tuning_cache & C, size_t reserved) { /* {{{ */
    unsigned int nr0 = U.mpi_split0(mesh).block_size_upper_bound();
    unsigned int nr1 = U.mpi_split1(mesh).block_size_upper_bound();
    unsigned int nr2 = U.mpi_split2(mesh).block_size_upper_bound();
    size_t min_my_ram = SIZE_MAX;
    lingen_substep_schedule S_lean;
    std::vector<lingen_substep_schedule> all_schedules;
    std::vector<lingen_substep_schedule::fft_type_t> allowed_ffts { lingen_substep_schedule::FFT_NONE };
#ifndef SELECT_MPFQ_LAYER_u64k1
    allowed_ffts.push_back(lingen_substep_schedule::FFT_FLINT);
#else
    allowed_ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
    allowed_ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
#endif

    for(lingen_substep_schedule::fft_type_t fft : allowed_ffts) {
        if (!U.fft_type_valid(fft)) continue;
        for(unsigned int shrink0 : all_splits_of(nr0)) {
            for(unsigned int shrink2 : all_splits_of(nr2)) {
                if (fft == lingen_substep_schedule::FFT_NONE)
                    if (shrink0 > 1 || shrink2 > 1) continue;
                unsigned int nrs0 = U.shrink_split0(mesh, shrink0).block_size_upper_bound();
                unsigned int nrs2 = U.shrink_split2(mesh, shrink2).block_size_upper_bound();
                /* first the splits with b0 == nrs0 */
                {
                    unsigned int b0 = nrs0;
                    for(unsigned int b1 : all_splits_of(nr1)) {
                        for(unsigned int b2 : all_splits_of(nrs2)) {
                            lingen_substep_schedule S;
                            S.fft_type = fft;
                            S.shrink0 = shrink0;
                            S.shrink2 = shrink2;
                            S.batch = {{ b0, b1, b2 }};
                            size_t my_ram = U.get_peak_ram(mesh, S);
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
                            S.fft_type = fft;
                            S.shrink0 = shrink0;
                            S.shrink2 = shrink2;
                            S.batch = {{ b0, b1, b2 }};
                            size_t my_ram = U.get_peak_ram(mesh, S);
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
            }
        }
    }

#ifdef SELECT_MPFQ_LAYER_u64k1
    os << fmt::sprintf("# %zu possible schedules to sort\n", all_schedules.size());
#endif

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
            << size_disp(U.get_peak_ram(mesh, S_lean), buf);
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
        /* This may print timing info to the output stream */
        U.get_call_time(os, P, mesh, S, C);
    }

    U.sort_schedules(os, all_schedules, P, mesh, C);

    std::map<std::string, lingen_substep_schedule> families;
    std::vector<lingen_substep_schedule> res;

    for(auto const & S : all_schedules) {
        std::string f = fmt::sprintf("%s%s;%s",
                mesh > 1 ? "MPI-" : "",
                op_mul_or_mp_base::op_name(U.op_type),
                S.fft_name());
        if (families.find(f) == families.end()) {
            families[f] = S;
            res.push_back(S);
        }
    }

    return res;
}
/* }}} */

struct lingen_tuner {
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;

    abdst_field ab; /* imported from the dims struct */
    cxx_mpz p;
    unsigned int m,n;
    size_t L;
    lingen_platform P;
    lingen_tuning_cache C;
    gmp_randstate_t rstate;
    const char * timing_cache_filename = NULL;
    const char * schedule_filename = NULL;
    std::ostream& os;

    /* stop measuring the time taken by the basecase when it is
     * more than this number times the time taken by the other
     * alternatives
     */
    double basecase_keep_until = 1.8;

    std::map<std::string, unsigned int> tuning_thresholds;
    /* length(E), length(E_left), length(E_right), number of occurrences
     */
    typedef std::tuple<size_t, size_t, size_t, unsigned int> weighted_call_t;
    std::vector<unsigned int> mesh_all;
    std::map<unsigned int, std::string> strat_name;


    struct output_info {/*{{{*/
        int quiet = 0;
        const char * tuning_log_filename = NULL;
        static void declare_usage(cxx_param_list & pl) {/*{{{*/
            param_list_decl_usage(pl, "tuning_log_filename",
                    "Output tuning log to this file\n");
            param_list_decl_usage(pl, "tuning_quiet",
                    "Silence tuning log\n");
        }/*}}}*/
        static void lookup_parameters(cxx_param_list & pl) {/*{{{*/
            lingen_platform::lookup_parameters(pl);
            param_list_lookup_string(pl, "tuning_quiet");
            param_list_lookup_string(pl, "tuning_log_filename");
        }/*}}}*/
        output_info(cxx_param_list & pl) {
            tuning_log_filename = param_list_lookup_string(pl, "tuning_log_filename");
            param_list_parse_int(pl, "tuning_quiet", &quiet);
        }
    };/*}}}*/
    static void declare_usage(cxx_param_list & pl) {/*{{{*/
        lingen_platform::declare_usage(pl);
        output_info::declare_usage(pl);
        param_list_decl_usage(pl, "tuning_schedule_filename",
                "Save (and re-load if it exists) tuning schedule from this file");
        param_list_decl_usage(pl, "tuning_timing_cache_filename",
                "Save (and re-load) timings for individual transforms in this file\n");
        param_list_decl_usage(pl, "basecase-keep-until",
                "When tuning, stop measuring basecase timing when it exceeds the time of the recursive algorithm (counting its leaf calls) by this factor\n");
        param_list_decl_usage(pl, "tuning_thresholds",
                "comma-separated list of threshols, given in the form <algorithm>:<threshold> value. Recognized values for <algorithm> are a subset of recursive,gfp_plain,flint,cantor,gf2x_plain. Thresholds are integers corresponding to the input size of E\n");
    }/*}}}*/
    static void lookup_parameters(cxx_param_list & pl) {/*{{{*/
        lingen_platform::lookup_parameters(pl);
        output_info::lookup_parameters(pl);
        param_list_lookup_string(pl, "tuning_schedule_filename");
        param_list_lookup_string(pl, "tuning_timing_cache_filename");
        param_list_lookup_string(pl, "basecase-keep-until");
        param_list_lookup_string(pl, "tuning_thresholds");
    }/*}}}*/
    lingen_tuner(std::ostream& os, bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl) :/*{{{*/
        ab(d.ab), 
        m(d.m), n(d.n), L(L), P(comm, pl), os(os)
    {
#ifdef SELECT_MPFQ_LAYER_u64k1
        mpz_set_ui(p, 2);
#else
        mpz_set (p, abfield_characteristic_srcptr(ab));
#endif
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, 1);

        param_list_parse_double(pl, "basecase-keep-until", &basecase_keep_until);

        schedule_filename = param_list_lookup_string(pl, "tuning_schedule_filename");
        /* only the leader will do the tuning, so only the leader cares
         * about loading/saving it...
         */
        timing_cache_filename = param_list_lookup_string(pl, "tuning_timing_cache_filename");
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        if (rank == 0)
            C.load(timing_cache_filename);

        const char * tmp = param_list_lookup_string(pl, "tuning_thresholds");
        if (tmp) {
            std::string tlist = tmp;
            for(size_t pos = 0 ; pos != string::npos ; ) {
                size_t next = tlist.find(',', pos);
                std::string tok;
                if (next == string::npos) {
                    tok = tlist.substr(pos);
                    pos = next;
                } else {
                    tok = tlist.substr(pos, next - pos);
                    pos = next + 1;
                }
                size_t colon = tok.find(':');
                if (colon == string::npos)
                    throw std::invalid_argument("tuning_thresholds is bad");
                std::string algorithm = tok.substr(0, colon);
                if (!(std::istringstream(tok.substr(colon + 1)) >> tuning_thresholds[algorithm])) {
                    throw std::invalid_argument("tuning_thresholds is bad");
                }
            }
        }
    }/*}}}*/
    ~lingen_tuner() {/*{{{*/
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        if (rank == 0)
            C.save(timing_cache_filename);
        gmp_randclear(rstate);
    }/*}}}*/
    std::tuple<size_t, double> mpi_threshold_comm_and_time() {/*{{{*/
        /* This is the time taken by gather() and scatter() right at the
         * threshold point. This total time is independent of the
         * threshold itself, in fact.
         */
        size_t data0 = m*(m+n)*(P.r*P.r-1);
        data0 = data0 * L;
#ifdef SELECT_MPFQ_LAYER_u64k1
        data0 = iceildiv(data0, ULONG_BITS) * sizeof(unsigned long);
#else
        data0 = abvec_elt_stride(ab, data0);
#endif

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
    std::vector<weighted_call_t> calls_and_weights_at_depth(int i) {/*{{{*/
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
    }/*}}}*/
    lingen_substep_characteristics substep(weighted_call_t const & cw, op_mul_or_mp_base::op_type_t op) { /* {{{ */
        size_t length_E;
        size_t length_E_left;
        size_t length_E_right;
        unsigned int weight;

        std::tie(length_E, length_E_left, length_E_right, weight) = cw;

        ASSERT_ALWAYS(length_E >= 2);
        ASSERT_ALWAYS(weight);

        size_t asize, bsize, csize;

        if (op == op_mul_or_mp_base::OP_MP) {
            csize = length_E_right;
            /* pi is the identity matrix for zero coefficients, but the
             * identity matrix already has length 1.
             */
            bsize = 1 + iceildiv(m * length_E_left, m+n);
            asize = csize + bsize - 1;
        } else {
            asize = 1 + iceildiv(m * length_E_left, m+n);
            bsize = 1 + iceildiv(m * length_E_right, m+n);
            csize = asize + bsize - 1;
        }

        ASSERT_ALWAYS(asize);
        ASSERT_ALWAYS(bsize);

        return lingen_substep_characteristics(
                ab, rstate, length_E,
                op,
                op == op_mul_or_mp_base::OP_MP ? m : (m+n),
                m+n, m+n,
                asize, bsize, csize);
    } /* }}} */
    bool recursion_makes_sense(size_t L) const {/*{{{*/
        return L >= 2;
    }/*}}}*/
    struct tuner_persistent_data {/*{{{*/
        typedef std::map<size_t, std::pair<unsigned int, double>, lingen_tuning_cache::coarse_compare> level_strategy_map;
        lingen_hints hints;
        lingen_hints const & stored_hints;
        level_strategy_map best;
        /* The "minimum mesh" field is only used by
         * tune_local_at_depth. This determines the "mesh" of calls
         * that will be tried. Options are:
         *
         * mesh=0 (if 0>=minimum_mesh) : basecase
         * mesh=1 (if 1>=minimum_mesh) : recursive, single-node
         * other                       : recursive, multi-node
         *
         * Currently, because scatter_mat and gather_mat are limited to
         * 1-n and n-1 conversions, we transition to single- to
         * multi-node all in one go.
         */
        int minimum_mesh = 0;
        double last_save = 0;
        size_t peak = 0;
        int ipeak = 0;
        size_t upper_threshold = 0;
        bool impose_hints;
        tuner_persistent_data(lingen_hints const & stored_hints) : stored_hints(stored_hints) {
            last_save = wct_seconds();
            impose_hints = !stored_hints.empty();
        }
    };/*}}}*/
    lingen_call_companion::mul_or_mp_times tune_local_at_depth_mp_or_mul(tuner_persistent_data & persist, weighted_call_t cw, int depth, unsigned int mesh, op_mul_or_mp_base::op_type_t op_type)/*{{{*/
    {
        lingen_call_companion::mul_or_mp_times U { op_type };
        lingen_hints const & stored_hints(persist.stored_hints);
        double & last_save(persist.last_save);
        size_t & peak(persist.peak);
        int & ipeak(persist.ipeak);

        /* For input length L, the reserved
         * storage at depth i is
         *   RMP'(i)  = [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[2\alpha\ell_i]
         *   RMUL'(i) = [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[4\alpha\ell_i]
         * with the notations \alpha=m/(m+n), \ell_i=L/2^(i+1), and
         * [] denotes ceiling.
         * The details of the computation are in the comments in
         * lingen.cpp
         */
        size_t base_E  = iceildiv(m,mesh)*iceildiv(m+n,mesh)*mpz_size(p)*sizeof(mp_limb_t);
        size_t base_pi = iceildiv(m+n,mesh)*iceildiv(m+n,mesh)*mpz_size(p)*sizeof(mp_limb_t);
        constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
        size_t reserved_base = base_E * iceildiv(L - (L >> depth), simd);
        size_t reserved_mp  = base_pi * iceildiv(iceildiv(m * iceildiv(L, 1<<depth), m+n), simd);
        size_t reserved_mul = base_pi * iceildiv(iceildiv(m * iceildiv(2*L, 1<<depth), m+n), simd);
        reserved_mp += reserved_base;
        reserved_mul += reserved_base;

        size_t reserved = op_type == op_mul_or_mp_base::OP_MP ? reserved_mp : reserved_mul;

        os << fmt::sprintf("# %s reserved storage = %s\n",
                op_mul_or_mp_base::op_name(op_type),
                size_disp(reserved));

        size_t L, Lleft, Lright;
        unsigned int weight;
        std::tie(L, Lleft, Lright, weight) = cw;
        lingen_call_companion::key K { depth, L };
        ASSERT_ALWAYS(weight);

        ASSERT_ALWAYS (recursion_makes_sense(L));
        auto step = substep(cw, op_type);
        bool print_here = true;
#if 0
        if (print_here)
            step.report_size_stats_human(os);
#endif

        std::vector<lingen_substep_schedule> SS;

        if (stored_hints.find(K) != stored_hints.end()) {
            lingen_substep_schedule S = stored_hints.at(K)[op_type].S;
            if (step.is_valid_schedule(P, mesh, S)) {
                os << "# Using imposed schedule " << S.serialize() << "\n";
                SS.push_back(S);
            } else {
                os << "# IGNORING invalid schedule " << S.serialize() << "\n";
            }
        }

        if (SS.empty()) {
            /* get the schedule by trying all possibilities */
            SS = optimize(os, step, P, mesh, C, reserved);
        }

        for(auto const & S : SS)
            step.get_and_report_call_time(os, P, mesh, S, C);

        lingen_substep_schedule S = SS.front();

        if (print_here)
            step.report_op_winner(os, mesh, S);

        if (wct_seconds() > last_save + 10) {
            int rank;
            MPI_Comm_rank(P.comm, &rank);
            if (rank == 0)
                C.save(timing_cache_filename);
            last_save = wct_seconds();
        }
        U = step.get_companion(os, P, mesh, S, C);
        U.reserved_ram = reserved;
        os << "#\n";

        size_t mm = U.ram_total();
        if (mm > peak) { ipeak = depth; peak = mm; }

        return U;
    }/*}}}*/
    void tune_local_at_depth(tuner_persistent_data & persist, int depth)/*{{{*/
    {
        tuner_persistent_data::level_strategy_map & best(persist.best);
        int & minimum_mesh(persist.minimum_mesh);
        lingen_hints & hints(persist.hints);
        lingen_hints const & stored_hints(persist.stored_hints);
        bool impose_hints(persist.impose_hints);

        auto cws = calls_and_weights_at_depth(depth);

        os << fmt::sprintf("####################### Measuring time at depth %d #######################\n", depth);

        ASSERT_ALWAYS(cws.size() <= 2);

        std::map<unsigned int, std::pair<unsigned int, double>> mesh_tt_weighted;

        lingen_call_companion U_typical;

        for(size_t idx = 0 ; idx < cws.size() ; idx++) {
            auto const & cw(cws[idx]);
            size_t L, Lleft, Lright;
            unsigned int weight;
            std::tie(L, Lleft, Lright, weight) = cw;
            if (!L) continue;
            double ratio = weight / (double) (1U << depth);
            os << fmt::sprintf("# input size %zu, %u times [%.1f%%]\n",
                    L, weight, 100*ratio);
        }
        os << "#\n";

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
                best[L] = { 0, 0 };
                mesh_tt_weighted[0].first += weight;
                continue;
            }

            lingen_call_companion::key K { depth, L };

            /* We **MUST** create hints[K], at this point. We will gain
             * _some_ insight from stored_hints[] and maybe from the
             * thresholds passed on the command line, but in any case we
             * will have to recompute the timings and the RAM usage */

            if (hints.find(K) == hints.end()) {
                unsigned int mesh = 0;
                bool forced = false;

                if (stored_hints.find(K) != stored_hints.end()) {
                    os << ("# Re-using stored schedule\n");
                    /* This is _not_ a complete set ! We still must
                     * compute the timings, the RAM, and so on. 
                     *
                     * And anyway, U.mp and U.mul are set by
                     * tune_local_at_depth_mp and
                     * tune_local_at_depth_mul, which read stored_hints
                     * again.
                     */
                    mesh = stored_hints.at(K).mesh;
                    forced = true;
                    if (mesh == UINT_MAX) mesh = P.r;
                    if (mesh > 1 && mesh != P.r) {
                        throw std::runtime_error(
                                fmt::sprintf(
                                    "stored schedule is invalid,"
                                    " we cannot (yet) run"
                                    " on a %d*%d grid"
                                    " a schedule meant"
                                    " for a %d*%d grid\n",
                                    P.r, P.r,
                                    mesh, mesh));
                    }

                    os << fmt::sprintf("# Forcing %s at this level\n",
                            strat_name[mesh]);
                } else if (recursion_makes_sense(L)) {
                    if (impose_hints) {
                        os << "# No stored schedule found in provided file,"
                            " computing new one\n";
                    }
                    std::string key_rec = "recursive";
                    std::string key_coll = "collective";
                    auto e = tuning_thresholds.end();
                    bool r = tuning_thresholds.find(key_rec)  != e;
                    bool c = tuning_thresholds.find(key_coll) != e;

                    size_t tr = r ? tuning_thresholds.at(key_rec)  : SIZE_MAX;
                    size_t tc = c ? tuning_thresholds.at(key_coll) : SIZE_MAX;
                    if (r && !c) tc = 0;
                    if (tc < tr) tc = tr;

                    if (r || c) {
                        forced = true;
                        mesh = L >= tc ? P.r : (L >= tr ? 1 : 0);

                        os << fmt::sprintf("# Forcing %s at this level,"
                                " since L=%zu"
                                ", tuning_threshold[%s]=%s"
                                ", tuning_threshold[%s]=%s"
                                "\n",
                                strat_name[mesh],
                                L,
                                key_rec,  r ? fmt::sprintf("%u", tr) : "undef",
                                key_coll, c ? fmt::sprintf("%u", tc) : "undef");
                    }
                }

                std::vector<unsigned int> mesh_try;
                std::map<unsigned int, lingen_call_companion> mesh_res;
                std::map<unsigned int, double> mesh_tt;
                std::map<unsigned int, double> mesh_tt_children;

                if (forced) {
                    mesh_try.push_back(mesh);
                } else {
                    if (minimum_mesh <= 0)
                        mesh_try.push_back(0);
                    if (recursion_makes_sense(L)) {
                        if (minimum_mesh <= 1)
                            mesh_try.push_back(1);
                        if (P.r > 1)
                            mesh_try.push_back(P.r);
                    }
                }
                for(auto mesh : mesh_try) {
                    lingen_call_companion U;
                    U.mesh = mesh;
                    if (mesh == 0) {
                        U.ttb = compute_and_report_basecase(L);
                        mesh_res[mesh] = U;
                        mesh_tt[mesh] = U.ttb;
                        mesh_tt_weighted[mesh].first += weight;
                        mesh_tt_weighted[mesh].second += mesh_tt[mesh] * weight;
                    } else {
                        U.mp = tune_local_at_depth_mp_or_mul(
                                persist, cw, depth, mesh,
                                op_mul_or_mp_base::OP_MP);
                        U.mul = tune_local_at_depth_mp_or_mul(
                                persist, cw, depth, mesh,
                                op_mul_or_mp_base::OP_MUL);
                        mesh_res[mesh] = U;
                        mesh_tt[mesh] = U.mp.tt.t + U.mul.tt.t;
                        mesh_tt_children[mesh] += best[Lleft].second;
                        mesh_tt_children[mesh] += best[Lright].second;
                        mesh_tt[mesh] += mesh_tt_children[mesh];
                        mesh_tt_weighted[mesh].first += weight;
                        mesh_tt_weighted[mesh].second += mesh_tt[mesh] * weight;
                    }
                }

                /* find the best mesh value */
                double ttbest = DBL_MAX;
                for(auto x : mesh_tt) {
                    if (x.second < ttbest) {
                        mesh = x.first;
                        ttbest = x.second;
                    }
                }

                lingen_call_companion U = mesh_res[mesh];
                U.total_ncalls = 0;
                U.complete = true;
                /* we no longer store ttb in the call companions which
                 * are intended for recursion */
                hints[K] = U;
                U_typical = U;

                /* discard mesh values that are less than the winner and
                 * appear to be slow enough that we don't think they'll
                 * ever catch up.
                 */
                for(auto x : mesh_tt) {
                    if (x.first == mesh || x.first > mesh)
                        continue;
                    if (x.second >= basecase_keep_until * mesh_tt[mesh]) {
                        os << fmt::sprintf("# Discarding %s from now on\n",
                                strat_name[x.first]);
                        minimum_mesh = x.first + 1;
                    }
                }
                /* At this point, we started using collective operations,
                 * which means that we don't want to go back. A priori.
                 * But it's quite difficult indeed because we might see
                 * some spurious results at small sizes. */
                // if (mesh > 1 && minimum_mesh <= 1) minimum_mesh = 2;

                best[L] = { mesh, mesh_tt[mesh] };
            }

            hints[K].total_ncalls += weight;
        }

        /* Now give a summary at this level */

        size_t L0 = std::get<0>(cws.front());
        size_t L1 = std::get<0>(cws.back());
        /* calls_and_weights_at_depth must return a sorted list */
        ASSERT_ALWAYS(L0 <= L1);
        unsigned int mesh0 = best[L0].first;
        unsigned int mesh1 = best[L1].first;

        for(auto mesh : mesh_all) {
            if (mesh_tt_weighted.find(mesh) != mesh_tt_weighted.end()) {
                double tt = mesh_tt_weighted[mesh].second;
                std::string rescaled;
                if (mesh_tt_weighted[mesh].first != (1U << depth)) {
                    double ratio = mesh_tt_weighted[mesh].first / (double) (1U << depth);
                    rescaled = fmt::sprintf("[rescaled from %.1f%%] ", 100*ratio);
                    tt /= ratio;
                }
                os << fmt::sprintf("# %s: %s%.2f [%.1fd]\n",
                        strat_name[mesh],
                        rescaled, tt, tt / 86400);
            }
        }

        double tt_total = best[L0].second * std::get<3>(cws.front())
                        + best[L1].second * std::get<3>(cws.back());

        if (mesh0 == mesh1) {
            os << fmt::sprintf("# BEST: %s: %.2f [%.1fd]\n",
                    strat_name[mesh0],
                    tt_total, tt_total / 86400);
        } else {
            os << fmt::sprintf("# BEST: mix of %s and %s: %.2f [%.1fd]\n",
                    strat_name[mesh0],
                    strat_name[mesh1],
                    tt_total, tt_total / 86400);
        }

        if (mesh0 || mesh1) {
            lingen_call_companion U = U_typical;
            if (U.mp.ram_total() > U.mul.ram_total()) {
                os << fmt::sprintf("#   (memory(MP): %s, incl %s reserved)\n",
                        size_disp(U.mp.ram_total()),
                        size_disp(U.mp.reserved_ram));
            } else {
                os << fmt::sprintf("#   (memory(MUL): %s, incl %s reserved)\n",
                        size_disp(U.mul.ram_total()),
                        size_disp(U.mul.reserved_ram));
            }

        }
    }/*}}}*/
    lingen_hints tune_local(lingen_hints & stored_hints) {/*{{{*/
        size_t N = m*n*L/(m+n);
        char buf[20];
        os << fmt::sprintf("# Measuring lingen data"
                " for N ~ %zu m=%u n=%u"
                " for a %zu-bit prime p,"
                " using a %u*%u grid of %u-thread nodes"
                " [max target RAM = %s]\n",
                N, m, n, mpz_sizeinbase(p, 2),
                P.r, P.r, P.T,
                size_disp(P.available_ram, buf));
#ifdef HAVE_OPENMP
        os << fmt::sprintf("# Note: non-cached basecase measurements"
                " are done using openmp as it is configured"
                " for the running code, that is, with %d threads\n",
                P.openmp_threads);
#endif
        
        mesh_all.push_back(0); strat_name[0] = "basecase";
        mesh_all.push_back(1); strat_name[1] = "recursive(single-node)";
        if (P.r > 1) {
            mesh_all.push_back(P.r);
            strat_name[P.r] = fmt::sprintf("recursive(%d*%d-nodes)", P.r, P.r);
        }

        int fl = log2(L) + 1;

        tuner_persistent_data persist(stored_hints);

        /* with basecase_keep_until == 0, then we never measure basecase */
        if (basecase_keep_until == 0)
            persist.minimum_mesh = 1;

        if (persist.impose_hints) {
            os << fmt::sprintf("# While we are doing timings here,"
                    " we'll take schedule decisions based on the hints"
                    " found in %s when they apply\n", schedule_filename);
        }

        for(int i = fl ; i>=0 ; i--)
            tune_local_at_depth(persist, i);

        tuner_persistent_data::level_strategy_map & best(persist.best);
        lingen_hints & hints(persist.hints);
        size_t peak(persist.peak);
        int ipeak(persist.ipeak);

        /*****************************************************************/
        /* make the mesh sizes monotonic. In truth, we can only do this
         * safely for mesh size 0 (basecase), since otherwise we would
         * ave an inconsistency in the sub-block sizes, which would make
         * the batch values invalid. */
        /* keys in the hint table are sorted as "top-level first" */
        std::map<unsigned int, lingen_call_companion::key> max_win_per_mesh;
        for(auto const & x : hints) {
            if (max_win_per_mesh.size() == mesh_all.size()) break;
            if (max_win_per_mesh.find(x.second.mesh) != max_win_per_mesh.end())
                max_win_per_mesh[x.second.mesh] = x.first;
        }
        for(auto & x : hints) {
            for(auto const & y : max_win_per_mesh) {
                if (y.first) continue;  // see above
                if (x.second.mesh > y.first && x.first < y.second) {
                    std::cout << fmt::format("## forcing %s at ({}) since it is known to win at ({})\n", strat_name[y.first],
                            y.first, y.second);
                    x.second.mesh = y.first;
                }
            }
        }

        os << ("################################# Total ##################################\n");
        if (!tuning_thresholds.empty()) {
            std::ostringstream ss;
            for(auto const & x : tuning_thresholds) {
                if (!ss.str().empty()) ss << ",";
                ss << x.first << ':' << x.second;
            }
            os << fmt::sprintf("# Using explicit tuning_thresholds=%s (from command-line)\n", ss.str());
        }
        size_t size_com0;
        double tt_com0;
        std::tie(size_com0, tt_com0) = mpi_threshold_comm_and_time();
        os << fmt::sprintf("# Communication time at lingen_mpi_threshold (%s): %.2f [%.1fd]\n", size_disp(size_com0, buf), tt_com0, tt_com0/86400);
        double time_best = best[L].second;
        time_best += tt_com0;
        os << fmt::sprintf("# Expected total time: %.2f [%.1fd], peak memory %s (at depth %d)\n", time_best, time_best / 86400, size_disp(peak, buf), ipeak);
        hints.ipeak=ipeak;
        hints.peak=peak;
        os << fmt::sprintf("(%u,%u,%u,%.1f,%1.f)\n",m,n,P.r,time_best,(double)peak/1024./1024./1024.);

        /* This one is strictly linear anyway */
        hints.tt_gather_per_unit = tt_com0 / 2 / L;
        hints.tt_scatter_per_unit = tt_com0 / 2 / L;

        return hints;
    }/*}}}*/
    lingen_hints tune() {/*{{{*/
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        lingen_hints hints;

        if (rank == 0) {
            lingen_hints stored_hints;
            if (schedule_filename) {
                std::ifstream is(schedule_filename);
                if (is && is >> stored_hints) {
                    /* This one _always_ goes to stdout */
                    std::cout << fmt::sprintf("# Read tuning schedule from %s\n", schedule_filename);
                } else {
                    std::cerr << fmt::sprintf("# Failed to read tuning schedule from %s\n", schedule_filename);
                }
            }
            hints = tune_local(stored_hints);
            if (schedule_filename && stored_hints.empty()) {
                std::ofstream os(schedule_filename);
                if (os && os << hints) {
                    /* This one _always_ goes to stdout */
                    std::cout << fmt::sprintf("# Written tuning schedule to %s\n", schedule_filename);
                } else {
                    std::cerr << fmt::sprintf("# Failed to write tuning schedule to %s\n", schedule_filename);
                }
            }
        }

        hints.share(0, P.comm);

        return hints;
    }/*}}}*/
};

lingen_hints lingen_tuning(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl)
{
    lingen_tuner::output_info O(pl);
    if (O.quiet) {
        std::ostringstream os;
        return lingen_tuner(os, d, L, comm, pl).tune();
    } else if (O.tuning_log_filename) {
        std::ofstream os(O.tuning_log_filename, std::ios_base::out);
        return lingen_tuner(os, d, L, comm, pl).tune();
    } else {
        return lingen_tuner(std::cout, d, L, comm, pl).tune();
    }
}

void lingen_tuning_decl_usage(cxx_param_list & pl)
{
    lingen_tuner::declare_usage(pl);
}

void lingen_tuning_lookup_parameters(cxx_param_list & pl)
{
    lingen_tuner::lookup_parameters(pl);
}

