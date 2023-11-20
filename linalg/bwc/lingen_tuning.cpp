#include "cado.h" // IWYU pragma: keep

#include <climits>                            // for UINT_MAX
#include <cstdint>                            // for SIZE_MAX
#include <cfloat>                              // for DBL_MAX
#include <cstdio>                              // for fputs, stderr
#include <cmath>                              // for log2

#include <algorithm>                           // for max, sort, find, unique
#include <array>                               // for array
#include <iostream>                            // for operator<<, basic_ostream
#include <fstream> // ofstream // IWYU pragma: keep
#include <map>                                 // for map, operator!=, map<>...
#include <memory>                              // for allocator_traits<>::va...
#include <stdexcept>                           // for invalid_argument, over...
#include <string>                              // for string, operator<<
#include <tuple>                               // for tie, get, tuple
#include <utility>                             // for pair, make_pair
#include <vector>                              // for vector

#include <gmp.h>                               // for mp_limb_t, mpz_size

#include "lingen_tuning.hpp"
#include "cxx_mpz.hpp"                         // for cxx_mpz
#include "fmt/core.h"                          // for check_format_string
#include "fmt/format.h"                        // for basic_buffer::append
#include "fmt/ostream.h"                       // for formatbuf<>::int_type
#include "arith-hard.hpp"                  // for abdst_field
#include "lingen_bw_dimensions.hpp"            // for bw_dimensions
#include "lingen_call_companion.hpp"           // for lingen_call_companion
#include "lingen_matpoly_select.hpp"           // for matpoly, matpoly::over...
#include "lingen_mul_substeps_base.hpp"        // for op_mul_or_mp_base, op_...
#include "lingen_platform.hpp"                 // for lingen_platform
#include "lingen_qcode_select.hpp"             // for test_basecase
#include "lingen_substep_characteristics.hpp"  // for lingen_substep_charact...
#include "lingen_substep_schedule.hpp"         // for lingen_substep_schedule
#include "lingen_tuning_cache.hpp"             // for lingen_tuning_cache::c...
#include "macros.h"                            // for iceildiv, ASSERT_ALWAYS
#include "misc.h"                              // for size_disp
#include "params.h"                            // for cxx_param_list, param_...
#include "subdivision.hpp"                     // for subdivision
#include "timing.h"                            // for wct_seconds, weighted_...

/* {{{ all_splits_of
 * Given n>=1 return the list of all integers k (1<=k<=n) such
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
#ifdef LINGEN_BINARY
    n /= 8;
#endif
    std::vector<unsigned int> res;
    for(unsigned int k = 1 ; k * k <= n ; k++) res.push_back(k);
    unsigned int j = res.size();
    if (res[j-1] * res[j-1] == n) j--;
    for( ; j-- ; ) res.push_back(iceildiv(n, res[j]));
    return res;
}
/* }}} */

std::vector<lingen_substep_schedule> optimize(
        std::ostream& os,
        lingen_substep_characteristics const & U,
        lingen_platform const & P,
        unsigned int mesh,
        std::vector<lingen_substep_schedule::fft_type_t> const & allowed_ffts,
        bool do_timings,
        lingen_tuning_cache & C,
        size_t reserved)
{ /* {{{ */
    unsigned int nr0 = U.mpi_split0(mesh).block_size_upper_bound();
    unsigned int nr1 = U.mpi_split1(mesh).block_size_upper_bound();
    unsigned int nr2 = U.mpi_split2(mesh).block_size_upper_bound();
    size_t min_my_ram = SIZE_MAX;
    lingen_substep_schedule S_lean;
    S_lean.fft_type = lingen_substep_schedule::FFT_NONE; /* placate compiler */
    std::vector<lingen_substep_schedule> all_schedules;

    unsigned int nvalid = 0;
    for(lingen_substep_schedule::fft_type_t fft : allowed_ffts) {
        if (!U.fft_type_valid(fft)) continue;
        nvalid++;
    }
    if (!nvalid) {
        throw std::invalid_argument("FFT choice restricted to zero valid FFTs, tuning_thresholds is probably wrong");
    }

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
                            if (reserved + my_ram <= P.available_ram || P.available_ram == 0) {
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
                            if (reserved + my_ram <= P.available_ram || P.available_ram == 0) {
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

    std::sort(all_schedules.begin(), all_schedules.end());

    all_schedules.erase(
            std::unique(all_schedules.begin(), all_schedules.end()),
            all_schedules.end());

    if (all_schedules.empty()) {
        char buf[20];
        std::ostringstream os;
        os << "Fatal error:"
            << " it is not possible to complete this calculation with only "
            << (P.available_ram ? size_disp(P.available_ram, buf) : "+infinity (?)")
            << " of memory for intermediate transforms.\n";
        os << "Based on the cost for input length "
            << U.input_length
            << ", we need at the very least "
            << size_disp(U.get_peak_ram(mesh, S_lean), buf);
        os  << ", plus "
            << size_disp(reserved, buf)
            << " for reserved memory at upper levels";
        os << " [with shrink="
            << fmt::format(FMT_STRING("({},{})"),
                    S_lean.shrink0, S_lean.shrink2)
            << ", batch="
            << fmt::format(FMT_STRING("({},{},{})"),
                    S_lean.batch[0], S_lean.batch[1], S_lean.batch[2])
            << "]\n";
        fputs(os.str().c_str(), stderr);
        throw std::overflow_error(os.str());
    }

    U.sort_schedules(os, all_schedules, P, mesh, C, do_timings);

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

    matpoly::arith_hard * ab; /* imported from the dims struct */
    cxx_mpz p;
    unsigned int m,n;
    size_t L;
    lingen_platform P;
    lingen_tuning_cache C;
    gmp_randstate_t rstate;
    const char * timing_cache_filename = NULL;
    const char * schedule_filename = NULL;

    /* stop measuring the time taken by the basecase when it is
     * more than this number times the time taken by the other
     * alternatives
     */
    double basecase_keep_until = 1.8;

    struct tuning_thresholds_t : public std::map<std::string, unsigned int> {/*{{{*/
        typedef std::map<std::string, unsigned int> super;
        static constexpr const char * recursive = "recursive";
        static constexpr const char * collective = "collective";
        static constexpr const char * ternary = "ternary";
        static constexpr const char * cantor = "cantor";
        static constexpr const char * flint = "flint";
        static constexpr const char * notiming = "notiming";
        static const std::vector<const char *> thresholds_verbs;
        static const std::vector<std::pair<lingen_substep_schedule::fft_type_t, const char *>> code_to_key;

        bool has(std::string const & key) const {
            return find(key) != end();
        }
        private:
        unsigned int& getref(std::string const & key) {
            return ((super&)(*this))[key];
        }
        public:
        unsigned int operator[](std::string const & key) const {
            return has(key) ? at(key) : UINT_MAX;
        }
        tuning_thresholds_t(cxx_param_list & pl, std::ostream& os, lingen_platform const & P) {/*{{{*/
            const char * tmp = param_list_lookup_string(pl, "tuning_thresholds");
            if (!tmp) return;
            std::string tlist = tmp;
            for(size_t pos = 0 ; pos != std::string::npos ; ) {
                size_t next = tlist.find(',', pos);
                std::string tok;
                if (next == std::string::npos) {
                    tok = tlist.substr(pos);
                    pos = next;
                } else {
                    tok = tlist.substr(pos, next - pos);
                    pos = next + 1;
                }
                auto error = [&tok](std::string const& reason) {
                    std::string base = fmt::format(
                            "tuning_thresholds is bad:"
                            " pair \"{}\" ", tok);
                    throw std::invalid_argument(base + reason);
                };

                size_t colon = tok.find(':');
                if (colon == std::string::npos)
                    error("has no colon");

                std::string algorithm = tok.substr(0, colon);
                if (std::find(thresholds_verbs.begin(), thresholds_verbs.end(), algorithm) == thresholds_verbs.end()) {
                    std::ostringstream os;
                    for(auto const & x : thresholds_verbs)
                        os << " " << x;
                    error(fmt::format(
                                "uses unrecognized key \"{}\""
                                " (recognized keys:{})",
                                algorithm, os.str()));
                }

                unsigned int & dst(getref(algorithm));
                if (!(std::istringstream(tok.substr(colon + 1)) >> dst))
                    error("has no understandable integer threshold");
            }
            if (has(collective) && P.r == 1) {
                const char * what = "interpreted as \"recursive\"";
                if (!has(recursive)) {
                    getref(recursive) = getref(collective);
                } else if (getref(collective) < getref(recursive)) {
                    getref(recursive) = getref(collective);
                } else {
                    what = "ignored";
                }
                erase(find(collective));
                os << "# Note: the tuning threshold \"collective\" is "
                    << what << " here, since we have a non-MPI run\n";
            }
        }/*}}}*/
    };/*}}}*/

    tuning_thresholds_t tuning_thresholds;

    /* length(E), length(E_left), length(E_right), number of occurrences */
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
        ab(&d.ab), 
        m(d.m), n(d.n), L(L), P(comm, pl),
        tuning_thresholds(pl, os, P)
    {
        mpz_set (p, ab->characteristic());
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
    }/*}}}*/
    // coverity[exn_spec_violation]
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
#ifdef LINGEN_BINARY
        data0 = iceildiv(data0, ULONG_BITS) * sizeof(unsigned long);
#else
        data0 = ab->vec_elt_stride(data0);
#endif

        // We must **NOT** divide by r*r, because the problem is
        // precisely caused by the fact that gather() and scatter() all
        // imply one contention point which is the central node.
        // data0 = data0 / (r*r);
        std::tuple<size_t, double> vv { 2 * data0, 2 * data0 / P.mpi_xput};
        return vv;
    }/*}}}*/
    double compute_and_report_basecase(std::ostream& os, size_t length, bool do_timings) { /*{{{*/
        double tt;

        lingen_tuning_cache::basecase_key K { mpz_sizeinbase(p, 2), m, n, length, P.openmp_threads };

        os << fmt::format(FMT_STRING("# basecase (@{}): "), length);

        if (!do_timings) {
            C[K] = { 0 };
            os << "not timed\n";
            return C[K];
        }

        if (!C.has(K)) {
            os << std::flush;
            tt = wct_seconds();
            test_basecase(ab, m, n, length, rstate);
            tt = wct_seconds() - tt;
            C[K] = { tt };
            os << fmt::sprintf("%.2f\n", tt);
        } else {
            os << fmt::sprintf("%.2f [from cache]\n", C[K]);
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

        /* at i == floor level, we have (1<<(i-1)) < L <= (1<<i), Q = 0, and:
         *
         * if v, then u=0 and v=L. we have v calls with size 1 and (1<<i)-v
         * with size 0 on the first read, but the calls with size 0 don't
         * really exist, so that in fact we're seeing leftovers of calls
         * with size 1 at the level above. So we actually have _at depth
         * i proper_ 2*v-(1<<i) calls with size 1, and that's all.
         *
         * if v==0, then u=1. we have all calls with size 1, code below
         * is fine
         */

        if (v && !((L-1)>>i)) {
            ASSERT_ALWAYS(Q == 0 && u == 0);
            weighted_call_t w1 { 2*Q + u + 1, Q + 1, Q + u, 2*v - (1 << i) };
            std::vector<weighted_call_t> res {{ w1 }};
            return res;
        } else if (v) {
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
        unsigned int minimum_mesh = 0;
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
    lingen_call_companion::mul_or_mp_times tune_local_at_depth_mp_or_mul(
            std::ostream& os,
            tuner_persistent_data & persist,
            weighted_call_t cw,
            int depth,
            unsigned int mesh,
            std::vector<lingen_substep_schedule::fft_type_t> const & allowed_ffts,
            bool do_timings,
            op_mul_or_mp_base::op_type_t op_type)/*{{{*/
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
         *
         * Note that here, we're using P.r and not the mesh parameter,
         * because what matters is what we're reserving at the _upper_
         * levels.  And most importantly, chances are that the largest
         * share of the reserved memory comes from the topmost levels
         * anyway, so we surmise that P.r is appropriate there.
         */
        size_t base_E  = iceildiv(m,P.r)*iceildiv(m+n,P.r)*mpz_size(p)*sizeof(mp_limb_t);
        size_t base_pi = iceildiv(m+n,P.r)*iceildiv(m+n,P.r)*mpz_size(p)*sizeof(mp_limb_t);
        constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
        size_t reserved_base = base_E * iceildiv(L - (L >> depth), simd);
        size_t reserved_mp  = base_pi * iceildiv(iceildiv(m * iceildiv(L, 1<<depth), m+n), simd);
        size_t reserved_mul = base_pi * iceildiv(iceildiv(m * iceildiv(2*L, 1<<depth), m+n), simd);
        reserved_mp += reserved_base;
        reserved_mul += reserved_base;

        size_t reserved = op_type == op_mul_or_mp_base::OP_MP ? reserved_mp : reserved_mul;

        os << fmt::format(FMT_STRING("# {} reserved storage = {}\n"),
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
            /* get the schedule by trying all possibilities for the
             * shrink and batch options. */
            SS = optimize(os, step, P, mesh,
                    allowed_ffts, do_timings, C, reserved);
        }

        for(auto const & S : SS)
            step.get_and_report_call_time(os, P, mesh, S, C, do_timings);

        lingen_substep_schedule S = SS.front();

        if (print_here)
            step.report_op_winner(os, mesh, S);

        if (do_timings) {
            if (wct_seconds() > last_save + 10) {
                int rank;
                MPI_Comm_rank(P.comm, &rank);
                if (rank == 0)
                    C.save(timing_cache_filename);
                last_save = wct_seconds();
            }
        }

        U = step.get_companion(os, P, mesh, S, C, reserved, do_timings);
        os << "#\n";

        size_t mm = U.ram_total();
        if (mm > peak) { ipeak = depth; peak = mm; }

        return U;
    }/*}}}*/
    struct configurations_to_test {/*{{{*/
        /* This is the list of dimensions or the MPI mesh that
         * we're going to try. For the moment, it's limited to
         * either {1*1, r*r}, or only one of them. In principle,
         * the framework could be extended to support further
         * mesh sizes, but that requires support in scatter_mat
         * and gather_mat to dispatch from m nodes to n nodes.
         */
        std::vector<unsigned int> mesh_sizes;
        std::vector<lingen_substep_schedule::fft_type_t> ffts;
        /* The various choices for the "shrink" and "batch" options are
         * not listed here, but are rather computed on the fly in the
         * optimize() function. This seems more convenient there, as it's
         * easier to do pruning based on the mesh sizes and such.
         */
        bool do_timings = true;
        configurations_to_test() {
            ffts = { lingen_substep_schedule::FFT_NONE };
#ifndef LINGEN_BINARY
            ffts.push_back(lingen_substep_schedule::FFT_FLINT);
#else
            ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
            ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
#endif
        }
        bool from_stored_hints(std::ostream& os, lingen_tuner & tuner, tuner_persistent_data & persist, lingen_call_companion::key const & K) {/*{{{*/
            lingen_hints const & stored_hints(persist.stored_hints);
            lingen_platform const & P(tuner.P);

            if (stored_hints.find(K) == stored_hints.end())
                return false;
            os << ("# Re-using stored schedule\n");
            /* This is _not_ a complete set ! We still must compute the
             * timings, the RAM, and so on.
             *
             * This is done via tune_local_at_depth_mp_or_mul, which
             * reads stored_hints a second time. In particular, it
             * fetches the info of the previously selected fft types for
             * both operations, and eventually sets U.mp and U.mul. This
             * implies that we do not need to restrict the list of
             * allowed ffts here.
             */
            unsigned int mesh = stored_hints.at(K).mesh;

            if (mesh == UINT_MAX) mesh = P.r;
            if (mesh > 1 && mesh != P.r) {
                throw std::runtime_error(
                        fmt::sprintf(
                            "stored schedule is invalid,"
                            " we cannot (yet) run on a %d*%d grid"
                            " a schedule meant for a %d*%d grid\n",
                            P.r, P.r,
                            mesh, mesh));
            }

            mesh_sizes = { mesh };
            os << fmt::format(FMT_STRING("# Forcing {} at depth {} L={}\n"),
                    tuner.strat_name[mesh], K.depth, K.L);

            return true;
        }/* }}} */
        std::string explain(tuning_thresholds_t const & T, std::string const & k) {
            if (T.has(k)) {
                return fmt::format(FMT_STRING(" tuning_threshold[{}]={}"), k, T[k]);
            } else {
                return fmt::format(FMT_STRING(" tuning_threshold[{}]=undef"), k);
            }
        }
        typedef tuning_thresholds_t T_t;
        bool from_thresholds_mesh_level(std::ostream& os, lingen_tuner & tuner, tuner_persistent_data & persist, lingen_call_companion::key const & K) {/*{{{*/
            T_t const & T(tuner.tuning_thresholds);
            lingen_platform const & P(tuner.P);
            unsigned int L = K.L;
            /* first decision: do we force recursion. There are various
             * reasons to do so, since we have various thresholds that
             * implicitly enable recursion */
            std::vector<std::pair<unsigned int, const char *>> all_rec;
            for(auto k : T_t::thresholds_verbs) {
                if (k == T_t::notiming) continue;
                if (T.has(k))
                    all_rec.push_back(std::make_pair(T[k], k));
            }
            /* If no threshold talks about being recursive, then it's
             * left to our decision. We keep all possible mesh sizes a
             * priori. However, we must pay attention to the fact that we
             * may know already that some mesh sizes should be discarded
             * (this can only happen if at an earlier level, the current mesh
             * size X was beaten by a larger one Y, which means that the
             * minimum mesh size was at most X for the previous size, and
             * is not at most Y, not more).
             */
            if (all_rec.empty()) {
                mesh_sizes.clear();
                if (persist.minimum_mesh <= 0) mesh_sizes.push_back(0);
                if (persist.minimum_mesh <= 1) mesh_sizes.push_back(1);
                if (P.r != 1 && persist.minimum_mesh <= P.r)
                    mesh_sizes.push_back(P.r);
                return true;
            }
            std::sort(all_rec.begin(), all_rec.end());

            mesh_sizes.clear();
            unsigned int next_mesh = 1;
            std::string explainer;
            for(auto const & tk : all_rec) {
                if (L < tk.first)
                    break;
                if (tk.second == T_t::collective)
                    next_mesh = P.r;

                mesh_sizes = { next_mesh };
                explainer = tk.second;
            }
            std::ostringstream explanation;
            bool done = false;
            if (mesh_sizes.empty()) {
                /* Given that we have identified at least one threshold
                 * that says "go recursive" and that we decided that we
                 * _don't_ go recursive, the conclusion is here that we
                 * want to do basecase only.
                 */
                mesh_sizes = { 0 };
                explanation << explain(T, all_rec.front().second);
                done = true;
            } else {
                explanation << explain(T, explainer);
                constexpr const char * ck = T_t::collective;
                if (!T.has(ck) && P.r > 1) {
                    mesh_sizes.push_back(P.r);
                    explanation << explain(T, ck);
                }
            }

            os << "# Testing only";
            for(auto mesh : mesh_sizes)
                os << " " << tuner.strat_name[mesh];
            os << fmt::format(FMT_STRING(" at depth {} L={} since"), K.depth, K.L)
                << explanation.str()
                << "\n";
            return done;
        }/*}}}*/
        bool from_thresholds_fft_level(std::ostream& os, lingen_tuner & tuner, lingen_call_companion::key const & K) {/*{{{*/
            T_t const & T(tuner.tuning_thresholds);
            /* second decision: should we restrict the set of ffts ? */
            /* In the following, we rely on the fact that the asymptotic
             * ordering is known in advance. If we add more cases, this
             * assumption will likely not hold anymore
             *
             * Note that by design, if we arrive here, we completed the
             * function find_configurations_mesh_level and got a true value,
             * meaning that L is larger than at least one of the recursive
             * thresholds. Therefore, checking against "recursive" or
             * "collective" doesn't make sense.
             */
            unsigned int L = K.L;
            ffts.clear();
            /* Okay, it's a bit of spaghetti, but the small and concise
             * loop would not be easier to understand.
             */
            std::ostringstream explanation;
#ifdef LINGEN_BINARY
            bool hc = T.has(T_t::cantor);
            bool ht = T.has(T_t::ternary);
            unsigned int tc = T[T_t::cantor];
            unsigned int tt = T[T_t::ternary];
            if (hc && ht) {
                if (tc >= tt) {
                    if (L >= tc) {
                        ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
                        explanation << explain(T, T_t::cantor);
                    } else if (L >= tt) {
                        ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
                        explanation << explain(T, T_t::ternary);
                        explanation << explain(T, T_t::cantor);
                    } else {
                        ffts.push_back(lingen_substep_schedule::FFT_NONE);
                        explanation << explain(T, T_t::recursive);
                        explanation << explain(T, T_t::ternary);
                    }
                } else if (tt >= tc) {
                    /* quite unlikely except for ridiculously small
                     * matrices perhaps */
                    if (L >= tt) {
                        ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
                        explanation << explain(T, T_t::ternary);
                    } else if (L >= tc) {
                        ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
                        explanation << explain(T, T_t::cantor);
                        explanation << explain(T, T_t::ternary);
                    } else {
                        ffts.push_back(lingen_substep_schedule::FFT_NONE);
                        explanation << explain(T, T_t::recursive);
                        explanation << explain(T, T_t::cantor);
                    }
                }
            } else if (hc) {
                if (L >= tc) {
                    ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
                    explanation << explain(T, T_t::cantor);
                } else {
                    ffts.push_back(lingen_substep_schedule::FFT_NONE);
                    ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
                    explanation << explain(T, T_t::recursive);
                    explanation << explain(T, T_t::cantor);
                }
            } else if (ht) {
                if (L >= tt) {
                    ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
                    ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
                    explanation << explain(T, T_t::ternary);
                } else {
                    ffts.push_back(lingen_substep_schedule::FFT_NONE);
                    explanation << explain(T, T_t::recursive);
                    explanation << explain(T, T_t::ternary);
                }
            } else {
                ffts.push_back(lingen_substep_schedule::FFT_NONE);
                ffts.push_back(lingen_substep_schedule::FFT_TERNARY);
                ffts.push_back(lingen_substep_schedule::FFT_CANTOR);
                explanation << explain(T, T_t::recursive);
            }
#else
            bool hf = T.has(T_t::flint);
            bool tf = T[T_t::flint];
            if (hf) {
                if (L >= tf) {
                    ffts.push_back(lingen_substep_schedule::FFT_FLINT);
                    explanation << explain(T, T_t::flint);
                } else {
                    ffts.push_back(lingen_substep_schedule::FFT_NONE);
                    explanation << explain(T, T_t::recursive);
                    explanation << explain(T, T_t::flint);
                }
            } else {
                ffts.push_back(lingen_substep_schedule::FFT_NONE);
                ffts.push_back(lingen_substep_schedule::FFT_FLINT);
                explanation << explain(T, T_t::recursive);
            }
#endif
            os << "# Testing only";
            for(auto fft : ffts)
                os << " " << lingen_substep_schedule::fft_name(fft);
            os << fmt::format(FMT_STRING(" at depth {} L={} since"), K.depth, K.L)
                << explanation.str()
                << "\n";
            return true;
        }/*}}}*/
    };/*}}}*/

    /* Determine which configurations will undergo testing. This covers
     * the single/collective, basecase/recursive choices, as well as the
     * fft choice. This does NOT cover the shrink and batch setting,
     * which are covered in tune_local_at_depth_mp_or_mul() and optimize()
     */
    configurations_to_test find_configurations_to_test(std::ostream & os,  tuner_persistent_data & persist, lingen_call_companion::key const & K) {/*{{{*/
        unsigned int L = K.L;
        configurations_to_test C;
        bool impose_hints(persist.impose_hints);

        if (C.from_stored_hints(os, *this, persist, K))
            return C;

        if (L >= tuning_thresholds[tuning_thresholds_t::notiming])
            C.do_timings = false;

        if (!recursion_makes_sense(L)) {
            C.mesh_sizes = { 0 };
            C.ffts.clear();
            return C;
        }

        if (impose_hints) {
            os << "# No stored schedule found in provided file,"
                " computing new one\n";
        }
        /* Try with the tuning_thresholds */

        if (C.from_thresholds_mesh_level(os, *this, persist, K)) {
            /* then the decision was trivial because lingen_thresholds
             * told us nothing. No need to bother with fft types either.
             */
            return C;
        }

        C.from_thresholds_fft_level(os, *this, K);
        return C;
    }/*}}}*/

    bool tune_local_at_depth(std::ostream& os, tuner_persistent_data & persist, int depth)/*{{{*/
    {
        tuner_persistent_data::level_strategy_map & best(persist.best);
        unsigned int & minimum_mesh(persist.minimum_mesh);
        lingen_hints & hints(persist.hints);

        auto cws = calls_and_weights_at_depth(depth);

        // output first goes to a temp string until we're sure that we're
        // doing timings.
        
        std::ostringstream os_pre;
        bool timed_something = false;

        os_pre << fmt::format(FMT_STRING("####################### Measuring time at depth {} #######################\n"), depth);

        ASSERT_ALWAYS(cws.size() <= 2);

        /* Given a mesh size (or 0 to mean "basecase"), this returns a
         * pair (number of calls, sum of the time over all calls)
         */
        std::map<unsigned int, std::pair<unsigned int, double>> mesh_tt_weighted;

        lingen_call_companion U_typical;

        for(size_t idx = 0 ; idx < cws.size() ; idx++) {
            auto const & cw(cws[idx]);
            size_t L, Lleft, Lright;
            unsigned int weight;
            std::tie(L, Lleft, Lright, weight) = cw;
            if (!L) continue;
            double ratio = weight / (double) (1U << depth);
            os_pre << fmt::sprintf("# input size %zu, %u times [%.1f%%]\n",
                    L, weight, 100*ratio);
        }
        os_pre << "#\n";

        bool do_all_timings = true;

        unsigned int ncalls_at_depth = 0;

        /* This goes to the "preliminary" buffer until we're sure that
         * there's timing to print!
         */
        std::ostream * p_talk = &os_pre;

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
                /* This is a non-call. It's only here in situations where
                 * at depth D, we have x calls of size 1, and then
                 * ((1<<D)-x) calls of size 0, which obviously correspond
                 * to nothing.
                 *
                 * One may question the meaning of the "weight" of these
                 * non-calls. I'd hazard the idea that 0 is most
                 * appropriate.
                 */
                best[L] = { 0, 0 };
                mesh_tt_weighted[0].first = 0; // += weight ?
                continue;
            }

            lingen_call_companion::key K { depth, L };

            /* We **MUST** create hints[K], at this point. We will gain
             * _some_ insight from stored_hints[] and maybe from the
             * thresholds passed on the command line, but in any case we
             * will have to recompute the timings and the RAM usage */

            if (hints.find(K) != hints.end()) {
                /* Then we only have to adjust, I think. At any rate, if
                 * we end up doing hints[K]=U here, we're going to
                 * overwrite the previously computed data, which sounds a
                 * bad idea. */
                hints[K].total_ncalls += weight;
                ncalls_at_depth += weight;
                continue;
            }

            configurations_to_test CF = find_configurations_to_test(os_pre, persist, K);

            std::map<unsigned int, lingen_call_companion> mesh_res;
            std::map<unsigned int, double> mesh_tt;
            std::map<unsigned int, double> mesh_tt_children;

            if (!CF.do_timings) do_all_timings = false;

            if (!CF.do_timings && CF.mesh_sizes.size() != 1)
                throw std::invalid_argument("Cannot have \"notiming\" in tuning_thresholds without specified thresholds for all choices");

            if (CF.do_timings && !timed_something) {
                os << os_pre.str();
                os_pre.str();
                p_talk = &os;
                timed_something = true;
            }

            for(auto mesh : CF.mesh_sizes) {
                lingen_call_companion U;
                U.mesh = mesh;
                if (mesh == 0) {
                    U.ttb = compute_and_report_basecase(*p_talk, L, CF.do_timings);
                    mesh_res[mesh] = U;
                    if (!CF.do_timings) continue;
                    mesh_tt[mesh] = U.ttb;
                    mesh_tt_weighted[mesh].first += weight;
                    mesh_tt_weighted[mesh].second += mesh_tt[mesh] * weight;
                } else {
                    U.mp = tune_local_at_depth_mp_or_mul(
                            *p_talk,
                            persist, cw, depth, mesh,
                            CF.ffts, CF.do_timings,
                            op_mul_or_mp_base::OP_MP);
                    U.mul = tune_local_at_depth_mp_or_mul(
                            *p_talk,
                            persist, cw, depth, mesh,
                            CF.ffts, CF.do_timings,
                            op_mul_or_mp_base::OP_MUL);
                    mesh_res[mesh] = U;
                    if (!CF.do_timings) continue;
                    mesh_tt[mesh] = U.mp.tt.t + U.mul.tt.t;
                    mesh_tt_children[mesh] += best[Lleft].second;
                    mesh_tt_children[mesh] += best[Lright].second;
                    mesh_tt[mesh] += mesh_tt_children[mesh];
                    mesh_tt_weighted[mesh].first += weight;
                    mesh_tt_weighted[mesh].second += mesh_tt[mesh] * weight;
                }
            }

            /* find the best mesh value */
            unsigned int meshbest = UINT_MAX;
            double ttbest = DBL_MAX;
            if (CF.do_timings) {
                for(auto x : mesh_tt) {
                    if (x.second < ttbest) {
                        meshbest = x.first;
                        ttbest = x.second;
                    }
                }
            } else {
                meshbest = CF.mesh_sizes.front();
            }

            lingen_call_companion U = mesh_res[meshbest];
            U.total_ncalls = 0;
            U.complete = true;
            /* we no longer store ttb in the call companions which
             * are intended for recursion */
            hints[K] = U;
            U_typical = U;

            if (CF.do_timings) {
                /* discard mesh values that are less than the winner and
                 * appear to be slow enough that we don't think they'll
                 * ever catch up.
                 */
                for(auto x : mesh_tt) {
                    if (x.first == meshbest || x.first > meshbest)
                        continue;
                    if (x.second >= basecase_keep_until * mesh_tt[meshbest]) {
                        (*p_talk) << fmt::sprintf("# Discarding %s from now on\n",
                                strat_name[x.first]);
                        minimum_mesh = x.first + 1;
                    }
                }
            }

            /* At this point, we started using collective operations,
             * which means that we don't want to go back. A priori.
             * But it's quite difficult indeed because we might see
             * some spurious results at small sizes. */
            // if (mesh > 1 && minimum_mesh <= 1) minimum_mesh = 2;

            best[L] = { meshbest, CF.do_timings ? mesh_tt[meshbest] : -1 };

            hints[K].total_ncalls += weight;
            ncalls_at_depth += weight;

        }

        /* Now give a summary at this level */

        size_t L0 = std::get<0>(cws.front());
        size_t L1 = std::get<0>(cws.back());
        /* calls_and_weights_at_depth must return a sorted list */
        ASSERT_ALWAYS(L0 <= L1);
        unsigned int mesh0 = best[L0].first;
        unsigned int mesh1 = best[L1].first;

        if (do_all_timings) {
            for(auto mesh : mesh_all) {
                if (mesh_tt_weighted.find(mesh) != mesh_tt_weighted.end()) {
                    double tt = mesh_tt_weighted[mesh].second;
                    // std::string rescaled;
                    if (mesh_tt_weighted[mesh].first != (1U << depth)) {
                        double ratio = mesh_tt_weighted[mesh].first / (double) (1U << depth);
                        // rescaled = fmt::sprintf("[rescaled from %.1f%%] ", 100*ratio);
                        tt /= ratio;
                    }
                    (*p_talk) << fmt::sprintf("# %s (%u calls): %.2f [%.1fd]\n",
                            strat_name[mesh],
                            ncalls_at_depth,
                            // rescaled,
                            tt, tt / 86400);
                }
            }

            double tt_total = best[L0].second * std::get<3>(cws.front())
                            + best[L1].second * std::get<3>(cws.back());

            if (mesh0 == mesh1) {
                (*p_talk) << fmt::sprintf("# BEST: %s: %.2f [%.1fd]\n",
                        strat_name[mesh0],
                        tt_total, tt_total / 86400);
            } else {
                (*p_talk) << fmt::sprintf("# BEST: mix of %s and %s: %.2f [%.1fd]\n",
                        strat_name[mesh0],
                        strat_name[mesh1],
                        tt_total, tt_total / 86400);
            }
        }

        if (mesh0 || mesh1) {
            lingen_call_companion U = U_typical;
            if (U.mp.ram_total() > U.mul.ram_total()) {
                (*p_talk) << fmt::sprintf("#   (memory(MP): %s, incl %s reserved)\n",
                        size_disp(U.mp.ram_total()),
                        size_disp(U.mp.reserved_ram));
            } else {
                (*p_talk) << fmt::sprintf("#   (memory(MUL): %s, incl %s reserved)\n",
                        size_disp(U.mul.ram_total()),
                        size_disp(U.mul.reserved_ram));
            }

        }
        return timed_something;
    }/*}}}*/
    lingen_hints tune_local(std::ostream& os, lingen_hints & stored_hints) {/*{{{*/
        size_t N = m*n*L/(m+n);
        char buf[20];

        /* same mechanism in tune_local_at_depth ; output goes to a
         * separate stringstream until we're sure that there's
         * interesting output somewhere.
         */
        std::ostringstream os_pre;
        bool timed_something = false;
        std::ostream * p_talk = &os_pre;

        (*p_talk) << fmt::sprintf("# Measuring lingen data"
                " for N ~ %zu m=%u n=%u"
                " for a %zu-bit prime p,"
                " using a %u*%u grid of %u-thread nodes"
                " [max target RAM = %s]\n",
                N, m, n, mpz_sizeinbase(p, 2),
                P.r, P.r, P.T,
                size_disp(P.available_ram, buf));
#ifdef HAVE_OPENMP
        (*p_talk) << fmt::sprintf("# Note: non-cached basecase measurements"
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

        int fl = log2(L-1) + 1; /* ceil(log_2(L)) */

        tuner_persistent_data persist(stored_hints);

        /* with basecase_keep_until == 0, then we never measure basecase */
        if (basecase_keep_until == 0)
            persist.minimum_mesh = 1;

        if (persist.impose_hints) {
            (*p_talk) << fmt::sprintf("# While we are doing timings here,"
                    " we'll take schedule decisions based on the hints"
                    " found in %s when they apply\n", schedule_filename);
        }

        for(int i = fl ; i>=0 ; i--) {
            bool timed_here = tune_local_at_depth((*p_talk), persist, i);
            if (timed_here && !timed_something) {
                os << os_pre.str();
                os_pre.str();
                p_talk = &os;
                timed_something = true;
            }
        }

        tuner_persistent_data::level_strategy_map & best(persist.best);
        lingen_hints & hints(persist.hints);
        size_t peak(persist.peak);
        int ipeak(persist.ipeak);

        /*****************************************************************/
        /* make the mesh sizes monotonic. In truth, we can only do this
         * safely for mesh size 0 (basecase), since otherwise we would
         * have an inconsistency in the sub-block sizes, which would make
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
                    (*p_talk) << fmt::format(FMT_STRING("## forcing %s at ({}) since it is known to win at ({})\n"), strat_name[y.first],
                            y.first, y.second);
                    x.second.mesh = y.first;
                }
            }
        }

        os << ("#########################  Summary of timing data  #######################\n");
        if (!tuning_thresholds.empty()) {
            std::ostringstream ss;
            for(auto const & x : tuning_thresholds) {
                if (!ss.str().empty()) ss << ",";
                ss << x.first << ':' << x.second;
            }
            os << fmt::sprintf("# Using explicit tuning_thresholds=%s (from command-line)\n", ss.str());
        }
        if (!timed_something) {
            os << "# note: no real timing data is reported, since the notiming\n"
               << "#       threshold has explicitly prevented all measurements\n";
        }

        size_t size_com0;
        double tt_com0;
        std::tie(size_com0, tt_com0) = mpi_threshold_comm_and_time();
        os << fmt::sprintf("# Communication time at lingen_mpi_threshold (%s): %.2f [%.1fd]\n", size_disp(size_com0, buf), tt_com0, tt_com0/86400);
        double time_best = best[L].second;
        if (time_best != -1) {
            time_best += tt_com0;
            os << fmt::sprintf("# Expected total time: %.2f [%.1fd], peak memory %s (at depth %d)\n", time_best, time_best / 86400, size_disp(peak, buf), ipeak);
        }
        hints.ipeak=ipeak;
        hints.peak=peak;
        os << fmt::sprintf("(%u,%u,%u,%.1f,%1.f)\n",m,n,P.r,time_best,(double)peak/1024./1024./1024.);

        /* This one is strictly linear anyway */
        hints.tt_gather_per_unit = tt_com0 / 2 / L;
        hints.tt_scatter_per_unit = tt_com0 / 2 / L;

        return hints;
    }/*}}}*/
    lingen_hints tune(std::ostream& os) {/*{{{*/
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
            hints = tune_local(os, stored_hints);
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

const std::vector<const char *>
lingen_tuner::tuning_thresholds_t::thresholds_verbs
{
    recursive,
    collective,
    // "collective" should allow finer grain, so as to allow
    // multiple mesh sizes.
    // "plain" makes no sense per se, as it's the base thing
    // to use as soon as we go recursive
#ifdef LINGEN_BINARY
    ternary,
    cantor,
#else
    flint,
#endif
    notiming,
};

const std::vector<std::pair<lingen_substep_schedule::fft_type_t, const char *>> lingen_tuner::tuning_thresholds_t::code_to_key
#ifdef LINGEN_BINARY
{
    { lingen_substep_schedule::FFT_NONE, recursive, },
    { lingen_substep_schedule::FFT_TERNARY, ternary, },
    { lingen_substep_schedule::FFT_CANTOR, cantor },
};
#else
{
    { lingen_substep_schedule::FFT_NONE, recursive, },
    { lingen_substep_schedule::FFT_FLINT, flint, },
};
#endif

/* Need this for conformance */
constexpr const char * lingen_tuner::tuning_thresholds_t::recursive;
constexpr const char * lingen_tuner::tuning_thresholds_t::collective;
constexpr const char * lingen_tuner::tuning_thresholds_t::ternary;
constexpr const char * lingen_tuner::tuning_thresholds_t::cantor;
constexpr const char * lingen_tuner::tuning_thresholds_t::flint;
constexpr const char * lingen_tuner::tuning_thresholds_t::notiming;


lingen_hints lingen_tuning(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl)
{
    lingen_tuner::output_info O(pl);
    if (O.quiet) {
        std::ostringstream os;
        return lingen_tuner(os, d, L, comm, pl).tune(os);
    } else if (O.tuning_log_filename) {
        std::ofstream os(O.tuning_log_filename, std::ios_base::out);
        return lingen_tuner(os, d, L, comm, pl).tune(os);
    } else {
        return lingen_tuner(std::cout, d, L, comm, pl).tune(std::cout);
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

