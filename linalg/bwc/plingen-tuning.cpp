#include "cado.h"

#include <cstddef>      /* see https://gcc.gnu.org/gcc-4.9/porting_to.html */
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include <float.h>
#ifdef  HAVE_SIGHUP
#include <signal.h>
#endif
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "mpfq_layer.h"
#include "lingen-polymat.hpp"
#include "lingen-matpoly.hpp"
// #include "lingen-bigpolymat.hpp" // 20150826: deleted.
#include "lingen-matpoly-ft.hpp"
#include "plingen.hpp"
#include "plingen-tuning.hpp"
#include "plingen-tuning-cache.hpp"
#include "lingen-platform.hpp"

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

template<typename>
struct lingen_substep_characteristics;
struct lingen_substep_schedule;
struct plingen_tuner;

size_t plingen_round_operand_size(size_t x, int bits) {/*{{{*/
    /* round x up to the next size that has all but its six most significant
     * bits set to 0.
     */
    if (x == 0) return x;
    x -= 1;
    size_t y = x >> bits;
    for(int i = 1 ; y ; y >>= i) x |= y;
    x += 1;
    return x;
}/*}}}*/

struct op_mul {/*{{{*/
    static constexpr const char * name = "MUL";
    typedef plingen_tuning_cache::mul_key key_type;
    static inline void fti_prepare(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc) {
        fft_get_transform_info_fppol(fti, p, n1, n2, nacc);
    }
    static inline void ift(matpoly & a, matpoly_ft & t, const struct fft_transform_info *)
    {
        ::ift(a, t);
    }
};/*}}}*/
struct op_mp {/*{{{*/
    static constexpr const char * name = "MP";
    typedef plingen_tuning_cache::mp_key key_type;
    static inline void fti_prepare(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc) {
        fft_get_transform_info_fppol_mp(fti, p, nmin, nmax, nacc);
    }
    static inline void ift(matpoly & c, matpoly_ft & tc, const struct fft_transform_info * fti)
    {
        mp_bitcnt_t cbits = fti->ks_coeff_bits;
        unsigned shift = MIN(fti->bits1 / cbits, fti->bits2 / cbits) - 1;
        ::ift_mp(c, tc, shift);
    }
};/*}}}*/

template<typename OP>
struct lingen_substep_characteristics {/*{{{*/
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;
    typedef plingen_tuning_cache tc_t;

    abdst_field ab;
    gmp_randstate_t & rstate;
    cxx_mpz p;

    /* length of the input (E) for the call under consideration ; this is
     * not the input length for the overall algorithm ! */
    size_t input_length;

    /* We're talking products (or any other bilinear operation that uses
     * transforms 1 by 1, e.g. MP) of matrices of size n0*n1 and n1*n2
     */
    unsigned int n0;
    unsigned int n1;
    unsigned int n2;

    /* length of the three operands.
     * A has dimension n0*n1, length asize
     * B has dimension n1*n2, length bsize
     * C has dimension n1*n2, length csize
     */
    size_t asize, bsize, csize;

    private:

    size_t transform_ram;
    struct fft_transform_info fti[1];

    public:

    lingen_substep_characteristics(abdst_field ab, gmp_randstate_t & rstate, size_t input_length, unsigned int n0, unsigned int n1, unsigned int n2, size_t asize, size_t bsize, size_t csize) :/*{{{*/
        ab(ab), rstate(rstate),
        input_length(input_length),
        n0(n0), n1(n1), n2(n2),
        asize(asize), bsize(bsize), csize(csize)
    {
        abfield_characteristic(ab, p);
        OP::fti_prepare(fti, p, asize, bsize, n1);
        size_t fft_alloc_sizes[3];
        fft_get_transform_allocs(fft_alloc_sizes, fti);
        transform_ram = fft_alloc_sizes[0];
    }/*}}}*/

    bool has_cached_time(tc_t const & C) const {
        typename OP::key_type K { mpz_sizeinbase(p, 2), asize, bsize };
        return C.has(K);
    }

    std::array<double, 4> get_ft_times(tc_t & C) const {/*{{{*/
        typename OP::key_type K { mpz_sizeinbase(p, 2), asize, bsize };

        if (C.has(K))
            return C[K];

        double tt_dft0;
        double tt_dft2;
        double tt_conv;
        double tt_ift;

        /* make all of these 1*1 matrices, just for timing purposes */
        matpoly a(ab, 1, 1, asize);
        matpoly b(ab, 1, 1, bsize);
        a.fill_random(asize, rstate);
        b.fill_random(bsize, rstate);

        matpoly c(ab, a.m, b.n, csize);

        matpoly_ft ta(ab, a.m, a.n, fti);
        matpoly_ft tb(ab, b.m, b.n, fti);
        matpoly_ft tc(ab, a.m, b.n, fti);

        double tt = 0;

        dft(ta, a);     /* warm up */
        tt = -wct_seconds();
        dft(ta, a);
        tt_dft0 = wct_seconds() + tt;

        dft(tb, b);     /* warm up */
        tt = -wct_seconds();
        dft(tb, b);
        tt_dft2 = wct_seconds() + tt;

        mul(tc, ta, tb); /* warm up */
        tt = -wct_seconds();
        mul(tc, ta, tb);
        tt_conv = wct_seconds() + tt;

        c.size = csize;
        ASSERT_ALWAYS(c.size <= c.alloc);
        OP::ift(c, tc, fti);    /* warm up */
        tt = -wct_seconds();
        OP::ift(c, tc, fti);
        tt_ift = wct_seconds() + tt;

        C[K] = { tt_dft0, tt_dft2, tt_conv, tt_ift };

        return C[K];
    }/*}}}*/

    size_t get_transform_ram() const { return transform_ram; }
    std::tuple<size_t, size_t, size_t> get_operand_ram() const {
        return {
            n0 * n1 * asize * mpz_size(p) * sizeof(mp_limb_t),
            n1 * n2 * bsize * mpz_size(p) * sizeof(mp_limb_t),
            n0 * n2 * csize * mpz_size(p) * sizeof(mp_limb_t) };
    }

    size_t get_peak_ram(pc_t const & P, sc_t const & S) const { /* {{{ */
        unsigned int nr0 = iceildiv(n0, P.r);
        unsigned int nr2 = iceildiv(n2, P.r);

        unsigned int nrs0 = iceildiv(nr0, S.shrink0);
        unsigned int nrs2 = iceildiv(nr2, S.shrink2);

        return get_transform_ram() * (S.batch * P.r * (nrs0 + nrs2) + nrs0*nrs2);
    }/*}}}*/

    private:
    /* This returns the communication time _per call_ */
    double get_comm_time(pc_t const & P, sc_t const & S) const { /* {{{ */
        /* TODO: maybe include some model of latency... */

        /* Each of the r*r nodes has local matrices size (at most)
         * nr0*nr1, nr1*nr2, and nr0*nr2. For the actual computations, we
         * care more about the shrinked submatrices, though */
        unsigned int nr0 = iceildiv(n0, P.r);
        unsigned int nr1 = iceildiv(n1, P.r);
        unsigned int nr2 = iceildiv(n2, P.r);

        /* The shrink parameters will divide the size of the local
         * matrices we consider by numbers shrink0 and shrink2. This
         * increases the time, and decreases the memory footprint */
        unsigned int nrs0 = iceildiv(nr0, S.shrink0);
        unsigned int nrs2 = iceildiv(nr2, S.shrink2);
        unsigned int ns0 = P.r * nrs0;
        unsigned int ns2 = P.r * nrs2;
        size_t comm = get_transform_ram() * (ns0+ns2) * nr1;
        comm *= S.shrink0 * S.shrink2;
        return comm;
    }/*}}}*/

    public:
    double get_call_time(pc_t const & P, sc_t const & S, tc_t & C) const { /* {{{ */
        auto ft = get_ft_times(C);
        double tt_dft0 = ft[0];
        double tt_dft2 = ft[1];
        double tt_conv = ft[2];
        double tt_ift = ft[3];

        double T_dft0, T_dft2, T_conv, T_ift, T_comm;


        /* Each of the r*r nodes has local matrices size (at most)
         * nr0*nr1, nr1*nr2, and nr0*nr2. For the actual computations, we
         * care more about the shrinked submatrices, though */
        unsigned int nr0 = iceildiv(n0, P.r);
        unsigned int nr1 = iceildiv(n1, P.r);
        unsigned int nr2 = iceildiv(n2, P.r);

        /* The shrink parameters will divide the size of the local
         * matrices we consider by numbers shrink0 and shrink2. This
         * increases the time, and decreases the memory footprint */
        unsigned int nrs0 = iceildiv(nr0, S.shrink0);
        unsigned int nrs2 = iceildiv(nr2, S.shrink2);
        // unsigned int ns0 = r * nrs0;
        // unsigned int ns2 = r * nrs2;

        unsigned int nr1b = iceildiv(nr1, S.batch);

        /* The 1-threaded timing will be obtained by setting T=batch=1.  */

        /* Locally, we'll add together r matrices of size nrs0*nrs2,
         * which will be computed as products of vectors of size
         * nrs0*batch and batch*nrs2. The operation will be repeated
         * nr1/batch * times.
         */

        /* These are just base values, we'll multiply them later on */
        T_dft0 = nr1b * tt_dft0;
        T_dft2 = nr1b * tt_dft2;
        T_conv = nr1b * tt_conv;
        T_ift  = tt_ift;

        /* First, we must decide on the scheduling order for the loop.
         * Either we compute, collect, and store all transforms on side
         * 0, and then use the index on side 0 as the inner loop (so that
         * the transforms on side 0 stay live much longer), or the
         * converse. This depends on which size has the larger number of
         * transforms.
         */

        if (nrs0 < nrs2) {
            /* then we want many live transforms on side 0. */
            T_dft0 *= iceildiv(nrs0 * S.batch, P.T);
            T_dft2 *= nrs2 * iceildiv(S.batch, P.T);
            T_conv *= nrs2 * P.r * S.batch * iceildiv(nrs0, P.T);
        } else {
            /* then we want many live transforms on side 2. */
            /* typical loop (on each node)
                allocate space for nrs0 * nrs2 transforms.
                for(k = 0 ; k < nr1 ; k += batch)
                    allocate space for batch * nrs2 * r transforms for b_(k..k+batch-1,j) (this takes into account the * r multiplier that is necessary for the received data).
                    compute batch * nrs2 transforms locally. This can be done in parallel
                    share batch * nrs2 transforms with r peers.
                    for(i = 0 ; i < nrs0 ; i++)
                        allocate space for batch * r transforms for a_(i,k..k+batch-1) (this takes into account the * r multiplier that is necessary for the received data).
                        compute batch transforms locally. This can be done in parallel
                        share batch transforms with r peers.
                        compute and accumulate r*batch*nrs2 convolution products into only nrs2 destination variables (therefore r*batch passes are necessary in order to avoid conflicts).
                        free batch * r transforms
                    free batch * nrs2 * r transforms
                compute nrs0 * nrs2 inverse transforms locally, in parallel
                free nrs0 * nrs2 transforms.
            */
            T_dft0 *= nrs0 * iceildiv(S.batch, P.T);
            T_dft2 *= iceildiv(nrs2 * S.batch, P.T);
            T_conv *= nrs0 * P.r * S.batch * iceildiv(nrs2, P.T);
        }

        T_ift *= iceildiv(nrs0*nrs2, P.T);

        T_dft0 *= S.shrink0 * S.shrink2;
        T_dft2 *= S.shrink0 * S.shrink2;
        T_conv *= S.shrink0 * S.shrink2;
        T_ift  *= S.shrink0 * S.shrink2;

        T_comm = get_comm_time(P, S) / P.mpi_xput;

        /* This is for _one_ call only (one node in the tree).
         * We must apply multiplier coefficients (typically 1<<depth) */
        return T_dft0 + T_dft2 + T_conv + T_ift + T_comm;
    }/*}}}*/

    double get_and_report_call_time(pc_t const & P, sc_t const & S, tc_t & C) const { /* {{{ */
        const char * step = OP::name;
        bool cached = has_cached_time(C);
        char buf[20];

        printf("# %s(@%zu) [shrink=(%u,%u) batch=%u] %s, ",
                step,
                input_length,
                S.shrink0,
                S.shrink2,
                S.batch,
                size_disp(get_peak_ram(P, S), buf));
        fflush(stdout);
        double tt = get_call_time(P, S, C);
        printf("%.2f%s\n",
                tt,
                cached ? " [from cache]" : "");
        return tt;
    }/*}}}*/

    void report_size_stats_human() const {/*{{{*/
        char buf[4][20];
        const char * step = OP::name;

        printf("# %s (per op): %s+%s+%s, transforms 3*%s\n",
                step,
                size_disp(asize*mpz_size(p)*sizeof(mp_limb_t), buf[0]),
                size_disp(bsize*mpz_size(p)*sizeof(mp_limb_t), buf[1]),
                size_disp(csize*mpz_size(p)*sizeof(mp_limb_t), buf[2]),
                size_disp(transform_ram, buf[3]));
        printf("# %s (total for %u*%u * %u*%u): %s, transforms %s\n",
                step,
                n0,n1,n1,n2,
                size_disp((n0*n1*asize+n1*n2*bsize+n0*n2*csize)*mpz_size(p)*sizeof(mp_limb_t), buf[0]),
                size_disp((n0*n1+n1*n2+n0*n2)*transform_ram, buf[1]));
    }/*}}}*/
};/*}}}*/

template<typename OP>
void optimize(lingen_substep_schedule & S, lingen_substep_characteristics<OP> const & U, lingen_platform const & P, size_t reserved) { /* {{{ */
        unsigned int nr0 = iceildiv(U.n0, P.r);
        unsigned int nr1 = iceildiv(U.n1, P.r);
        unsigned int nr2 = iceildiv(U.n2, P.r);

        lingen_substep_schedule res = S;

        for( ; S.batch < nr1 && (reserved + U.get_peak_ram(P, S)) <= P.available_ram ; ) {
            S = res;
            res.batch++;
        }

        for( ; (reserved + U.get_peak_ram(P, S)) > P.available_ram ; ) {
            unsigned int nrs0 = iceildiv(nr0, S.shrink0);
            unsigned int nrs2 = iceildiv(nr2, S.shrink2);
            if (nrs0 < nrs2 && S.shrink2 < nr2) {
                S.shrink2++;
            } else if (S.shrink0 < nr0) {
                S.shrink0++;
            } else {
                char buf[20];
                std::ostringstream os;
                os << "Fatal error:"
                    << " it is not possible to complete this calculation with only "
                    << size_disp(P.available_ram, buf)
                    << " of memory for intermediate transforms.\n";
                os << "Based on the cost for input length "
                    << U.input_length
                    << ", we need at the very least "
                    << size_disp(U.get_peak_ram(P, S), buf);
                os  << ", plus "
                    << size_disp(reserved, buf)
                    << " for reserved memory at upper levels";
                os << " [with shrink=(" << S.shrink0 << "," << S.shrink2
                    << "), batch=" << S.batch << "]\n";
                fputs(os.str().c_str(), stderr);
                exit(EXIT_FAILURE);
            }
        }
    }
    /* }}} */

struct plingen_tuner {
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;

    /* imported from the dims struct */
    abdst_field ab;
    unsigned int m,n;

    size_t L;

    lingen_platform P;

    plingen_tuning_cache C;

    cxx_mpz p;

    gmp_randstate_t rstate;

    const char * timing_cache_filename = NULL;

    /* stop measuring the time taken by the basecase when it is
     * more than this number times the time taken by the other
     * alternatives
     */
    double basecase_keep_until = 1.8;

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

    plingen_tuner(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl) :
        ab(d.ab), m(d.m), n(d.n), L(L), P(comm, pl)
    {
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, 1);
        abfield_characteristic(ab, p);

        param_list_parse_double(pl, "basecase-keep-until", &basecase_keep_until);

        timing_cache_filename = param_list_lookup_string(pl, "tuning_timing_cache_filename");
        C.load(timing_cache_filename);
    }

    ~plingen_tuner() {
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
        return { 2 * data0, 2 * data0 / P.mpi_xput};
    }/*}}}*/

    double compute_and_report_basecase(size_t length) { /*{{{*/
        extern void test_basecase(abdst_field ab, unsigned int m, unsigned int n, size_t L, gmp_randstate_t rstate);

        double tt;

        plingen_tuning_cache::basecase_key K { mpz_sizeinbase(p, 2), m, n, length, P.openmp_threads };

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
            std::vector<weighted_call_t> res
            {{
                 { 2*Q + u,     Q + u, Q,     (1 << i) - v },
                 { 2*Q + u + 1, Q + 1, Q + u, v },
            }};
            return res;
        } else {
            std::vector<weighted_call_t> res
            {{
                 { 2*Q + u,     Q + u, Q,     (1 << i) - v },
            }};
            return res;
        }
    }

    lingen_substep_characteristics<op_mp> mp_substep(weighted_call_t const & cw) { /* {{{ */
        size_t length_E;
        size_t length_E_left;
        size_t length_E_right;
        unsigned int weight;

        std::tie(length_E, length_E_left, length_E_right, weight) = cw;

        ASSERT_ALWAYS(length_E >= 2);
        ASSERT_ALWAYS(weight);

        size_t csize = length_E_right;
        size_t bsize = iceildiv(m * length_E_left, m+n);
        size_t asize = csize + bsize - 1;

        ASSERT_ALWAYS(asize);
        ASSERT_ALWAYS(bsize);

        return lingen_substep_characteristics<op_mp>(
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
            lingen_substep_schedule S;
            optimize(S, step, P, reserved);
            if (print && !printed_mem_once++) {
                step.report_size_stats_human();
                step.get_and_report_call_time(P, S, C);
            } else {
                step.get_call_time(P, S, C);
            }
            schedules_mp[L] = S;
        }
    } /* }}} */
    lingen_substep_characteristics<op_mul> mul_substep(weighted_call_t const & cw) { /* {{{ */
        size_t length_E;
        size_t length_E_left;
        size_t length_E_right;
        unsigned int weight;

        std::tie(length_E, length_E_left, length_E_right, weight) = cw;

        ASSERT_ALWAYS(length_E >= 2);
        ASSERT_ALWAYS(weight);

        size_t asize = iceildiv(m * length_E_left, m+n);
        size_t bsize = iceildiv(m * length_E_right, m+n);
        size_t csize = asize + bsize - 1;

        ASSERT_ALWAYS(asize);
        ASSERT_ALWAYS(bsize);

        return lingen_substep_characteristics<op_mul>(
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

                lingen_substep_schedule S;
                optimize(S, step, pp, reserved);
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

    lingen_hints_t tune_local() {
        size_t N = m*n*L/(m+n);
        char buf[20];
        printf("# Measuring lingen data for N ~ %zu m=%u n=%u for a %zu-bit prime p, using a %u*%u grid of %u-thread nodes [max target RAM = %s]\n",
                N, m, n, mpz_sizeinbase(p, 2),
                P.r, P.r, P.T,
                size_disp(P.available_ram, buf));
#ifdef HAVE_OPENMP
        printf("# Note: non-cached basecase measurements are done using openmp as it is configured for the running code, that is, with %d threads\n", P.openmp_threads);
#endif

        lingen_hints_t hints;

        int fl = log2(L) + 1;

        /* with basecase_keep_until == 0, then we never measure basecase */
        bool basecase_eliminated = basecase_keep_until == 0;
        std::map<size_t, std::tuple<bool, std::array<double, 3> >, plingen_tuning_cache::coarse_compare> best;
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
             * plingen.cpp
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

            for(size_t idx = 0 ; idx < cws.size() ; idx++) {
                auto const & cw(cws[idx]);
                size_t L, Lleft, Lright;
                unsigned int weight;
                std::tie(L, Lleft, Lright, weight) = cw;

                ASSERT_ALWAYS(weight);

                if (!L) {
                    best[L] = { false, { 0, 0, 0 }};
                    continue;
                }

                if (best.find(L) == best.end()) {
                    double ttb = DBL_MAX;
                    double ttr = DBL_MAX;
                    double ttrchildren = DBL_MAX;

                    lingen_call_companion::key K { i, L };
                    lingen_call_companion U;
                    U.weight = weight;

                    if (!recursion_makes_sense(L) || !basecase_eliminated)
                        ttb = compute_and_report_basecase(L);

                    if (recursion_makes_sense(L)) {
                        auto MUL = mul_substep(cw);
                        auto MP  = mp_substep(cw);

                        U.mp.S = schedules_mp[L];
                        U.mul.S = schedules_mul[L];

                        U.mp.tt = MP.get_call_time(P, U.mp.S, C);
                        U.mul.tt = MUL.get_call_time(P, U.mul.S, C);

                        ttr = U.mp.tt + U.mul.tt;
                        ttrchildren = 0;
                        ttrchildren += std::get<1>(best[Lleft])[std::get<0>(best[Lleft])];
                        ttrchildren += std::get<1>(best[Lright])[std::get<0>(best[Lright])];

                        size_t m;
                        m = U.mp.ram = MP.get_peak_ram(P, schedules_mp[L]);
                        m += U.mp.reserved_ram = reserved_mp;
                        if (m > ram_mp) ram_mp = m;
                        if (m > peak) { ipeak = i; peak = m; }
                        m = U.mul.ram = MUL.get_peak_ram(P, schedules_mul[L]);
                        m += U.mul.reserved_ram = reserved_mul;
                        if (m > ram_mul) ram_mul = m;
                        if (m > peak) { ipeak = i; peak = m; }
                    }

                    if (ttb >= basecase_keep_until * (ttr + ttrchildren))
                        basecase_eliminated = true;

                    bool rwin = ttb >= (ttr + ttrchildren);
                    best[L] = { rwin, {ttb, ttr + ttrchildren, ttr} };

                    U.recurse = rwin;
                    /* See comment in compute_schedules_for_mul.
                     * Presently we don't identify cases where
                     * lingen_threshold makes sense at all */
                    U.go_mpi = rwin;
                    U.ttb = ttb;

                    hints[K] = U;
                }

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
            size_t L0r = plingen_round_operand_size(L0);
            size_t L1r = plingen_round_operand_size(L1);
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
                printf("#%*s basecase(threshold>%zu): %.2f [%.1fd]%s\n", pad, msg, L1,
                        time_b, time_b / 86400, isbest);
                msg = msg2;
            }
            if (!approx_same && recursion_makes_sense(L1)) {
                const char * isbest = (rec0 && !rec1) ? strbest : "";
                std::ostringstream os2;
                os2 << " mixed(threshold=" << L1 << "): ";
                std::string ss2 = os2.str();
                printf("#%*s%s%.2f [%.1fd] (self: %.2f [%.1fd])%s\n", pad, msg, ss2.c_str(),
                        time_m, time_m / 86400, time_m_self, time_m_self / 86400, isbest);
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
            } else if (rec0 && !rec1) {
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
        printf("# Automatically tuned lingen_mpi_threshold=%zu\n", upper_threshold);
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
    lingen_hints_t tune() {
        int rank;
        MPI_Comm_rank(P.comm, &rank);
        lingen_hints_t hints;

        if (rank == 0)
            hints = tune_local();

        hints.share(0, P.comm);

        return hints;
    }
};

template<typename T>
    typename std::enable_if<std::is_trivially_copyable<typename T::mapped_type>::value && std::is_trivially_copyable<typename T::mapped_type>::value, void>::type
share(T & m, int root, MPI_Comm comm)
{
    typedef typename T::key_type K;
    typedef typename T::mapped_type M;
    typedef std::pair<K, M> V;

    int rank;
    MPI_Comm_rank(comm, &rank);
    std::vector<V> serial;
    serial.insert(serial.end(), m.begin(), m.end());
    unsigned long ns = serial.size();
    MPI_Bcast(&ns, 1, MPI_UNSIGNED_LONG, root, comm);
    if (rank) serial.assign(ns, V());
    MPI_Bcast(&serial.front(), ns * sizeof(V), MPI_BYTE, 0, comm);
    if (rank) m.insert(serial.begin(), serial.end());
}
void lingen_hints_t::share(int root, MPI_Comm comm)
{
    ::share((super&) *this, root, comm);
    MPI_Bcast(&tt_gather_per_unit, 1, MPI_DOUBLE, root, comm);
    MPI_Bcast(&tt_scatter_per_unit, 1, MPI_DOUBLE, root, comm);
}

lingen_hints_t plingen_tuning(bw_dimensions & d, size_t L, MPI_Comm comm, cxx_param_list & pl)
{
    return plingen_tuner(d, L, comm, pl).tune();
}

void plingen_tuning_decl_usage(cxx_param_list & pl)
{
    plingen_tuner::declare_usage(pl);
}

void plingen_tuning_lookup_parameters(cxx_param_list & pl)
{
    plingen_tuner::lookup_parameters(pl);
}

