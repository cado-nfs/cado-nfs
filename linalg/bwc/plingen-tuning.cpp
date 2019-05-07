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
#include "lingen-polymat.h"
#include "lingen-matpoly.h"
// #include "lingen-bigpolymat.h" // 20150826: deleted.
#include "lingen-matpoly-ft.h"
#include "plingen.h"
#include "plingen-tuning.hpp"
#include "plingen-tuning-cache.hpp"

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
struct lingen_platform_characteristics;
struct lingen_substep_scheduling_characteristics;
struct plingen_tuner;

struct op_mul {/*{{{*/
    static constexpr const char * name = "MUL";
    typedef plingen_tuning_cache::mul_key key_type;
    static inline void fti_prepare(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc) {
        fft_get_transform_info_fppol(fti, p, n1, n2, nacc);
    }
    static inline double ift(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, const struct fft_transform_info * fti, int draft)
    {
        return matpoly_ft_ift(ab, a, t, fti, draft);
    }
};/*}}}*/
struct op_mp {/*{{{*/
    static constexpr const char * name = "MP";
    typedef plingen_tuning_cache::mp_key key_type;
    static inline void fti_prepare(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc) {
        fft_get_transform_info_fppol_mp(fti, p, nmin, nmax, nacc);
    }
    static inline double ift(abdst_field ab, matpoly_ptr c, matpoly_ft_ptr tc, const struct fft_transform_info * fti, int draft)
    {
        mp_bitcnt_t cbits = fti->ks_coeff_bits;
        unsigned shift = MIN(fti->bits1 / cbits, fti->bits2 / cbits) - 1;
        return matpoly_ft_ift_mp(ab, c, tc, shift, fti, draft);
    }
};/*}}}*/

struct lingen_platform_characteristics {/*{{{*/
    /* input characteristics -- the ones we have to live with */

    /* We give timings for a run on r*r nodes, with T threads per node */
    /* Note that all of this can also be auto-detected, a priori */
    unsigned int r;
    unsigned int T;     /* **PHYSICAL** cores, or we say rubbish  */
    int openmp_threads;

    size_t available_ram;

    /* Assume we output something like one gigabyte per second. This is
     * rather conservative for HPC networks */
    static constexpr double mpi_xput = 1e9;

    static void lookup_parameters(cxx_param_list & pl) {
        param_list_lookup_string(pl, "max_ram");
        param_list_lookup_string(pl, "tuning_T");
        param_list_lookup_string(pl, "tuning_r");
        param_list_lookup_string(pl, "mpi");
        param_list_lookup_string(pl, "thr");
    }

    static void declare_usage(cxx_param_list & pl) {
        /* TODO: this shall supersede mpi= and thr= that are currently
         * parsed from within plingen.cpp */
        param_list_decl_usage(pl, "max_ram",
                "Maximum local memory to be used for transforms and matrices, in GB");
        param_list_decl_usage(pl, "tuning_T",
                "For --tune only: target number of threads (if for different platform)");
        param_list_decl_usage(pl, "tuning_r",
                "For --tune only: size of the mpi grid (the grid would be r times r)");
    }

    lingen_platform_characteristics(cxx_param_list & pl) {/*{{{*/

        int mpi[2];
        int thr[2];

        param_list_parse_intxint(pl, "mpi", mpi);
        param_list_parse_intxint(pl, "thr", thr);

        T = thr[0] * thr[1];
        r = mpi[0];

        double dtmp = 1;
        param_list_parse_double(pl, "max_ram", &dtmp);
        available_ram = dtmp * (1 << 30);

        param_list_parse_uint(pl, "tuning_T", &T);
        param_list_parse_uint(pl, "tuning_r", &r);

        openmp_threads = 1;
#ifdef HAVE_OPENMP
        openmp_threads = omp_get_max_threads();
#endif
    }/*}}}*/
};/*}}}*/

struct lingen_substep_scheduling_characteristics {/*{{{*/
    typedef lingen_platform_characteristics pc_t;
    typedef lingen_substep_scheduling_characteristics sc_t;

    /* output characteristics -- the ones we have to choose */

    /* batch, shrink0, shrink2. Used to be static parameters, now they're
     * dynamic.
     *
     * batch: with batch=b, schedule b times more simultaneous
     *  transforms, so that we can have more parallelism. Costs b times
     *  more RAM
     * shrink0: divide the number of transforms on dimension 0 (m/n for
     *  MP, (m+n)/r for MUL) so as to save memory
     * shrink2: divide the number of transforms on dimension 2 (m+n)/r
     *  for both MP and MUL) so as to save memory
     */

    /* shrink0 and shrink2 are important. Here, we restrict our operation
     * to shrink0*shrink2 multiplications of matrices of size
     * (n0/shrink0) times (n2/shrink2).
     */
    unsigned int shrink0 = 1, shrink2 = 1;

    /* batch is not used much, because it goes in the direction of
     * spending _more_ memory, which is rarely what we're interested in.
     * What this does is that several steps of the outer loop (off size
     * n1/r) are done simultaneously: more transforms are kept live.
     */
    unsigned int batch = 1;

    template<typename OP>
    void optimize(lingen_substep_characteristics<OP> const & U, pc_t const & P) { /* {{{ */
        unsigned int nr0 = iceildiv(U.n0, P.r);
        unsigned int nr1 = iceildiv(U.n1, P.r);
        unsigned int nr2 = iceildiv(U.n2, P.r);

        sc_t res = *this;

        for( ; batch < nr1 && U.get_peak_ram(P, *this) <= P.available_ram ; ) {
            *this = res;
            res.batch++;
        }

        for( ; U.get_peak_ram(P, *this) > P.available_ram ; ) {
            unsigned int nrs0 = iceildiv(nr0, shrink0);
            unsigned int nrs2 = iceildiv(nr2, shrink2);
            if (nrs0 < nrs2 && shrink2 < nr2) {
                shrink2++;
            } else if (shrink0 < nr0) {
                shrink0++;
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
                    << size_disp(U.get_peak_ram(P, *this), buf);
                os << " [with shrink=(" << shrink0 << "," << shrink2
                    << "), batch=" << batch << "]\n";
                fputs(os.str().c_str(), stderr);
                exit(EXIT_FAILURE);
            }
        }
    }
    /* }}} */
};/*}}}*/

template<typename OP>
struct lingen_substep_characteristics {/*{{{*/
    typedef lingen_platform_characteristics pc_t;
    typedef lingen_substep_scheduling_characteristics sc_t;
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

    /* it's easier to compute it at ctor time (TODO) */
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
        matpoly a, b, c;
        matpoly_init(ab, a, 1, 1, asize);
        matpoly_init(ab, b, 1, 1, bsize);
        matpoly_init(ab, c, 1, 1, csize);
        matpoly_fill_random(ab, a, asize, rstate);
        matpoly_fill_random(ab, b, bsize, rstate);

        matpoly_ft tc, ta, tb;
        matpoly_clear(ab, c);
        matpoly_init(ab, c, a->m, b->n, csize);

        matpoly_ft_init(ab, ta, a->m, a->n, fti);
        matpoly_ft_init(ab, tb, b->m, b->n, fti);
        matpoly_ft_init(ab, tc, a->m, b->n, fti);

        double tt = 0;

        tt = -wct_seconds();
        matpoly_ft_dft(ab, ta, a, fti, 0);
        tt_dft0 = wct_seconds() + tt;

        tt = -wct_seconds();
        matpoly_ft_dft(ab, tb, b, fti, 0);
        tt_dft2 = wct_seconds() + tt;

        tt = -wct_seconds();
        matpoly_ft_mul(ab, tc, ta, tb, fti, 0);
        tt_conv = wct_seconds() + tt;

        tt = -wct_seconds();
        c->size = csize;
        ASSERT_ALWAYS(c->size <= c->alloc);
        OP::ift(ab, c, tc, fti, 0);
        tt_ift = wct_seconds() + tt;

        matpoly_ft_clear(ab, ta, fti);
        matpoly_ft_clear(ab, tb, fti);
        matpoly_ft_clear(ab, tc, fti);

        C[K] = { tt_dft0, tt_dft2, tt_conv, tt_ift };

        matpoly_clear(ab, a);
        matpoly_clear(ab, b);
        matpoly_clear(ab, c);

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
        double tt=get_call_time(P, S, C);
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


/* This object is passed as a companion info to a call
 * of bw_biglingen_recursive ; it is computed by the code in this file,
 * but once tuning is over, it is essentially fixed.
 */
#if 0
struct matrix_product_schedule {

    /* recursion depth */
    unsigned int depth;

    matrix_product_schedule_with_transforms(size_t input_length, const char * step, tcache_key const & K)
        : input_length(input_length), step(step), K(K) {}

#if 0
    void report_time_stats_human(double reference=0) const {/*{{{*/
        double T_dft0, T_dft2, T_conv, T_ift, T_comm;
        std::tie(T_dft0, T_dft2, T_conv, T_ift, T_comm) = global;
        double T_total = T_dft0 + T_dft2 + T_conv + T_ift + T_comm;
        char buf[80] = { '\0' };

        if (reference) {
            snprintf(buf, sizeof(buf), " [efficiency %.1f%%]", 100.0 * reference / (r * r * T * T_total));
        }

        printf("# total %s %u-threaded time on %u nodes [shrink=(%u,%u) batch=%u]:"
                " %.1f + %.1f + %.1f + %.1f + %.1f = %.1f [%.1f days] %s\n",
                step,
                T,
                r*r,
                shrink0,
                shrink2,
                batch,
                T_dft0, T_dft2, T_conv, T_ift, T_comm,
                T_total, T_total/86400, buf);
    }/*}}}*/
#endif

    void report_comm_and_peakmem_human() const {/*{{{*/
        char buf[20];
        size_t comm = compute_derived_comm();
        printf("# %s comm per node %s\n", step, size_disp(comm, buf));
        size_t peakram = compute_derived_peak_ram();
        printf("# %s peak memory [shrink=(%u,%u) batch=%u]: %s\n",
                step, shrink0, shrink2, batch, size_disp(peakram, buf));
    }/*}}}*/
#if 0
    void report_stats_parsable() const {/*{{{*/
        /* This data will go to a timing cache file, so that we can
         * easily rerun the program.
         */
        double tt_dft0;
        double tt_dft2;
        double tt_conv;
        double tt_ift;
        double dummy_tt_comm; /* zero */
        char buf[20];

        std::tie(tt_dft0, tt_dft2, tt_conv, tt_ift, dummy_tt_comm) = tcache[K];

        size_t peakram = compute_derived_peak_ram();

        printf("#;depth;recursive_input_length;step;step_op1_length;step_op2_length;size_one_transform;tt_dft0;tt_dft2;tt_conv;tt_ift\n");
        printf(";%d;%zu;%s;%zu;%zu;%s;%.2f;%.2f;%.2f;%.2f\n",
                depth,
                input_length,
                step, asize, bsize, 
                size_disp(peakram, buf),
                tt_dft0, tt_dft2, tt_conv, tt_ift);
    }/*}}}*/
#endif
};
#endif

struct plingen_tuner {
    typedef lingen_platform_characteristics pc_t;
    typedef lingen_substep_scheduling_characteristics sc_t;

    /* imported from the dims struct */
    abdst_field ab;
    unsigned int m,n;

    size_t N;
    size_t L;

    lingen_platform_characteristics P;

    plingen_tuning_cache C;

    cxx_mpz p;

    gmp_randstate_t rstate;

    const char * timing_cache_filename = NULL;

    /* stop measuring the time taken by the basecase when it is
     * more than this number times the time taken by the other operations
     */
    unsigned int basecase_keep_until = 4;

    std::map<size_t, lingen_substep_scheduling_characteristics> schedules_mp, schedules_mul;

    static void declare_usage(cxx_param_list & pl) {/*{{{*/
        lingen_platform_characteristics::declare_usage(pl);
        param_list_decl_usage(pl, "tuning_N",
                "For --tune only: target size for the corresponding matrix");
        param_list_decl_usage(pl, "tuning_timing_cache_filename",
                "For --tune only: save (and re-load) timings for individual transforms in this file\n");
        param_list_decl_usage(pl, "basecase-keep-until",
                "For --tune only: stop measuring basecase timing when it exceeds the time of the recursive algorithm (counting its leaf calls) by this factor\n");
    }/*}}}*/

    static void lookup_parameters(cxx_param_list & pl) {/*{{{*/
        lingen_platform_characteristics::lookup_parameters(pl);
        param_list_lookup_string(pl, "tuning_N");
        param_list_lookup_string(pl, "tuning_timing_cache_filename");
        param_list_lookup_string(pl, "basecase-keep-until");
    }/*}}}*/

    plingen_tuner(dims * d, cxx_param_list & pl) :
        ab(d->ab), m(d->m), n(d->n), P(pl)
    {
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, 1);
        abfield_characteristic(ab, p);

        param_list_parse_uint(pl, "basecase-keep-until", &basecase_keep_until);

        /* TODO */
        N = 37000000;
        param_list_parse_size_t(pl, "tuning_N", &N);
        L = N/n + N/m;

        char buf[20];
        printf("# Measuring lingen data for N=%zu m=%u n=%u for a %zu-bit prime p, using a %u*%u grid of %u-thread nodes [max target RAM = %s]\n",
                N, m, n, mpz_sizeinbase(p, 2),
                P.r, P.r, P.T,
                size_disp(P.available_ram, buf));
#ifdef HAVE_OPENMP
        printf("# Note: non-cached basecase measurements are done using openmp as it is configured for the running code, that is, with %d threads\n", P.openmp_threads);
#endif

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
    void compute_schedules_for_mp(int i) { /* {{{ */
        int printed_mem_once=0;
        for(auto const & cw : calls_and_weights_at_depth(i)) {
            size_t L = std::get<0>(cw);
            if (!recursion_makes_sense(L)) continue;
            auto step = mp_substep(cw);
            if (!printed_mem_once++)
                step.report_size_stats_human();
            lingen_substep_scheduling_characteristics S;
            S.optimize(step, P);
            step.get_and_report_call_time(P, S, C);
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

    void compute_schedules_for_mul(int i) { /* {{{ */
        int printed_mem_once=0;
        for(auto const & cw : calls_and_weights_at_depth(i)) {
            size_t L = std::get<0>(cw);
            if (!recursion_makes_sense(L)) continue;
            auto step = mul_substep(cw);
            if (!printed_mem_once++)
                step.report_size_stats_human();
            lingen_substep_scheduling_characteristics S;
            S.optimize(step, P);
            step.get_and_report_call_time(P, S, C);
            schedules_mul[L] = S;
        }
    } /* }}} */

    void tune() {
        int fl = log2(L) + 1;

        /* with basecase_keep_until == 0, then we never measure basecase */
        bool basecase_eliminated = !basecase_keep_until;
        std::map<size_t, std::tuple<bool, std::array<double, 3> >, plingen_tuning_cache::coarse_compare> best;
        size_t upper_threshold = SIZE_MAX;

        for(int i = fl ; i>=0 ; i--) {
            auto cws = calls_and_weights_at_depth(i);

            // printf("########## Measuring time at depth %d ##########\n", i);

            /* calls_and_weights_at_depth must return a sorted list */
            // bool same = L0 == L1;
            // bool approx_same = L0r == L1r;
            // unsigned int weight0 = std::get<3>(cws.first());
            // unsigned int weight1 = (1u << i) - weight0;

            double time_b = 0;
            double time_r = 0;
            double time_m = 0;
            double time_r_self = 0;
            double time_m_self = 0;

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

                    if (!recursion_makes_sense(L) || !basecase_eliminated)
                        ttb = compute_and_report_basecase(L);

                    if (recursion_makes_sense(L)) {
                        auto MUL = mul_substep(cw);
                        auto MP  = mp_substep(cw);
                        ttr = MP.get_call_time(P, schedules_mp[L], C)
                            + MUL.get_call_time(P, schedules_mul[L], C);
                        ttrchildren = 0;
                        ttrchildren += std::get<1>(best[Lleft])[std::get<0>(best[Lleft])];
                        ttrchildren += std::get<1>(best[Lright])[std::get<0>(best[Lright])];
                    }

                    if (ttb >= basecase_keep_until * (ttr + ttrchildren))
                        basecase_eliminated = true;

                    bool rwin = ttb >= (ttr + ttrchildren);
                    best[L] = { rwin, {ttb, ttr + ttrchildren, ttr} };
                }

                time_b += std::get<1>(best[L])[0] * weight;
                time_r += std::get<1>(best[L])[1] * weight;
                time_r_self += std::get<1>(best[L])[2] * weight;
                time_m += std::get<1>(best[L])[idx] * weight;
                time_m_self += std::get<1>(best[L])[2*idx] * weight;
            }

            size_t L0 = std::get<0>(cws.front());
            size_t L1 = std::get<0>(cws.back());
            // ASSERT_ALWAYS(L0 <= L1);
            size_t L0r = plingen_tuning_cache::round_operand_size(L0);
            size_t L1r = plingen_tuning_cache::round_operand_size(L1);

            bool rec0 = std::get<0>(best[L0]);
            bool rec1 = std::get<0>(best[L1]);

            std::ostringstream os;
            os << "Time at depth " << i << ":";
            std::string oss = os.str();
            int pad = oss.size();
            const char * msg = oss.c_str();
            const char * msg2 = " ";
            const char * strbest = " [BEST]";
            if (basecase_eliminated || !recursion_makes_sense(L1))
                strbest="";
            if (!basecase_eliminated) {
                const char * isbest = (!rec0 && !rec1) ? strbest : "";
                printf("# %*s basecase(threshold>%zu): %.2f [%.1fd]%s\n", pad, msg, L1, time_b, time_b / 86400, isbest);
                msg = msg2;
            }
            if (L1r != L0r && recursion_makes_sense(L1)) {
                const char * isbest = (rec0 && !rec1) ? strbest : "";
                printf("# %*s mixed(threshold=%zu): %.2f [%.1fd] (self: %.2f [%.1fd])%s\n", pad, msg, L1, time_m, time_m / 86400, time_m_self, time_m_self / 86400, isbest);
                msg = msg2;
            }
            if (recursion_makes_sense(L0)) {
                const char * isbest = rec0 ? strbest : "";
                printf("# %*s recursive(threshold<=%zu): %.2f [%.1fd] (self: %.2f [%.1fd])%s\n", pad, msg, L0, time_r, time_r / 86400, time_r_self, time_r_self / 86400, isbest);
                msg = msg2;
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
        printf("# Automatically tuned lingen_mpi_threshold=%zu\n", upper_threshold);
        size_t size_com0;
        double tt_com0;
        std::tie(size_com0, tt_com0) = mpi_threshold_comm_and_time();
        char buf[20];
        printf("# Communication time at lingen_mpi_threshold (%s): %.2f [%.1fd]\n", size_disp(size_com0, buf), tt_com0, tt_com0/86400);
        double time_best = std::get<1>(best[L])[std::get<0>(best[L])];
        time_best += tt_com0;
        printf("# Expected total time: %.2f [%.1fd]\n", time_best, time_best / 86400);
    }
};

void plingen_tune_full(dims * d, cxx_param_list & pl)/*{{{*/
{
    plingen_tuner tuner(d, pl);

    tuner.tune();

#if 0
    size_t L = N/n + N/m;

    /* with basecase_keep_until == 0, then we never measure basecase */
    bool basecase_eliminated = !basecase_keep_until;

    size_t suggest_threshold = 0;

    double tt_com0 = tuner.compute_and_report_mpi_threshold_comm();

    double tt_total = tt_com0;
    double tt_total_nocomm = 0;
    size_t peakpeak=0;


    /* At depth i=log2(L) we run with input length L/2^i < 1, which means
     * sometimes a matrix, sometimes not. For sure, it makes no sense to
     * go _below_ that and use recursive calls at this depth */

    int floor = log2(L) + 1;

    // std::vector<matrix_product_schedule_with_transforms> schedules_mp;
    // std::vector<matrix_product_schedule_with_transforms> schedules_mul;
    // std::vector<std::pair<size_t, size_t>> aux_mems;

    for(int i = floor ; i>=0 ; i--) {

        printf("########## Measuring time at depth %d ##########\n", i);

        char buf[20];

        /* During recursive calls proper, we need:
         *  - storage for our input and our output, that need often not
         *    live simultaneously. For simplicity we'll count both, but
         *    this overlap is clearly a gross upper bound.
         *    our input is 1/r^2*m*(m+n)*\ell  (with \ell = L/2^i)
         *    our output is 1/r^2*(m+n)*(m+n)*(m/(m+n))*\ell , i.e. the
         *    same
         *  - storage *at all levels above the current one* (i.e. with
         *    lower depth) for the data that is still live and need to
         *    exist until after we return. This count is maximized in the
         *    rightmost branch, where pi_left at all levels must be kept.
         *    pi_left at depth 0 is
         *                  1/r^2*(m+n)*(m+n)*(m/(m+n))*L/2
         *    so the cumulated cost above is twice that minus our
         *    pi_right (i.e. L/2 above should instead be counted as
         *    L-L/2^i).
         */
        size_t c0 = iceildiv(m,r)*iceildiv(m+n,r)*mpz_size(tuner.p);
        size_t mem_our_input_and_output = 2 * c0 * iceildiv(L, 1 << i);
        size_t mem_upper_layers = c0 * (L - (L >> i));

        double tt_basecase_total = 0;
        double tt_mp_total = 0;
        double tt_mul_total = 0;
        size_t peak_mp = 0;
        size_t peak_mul = 0;
        aux_mems.push_back(std::make_pair(mem_our_input_and_output, mem_upper_layers));

        /* Arrange so that transforms + reserved storage does not exceed
         * our ram budget.
         */
        tuner.max_transform_ram = max_ram - mem_our_input_and_output - mem_upper_layers;
        if (!basecase_eliminated)
            tt_basecase_total = tuner.compute_and_report_basecase(i);

        if (i < floor) {
            try {
                tuner.compute_schedules_for_mp(i);
                tuner.compute_schedules_for_mul(i);
                tt_mp_total  = schedules_mp.back() .compute_derived_timings();
                peak_mp      = schedules_mp.back() .compute_derived_peak_ram();
                tt_mul_total = schedules_mul.back().compute_derived_timings();
                peak_mul     = schedules_mul.back().compute_derived_peak_ram();
            } catch (std::string & e) {
                fputs("\n", stderr);
                fprintf(stderr, "##################################\n");
                fputs(e.c_str(), stderr);
                fprintf(stderr, "##################################\n");
                fputs("\n", stderr);
                exit(EXIT_FAILURE);
            }
        } else {
            /* We have no other option than to use the basecase
             * algorithm.
             */
            tt_mp_total = DBL_MAX;
            tt_mul_total = DBL_MAX;
        }

        double tt_ontop = tt_mp_total + tt_mul_total;

        if (i == floor) {
            suggest_threshold = 0;
            tt_total_nocomm = tt_basecase_total;
        } else if (!basecase_eliminated) {
            printf("# cost of doing depth %d by a recursive call:\n"
                    "#     cost at depth %d: %.1f [%.1f days]\n"
                    "#     cost at depth %d (MP+MUL): %.1f [%.1f days]\n"
                    "#     total: %.1f [%.1f days]\n",
                    i,
                    i+1, tt_total_nocomm, tt_total_nocomm / 86400,
                    i, tt_ontop, tt_ontop / 86400,
                    (tt_total_nocomm + tt_ontop),
                    (tt_total_nocomm + tt_ontop) / 86400
                    );
            if (tt_total_nocomm != 0 && tt_total_nocomm + tt_ontop < tt_basecase_total) {
                if (suggest_threshold == 0) {
                    suggest_threshold = L>>i;
                    printf("# suggest lingen_threshold = %zu\n", L>>i);
                }
                if (basecase_keep_until * (tt_total_nocomm + tt_ontop) < tt_basecase_total) {
                    printf("# no longer counting basecase cost\n");
                    basecase_eliminated = true;
                }
                tt_total_nocomm += tt_ontop;
            } else {
                printf("# basecase still seems to win for depth %d, input length <= %zu\n", i, iceildiv(L, 1<<i));
                suggest_threshold = 0;
                tt_total_nocomm = tt_basecase_total;
            }
        } else {
            printf( "# added cost at depth %d (MP+MUL): %.1f [%.1f days]\n",
                    i, tt_ontop, tt_ontop / 86400);
            tt_total_nocomm += tt_ontop;
        }

        size_t tmem = std::max(peak_mp, peak_mul);
        size_t rmem = mem_our_input_and_output + mem_upper_layers + tmem;

        printf("# Total reserved memory per node at depth %d: %s\n",
                i,
                size_disp(rmem, buf));
        printf("#     operands at the levels above: %s\n",
                size_disp(mem_upper_layers, buf));
        printf("#     operands at this level: %s\n",
                size_disp(mem_our_input_and_output, buf));
        printf("#     transforms at this level: %s\n",
                size_disp(tmem, buf));

        tt_total = tt_total_nocomm + tt_com0;
        printf("# Cumulated time so far %.1f [%.1f days],"
                " peak RAM = %s, incl. %.1f days comm"
                " at lingen_mpi_threshold=%zu\n",
                tt_total, tt_total / 86400,
                size_disp(rmem, buf),
                tt_com0 / 86400,
                suggest_threshold);
        peakpeak = std::max(peakpeak, rmem);
    }
    {
        printf("#### Summary\n");
        unsigned int i;
        char buf[20];
        for(i = 0 ; i < schedules_mp.size() ; ++i) {
            /* should be const, really, but not until we get rid of
             * matrix_product_schedule_with_transforms::per_matrix */
            if ((L >> i) < suggest_threshold) break;
            auto & Smp(*(schedules_mp.end()-i-1));
            auto & Smul(*(schedules_mul.end()-i-1));
            auto const & MM(*(aux_mems.end()-i-1));
            printf("# depth=%u, input size %zu..%zu:\n", i, L >> i, iceildiv(L,1<<i));
            printf("#    operands %s\n", size_disp(MM.first, buf));
            printf("#    operands above %s\n", size_disp(MM.second, buf));
            printf("#    MP[shrink=(%u,%u) batch=%u] time=%.1f [%.1f days] RAM=%s\n",
                    Smp.shrink0,
                    Smp.shrink2,
                    Smp.batch,
                    Smp.compute_derived_timings(),
                    Smp.compute_derived_timings() / 86400,
                    size_disp(Smp.compute_derived_peak_ram(), buf));
            printf("#    MUL[shrink=(%u,%u) batch=%u] time=%.1f [%.1f days] RAM=%s\n",
                    Smul.shrink0,
                    Smul.shrink2,
                    Smul.batch,
                    Smul.compute_derived_timings(),
                    Smul.compute_derived_timings() / 86400,
                    size_disp(Smul.compute_derived_peak_ram(), buf));
        }
        {
            printf("# depth=%u: communication time=%.1f [%.1f days]\n",
                    i, tt_com0, tt_com0/86400);
            double ttb = tuner.compute_and_report_basecase(i);
            printf("# depth=%u: basecase time=%.1f [%.1f days]\n",
                    i, ttb, ttb/86400);
        }
        printf("# Cumulated time %.1f [%.1f days],"
                " peak RAM = %s, incl. %.1f days comm"
                " at lingen_mpi_threshold=%zu\n",
                tt_total, tt_total / 86400,
                size_disp(peakpeak, buf),
                tt_com0 / 86400,
                suggest_threshold);
    }

    printf("(%u,%u,%u,%.1f,%1.f)\n",m,n,r,tt_total,(double)peakpeak/1024./1024./1024.);
#endif
}/*}}}*/

void plingen_tuning(dims * d, MPI_Comm comm, cxx_param_list & pl)
{
    int rank;

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        plingen_tune_full(d, pl);
    }

    return;
}

void plingen_tuning_decl_usage(cxx_param_list & pl)
{
    plingen_tuner::declare_usage(pl);
}

void plingen_tuning_lookup_parameters(cxx_param_list & pl)
{
    plingen_tuner::lookup_parameters(pl);
}

