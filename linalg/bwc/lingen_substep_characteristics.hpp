#ifndef LINGEN_SUBSTEP_CHARACTERISTICS_HPP_
#define LINGEN_SUBSTEP_CHARACTERISTICS_HPP_

#include <iostream>

#include "cxx_mpz.hpp"
#include "fmt/format.h"
#include "lingen_platform.hpp"
#include "lingen_substep_schedule.hpp"
#include "lingen_tuning_cache.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_mul_substeps.hpp"
#include "lingen_call_companion.hpp"

/* This file is an offspring from lingen_tuning.cpp - some of it is not
 * really polished to the expected level of a usable interface. It
 * happens to be also used by the binary
 * time_matpoly_ft_parallel_${gfp_layer}
 */
template<typename OP>
struct lingen_substep_characteristics {
    typedef lingen_substep_characteristics<OP> ch_t;
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;
    typedef lingen_tuning_cache tc_t;

    abdst_field ab;
    cxx_mpz p;
    gmp_randstate_t & rstate;

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
    OP op;
    typedef typename lingen_tuning_cache_key<OP>::key_type cache_key_type;
    typedef typename lingen_tuning_cache_key<OP>::value_type cache_value_type;

    public:

    lingen_substep_characteristics(abdst_field ab, gmp_randstate_t & rstate, size_t input_length, unsigned int n0, unsigned int n1, unsigned int n2, size_t asize, size_t bsize, size_t csize) :/*{{{*/
        ab(ab),
        rstate(rstate),
        input_length(input_length),
        n0(n0), n1(n1), n2(n2),
        asize(asize), bsize(bsize), csize(csize)
    {
#ifdef SELECT_MPFQ_LAYER_u64k1
        mpz_set_ui(p, 2);
        op = OP(asize, bsize);
#else
        mpz_set(p, abfield_characteristic_srcptr(ab));
        op = OP(p, asize, bsize, n1);
#endif
        transform_ram = op.get_transform_ram();
    }/*}}}*/

    bool has_cached_time(tc_t const & C) const {/*{{{*/
        cache_key_type K { mpz_sizeinbase(p, 2), asize, bsize };
        return C.has(K);
    }/*}}}*/

    private:
    static inline int max_threads() {
#ifdef HAVE_OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

    /* tl;dr: array of doubles is input, weighted_double is output.
     *
     * A parallelizable timing is created from an array of
     * micro-measurements that tell how long (WCT) it takes to do k times
     * the same operation in parallel, for k less than some bound
     * (typically the number of threads, but it can be less if RAM
     * allows only so many operations in parallel).
     *
     * We _use_ a parallelizable timing in the form of the parent
     * weighted_double, essentially. This carries the info on the
     * _amortized time_ (t) it takes to do some number (n) of operations.
     * Originally t is initialized to what the array of concurrent
     * timings says of what happens for just one thread.
     *
     * The gist of it is that we have the .parallelize() method. It
     * changes the underlying weighted_double to reflect the amortized
     * timing, again, but when some number of operations are done in
     * parallel.
     *
     * This is in contrast with operator*(), which counts more operations
     * to do but not in parallel.
     */
    class parallelizable_timing : public weighted_double/*{{{*/
    {
        std::vector<double> concurrent_timings;
        bool counted_as_parallel = false;
        unsigned int max_parallel() const { return concurrent_timings.size() - 1; }
        public:
        struct cannot : public std::exception {
            const char * what() const noexcept override { return "cannot do this level of parallelism"; }
        };
        private:
        double get_concurrent_time(unsigned int k) {
            /* How long does it take (WCT) to do k operations in parallel ? */
            /* Ideally, this should be almost constant as long as k is
             * less than the number of threads, and then it should reach
             * linear behaviour (but we're restricting to small k
             * anyway).
             * In practice, "constant" isn't really constant because
             * there's some overhead. After all, all cores and
             * hyperthreads are competing for cache, for example. So we
             * allow a linear regression.
             */
            ASSERT_ALWAYS(k > 0);
            if (k > max_parallel()) throw cannot();
            unsigned int kl=1, kh=1;
            double tl=t, th=t;
            for(unsigned int i = 1 ; i <= k ; i++) {
                if (concurrent_timings[i] < 0) continue;
                kl = i;
                tl = concurrent_timings[i];
            }
            for(unsigned int i = max_threads() ; i >= k ; i--) {
                if (concurrent_timings[i] < 0) continue;
                kh = i;
                th = concurrent_timings[i];
            }
            if (kl == kh) return tl;
            ASSERT_ALWAYS(tl >= 0);
            ASSERT_ALWAYS(th >= 0);
            return tl + (th-tl)*(k-kl)/(kh-kl);
        }
        public:
        parallelizable_timing(std::vector<double> const & c) : concurrent_timings(c) {
            // ASSERT_ALWAYS(c.size() == (size_t) max_threads() + 1);
            this->n = 1;
            this->t = concurrent_timings[1];
        }
        /* This ctor is only for something that we will _never_
         * parallelize, so that we have no use for an array of
         * concurrent_timings.
         * In other cases, we insist on providing the
         * concurrent_timings array, via the other ctor.
         */
        parallelizable_timing(double d)
            : weighted_double { 1, d }
            , concurrent_timings { d } {}
        parallelizable_timing() : weighted_double { 0, 0 } {}

        parallelizable_timing& operator*=(unsigned int x) {
            n *= x;
            return *this;
        }
        parallelizable_timing operator*(unsigned int x) const { 
            parallelizable_timing X = *this;
            X *= x;
            return X;
        }
        parallelizable_timing& parallelize(unsigned int k, unsigned int T) {
            ASSERT_ALWAYS(!counted_as_parallel);
            counted_as_parallel = true;
            /* We suppose that on input, *this reprensents n operations
             * that sequentially take time t, and that this time
             * represents something single-threaded.  This function
             * replaces the inner operation by the fact of doing k times
             * this same thing, in parallel. The max level of expected
             * parallelism is given by T. The field t is replaced then
             * by the _average time_ that these operation take (that is,
             * 1/k times the cumulated time it takes to do k operations).
             * If perfect parallelism is achieved (all k operations done
             * up to T times in parallel, with linear scaling), this
             * effectively does
             *  this->n *= k
             *  this->t *= iceildiv(k, T) / (double) k;
             * However, even with T <= max_threads(), this
             * process is likely to give somewhat different results
             * because not everything achieves perfect parallelism at the
             * CPU level.
             *
             * Note that T should really be max_threads().
             */
            if (T > max_parallel()) throw cannot();
            unsigned int qk = k / T;
            unsigned int rk = k - qk * T;
            /* get_concurrent_time returns the wall-clock time needed to
             * do T dfts in parallel
             */
            double tt = qk * get_concurrent_time(T);
            if (rk) tt += get_concurrent_time(rk);
            n *= k;
            t = tt / k;
            return *this;
        }
        operator double() const { return n * t; }
        double operator+(parallelizable_timing const & b) const { return *this + (double) b; }
        double operator+(double y) const { return (double) *this + y; }
    };/*}}}*/


    /* We can measure timings in one of three manners here.
     *  - "cache-hot" do repeatedly the same thing, on identical data.
     *  - "cache-cold", quick: do the calculation only once (malloc,
     *  compute, free).
     *  - "cache-cold", thorough: do the calculation on matrices of
     *  different inputs.
     *
     * None is a true winner. "cache-hot" is easily optimistic, and both
     * "cache-cold" options are somewhat pessimistic.
     */

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

    public:
    typedef double (lingen_substep_characteristics<OP>::*raw_timer_member)(unsigned int, unsigned int, unsigned int) const;

    /* the three functions below return the WCT needed to do the
     * operation named (name) doing it with (nparallel) threads working
     * at the same time (at our level -- we don't preclude the use of
     * more omp threads at the level below). The computation is repeated
     * (n) times, and we work with allocated space that is (k) times
     * bigger than what is necessary to accomodate the (nparallel)
     * parallel calls.
     */

    struct microbench_dft { /*{{{*/
        pc_t P;
        ch_t const & U;
        static constexpr const char * name = "dft";
        microbench_dft(pc_t const& P, ch_t const & U) : P(P), U(U) {}
        unsigned int max_parallel() {
            size_t R = U.asize * mpz_size(U.p) * sizeof(mp_limb_t);
            R += U.get_transform_ram();
            unsigned int n = P.available_ram / R;
            if (n == 0) throw std::overflow_error("not enough RAM");
            return n;
        }
        double operator()(unsigned int nparallel) {
            matpoly a(U.ab, nparallel, 1, U.asize);
            a.fill_random(U.asize, U.rstate);
            matpoly_ft<typename OP::FFT> ta(a.m, a.n, U.op.fti);
            double tt = -wct_seconds();
            matpoly_ft<typename OP::FFT>::dft(ta, a);
            return tt + wct_seconds();
        }
    };/*}}}*/
    struct microbench_ift { /*{{{*/
        pc_t P;
        ch_t const & U;
        static constexpr const char * name = "ift";
        microbench_ift(pc_t const& P, ch_t const & U) : P(P), U(U) {}
        unsigned int max_parallel() {
            size_t R = U.csize * mpz_size(U.p) * sizeof(mp_limb_t);
            R += U.get_transform_ram();
            unsigned int n = P.available_ram / R;
            if (n == 0) throw std::overflow_error("not enough RAM");
            return n;
        }
        double operator()(unsigned int nparallel) {
            matpoly c(U.ab, nparallel, 1, U.csize);
            matpoly_ft<typename OP::FFT> tc(c.m, c.n, U.op.fti);
            /* This is important, since otherwise the inverse transform won't
             * work */
            c.set_size(U.csize);
            tc.zero(); /* would be .fill_random(rstate) if we had it */
            double tt = -wct_seconds();
            matpoly_ft<typename OP::FFT>::ift(c, tc);
            return tt + wct_seconds();
        }
    };/*}}}*/
    struct microbench_conv { /*{{{*/
        pc_t P;
        ch_t const & U;
        static constexpr const char * name = "conv";
        microbench_conv(pc_t const& P, ch_t const & U) : P(P), U(U) {}
        unsigned int max_parallel() {
            /* for k = nparallel, we need 2k+1 transforms in ram */
            unsigned int n = (P.available_ram / U.get_transform_ram() - 1) / 2;
            if (n == 0) throw std::overflow_error("not enough RAM");
            return n;
        }
        double operator()(unsigned int nparallel) {
            matpoly c(U.ab, nparallel, 1, U.csize);
            matpoly_ft<typename OP::FFT> tc(c.m, c.n, U.op.fti);
            /* This is important, since otherwise the inverse transform won't
             * work */
            c.set_size(U.csize);
            tc.zero(); /* would be .fill_random(rstate) if we had it */
            double tt = -wct_seconds();
            matpoly_ft<typename OP::FFT>::ift(c, tc);
            return tt + wct_seconds();
        }
    };/*}}}*/

    private:
    template<typename T>
    void complete_tvec(std::ostream& os, std::vector<double> & tvec, T & F, unsigned int kl, unsigned int kh) const/*{{{*/
    {
        /* The ideal expected behaviour is that this should be constant.
         * But it's not, since flint itself uses openmp. We're doing
         * better, since we're coarser grain. But still, as timing goes,
         * this means that we should be prepared to have a roughly linear
         * curve.
         */
        if (kh <= kl + 1) return;
        double tl = tvec[kl];
        double th = tvec[kh];
        if (th - tl <= 0.1 * tl) return;
        unsigned int k = (kl + kh) / 2;
        double tk = F(k);
        os << fmt::format(" {}:{:.2g}", k, tk);
        os.flush();
        tvec[k] = tk;
        double linfit_tk = tl + (th-tl)*(k-kl)/(kh-kl);
        if ((tk - linfit_tk) <= 0.1*linfit_tk) return;
        complete_tvec(os, tvec, F, kl, k);
        complete_tvec(os, tvec, F, k, kh);
    }/*}}}*/

    public:
    template<typename T>
    std::vector<double> fill_tvec(T F) const {/*{{{*/
        /* can't make my mind as to whether I should just output that to
         * the terminal, or stow it in a string first...
         */
        std::ostream& os(std::cout);

        os << "# wct for " << F.name << " by nthreads:";

        unsigned int TMAX = max_threads();
        if (F.max_parallel() < TMAX) {
            TMAX = F.max_parallel();
            os << fmt::format(" [capped to {}]", TMAX);
        }

        std::vector<double> tvec(TMAX + 1, -1);

        unsigned int kl = 1;
        double tl = F(kl);
        tvec[kl] = tl;
        os << fmt::format(" {}:{:.2g}", kl, tl);
        os.flush();


        unsigned int kh = TMAX;
        if (kh > 1) {
            double th = F(kh);
            tvec[kh] = th;
            os << fmt::format(" {}:{:.2g}", kh, th);
            os.flush();
        }

        complete_tvec(os, tvec, F, kl, kh);

        os << std::endl;
        return tvec;
    }/*}}}*/

    private:
    std::array<parallelizable_timing, 4> get_ft_times(pc_t const& P, tc_t & C) const {/*{{{*/
        cache_key_type K { mpz_sizeinbase(p, 2), asize, bsize };
        
        unsigned int TMAX = max_threads();

        if (C.has(K)) {
            std::array<parallelizable_timing, 4> res;
            for(unsigned int i = 0 ; i < 4 ; i++) {
                std::vector<double> tvec(TMAX + 1, -1);
                for(auto const & x : C[K][i])
                    if (x.first < tvec.size())
                        tvec[x.first] = x.second;
                res[i] = parallelizable_timing(tvec);
            }
            return res;
        }

        std::vector<double> tvec_dft = fill_tvec(microbench_dft(P, *this));
        std::vector<double> tvec_ift = fill_tvec(microbench_ift(P, *this));
        std::vector<double> tvec_conv= fill_tvec(microbench_conv(P, *this));

        cache_value_type cached_res;
        for(unsigned int i = 1 ; i <= TMAX ; i++) {
            if (tvec_dft[i] >= 0) cached_res[0].push_back( { i, tvec_dft[i] } );
            if (tvec_dft[i] >= 0) cached_res[1].push_back( { i, tvec_dft[i] } );
            if (tvec_conv[i] >= 0) cached_res[2].push_back( { i, tvec_conv[i] } );
            if (tvec_ift[i] >= 0) cached_res[3].push_back( { i, tvec_ift[i] } );
        }

        C[K] = cached_res;

        parallelizable_timing tt_dft0(tvec_dft);
        parallelizable_timing tt_dft2(tvec_dft);
        parallelizable_timing tt_ift(tvec_ift);
        parallelizable_timing tt_conv(tvec_conv);

        return {{ tt_dft0, tt_dft2, tt_conv, tt_ift }};
    }/*}}}*/

    public:
    size_t get_transform_ram() const { return transform_ram; }
    /* {{{ std::tuple<size_t, size_t, size_t> get_operand_ram() const {
        return {
            n0 * n1 * asize * mpz_size(p) * sizeof(mp_limb_t),
            n1 * n2 * bsize * mpz_size(p) * sizeof(mp_limb_t),
            n0 * n2 * csize * mpz_size(p) * sizeof(mp_limb_t) };
    }
    *//*}}}*/
    int mesh_inner_size(pc_t const & P) const {/*{{{*/
        return P.r;
    }/*}}}*/
    subdivision mpi_split0(pc_t const & P) const {/*{{{*/
        return subdivision(n0, mesh_inner_size(P));
    }/*}}}*/
    subdivision mpi_split1(pc_t const & P) const {/*{{{*/
        return subdivision(n1, mesh_inner_size(P));
    }/*}}}*/
    subdivision mpi_split2(pc_t const & P) const {/*{{{*/
        return subdivision(n2, mesh_inner_size(P));
    }/*}}}*/
    subdivision shrink_split0(pc_t const & P, unsigned int shrink0) const {/*{{{*/
        return subdivision(mpi_split0(P).block_size_upper_bound(), shrink0);
    }/*}}}*/
    subdivision shrink_split2(pc_t const & P, unsigned int shrink2) const {/*{{{*/
        return subdivision(mpi_split2(P).block_size_upper_bound(), shrink2);
    }/*}}}*/
    subdivision shrink_split0(pc_t const & P, sc_t const & S) const {/*{{{*/
        return shrink_split0(P, S.shrink0);
    }/*}}}*/
    subdivision shrink_split2(pc_t const & P, sc_t const & S) const {/*{{{*/
        return shrink_split2(P, S.shrink2);
    }/*}}}*/
    size_t get_peak_ram(pc_t const & P, sc_t const & S) const { /* {{{ */
        unsigned int nrs0 = shrink_split0(P, S).block_size_upper_bound();
        unsigned int nrs2 = shrink_split2(P, S).block_size_upper_bound();
        return get_transform_ram() * (S.batch[1] * P.r * (S.batch[0] + S.batch[2]) + nrs0*nrs2);
    }/*}}}*/

    private:

    std::array<parallelizable_timing, 6> get_call_time_backend(pc_t const & P, sc_t const & S, tc_t & C) const { /* {{{ */
        auto ft = get_ft_times(P, C);
        /* These are just base values, we'll multiply them later on */
        parallelizable_timing T_dft0 = ft[0];
        parallelizable_timing T_dft2 = ft[1];
        parallelizable_timing T_conv = ft[2];
        parallelizable_timing T_ift = ft[3];

        /* Each of the r*r nodes has local matrices size (at most)
         * nr0*nr1, nr1*nr2, and nr0*nr2. For the actual computations, we
         * care more about the shrunk submatrices, though */
        unsigned int nr1 = mpi_split1(P).block_size_upper_bound();

        /* The shrink parameters will divide the size of the local
         * matrices we consider by numbers shrink0 and shrink2. This
         * increases the time, and decreases the memory footprint */
        unsigned int nrs0 = shrink_split0(P, S).block_size_upper_bound();
        unsigned int nrs2 = shrink_split2(P, S).block_size_upper_bound();
        unsigned int ns0 = nrs0 * mesh_inner_size(P);
        unsigned int ns2 = nrs2 * mesh_inner_size(P);

        /* The 1-threaded timing will be obtained by setting T=batch=1.  */

        /* Locally, we'll add together r matrices of size nrs0*nrs2,
         * which will be computed as products of vectors of size
         * batch[0]*batch[1] and batch[1]*batch[2]. (with wither
         * batch[0]==nrs0 or batch[2]==nrs2). The operation will be
         * repeated nrs0/batch[0] * nr1/batch[1] * nrs2/batch[2] times
         * (not all multipliers apply to everyone).
         */

        T_dft0.parallelize(S.batch[0] * S.batch[1], P.T);
        T_dft0 *= iceildiv(nrs0, S.batch[0]);
        T_dft0 *= iceildiv(nr1, S.batch[1]);

        T_dft2.parallelize(S.batch[1] * S.batch[2], P.T);
        T_dft2 *= iceildiv(nrs2, S.batch[2]);
        T_dft2 *= iceildiv(nr1, S.batch[1]);

        T_conv.parallelize(S.batch[0] * S.batch[2], P.T);
        T_conv *= mesh_inner_size(P) * S.batch[1];
        T_conv *= iceildiv(nr1, S.batch[1]) * S.batch[1];
        T_conv *= iceildiv(nrs0, S.batch[0]);
        T_conv *= iceildiv(nrs2, S.batch[2]);

        T_ift.parallelize(S.batch[0] * S.batch[2], P.T);
        T_ift *= iceildiv(nrs0, S.batch[0]);
        T_ift *= iceildiv(nrs2, S.batch[2]);

        /* Now the communication time _per call_ */
        /* TODO: maybe include some model of latency... */
        parallelizable_timing T = get_transform_ram() / P.mpi_xput;
        T *= S.batch[1] * iceildiv(nr1, S.batch[1]);
        parallelizable_timing T_comm0 = T * ns0;
        parallelizable_timing T_comm2 = T * ns2;

        T_dft0 *= S.shrink0 * S.shrink2;
        T_dft2 *= S.shrink0 * S.shrink2;
        T_conv *= S.shrink0 * S.shrink2;
        T_ift  *= S.shrink0 * S.shrink2;
        T_comm0 *= S.shrink0 * S.shrink2;
        T_comm2 *= S.shrink0 * S.shrink2;

        /* This is for _one_ call only (one node in the tree).
         * We must apply multiplier coefficients (typically 1<<depth) */
        return {{ T_dft0, T_dft2, T_conv, T_ift, T_comm0, T_comm2 }};
    }/*}}}*/

    public:


    lingen_call_companion::mul_or_mp_times get_companion(pc_t const & P, sc_t const & S, tc_t & C) const { /* {{{ */
        lingen_call_companion::mul_or_mp_times D;
        D.S = S;
        auto A = get_call_time_backend(P, S, C); 
        double tt = 0;
        for(auto const & a : A) tt = tt + a;
        D.tt = { 1, tt };
        D.t_dft_A = A[0];
        D.t_dft_B = A[1];
        D.t_conv = A[2];
        D.t_ift_C = A[3];
        D.t_dft_A_comm = A[4];
        D.t_dft_B_comm = A[5];
        D.per_transform_ram = get_transform_ram();
        D.ram = get_peak_ram(P, S);
        D.asize = asize;
        D.bsize = bsize;
        D.csize = csize;
        return D;
    }/*}}}*/
    double get_call_time(pc_t const & P, sc_t const & S, tc_t & C) const {/*{{{*/
        auto A = get_call_time_backend(P, S, C);
        double tt = 0;
        for(auto const & a : A)
            tt = tt + a;
        return tt;
    }/*}}}*/

    private:
    struct schedule_sorter {/*{{{*/
        lingen_substep_characteristics const & U;
        lingen_platform const & P;
        lingen_tuning_cache & C;
        schedule_sorter(lingen_substep_characteristics const & U, lingen_platform const & P, lingen_tuning_cache & C) :
            U(U), P(P), C(C)
        {}
        bool operator()(lingen_substep_schedule const & a, lingen_substep_schedule const & b) const {/*{{{*/
            double ta = U.get_call_time(P, a, C);
            double tb = U.get_call_time(P, b, C);
            if (ta != tb) return ta < tb;
#if 0
            /* Doesn't make much sense, since timing comparisons will
             * almost certainly always return unequal. But well, who
             * knows...
             */
            size_t za = U.get_peak_ram(P, a);
            size_t zb = U.get_peak_ram(P, b);
            return za < zb;
#endif
            /* Favor more batching instead */
            unsigned int Ba = a.batch[1] * std::min(a.batch[0], a.batch[2]);
            unsigned int Bb = b.batch[1] * std::min(b.batch[0], b.batch[2]);
            if (Ba != Bb)
                return Ba > Bb;
            return false;
        }/*}}}*/
    };/*}}}*/

    public:
    void sort_schedules(std::vector<lingen_substep_schedule>& schedules, lingen_platform const & P, lingen_tuning_cache & C) const {/*{{{*/
        sort(schedules.begin(), schedules.end(), schedule_sorter(*this, P, C));
    }/*}}}*/

    double get_and_report_call_time(pc_t const & P, sc_t const & S, tc_t & C) const { /* {{{ */
        const char * step = OP::name;
        bool cached = has_cached_time(C);
        char buf[20];

        printf("# %s(@%zu) [shrink=(%u,%u) batch=(%u,%u,%u)] %s, ",
                step,
                input_length,
                S.shrink0,
                S.shrink2,
                S.batch[0],
                S.batch[1],
                S.batch[2],
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
};

#endif	/* LINGEN_SUBSTEP_CHARACTERISTICS_HPP_ */
