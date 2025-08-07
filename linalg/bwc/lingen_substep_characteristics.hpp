#ifndef CADO_LINGEN_SUBSTEP_CHARACTERISTICS_HPP
#define CADO_LINGEN_SUBSTEP_CHARACTERISTICS_HPP

#include <cstddef>

#include <algorithm>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>
#include <stdexcept>

#include "fmt/format.h"
#include <gmp.h>

#include "gmp_aux.h"
#include "cxx_mpz.hpp"
#include "subdivision.hpp"
#include "lingen_platform.hpp"
#include "lingen_substep_schedule.hpp"
#include "lingen_tuning_cache.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_matpoly_select.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_mul_substeps_base.hpp"
#include "lingen_mul_substeps.hpp"
#include "lingen_call_companion.hpp"
#include "macros.h"
#include "timing.h"


/* This file is an offspring from lingen_tuning.cpp - some of it is not
 * really polished to the expected level of a usable interface. It
 * happens to be also used by the binary
 * time_matpoly_ft_parallel_${gfp_layer}
 */

template<typename OP, bool = !std::is_same<typename OP::FFT, void>::value>
struct microbench_dft;
template<typename OP, bool = !std::is_same<typename OP::FFT, void>::value>
struct microbench_ift;
template<typename OP, bool = !std::is_same<typename OP::FFT, void>::value>
struct microbench_conv;

template<bool is_binary>
struct lingen_substep_characteristics {
    typedef lingen_substep_characteristics ch_t;
    typedef lingen_platform pc_t;
    typedef lingen_substep_schedule sc_t;
    typedef lingen_tuning_cache tc_t;

    typename matpoly<is_binary>::arith_hard * ab;
    cxx_mpz p;
    cxx_gmp_randstate & rstate;

    /* length of the input (E) for the call under consideration ; this is
     * not the input length for the overall algorithm !
     *
     * We only use it for printing.
     * */
    size_t input_length;

    op_mul_or_mp_base::op_type_t op_type;

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

    public:

    bool is_valid_schedule(lingen_platform, unsigned int mesh, lingen_substep_schedule const & S) const
    {
        const unsigned int shrink0 = S.shrink0;
        const unsigned int shrink2 = S.shrink2;
        const unsigned int r = mesh;
        const subdivision mpi_split0(n0, r);
        const subdivision mpi_split1(n1, r);
        const subdivision mpi_split2(n2, r);
        const subdivision shrink_split0(mpi_split0.block_size_upper_bound(), shrink0);
        const subdivision shrink_split2(mpi_split2.block_size_upper_bound(), shrink2);
        const unsigned int nrs0(shrink_split0.block_size_upper_bound());
        const unsigned int nrs2(shrink_split2.block_size_upper_bound());
        const unsigned int nr1 = mpi_split1.block_size_upper_bound();
        const unsigned int b0 = S.batch[0];
        const unsigned int b1 = S.batch[1];
        const unsigned int b2 = S.batch[2];
        const subdivision loop1 = subdivision::by_block_size(nr1, b1);
        for(unsigned int round = 0 ; round < shrink0 * shrink2 ; round++) {
            const unsigned round0 = round % shrink0;
            const unsigned round2 = round / shrink0;
            unsigned int i0, i1;
            unsigned int j0, j1;
            std::tie(i0, i1) = shrink_split0.nth_block(round0);
            std::tie(j0, j1) = shrink_split2.nth_block(round2);
            ASSERT_ALWAYS((i1 - i0) <= nrs0);
            ASSERT_ALWAYS((j1 - j0) <= nrs2);
            const subdivision loop0 = subdivision::by_block_size(i1 - i0, b0);
            const subdivision loop2 = subdivision::by_block_size(j1 - j0, b2);

            if (loop0.nblocks() != 1 && loop2.nblocks() != 1)
                return false;
            const bool process_blocks_row_major = b0 == nrs0;

            for(unsigned int iloop1 = 0 ; iloop1 < loop1.nblocks() ; iloop1++) {
                if (process_blocks_row_major) {
                    if (loop0.nblocks() != 1) return false;
                    if (nrs0 != b0) return false;
                } else {
                    if (loop2.nblocks() != 1) return false;
                    if (nrs2 != b2) return false;
                }
            }
        }
        return true;
    }

    lingen_substep_characteristics(typename matpoly<is_binary>::arith_hard * ab, cxx_gmp_randstate & rstate, size_t input_length, op_mul_or_mp_base::op_type_t op_type, unsigned int n0, unsigned int n1, unsigned int n2, size_t asize, size_t bsize, size_t csize) :/*{{{*/
        ab(ab),
        rstate(rstate),
        input_length(input_length),
        op_type(op_type),
        n0(n0), n1(n1), n2(n2),
        asize(asize), bsize(bsize), csize(csize)
    {
        mpz_set(p, ab->characteristic());
        // op = OP(asize, bsize);
        // op = OP(p, asize, bsize, n1);
        // fft_alloc_sizes = op.fti.get_alloc_sizes();
    }/*}}}*/

    static int max_threads() {
#ifdef HAVE_OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

    private:
    bool has_cached_time(tc_t const & C, lingen_substep_schedule::fft_type_t fft_type) const {/*{{{*/
        lingen_tuning_cache::mul_or_mp_key K { op_type, fft_type, mpz_sizeinbase(p, 2), asize, bsize };
        return C.has(K);
    }/*}}}*/

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
            ASSERT_ALWAYS(k <= max_parallel());
            unsigned int kl=1, kh=1;
            double tl=t, th=t;
            for(unsigned int i = 1 ; i <= k ; i++) {
                if (concurrent_timings[i] < 0) continue;
                kl = i;
                tl = concurrent_timings[i];
            }
            for(unsigned int i = max_parallel() ; i >= k ; i--) {
                if (concurrent_timings[i] < 0) continue;
                kh = i;
                th = concurrent_timings[i];
            }
            if (kl == kh) return tl;
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
            , concurrent_timings { 0, d }
        {
        }

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
            T = std::min(T, max_parallel());
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

// #define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

    public:
    // typedef double (lingen_substep_characteristics::*raw_timer_member)(unsigned int, unsigned int, unsigned int) const;

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
        const double tl = tvec[kl];
        const double th = tvec[kh];
        if (th - tl <= 0.1 * tl) return;
        unsigned int k = (kl + kh) / 2;
        double tk = F(k);
        os << fmt::format(" {}:{:.2g}", k, tk);
        os.flush();
        tvec[k] = tk;
        const double linfit_tk = tl + (th-tl)*(k-kl)/(kh-kl);
        if ((tk - linfit_tk) <= 0.1*linfit_tk) return;
        complete_tvec(os, tvec, F, kl, k);
        complete_tvec(os, tvec, F, k, kh);
    }/*}}}*/

    public:
    template<typename T>
    std::vector<double> fill_tvec(std::ostream& os, T const & F) const {/*{{{*/
        /* can't make my mind as to whether I should just output that to
         * the terminal, or stow it in a string first...
         */

        constexpr const char * Fname = T::name;

        unsigned int TMAX = max_threads();

        std::vector<double> tvec(TMAX + 1, -1);

        if (!Fname) {
            tvec[1] = tvec[TMAX] = 0;
            return tvec;
        }

        os << fmt::format("# {}{};{} (@{}) wct for {} by nthreads:",
                F.mesh > 1 ? "MPI-" : "",
                F.op.op_name(),
                lingen_substep_schedule::fft_name(encode_fft_type<typename T::OP::FFT>),
                input_length,
                Fname
                );

        if (F.max_parallel() < TMAX) {
            TMAX = F.max_parallel();
            os << fmt::format(" [capped to {}]", TMAX);
        }

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

        os << "\n";
        return tvec;
    }/*}}}*/

    private:
    template<typename T>
    parallelizable_timing get_ft_time_from_cache_or_recompute(std::ostream& os, T const & F, tc_t::mul_or_mp_value::value_type & store) const
    {
        if (!store.empty()) {
            /* what do we have in cache, to start with ? We need to know how
             * many threads max were considered at the time when the tuning
             * was made.
             */
            unsigned int th_cache = 0;
            for(auto const & x : store)
                th_cache = std::max(th_cache, x.first);

            /* How many threads do we support max on the current platform ?
             * This is not necessarily equal.
             */
            if (F.max_parallel() == th_cache) {
                std::vector<double> tvec(th_cache + 1, -1);
                for(auto const & x : store) {
                    tvec[x.first] = x.second;
                }
                return parallelizable_timing(tvec);
            }
            os << fmt::format("# ignoring cached entry,"
                    " computed for up to {} threads"
                    " (here {} max)\n", th_cache, F.max_parallel());
            store.clear();
        }

        std::vector<double> tvec = fill_tvec(os, F);
        for(unsigned int i = 1 ; i < tvec.size() ; i++) {
            if (tvec[i] >= 0) store.emplace_back(i, tvec[i]);
        }
        return parallelizable_timing(tvec);
    }

    template<typename OP>
    std::array<parallelizable_timing, 4> get_ft_times(std::ostream& os, OP const & op, pc_t const& P, unsigned int mesh, lingen_substep_schedule const & S, tc_t & C) const {/*{{{*/
        lingen_tuning_cache::mul_or_mp_key K { op_type, S.fft_type, mpz_sizeinbase(p, 2), asize, bsize };
        typedef lingen_tuning_cache::mul_or_mp_value cache_value_type;
        
        std::array<parallelizable_timing, 4> res;
        cache_value_type cached_res;

        if (C.has(K))
            cached_res = C[K];

        {
            microbench_dft<OP> F(op, P, mesh, *this);
            res[0] = get_ft_time_from_cache_or_recompute(os, F, cached_res[0]);
            res[1] = res[0];
            cached_res[1] = cached_res[0];
        }

        {
            microbench_ift<OP> F(op, P, mesh, *this);
            res[2] = get_ft_time_from_cache_or_recompute(os, F, cached_res[2]);
        }

        {
            microbench_conv<OP> F(op, P, mesh, *this);
            res[3] = get_ft_time_from_cache_or_recompute(os, F, cached_res[3]);
        }

        /* A priori this is either all-fresh or all-old. But nothing in
         * the code checks that the situation is indeed the same for the
         * four items.
         */
        C[K] = cached_res;

        return res;
    }/*}}}*/

    public:
    /* {{{ std::tuple<size_t, size_t, size_t> get_operand_ram() const {
        return {
            n0 * n1 * asize * mpz_size(p) * sizeof(mp_limb_t),
            n1 * n2 * bsize * mpz_size(p) * sizeof(mp_limb_t),
            n0 * n2 * csize * mpz_size(p) * sizeof(mp_limb_t) };
    }
    *//*}}}*/
#if 0
    int mesh_inner_size(pc_t const & P) const {/*{{{*/
        return P.r;
    }/*}}}*/
#endif
    subdivision mpi_split0(unsigned int mesh) const {/*{{{*/
        return { n0, mesh };
    }/*}}}*/
    subdivision mpi_split1(unsigned int mesh) const {/*{{{*/
        return { n1, mesh };
    }/*}}}*/
    subdivision mpi_split2(unsigned int mesh) const {/*{{{*/
        return { n2, mesh };
    }/*}}}*/
    subdivision shrink_split0(unsigned int mesh, unsigned int shrink0) const {/*{{{*/
        unsigned int nr0 = mpi_split0(mesh).block_size_upper_bound();
        return { nr0, shrink0 };
    }/*}}}*/
    subdivision shrink_split2(unsigned int mesh, unsigned int shrink2) const {/*{{{*/
        unsigned int nr2 = mpi_split2(mesh).block_size_upper_bound();
        return { nr2, shrink2 };
    }/*}}}*/
    subdivision shrink_split0(unsigned int mesh, sc_t const & S) const {/*{{{*/
        return shrink_split0(mesh, S.shrink0);
    }/*}}}*/
    subdivision shrink_split2(unsigned int mesh, sc_t const & S) const {/*{{{*/
        return shrink_split2(mesh, S.shrink2);
    }/*}}}*/
    std::array<std::array<unsigned int, 3>, 2> get_peak_ram_multipliers(unsigned int mesh, sc_t const & S) const { /* {{{ */
        unsigned int nrs0 = shrink_split0(mesh, S).block_size_upper_bound();
        unsigned int nrs2 = shrink_split2(mesh, S).block_size_upper_bound();

        const unsigned int b0 = S.batch[0];
        const unsigned int b1 = S.batch[1];
        const unsigned int b2 = S.batch[2];
        unsigned int mul0 = 0;
        mul0 += b1 * mesh * b0;
        mul0 += b1 * mesh * b2;
        mul0 += nrs0*nrs2;
        unsigned int mul1 = std::max(b0 * std::max(b1, b2), nrs0 * nrs2);
        unsigned int mul12 = b0 * b2;

        /* *IF* the pragma omp parallel statements go with an appropriate
         * num_threads() clause, then yes, it makes sense to do this
         * min(). Otherwise, the multipliers will be max_threads() in all
         * cases.
         */
        mul1 = std::min(mul1, (unsigned int) max_threads());
        mul12 = std::min(mul12, (unsigned int) max_threads());

        /*
        size_t live_transforms = fft_alloc_sizes[0] * mul0;
        size_t maxtemp_ft = std::min(nrs0 * nrs2,
                (unsigned int) max_threads()) * fft_alloc_sizes[1];
        size_t maxtemp_addmul = std::min(S.batch[0] * S.batch[2],
                (unsigned int) max_threads())
            * (fft_alloc_sizes[1] + fft_alloc_sizes[2]);
        return live_transforms + std::max(maxtemp_ft, maxtemp_addmul);
        */
        return {{
            // if (mul1 * op.get_alloc_sizes()[1] < mul12 *
            // (op.get_alloc_sizes()[1] + op.get_alloc_sizes()[2])) then
            // the largest amount of ram is with the first of the two
            // options below.
            {{ mul0, mul12, mul12 }},
            {{ mul0, mul1, 0 }}
        }};
    }
    /*}}}*/
    size_t get_peak_ram(op_mul_or_mp_base const & op, unsigned int mesh, sc_t const & S) const { /* {{{ */
        auto multipliers = get_peak_ram_multipliers(mesh, S);
        size_t rpeak = 0;
        for(auto const & M : multipliers) {
            size_t r = 0;
            for(unsigned int i = 0 ; i < 3 ; i++)
                r += M[i] * op.get_alloc_sizes()[i];
            rpeak = std::max(rpeak, r);
        }
        return rpeak;
    }
    /*}}}*/
    size_t get_peak_ram(unsigned int mesh, sc_t const & S) const { /* {{{ */
        std::shared_ptr<op_mul_or_mp_base> op = instantiate(S.fft_type);
        return get_peak_ram(*op, mesh, S);
    }
    /*}}}*/

    bool fft_type_valid(lingen_substep_schedule::fft_type_t fft_type) const;


    private:

    struct call_time_digest {
        parallelizable_timing T_dft0;
        parallelizable_timing T_dft2;
        parallelizable_timing T_conv;
        parallelizable_timing T_ift;
        parallelizable_timing T_comm0;
        parallelizable_timing T_comm2;
        std::array<size_t, 3> fft_alloc_sizes;
        double total() const {
            double tt = 0;
            tt += T_dft0;
            tt += T_dft2;
            tt += T_conv;
            tt += T_ift;
            tt += T_comm0;
            tt += T_comm2;
            return tt;
        }
    };

    template<typename OP>
    call_time_digest get_call_time_backend(
            std::ostream& os,
            OP const & op,
            pc_t const & P,
            unsigned int mesh,
            sc_t const & S,
            tc_t & C,
            bool do_timings) const
    { /* {{{ */
        ASSERT_FOR_STATIC_ANALYZER(mesh != 0);

        /* XXX Any change here must also be reflected in the mp_or_mul
         * structure in lingen_matpoly_bigmatpoly_ft_common.hpp
         */
        parallelizable_timing T_dft0 { -1 };
        parallelizable_timing T_dft2 { -1 };
        parallelizable_timing T_ift  { -1 };
        parallelizable_timing T_conv { -1 };

        if (do_timings) {
            auto ft = get_ft_times(os, op, P, mesh, S, C);
            /* These are just base values, we'll multiply them later on */
            T_dft0 = ft[0];
            T_dft2 = ft[1];
            T_ift = ft[2];
            T_conv = ft[3];
        }

        /* Each of the r*r nodes has local matrices size (at most)
         * nr0*nr1, nr1*nr2, and nr0*nr2. For the actual computations, we
         * care more about the shrunk submatrices, though */
        unsigned int nr1 = mpi_split1(mesh).block_size_upper_bound();

        /* The shrink parameters will divide the size of the local
         * matrices we consider by numbers shrink0 and shrink2. This
         * increases the time, and decreases the memory footprint */
        unsigned int nrs0 = shrink_split0(mesh, S).block_size_upper_bound();
        unsigned int nrs2 = shrink_split2(mesh, S).block_size_upper_bound();

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
        T_conv *= mesh;
        T_conv *= iceildiv(nr1, S.batch[1]) * S.batch[1];
        T_conv *= iceildiv(nrs0, S.batch[0]);
        T_conv *= iceildiv(nrs2, S.batch[2]);

        T_ift.parallelize(nrs0 * nrs2, P.T);

        /* Now the communication time _per call_ */
        /* TODO: maybe include some model of latency... */
        /* FIXME: This is wrong with the non-caching approach */
        parallelizable_timing T = op.get_alloc_sizes()[0] / P.mpi_xput;
        T *= S.batch[1] * iceildiv(nr1, S.batch[1]);
        /* allgather over n nodes */
        T *= mesh - 1;
        parallelizable_timing T_comm0 = T;
        parallelizable_timing T_comm2 = T;
        T_comm0 *= S.batch[0] * iceildiv(nrs0, S.batch[0]);
        T_comm2 *= S.batch[2] * iceildiv(nrs2, S.batch[2]);

        call_time_digest ret;

        ret.T_dft0 = T_dft0 *= S.shrink0 * S.shrink2;
        ret.T_dft2 = T_dft2 *= S.shrink0 * S.shrink2;
        ret.T_conv = T_conv *= S.shrink0 * S.shrink2;
        ret.T_ift  = T_ift  *= S.shrink0 * S.shrink2;
        ret.T_comm0 = T_comm0 *= S.shrink0 * S.shrink2;
        ret.T_comm2 = T_comm2 *= S.shrink0 * S.shrink2;
        ret.fft_alloc_sizes = op.get_alloc_sizes();

        /* This is for _one_ call only (one node in the tree).
         * We must apply multiplier coefficients (typically 1<<depth) */
        return ret;
    }/*}}}*/

    template<typename fft_type> friend struct matpoly_checker_ft;

    std::shared_ptr<op_mul_or_mp_base> instantiate(lingen_substep_schedule::fft_type_t fft_type) const;

    call_time_digest get_call_time_backend(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, bool do_timings) const;

    public:

    lingen_call_companion::mul_or_mp_times get_companion(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, size_t reserved_ram, bool do_timings) const { /* {{{ */
        lingen_call_companion::mul_or_mp_times D { op_type };
        D.S = S;
        auto A = get_call_time_backend(os, P, mesh, S, C, do_timings); 
        double tt = A.total();
        D.tt = { 1, tt };
        D.t_dft_A = A.T_dft0;
        D.t_dft_B = A.T_dft2;
        D.t_conv = A.T_conv;
        D.t_ift_C = A.T_ift;
        D.t_dft_A_comm = A.T_comm0;
        D.t_dft_B_comm = A.T_comm2;
        D.fft_alloc_sizes = A.fft_alloc_sizes;
        D.peak_ram_multipliers = get_peak_ram_multipliers(mesh, S);
        D.asize = asize;
        D.bsize = bsize;
        D.csize = csize;
        D.reserved_ram = reserved_ram;
        return D;
    }/*}}}*/
    double get_call_time(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, bool do_timings) const {/*{{{*/
        return get_call_time_backend(os, P, mesh, S, C, do_timings).total();
    }/*}}}*/

    private:
    struct sortop {
        bool operator()(std::tuple<double, unsigned int, lingen_substep_schedule> const & a, std::tuple<double, unsigned int, lingen_substep_schedule> const & b) const {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) < std::get<0>(b);
            else
                return std::get<1>(a) > std::get<1>(b);
        }
    };

    public:
    void sort_schedules(
            std::ostream& os,
            std::vector<lingen_substep_schedule>& schedules,
            lingen_platform const & P,
            unsigned int mesh,
            lingen_tuning_cache & C,
            bool do_timings) const
    {/*{{{*/
        std::vector<std::tuple<double, unsigned int, lingen_substep_schedule>> precomp;
        for(auto & S : schedules) {
            /* This should ensure that all timings are obtained from cache */
            /* This may print timing info to the output stream */
            double t = do_timings ? get_call_time(os, P, mesh, S, C, do_timings) : 0;
            unsigned int B = S.batch[1] * std::min(S.batch[0], S.batch[2]);
            precomp.emplace_back(t, B, S);
        }
        sort(precomp.begin(), precomp.end(), sortop());
        schedules.clear();
        for(auto & P : precomp) {
            schedules.emplace_back(std::get<2>(P));
        }
    }/*}}}*/

    double get_and_report_call_time(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, bool do_timings) const { /* {{{ */
        bool cached = has_cached_time(C, S.fft_type);
        std::shared_ptr<op_mul_or_mp_base> op = instantiate(S.fft_type);
        os << fmt::format("# {}{};{} (@{}) [shrink=({},{}) batch=({},{},{})] {}, ",
                mesh > 1 ? "MPI-" : "",
                op_mul_or_mp_base::op_name(op_type),
                S.fft_name(),
                input_length,
                S.shrink0,
                S.shrink2,
                S.batch[0],
                S.batch[1],
                S.batch[2],
                size_disp(get_peak_ram(*op, mesh, S)));
        os << std::flush;
        double tt = get_call_time(os, P, mesh, S, C, do_timings);
        if (do_timings) {
            os << fmt::format("{:.2f}{}\n",
                    tt,
                    cached ? " [from cache]" : "");
        } else {
            os << "not timed\n";
            tt = -1;
        }
        return tt;
    }/*}}}*/
    void report_op_winner(std::ostream& os, unsigned int mesh, sc_t const & S) const
    {
        std::shared_ptr<op_mul_or_mp_base> op = instantiate(S.fft_type);
        os << fmt::format("# {}{};{} wins : {}\n",
                mesh > 1 ? "MPI-" : "",
                op_mul_or_mp_base::op_name(op_type),
                S.fft_name(),
                op->explain());
    }

    void report_size_stats_human(std::ostream& os, op_mul_or_mp_base const & op) const {/*{{{*/
        os << fmt::format("# {} (per op): {}+{}+{}, transforms 3*{}\n",
                op_mul_or_mp_base::op_name(op_type),
                size_disp(asize*mpz_size(p)*sizeof(mp_limb_t)),
                size_disp(bsize*mpz_size(p)*sizeof(mp_limb_t)),
                size_disp(csize*mpz_size(p)*sizeof(mp_limb_t)),
                size_disp(op.get_alloc_sizes()[0]));
        os << fmt::format("# {} (total for {}*{} * {}*{}): {}, transforms {}\n",
                op_mul_or_mp_base::op_name(op_type),
                n0,n1,n1,n2,
                size_disp((n0*n1*asize+n1*n2*bsize+n0*n2*csize)*mpz_size(p)*sizeof(mp_limb_t)),
                size_disp((n0*n1+n1*n2+n0*n2)*op.get_alloc_sizes()[0]));
    }/*}}}*/
};

#ifdef LINGEN_BINARY
template<>
inline std::shared_ptr<op_mul_or_mp_base>
lingen_substep_characteristics<true>::instantiate(lingen_substep_schedule::fft_type_t fft_type) const
{
    /* must create the op type from the lingen_substep_schedule
     * (a.k.a. sc_t) type */ 
    switch(op_type) {
        case op_mul_or_mp_base::OP_MP:
            switch(fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return std::make_shared<op_mp<true, void>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_CANTOR:
                    return std::make_shared<op_mp<true, gf2x_cantor_fft_info>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_TERNARY:
                    return std::make_shared<op_mp<true, gf2x_ternary_fft_info>>(p, asize, bsize, n1);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        case op_mul_or_mp_base::OP_MUL:
            switch(fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return std::make_shared<op_mul<true, void>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_CANTOR:
                    return std::make_shared<op_mul<true, gf2x_cantor_fft_info>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_TERNARY:
                    return std::make_shared<op_mul<true, gf2x_ternary_fft_info>>(p, asize, bsize, n1);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        default:
            throw std::runtime_error("invalid data (op)");
    }
}
template<>
inline auto lingen_substep_characteristics<true>::get_call_time_backend(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, bool do_timings) const -> call_time_digest { /* {{{ */
    switch(op_type) {
        case op_mul_or_mp_base::OP_MP:
            switch(S.fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return get_call_time_backend(os, op_mp<true, void>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_CANTOR:
                    return get_call_time_backend(os, op_mp<true, gf2x_cantor_fft_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_TERNARY:
                    return get_call_time_backend(os, op_mp<true, gf2x_ternary_fft_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        case op_mul_or_mp_base::OP_MUL:
            switch(S.fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return get_call_time_backend(os, op_mul<true, void>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_CANTOR:
                    return get_call_time_backend(os, op_mul<true, gf2x_cantor_fft_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_TERNARY:
                    return get_call_time_backend(os, op_mul<true, gf2x_ternary_fft_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        default:
            throw std::runtime_error("invalid data (op)");
    }
}/*}}}*/

template<>
inline bool lingen_substep_characteristics<true>::fft_type_valid(lingen_substep_schedule::fft_type_t fft_type) const {
    try {
        std::shared_ptr<op_mul_or_mp_base> op = instantiate(fft_type);
        if (op_type == op_mul_or_mp_base::OP_MUL) {
            auto x = dynamic_cast<op_mul<true, gf2x_ternary_fft_info> const *>(op.get());
            if (x)
                return x->fti.K != 0;
        }
        if (op_type == op_mul_or_mp_base::OP_MP) {
            auto x = dynamic_cast<op_mp<true, gf2x_ternary_fft_info> const *>(op.get());
            if (x)
                return x->fti.K != 0;
        }
    } catch (std::exception const &) {
        return false;
    }
    return true;
}

#endif

#ifndef LINGEN_BINARY
template<>
inline std::shared_ptr<op_mul_or_mp_base>
lingen_substep_characteristics<false>::instantiate(lingen_substep_schedule::fft_type_t fft_type) const
{
    /* must create the op type from the lingen_substep_schedule
     * (a.k.a. sc_t) type */ 
    switch(op_type) {
        case op_mul_or_mp_base::OP_MP:
            switch(fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return std::make_shared<op_mp<false, void>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_FLINT:
                    return std::make_shared<op_mp<false, fft_transform_info>>(p, asize, bsize, n1);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        case op_mul_or_mp_base::OP_MUL:
            switch(fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return std::make_shared<op_mul<false, void>>(p, asize, bsize, n1);
                case lingen_substep_schedule::FFT_FLINT:
                    return std::make_shared<op_mul<false, fft_transform_info>>(p, asize, bsize, n1);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        default:
            throw std::runtime_error("invalid data (op)");
    }
}

template<>
inline auto lingen_substep_characteristics<false>::get_call_time_backend(std::ostream& os, pc_t const & P, unsigned int mesh, sc_t const & S, tc_t & C, bool do_timings) const -> call_time_digest { /* {{{ */
    switch(op_type) {
        case op_mul_or_mp_base::OP_MP:
            switch(S.fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return get_call_time_backend(os, op_mp<false, void>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_FLINT:
                    return get_call_time_backend(os, op_mp<false, fft_transform_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        case op_mul_or_mp_base::OP_MUL:
            switch(S.fft_type) {
                case lingen_substep_schedule::FFT_NONE:
                    return get_call_time_backend(os, op_mul<false, void>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                case lingen_substep_schedule::FFT_FLINT:
                    return get_call_time_backend(os, op_mul<false, fft_transform_info>(p, asize, bsize, n1), P, mesh, S, C, do_timings);
                default:
                    throw std::runtime_error("invalid data (fft_type)");
            }
        default:
            throw std::runtime_error("invalid data (op)");
    }
}/*}}}*/

template<>
inline bool lingen_substep_characteristics<false>::fft_type_valid(lingen_substep_schedule::fft_type_t fft_type) const {
    try {
        std::shared_ptr<op_mul_or_mp_base> op = instantiate(fft_type);
        /* eh ? Don't we need to have some testing here ? XXX */
    } catch (std::exception const &) {
        return false;
    }
    return true;
}

#endif


/* the three structs below return the WCT needed to do the
 * operation named (name) doing it with (nparallel) threads working
 * at the same time (at our level -- we don't preclude the use of
 * more omp threads at the level below). The computation is repeated
 * (n) times, and we work with allocated space that is (k) times
 * bigger than what is necessary to accomodate the (nparallel)
 * parallel calls.
 */

template<typename OP_T>
struct microbench_dft<OP_T, true> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = "dft";
    microbench_dft(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh)
        , U(U)
    {}
    unsigned int max_parallel() const {
        size_t R = 0;
        /* storage for one input coeff and one output transform */
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        R += iceildiv(U.asize, simd) * mpz_size(U.p) * sizeof(mp_limb_t);
        R += op.get_alloc_sizes()[0];
        /* plus the temp memory for the dft operation */
        R += op.get_alloc_sizes()[1];
        size_t n;
        if (P.available_ram) {
            n = P.available_ram / R;
            /* TODO: we probably want to ask P about max_threads. */
            n = std::min(n, (size_t) U.max_threads());
            if (n == 0) throw std::runtime_error("not enough RAM");
        } else {
            n = U.max_threads();
        }
        return n;
    }
    double operator()(unsigned int nparallel) const
    {
        matpoly<is_binary> a(U.ab, nparallel, 1, U.asize);
        a.zero_pad(U.asize);
        a.fill_random(0, U.asize, U.rstate);
        matpoly_ft<typename OP::FFT> ta(a.m, a.n, op.fti);
        double tt = -wct_seconds();
        matpoly_ft<typename OP::FFT>::dft(ta, a);
        return tt + wct_seconds();
    }
};/*}}}*/
template<typename OP_T>
struct microbench_ift<OP_T, true> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = "ift";
    microbench_ift(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh), U(U)
    {}
    unsigned int max_parallel() const {
        size_t R = 0;
        /* storage for one input transform and one output coeff */
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        R += iceildiv(U.csize, simd) * mpz_size(U.p) * sizeof(mp_limb_t);
        R += op.get_alloc_sizes()[0];
        /* plus the temp memory for the ift operation */
        R += op.get_alloc_sizes()[1];
        size_t n;
        if (P.available_ram) {
            n = P.available_ram / R;
            /* TODO: we probably want to ask P about max_threads. */
            n = std::min(n, (size_t) U.max_threads());
            if (n == 0) throw std::runtime_error("not enough RAM");
        } else {
            n = U.max_threads();
        }
        return n;
    }
    double operator()(unsigned int nparallel) const {
        matpoly<is_binary> c(U.ab, nparallel, 1, U.csize);
        matpoly_ft<typename OP::FFT> tc(c.m, c.n, op.fti);
        /* This is important, since otherwise the inverse transform won't
         * work */
        c.set_size(U.csize);
        tc.zero(); /* would be .fill_random(rstate) if we had it */
        const double tt = -wct_seconds();
        matpoly_ft<typename OP::FFT>::ift(c, tc);
        return tt + wct_seconds();
    }
};/*}}}*/
template<typename OP_T>
struct microbench_conv<OP_T, true> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = "conv";
    microbench_conv(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh)
        , U(U)
    {}
    unsigned int max_parallel() const {
        /* for k = nparallel, we need 2k+1 transforms in ram */
        size_t R = 0;
        /* storage for two input transform */
        R += 2 * op.get_alloc_sizes()[0];
        /* plus the temp memory for the addmul operation */
        R += op.get_alloc_sizes()[1];
        R += op.get_alloc_sizes()[2];
        /* take into account the +1 transform */
        size_t n;
        if (P.available_ram) {
            n = (P.available_ram - op.get_alloc_sizes()[0]) / R;
            /* TODO: we probably want to ask P about max_threads. */
            n = std::min(n, (size_t) U.max_threads());
            if (n == 0 || P.available_ram < op.get_alloc_sizes()[0])
                throw std::runtime_error("not enough RAM");
        } else {
            n = U.max_threads();
        }
        return n;
    }
    double operator()(unsigned int nparallel) const {
        matpoly_ft<typename OP::FFT> tc( nparallel, 1, op.fti);
        matpoly_ft<typename OP::FFT> tc0(nparallel, 1, op.fti);
        matpoly_ft<typename OP::FFT> tc1(1, 1, op.fti);
        tc0.zero(); /* would be .fill_random(rstate) if we had it */
        tc1.zero(); /* would be .fill_random(rstate) if we had it */
        double tt = -wct_seconds();
        matpoly_ft<typename OP::FFT>::addcompose(tc, tc0, tc1);
        return tt + wct_seconds();
    }
};/*}}}*/

template<typename OP_T>
struct microbench_dft<OP_T, false> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = nullptr;
    microbench_dft(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh)
        , U(U)
    {}
    unsigned int max_parallel() const { return U.max_threads(); }
    double operator()(unsigned int) const { return 0; }
};/*}}}*/
template<typename OP_T>
struct microbench_ift<OP_T, false> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = nullptr;
    microbench_ift(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh)
        , U(U)
    {}
    unsigned int max_parallel() const { return U.max_threads(); }
    double operator()(unsigned int) const { return 0; }
};/*}}}*/
template<typename OP_T>
struct microbench_conv<OP_T, false> { /*{{{*/
    static constexpr bool is_binary = OP_T::is_binary;
    typedef lingen_substep_characteristics<is_binary> ch_t;
    typedef lingen_platform pc_t;
    typedef OP_T OP;
    OP op;
    pc_t P;
    unsigned int mesh;
    ch_t const & U;
    static constexpr const char * name = "product";
    microbench_conv(OP op, pc_t const& P, unsigned int mesh, ch_t const & U)
        : op(std::move(op))
        , P(P)
        , mesh(mesh)
        , U(U)
    {}
    unsigned int max_parallel() const { return U.max_threads(); }
    double operator()(unsigned int nparallel) const {
        constexpr const unsigned int splitwidth = is_binary ? 64 : 1;
        matpoly<is_binary> a(U.ab, nparallel, splitwidth, U.asize);
        matpoly<is_binary> b(U.ab, splitwidth, 1, U.asize);
        a.zero_pad(U.asize);
        a.fill_random(0, U.asize, U.rstate);
        b.zero_pad(U.bsize);
        b.fill_random(0, U.bsize, U.rstate);
        matpoly<is_binary> c(U.ab, nparallel, 1, U.csize);
        c.zero_pad(U.csize);
        const double tt = -wct_seconds();
        OP::addcompose(c, a, b);
        return (tt + wct_seconds()) / splitwidth;
    }
};/*}}}*/

#endif	/* LINGEN_SUBSTEP_CHARACTERISTICS_HPP_ */
