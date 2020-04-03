#ifndef LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_
#define LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_

/* The code here is inserted by the compiler inside lingen_matpoly_ft.cpp
 * and lingen_bigmatpoly_ft.cpp, as a template function. The class
 * OP_CTX<T> is instantiated as either OP_CTX<matpoly> or
 * OP_CTX<bigmatpoly>, with specificities in each file (mostly dealing
 * with the extra MPI definitions used in bigmatpoly, as well as the
 * back-end of the mpi operations.
 */

/* For consistency and readability, we add a few includes, although
 * they're evidently already included in the two main consumers of this
 * file. */

/* Use matpoly as imported from lingen_matpoly_ft.hpp */
#include "lingen_matpoly_ft.hpp"
#include "lingen_mul_substeps.hpp"
#include "subdivision.hpp"
#include "logline.h"

template<typename T>
struct OP_CTX_base {
    T & c;
    T const & a;
    T const & b;
    OP_CTX_base(T & c, T const & a, T const & b) : c(c), a(a), b(b) {}
};

template<typename T, typename fft_type> struct OP_CTX;

template<typename OP_CTX_T, typename OP_T> struct mp_or_mul { 
    typedef typename OP_CTX_T::T T;
    typedef typename OP_CTX_T::FFT fft_type;
    OP_CTX_T & CTX;
    OP_T & OP;
    T & c;
    T const & a;
    T const & b;
    const lingen_call_companion::mul_or_mp_times * M;
    unsigned int shrink0;
    unsigned int shrink2;
    subdivision mpi_split0;
    subdivision mpi_split1;
    subdivision mpi_split2;
    subdivision shrink_split0;
    subdivision shrink_split2;
    unsigned int nrs0, nrs2;
    unsigned int b0, b1, b2;
    /* Declare ta, tb, tc early on so that we don't malloc/free n times.  */
    matpoly_ft<fft_type> ta, tb, tc;
    mp_or_mul(OP_CTX_T & CTX, OP_T & OP,
            const lingen_call_companion::mul_or_mp_times * M)
        : CTX(CTX)
        , OP(OP)
        , c(CTX.c)
        , a(CTX.a)
        , b(CTX.b)
        , M(M)
        , shrink0(M ? M->S.shrink0 : 1)
        , shrink2(M ? M->S.shrink2 : 1)
        , mpi_split0(a.m, CTX.mesh_inner_size())
        , mpi_split1(a.n, CTX.mesh_inner_size())
        , mpi_split2(b.n, CTX.mesh_inner_size())
        , shrink_split0(mpi_split0.block_size_upper_bound(), shrink0)
        , shrink_split2(mpi_split2.block_size_upper_bound(), shrink2)
        /* first, upper bounds on output block dimensions */
        , nrs0(shrink_split0.block_size_upper_bound())
        , nrs2(shrink_split2.block_size_upper_bound())
        /* By default we use full batching, which costs some memory ! */
        , b0(M ? M->S.batch[0] : nrs0)
        , b1(M ? M->S.batch[1] : a.n)
        , b2(M ? M->S.batch[2] : nrs2)
        , ta(b0, CTX.mesh_inner_size() * b1, OP.fti)
        , tb(CTX.mesh_inner_size() * b1, b2, OP.fti)
        , tc(nrs0, nrs2, OP.fti)
    {
        CTX.mesh_checks();

        if (!M) return;

        /* The smallstep "MP" or "MUL" has already been planned since the
         * first entry in the recursive function in lingen.cpp -- here
         * we're only beginning the planning of the small steps. This
         * used to be done together with the planning of MP and MUL
         * themselves, but we prefer to do that closer to the code.
         *
         * XXX Note that any changes to the control flow below, in
         * operator()() and the other functions, must be reflected in
         * lingen_substep_characteristics::get_call_time_backend
         */
        begin_plan_smallstep_microsteps(M->step_name());
        plan_smallstep("dft_A", M->t_dft_A);
        if (CTX.uses_mpi) {
            begin_plan_smallstep("dft_A_comm", M->t_dft_A_comm);
            plan_smallstep("export", M->t_dft_A_comm);
            plan_smallstep("comm", M->t_dft_A_comm);
            plan_smallstep("import", M->t_dft_A_comm);
            end_plan_smallstep();
        }
        plan_smallstep("dft_B", M->t_dft_B);
        if (CTX.uses_mpi) {
            begin_plan_smallstep("dft_B_comm", M->t_dft_B_comm);
            plan_smallstep("export", M->t_dft_B_comm);
            plan_smallstep("comm", M->t_dft_B_comm);
            plan_smallstep("import", M->t_dft_B_comm);
            end_plan_smallstep();
        }
        plan_smallstep("addmul", M->t_conv);
        plan_smallstep("ift_C", M->t_ift_C);
        end_plan_smallstep();
    }
    template<typename... Args>
    inline void begin_smallstep(Args&& ...args) {
        if (M) CTX.stats.begin_smallstep(args...);
    }
    template<typename... Args>
    inline void skip_smallstep(Args&& ...args) {
        if (M) CTX.stats.skip_smallstep(args...);
    }
    template<typename... Args>
    inline void end_smallstep(Args&& ...args) {
        if (M) CTX.stats.end_smallstep(args...);
    }
    template<typename... Args>
    inline void plan_smallstep(Args&& ...args) {
        if (M) CTX.stats.plan_smallstep(args...);
    }
    template<typename... Args>
    inline void begin_plan_smallstep_microsteps(Args&& ...args) {
        if (M) CTX.stats.begin_plan_smallstep_microsteps(args...);
    }
    template<typename... Args>
    inline void begin_plan_smallstep(Args&& ...args) {
        if (M) CTX.stats.begin_plan_smallstep(args...);
    }
    template<typename... Args>
    inline void end_plan_smallstep(Args&& ...args) {
        if (M) CTX.stats.end_plan_smallstep(args...);
    }
    inline bool local_smallsteps_done(bool compulsory = false) {
        return M ? CTX.stats.local_smallsteps_done(compulsory) : true;
    }


    /* loop0 and loop2 depend on the exact output block. At most we're
     * iterating on, respectively, ceil(ceil(n0/r)/shrink0), and
     * ceil(ceil(n2/r)/shrink2). However this depends on the exact block
     * indices within shrink_split0 and shrink_split2.
     */
    subdivision loop0;
    subdivision loop1;
    subdivision loop2;

    /* Based on shrink0 and shrink2, we process output in blocks of size
     * nrs0*nrs2, for a total number of rounds equal to shrink0*shrink2.
     * No allocation carries over from one of these rounds to the next.
     *
     * For each of these rounds, we allocate space for the input and
     * output transforms, but the order of processing varies. Entries in
     * the output block are processed as blocks of size b0*b2. For each
     * these, b1 pairs of input transforms are computed  at the same time
     * (that is, b1*(b0+b2)), and it total we collect
     * CTX.mesh_inner_size() as many from the MPI peers. The order in
     * which the blocks of size b0*b2 are processed to cover the range of
     * size of nrs0 * nrs2 output values controls the number of input
     * transforms we have to compute in total, namely:
     *  ceil(nrs0/b0)*(b0+b2)*b1 if we process blocks row-major.
     *  ceil(nrs2/b2)*(b0+b2)*b1 if we process blocks col-major.
     *
     * Note though that this is a rather exaggerated notion: we *must*
     * have either b0==nrs0 or n2==nrs2, since otherwise we'd be better
     * off having different shrink0 / shrink2 values. (say if b0<nrs0 and
     * b2==1, then shrink0 should be increased so that b0 == nrs0).
     * Therefore, if we heed this adjustment, the number of transforms
     * that are computed to process a block of size b0*b2 in the output
     * is always (b0+b2)*b1, with the processing order being determined
     * by which of b0==nrs0 or b2==nrs2 holds.
     */
    void dft_A_for_block(unsigned int i0, unsigned int iloop0, unsigned int iloop1)/*{{{*/
    {
        unsigned int aj = CTX.a_jrank();
        unsigned int ak0mpi, ak1mpi;
        std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(aj);

        unsigned int kk0,kk1;
        std::tie(kk0, kk1) = loop1.nth_block(iloop1);
        ASSERT_ALWAYS((kk1 - kk0) <= b1);
        /* XXX ak0 and co, and esp. ak1-ak0, are *NOT* identical across
         * mpi jobs. All that we have is ak1-ak0 <= b1 and bk1-bk0 <= b1.
         * In the non-mpi case, ak0==bk0 and ak1==bk1 */
        unsigned int ak0 = ak0mpi + kk0;
        unsigned int ak1 = std::min(ak1mpi, ak0 + b1);

        unsigned int ii0, ii1;
        std::tie(ii0, ii1) = loop0.nth_block(iloop0);

        submatrix_range Ra  (i0 + ii0, ak0-ak0mpi, ii1 - ii0, ak1-ak0);
        submatrix_range Rat (     0,   aj * b1,    ii1 - ii0, ak1-ak0);
        submatrix_range Ratx(     0,   aj * b1,    ii1 - ii0, b1);

        begin_smallstep("dft_A", b0 * b1);
        ta.zero();  // for safety because of rounding.
        matpoly_ft<fft_type>::dft(ta.view(Rat), CTX.a_local().view(Ra));
        end_smallstep();

        if (!CTX.uses_mpi) return;

        const unsigned int r = CTX.mesh_inner_size();
        submatrix_range Ratxx(0,   0,          ii1 - ii0, r*b1);

        // allgather ta among r nodes.
        begin_smallstep("dft_A_comm", (r-1) * b0 * b1);
        begin_smallstep("export", (r-1) * b0 * b1);
        ta.view(Ratx).to_export();
        end_smallstep();
        begin_smallstep("comm", (r-1) * b0 * b1);
        /* The data isn't contiguous, so we have to do
         * several allgather operations.  */
        for(unsigned int i = 0 ; i < ii1 - ii0 ; i++) {
            CTX.a_allgather(ta.part(i, 0), b1);
        }
        end_smallstep();
        begin_smallstep("import", (r-1) * b0 * b1);
        /* we imported data from everywhere, we can now put
         * it back to local format */
        ta.view(Ratxx).to_import();
        end_smallstep();
        end_smallstep();
    }/*}}}*/

    void dft_B_for_block(unsigned int j0, unsigned int iloop1, unsigned int iloop2)/*{{{*/
    {
        unsigned int bi = CTX.b_irank();
        unsigned int bk0mpi, bk1mpi;
        std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(bi);

        unsigned int kk0,kk1;
        std::tie(kk0, kk1) = loop1.nth_block(iloop1);
        ASSERT_ALWAYS((kk1 - kk0) <= b1);
        /* XXX ak0 and co, and esp. ak1-ak0, are *NOT* identical across
         * mpi jobs. All that we have is ak1-ak0 <= b1 and bk1-bk0 <= b1.
         * In the non-mpi case, ak0==bk0 and ak1==bk1 */
        unsigned int bk0 = bk0mpi + kk0;
        unsigned int bk1 = std::min(bk1mpi, bk0 + b1);

        unsigned int jj0, jj1;
        std::tie(jj0, jj1) = loop2.nth_block(iloop2);


        submatrix_range Rb   (bk0-bk0mpi, j0 + jj0, bk1-bk0, jj1 - jj0);
        submatrix_range Rbt  (bi * b1,      0,      bk1-bk0, jj1 - jj0);
        submatrix_range Rbtx (bi * b1,      0,      b1,      jj1 - jj0);

        begin_smallstep("dft_B", b1 * b2);
        tb.zero();
        matpoly_ft<fft_type>::dft(tb.view(Rbt), CTX.b_local().view(Rb));
        end_smallstep();

        if (!CTX.uses_mpi) return;

        const unsigned int r = CTX.mesh_inner_size();
        submatrix_range Rbtxx(0,            0, r*b1,    jj1 - jj0);

        // allgather tb among r nodes
        begin_smallstep("dft_B_comm", (r-1) * b1 * b2);
        begin_smallstep("export", (r-1) * b1 * b2);
        tb.view(Rbtx).to_export();
        end_smallstep();
        begin_smallstep("comm", (r-1) * b1 * b2);
        CTX.b_allgather(tb.data, b1 * b2);
        end_smallstep();
        begin_smallstep("import", (r-1) * b1 * b2);
        tb.view(Rbtxx).to_import();
        end_smallstep();
        end_smallstep();
    }/*}}}*/

    void addmul_for_block(unsigned int iloop0, unsigned int iloop2)/*{{{*/
    {
        const unsigned int r = CTX.mesh_inner_size();
        unsigned int ii0, ii1;
        unsigned int jj0, jj1;
        std::tie(ii0, ii1) = loop0.nth_block(iloop0);
        std::tie(jj0, jj1) = loop2.nth_block(iloop2);

        begin_smallstep("addmul", b0 * b1 * b2 * r);

        // rounding might surprise us.
        submatrix_range Ratxx(0,   0,   ii1 - ii0, r*b1);
        submatrix_range Rbtxx(0,   0,   r*b1,      jj1 - jj0);
        submatrix_range Rct  (ii0, jj0, ii1 - ii0, jj1 - jj0);

        /* we can't do much progress tracking, since we delve into openmp
         * almost immediately anyway */
        matpoly_ft<fft_type>::addcompose(tc.view(Rct), ta.view(Ratxx), tb.view(Rbtxx));

        end_smallstep();
    }/*}}}*/

    void operator()() {
        if (M) begin_smallstep(M->step_name());

        CTX.alloc_c_if_needed(OP.csize);

        const unsigned int r = CTX.mesh_inner_size();


        /* The order in which we do the transforms is not really our main
         * concern at this point. If sharing makes sense, then probably
         * shrink0 and shrink2 do not. So they're serving opposite purposes.
         */
        const unsigned int nr1 = mpi_split1.block_size_upper_bound();
        loop1 = subdivision::by_block_size(nr1, b1);

        // unsigned int imax = mpi_split0.nth_block_size(CTX.a_irank());
        // unsigned int jmax = mpi_split2.nth_block_size(CTX.b_jrank());

        /* We must both count the number of transforms we really have to deal
         * with, as well as the theoretical upper bound, because the latter
         * was used to count the theoretical time.
         *
         * For the upper bounds, ak1-ak0 and bk1-bk0 are always replaced by
         * batch.
         */

        /* In the non-mpi case, mpi_split1 has one chunk only, 
         * rank==0, so that ak0mpi=bk0mpi=0 and ak1mpi=bk1mpi=a.n
         */
        unsigned int aj = CTX.a_jrank();
        unsigned int bi = CTX.b_jrank();
        unsigned int ak0mpi, ak1mpi;
        unsigned int bk0mpi, bk1mpi;
        std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(aj);
        std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(bi);

        for(unsigned int round = 0 ; round < shrink0 * shrink2 ; round++) {
            unsigned round0 = round % shrink0;
            unsigned round2 = round / shrink0;

            /* Prepare the processing of the small blocks of size b0*b2
             */
            unsigned int i0, i1;
            unsigned int j0, j1;
            std::tie(i0, i1) = shrink_split0.nth_block(round0);
            std::tie(j0, j1) = shrink_split2.nth_block(round2);
            // TODO: we only have true data for [xi0, xi1[, not [i0, i1[
            // unsigned int xi0 = std::min(i0, imax);
            // unsigned int xi1 = std::min(i1, imax);
            // TODO: we only have true data for [xj0, xj1[, not [j0, j1[
            // unsigned int xj0 = std::min(j0, jmax);
            // unsigned int xj1 = std::min(j1, jmax);
            /* Note that i0, i1, j0, j1 are equal for all mpi jobs. Therefore
             * the three loops below are synchronous for all jobs. */
            ASSERT_ALWAYS((i1 - i0) <= nrs0);
            ASSERT_ALWAYS((j1 - j0) <= nrs2);
            loop0 = subdivision::by_block_size(i1 - i0, b0);
            loop2 = subdivision::by_block_size(j1 - j0, b2);

            ASSERT_ALWAYS(loop0.nblocks() == 1 || loop2.nblocks() == 1);
            bool process_blocks_row_major = b0 == nrs0;

            /* Now do a subblock */
            tc.zero();

            for(unsigned int iloop1 = 0 ; iloop1 < loop1.nblocks() ; iloop1++) {
                if (process_blocks_row_major) {
                    ASSERT_ALWAYS(loop0.nblocks() == 1);
                    ASSERT_ALWAYS(nrs0 == b0);
                    unsigned int iloop0 = 0;
                    logline_printf(1, "dft_A (%u*%u)\n", b0, b1);
                    dft_A_for_block(i0, iloop0, iloop1);
                    for(unsigned int iloop2 = 0 ; iloop2 < loop2.nblocks() ; iloop2++) {
                        logline_printf(1, "dft_B (%u*%u)\n", b1, b2);
                        dft_B_for_block(j0, iloop1, iloop2);

                        logline_printf(1, "addmul\n");
                        addmul_for_block(iloop0, iloop2);
                    }
                    /* adjust counts */
                    unsigned int e2 = (iceildiv(nrs2, b2) - loop2.nblocks());
                    unsigned int x2 = e2 * b2;
                    skip_smallstep("dft_B", b1 * x2);
                    if (CTX.uses_mpi) {
                        begin_smallstep("dft_B_comm", (r-1) * b1 * x2);
                        skip_smallstep("export", (r-1) * b1 * x2);
                        skip_smallstep("comm", (r-1) * b1 * x2);
                        skip_smallstep("import", (r-1) * b1 * x2);
                        end_smallstep();
                    }
                    skip_smallstep("addmul", b0 * b1 * x2 * r);
                } else {
                    ASSERT_ALWAYS(loop2.nblocks() == 1);
                    ASSERT_ALWAYS(nrs2 == b2);
                    unsigned int iloop2 = 0;
                    logline_printf(1, "dft_B (%u*%u)\n", b1, b2);
                    dft_B_for_block(j0, iloop1, iloop2);
                    for(unsigned int iloop0 = 0 ; iloop0 < loop0.nblocks() ; iloop0++) {
                        logline_printf(1, "dft_A (%u*%u)\n", b0, b1);
                        dft_A_for_block(i0, iloop0, iloop1);

                        logline_printf(1, "addmul\n");
                        addmul_for_block(iloop0, iloop2);
                    }
                    /* adjust counts */
                    unsigned int e0 = (iceildiv(nrs0, b0) - loop0.nblocks());
                    unsigned int x0 = e0 * b0;
                    skip_smallstep("dft_A", x0 * b1);
                    if (CTX.uses_mpi) {
                        begin_smallstep("dft_A_comm", (r-1) * x0 * b1);
                        skip_smallstep("export", (r-1) * x0 * b1);
                        skip_smallstep("comm", (r-1) * x0 * b1);
                        skip_smallstep("import", (r-1) * x0 * b1);
                        end_smallstep();
                    }
                    skip_smallstep("addmul", x0 * b1 * b2 * r);
                }
            }

            CTX.c.set_size(OP.csize);
            /* tc.get_size()() <= nrs0 * nrs2 */
            begin_smallstep("ift_C", nrs0 * nrs2);
            ASSERT_ALWAYS(CTX.c_local().get_size() <= CTX.c_local().capacity());
            submatrix_range Rc(i0, j0, i1-i0, j1-j0);
            submatrix_range Rct(0,  0, i1-i0, j1-j0);
            logline_printf(1, "ift_C (%u*%u)\n", nrs0, nrs2);
            matpoly_ft<fft_type>::ift(CTX.c_local().view(Rc), tc.view(Rct));
            end_smallstep();
        }
        /* make it compulsory so that we gain some error reporting */
        ASSERT_ALWAYS(local_smallsteps_done(true));

        if (M) end_smallstep();
    }
};

#endif	/* LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_ */
