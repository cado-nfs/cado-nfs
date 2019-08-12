#ifndef LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_
#define LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_

/* The code here is inserted by the compiler inside lingen_matpoly_ft.cpp
 * and lingen_bigmatpoly_ft.cpp, as a template function. The class
 * OP_CTX<T> is instantiated as either OP_CTX<matpoly> or
 * OP_CTX<bigmatpoly>, with specificities in each file (mostly dealing
 * with the extra MPI definitions used in bigmatpoly, as well as the
 * back-end of the mpi operations.
 */

/* For consistency and readability, we add these includes, although
 * they're evidently already included in the two main consumers of this
 * file.
 */
#include "lingen_matpoly.hpp"
#include "lingen_matpoly_ft.hpp"
#include "flint-fft/transform_interface.h"

template<typename T>
struct OP_CTX_base {
    T & c;
    T const & a;
    T const & b;
    const struct fft_transform_info * fti;
    OP_CTX_base(T & c, T const & a, T const & b, const struct fft_transform_info * fti) : c(c), a(a), b(b), fti(fti) {}
};

template<typename T> struct OP_CTX;


/* middle product and multiplication are really the same thing, so better
 * avoid code duplication */

template<typename T>
struct op_mul {/*{{{*/
    size_t csize;
    static const char * name() { return "MUL"; }
    op_mul(T const & a, T const & b, unsigned int adj, fft_transform_info * fti)
    {
        csize = a.size + b.size; csize -= (csize > 0);
        fft_get_transform_info_fppol(fti, a.ab->p, a.size, b.size, a.n);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t)
    {
        ::ift(a, t);
    }
};/*}}}*/
template<typename T>
struct op_mp {/*{{{*/
    size_t csize;
    static const char * name() { return "MP"; }
    unsigned int shift;
    op_mp(T const & a, T const & b, unsigned int adj, fft_transform_info * fti)
    {
        csize = MAX(a.size, b.size) - MIN(a.size, b.size) + 1;
        shift = MIN(a.size, b.size) - 1;
        fft_get_transform_info_fppol_mp(fti, a.ab->p, MIN(a.size, b.size), MAX(a.size, b.size), a.n);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t)
    {
        ::ift_mp(a, t, shift);
    }
};/*}}}*/

template<typename OP_CTX_T, typename OP_T>
static void mp_or_mul(OP_CTX_T & CTX, OP_T & OP, const struct fft_transform_info * fti, const struct lingen_substep_schedule * S)/*{{{*/
{
    typedef typename OP_CTX_T::T T;
    T & c(CTX.c);
    T const & a(CTX.a);
    T const & b(CTX.b);

    CTX.alloc_c_if_needed(OP.csize);

    CTX.mesh_checks();
    const unsigned int r = CTX.mesh_size();

    unsigned int batch = S ? S->batch : a.n;
    unsigned int shrink0 = S ? S->shrink0 : 1;
    unsigned int shrink2 = S ? S->shrink2 : 1;

    subdivision mpi_split0(a.m, r);
    subdivision mpi_split1(a.n, r);
    subdivision mpi_split2(b.n, r);
    subdivision shrink_split0(mpi_split0.block_size_upper_bound(), shrink0);
    subdivision shrink_split2(mpi_split2.block_size_upper_bound(), shrink2);
    /* The order in which we do the transforms is not really our main
     * concern at this point. If sharing makes sense, then probably
     * shrink0 and shrink2 do not. So they're serving opposite purposes.
     */
    /* Declare ta, tb, tc early on so that we don't malloc/free n times.
     */
    matpoly_ft ta, tb, tc;
    const unsigned int nr1 = mpi_split1.block_size_upper_bound();

    /* first, upper bounds on nrs0 and nrs2 */
    unsigned int nrs0 = shrink_split0.block_size_upper_bound();
    unsigned int nrs2 = shrink_split2.block_size_upper_bound();
    tc = matpoly_ft (c.ab, nrs0, nrs2, fti);

    /* We must decide on an ordering beforehand. We cannot do this
     * dynamically because of rounding issues: e.g. for 13=7+6, we
     * will do both 7*6 and 6*7 in the inner loops.
     */
    bool inner_is_row_major = nrs0 < nrs2;
    if (inner_is_row_major) {
        ta = matpoly_ft(a.ab, nrs0, r * batch, fti);
        tb = matpoly_ft(a.ab, r * batch, 1, fti);
    } else {
        ta = matpoly_ft(a.ab, 1, r * batch, fti);
        tb = matpoly_ft(a.ab, r * batch, nrs2, fti);
    }

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
    unsigned int ak0mpi, ak1mpi;
    unsigned int bk0mpi, bk1mpi;
    std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(CTX.a_jrank());
    std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(CTX.b_irank());

    for(unsigned int round0 = 0 ; round0 < shrink0 ; round0++) {
        unsigned int i0, i1;
        std::tie(i0, i1) = shrink_split0.nth_block(round0);
        // TODO: we only have true data for [xi0, xi1[, not [i0, i1[
        // unsigned int xi0 = std::min(i0, imax);
        // unsigned int xi1 = std::min(i1, imax);

        for(unsigned int round2 = 0 ; round2 < shrink2 ; round2++) {
            unsigned int j0, j1;
            std::tie(j0, j1) = shrink_split2.nth_block(round2);
            // TODO: we only have true data for [xj0, xj1[, not [j0, j1[
            // unsigned int xj0 = std::min(j0, jmax);
            // unsigned int xj1 = std::min(j1, jmax);

            submatrix_range Rc(i0, j0, i1-i0, j1-j0);
            submatrix_range Rct(0, 0, i1-i0, j1-j0);

            /* Now do a subblock */
            tc.zero();
            for(unsigned int k = 0 ; k < nr1 ; k += batch) {
                /* In the non-mpi case, ak0==bk0 and ak1==bk1 */
                unsigned int ak0 = ak0mpi + k;
                unsigned int ak1 = MIN(ak1mpi, ak0 + batch);
                unsigned int bk0 = bk0mpi + k;
                unsigned int bk1 = MIN(bk1mpi, bk0 + batch);
                if (inner_is_row_major) {
                    submatrix_range Ra(i0, ak0-ak0mpi, i1-i0, ak1-ak0);
                    submatrix_range Rat(0, CTX.a_jrank()*batch, i1-i0, ak1-ak0);
                    submatrix_range Ratx(0, CTX.a_jrank()*batch, i1-i0, batch);
                    submatrix_range Rbt(CTX.b_irank()*batch, 0, bk1-bk0, 1);
                    submatrix_range Rbtx(CTX.b_irank()*batch, 0, batch, 1);
                    matpoly_ft::view_t ta_loc = ta.view(Rat);
                    matpoly_ft::view_t tb_loc = tb.view(Rbt);

                    /* the true count is ta_loc.size(), but 
                     * upper bounds are:
                     *  i1-i0 <= nrs0
                     *  ak1-ak0 <= batch
                     */
                    CTX.begin_smallstep("dft_A", nrs0 * batch);
                    ta.zero();  // for safety because of rounding.
                    dft(ta_loc, CTX.a_local().view(Ra));
                    CTX.end_smallstep();

                    if (CTX.uses_mpi) {
                        // allgather ta among r nodes.
                        /* r * ta_loc.size() <= r * nrs0 * batch */
                        CTX.begin_smallstep("dft_A_comm", r * nrs0 * batch);
                        CTX.begin_smallstep("export", r * nrs0 * batch);
                        to_export(ta.view(Ratx));
                        CTX.end_smallstep();
                        CTX.begin_smallstep("comm", r * nrs0 * batch);
                        /* The data isn't contiguous, so we have to do
                         * several allgather operations.  */
                        for(unsigned int i = i0 ; i < i1 ; i++) {
                            CTX.do_allgather(ta.part(i - i0, 0), batch);
                        }
                        CTX.end_smallstep();
                        CTX.begin_smallstep("import", r * nrs0 * batch);
                        /* we imported data from everywhere, we can now put
                         * it back to local format */
                        to_import(ta.view(submatrix_range(0, 0, i1-i0, ta.ncols())));
                        CTX.end_smallstep();
                        CTX.end_smallstep();
                    }

                    unsigned int j;
                    for(j = j0 ; j < j1 ; j++) {
                        submatrix_range Rb(bk0-bk0mpi, j, bk1-bk0, 1);
                        /* tb_loc.size = bk1-bk0 <= batch */
                        CTX.begin_smallstep("dft_B", batch);
                        tb.zero();
                        dft(tb_loc, CTX.b_local().view(Rb));
                        CTX.end_smallstep();

                        if (CTX.uses_mpi) {
                            // allgather tb among r nodes
                            /* r * tb_loc.size() <= r * batch */
                            CTX.begin_smallstep("dft_B_comm", r * batch);
                            CTX.begin_smallstep("export", r * batch);
                            to_export(tb.view(Rbtx));
                            CTX.end_smallstep();
                            CTX.begin_smallstep("comm", r * batch);
                            CTX.do_allgather(tb.data, batch);
                            CTX.end_smallstep();
                            CTX.begin_smallstep("import", r * batch);
                            to_import(tb);
                            CTX.end_smallstep();
                            CTX.end_smallstep();
                        }

                        /* (i1-i0) * tb.nrows() <= nrs0 * tb.nrows() */
                        CTX.begin_smallstep("addmul", nrs0 * tb.nrows());
                        // rounding might surprise us.
                        addmul(tc.view(submatrix_range(0, j-j0, i1-i0, 1)), 
                                ta.view(submatrix_range(0, 0, i1-i0, ta.ncols())), 
                                tb.view(submatrix_range(0, 0, tb.nrows(), 1)));
                        CTX.end_smallstep();
                    }
                    for( ; j < j0 + nrs2 ; j++) {
                        CTX.skip_smallstep("dft_B", batch);
                        if (CTX.uses_mpi) {
                            CTX.begin_smallstep("dft_B_comm", r * batch);
                            CTX.skip_smallstep("export", r * batch);
                            CTX.skip_smallstep("comm", r * batch);
                            CTX.skip_smallstep("import", r * batch);
                            CTX.end_smallstep();
                        }
                        CTX.skip_smallstep("addmul", nrs0 * tb.nrows());
                    }
                } else {
                    submatrix_range Rb(bk0-bk0mpi, j0, bk1-bk0, j1-j0);
                    submatrix_range Rbt(CTX.b_irank()*batch, 0, bk1-bk0, j1-j0);
                    submatrix_range Rbtx(CTX.b_irank()*batch, 0, batch, j1-j0);
                    submatrix_range Rat(0, CTX.a_jrank()*batch, 1, ak1-ak0);
                    submatrix_range Ratx(0, CTX.a_jrank()*batch, 1, batch);
                    matpoly_ft::view_t ta_loc = ta.view(Rat);
                    matpoly_ft::view_t tb_loc = tb.view(Rbt);

                    /* the true count is tb_loc.size(), but 
                     * upper bounds are:
                     *  j1-j0 <= nrs2
                     *  bk1-bk0 <= batch
                     */
                    CTX.begin_smallstep("dft_B", nrs2 * batch);
                    tb.zero();
                    dft(tb_loc, CTX.b_local().view(Rb));
                    CTX.end_smallstep();

                    if (CTX.uses_mpi) {
                        // allgather tb among r nodes
                        /* r * tb_loc.size() <= r * nrs2 * batch */
                        CTX.begin_smallstep("dft_B_comm", r * nrs2 * batch);
                        CTX.begin_smallstep("export", r * nrs2 * batch);
                        to_export(tb.view(Rbtx));
                        CTX.end_smallstep();
                        CTX.begin_smallstep("comm", r * nrs2 * batch);
                        CTX.do_allgather(tb.data, nrs2 * batch);
                        CTX.end_smallstep();
                        CTX.begin_smallstep("import", r * nrs2 * batch);
                        to_import(tb.view(submatrix_range(0, 0, tb.nrows(), j1-j0)));
                        CTX.end_smallstep();
                        CTX.end_smallstep();
                    }

                    unsigned int i;
                    for(i = i0 ; i < i1 ; i++) {
                        /* ta_loc.size() = ak1-ak0 <= batch */
                        CTX.begin_smallstep("dft_A", batch);
                        ta.zero();
                        submatrix_range Ra(i, ak0-ak0mpi, 1, ak1-ak0);
                        dft(ta_loc, CTX.a_local().view(Ra));
                        CTX.end_smallstep();

                        if (CTX.uses_mpi) {
                            // allgather ta among r nodes
                            /* r * ta_loc.size() <= r * batch */
                            CTX.begin_smallstep("dft_A_comm", r * batch);
                            CTX.begin_smallstep("export", r * batch);
                            to_export(ta.view(Ratx));
                            CTX.end_smallstep();
                            CTX.begin_smallstep("comm", r * batch);
                            CTX.do_allgather(ta.data, batch);
                            CTX.end_smallstep();
                            CTX.begin_smallstep("import", r * batch);
                            to_import(ta);
                            CTX.end_smallstep();
                            CTX.end_smallstep();
                        }

                        /* (j1-j0) * ta.ncols() <= nrs2 * ta.ncols() */
                        CTX.begin_smallstep("addmul", ta.ncols() * nrs2);
                        addmul(tc.view(submatrix_range(i-i0, 0, 1, j1-j0)), 
                                ta.view(submatrix_range(0, 0, 1, ta.ncols())), 
                                tb.view(submatrix_range(0, 0, tb.nrows(), j1-j0)));
                        CTX.end_smallstep();
                    }
                    for( ; i < i0 + nrs0 ; i++) {
                        CTX.skip_smallstep("dft_A", batch);
                        if (CTX.uses_mpi) {
                            CTX.begin_smallstep("dft_A_comm", r * batch);
                            CTX.skip_smallstep("export", r * batch);
                            CTX.skip_smallstep("comm", r * batch);
                            CTX.skip_smallstep("import", r * batch);
                            CTX.end_smallstep();
                        }
                        CTX.skip_smallstep("addmul", ta.ncols() * nrs2);
                    }
                }
            }
            /* c.size and c_local.size are different fields in the mpi
             * case, and must be kept in sync */
            c.size = OP.csize;
            CTX.c_local().size = OP.csize;
            /* tc.size() <= nrs0 * nrs2 */
            CTX.begin_smallstep("ift_C", nrs0 * nrs2);
            ASSERT_ALWAYS(CTX.c_local().size <= CTX.c_local().alloc);
            OP.ift(CTX.c_local().view(Rc), tc.view(Rct));
            CTX.end_smallstep();
        }
    }
    ASSERT_ALWAYS(CTX.local_smallsteps_done());
}/*}}}*/

#endif	/* LINGEN_MATPOLY_BIGMATPOLY_FT_COMMON_HPP_ */
