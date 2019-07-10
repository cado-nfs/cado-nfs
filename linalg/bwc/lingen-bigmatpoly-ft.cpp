#include "cado.h"
#include "mpfq_layer.h"
#include <cstdlib>
#include "portability.h"
#include "macros.h"
#include "lingen-matpoly-ft.hpp"
#include "lingen-bigmatpoly-ft.hpp"
#include "logline.h"
#include "fmt/format.h"

/* This is the interface for doing products of polynomial matrices by
 * caching transforms, and transferring them over MPI. There are
 * currently two pitfalls.
 *
 * 1/ we have no streaming of computations and communications. This means
 * two things. First, communication/computation overlap can be
 * beneficial: while transforms of rows on the left term of a product are
 * computed (assuming transforms of the right term are ready), it is
 * possible to communicate the already computed row transforms. A second
 * point is that if we manage to get these computed soon, then there is
 * no need to keep them in memory. Such an approach would make it
 * possible to reduce the overall memory footprint. The place to fix this
 * is to modify the control flow in bigmatpoly_mul_caching_adj
 *
 * 2/ we are transferring full-length transforms, while it would 
 * be sufficient to transfer the truncated transforms of course.
 */


#if 0
/* {{{  init/zero/clear interface for bigmatpoly_ft */

/* This completes the initialization process. This is _not_ a collective
 * operation */
void bigmatpoly_ft::finish_init(abdst_field ab, unsigned int m, unsigned int n, const struct fft_transform_info * fti)
{
    this->ab = ab;
    ASSERT_ALWAYS(check_pre_init());
    this->m = m;
    this->n = n;
    this->fti = fti;
    /* XXX we might want to allocate a matrix of size
     * ceiling(m/m1)*ceiling(n/n1) no matter what
     */
    m0 = subdivision(m, m1).nth_block_size(irank());;
    n0 = subdivision(n, n1).nth_block_size(jrank());;
    /* For matpoly_ft, init or finish_init are identical */
    my_cell() = matpoly_ft(ab, m0, n0, fti);
}
bigmatpoly_ft::bigmatpoly_ft(bigmatpoly_model const & model)
    : bigmatpoly_model(model), ab(NULL), m(0), n(0)
{
    m0 = n0 = 0;
}

/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigmatpoly_ft_finish_init
 */
bigmatpoly_ft::bigmatpoly_ft(abdst_field ab, bigmatpoly_model const & model, unsigned int m, unsigned int n, const struct fft_transform_info * fti)
    : bigmatpoly_model(model)
    , ab(ab)
    , m(m)
    , n(n)
    , fti(fti)
{
    ASSERT_ALWAYS(m % m1 == 0);
    ASSERT_ALWAYS(n % n1 == 0);

    cells.reserve(m1*n1);
    for(unsigned int k = m1*n1; k--;) cells.emplace_back();

    /* Either none or all must be non-zero */
    ASSERT_ALWAYS((!m||!n) ^ (m&&n));
    if (!m)
        return;

    /* Allocate only our own cell. If we happen to appear as a left
     * multiplicand someday, we'll need recipient cells for our whole
     * row. If we happen to appear as a right multiplicand, we'll need
     * recipient cells for our whole columns. But there is no reason for
     * both situations (let alone any) to happen. Thus this allocation of
     * recipient cells is dynamic */
    finish_init(ab, m, n, fti);
}

/* This checks the validity of the m,n,size fields. If the structure
 * corresponds to a lazy allocation state, this function returns 1. If
 * the structure corresponds to something fully initializaed already,
 * this function returns 0. If the fields are filled in an inconsistent
 * manner, abort() is called.
 */
int bigmatpoly_ft::check_pre_init() const
{
    if (m && n) return 0;
    if (!m && !n) return 1;
    abort();
    return 0;
}

#if 0
/* We are a left multiplicand. This is a no-op if space for our row has
 * already been allocated */
void bigmatpoly_ft::provision_row()
{
    unsigned int qn = n / n1;
    unsigned int rn = n % n1;
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly_ft & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them = matpoly_ft(ab, m0, qn + (j < rn), fti);
    }
}

/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
void bigmatpoly_ft::provision_col()
{
    unsigned int qm = m / m1;
    unsigned int rm = m % m1;
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly_ft & them = cell(i, jrank());
        if (them.check_pre_init())
            them = matpoly_ft(ab, qm + (i < rm), n0, fti);
    }
}
#endif

/* This zeroes out _our_ cell */
void bigmatpoly_ft::zero()
{
    my_cell().zero();
}

#if 0
/* okay, it's ugly */
void bigmatpoly_ft_swap(bigmatpoly_ft_ptr a, bigmatpoly_ft_ptr b)
{
    bigmatpoly_ft x;
    memcpy(x, a, sizeof(bigmatpoly_ft));
    memcpy(a, b, sizeof(bigmatpoly_ft));
    memcpy(b, x, sizeof(bigmatpoly_ft));
}
#endif

/* }}} */

/* {{{ allgather operations */
void bigmatpoly_ft::allgather_row()
{
    abort();
    provision_row();
    
    /* Each node makes his cell exportable */
    my_cell().to_export();

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    MPI_Datatype mpi_ft;
    MPI_Type_contiguous(tsize, MPI_BYTE, &mpi_ft);
    MPI_Type_commit(&mpi_ft);

    /* TODO: only transfer up to the truncated length ? */
    for(unsigned int k = 0 ; k < n1 ; k++) {
        matpoly_ft & x = cell(irank(), k);
        MPI_Bcast(x.data, x.m * x.n, mpi_ft, k, com[1]);
    }

    MPI_Type_free(&mpi_ft);

    /* and now all nodes on the row import the cells from their friends */
    for(unsigned int j = 0 ; j < n1 ; j++)
        cell(irank(), j).to_import();
}

void bigmatpoly_ft::allgather_col()
{
    abort();
    provision_col();

    /* Each node makes his cell exportable */
    my_cell().to_export();

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    MPI_Datatype mpi_ft;
    MPI_Type_contiguous(tsize, MPI_BYTE, &mpi_ft);
    MPI_Type_commit(&mpi_ft);

    /* TODO: only transfer up to the truncated length ? */
    for(unsigned int k = 0 ; k < m1 ; k++) {
        matpoly_ft & x = cell(k, jrank());
        MPI_Bcast(x.data, x.m * x.n, mpi_ft, k, com[2]);
    }
    MPI_Type_free(&mpi_ft);

    /* and now all nodes on the row import the cells from their friends */
    for(unsigned int i = 0 ; i < m1 ; i++)
        cell(i, jrank()).to_import();
}
/* }}} */

void bigmatpoly_ft::mul(bigmatpoly_ft & a, bigmatpoly_ft & b)/*{{{*/
{
    abort();
    ASSERT_ALWAYS(a.n == b.m);
    ASSERT_ALWAYS(a.n1 == b.m1);
    if (check_pre_init()) {
        finish_init(a.ab, a.m, b.n, fti);
    }
    ASSERT_ALWAYS(m);
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);
    a.allgather_row();
    a.allgather_col();
    // ASSERT_ALWAYS(alloc >= size);

    matpoly_ft & lc = my_cell();
    lc.zero();
    ASSERT_ALWAYS(a.n == b.m);

    for(unsigned int k = 0 ; k < a.n1 ; k++)
        addmul(lc, a.cell(irank(), k), b.cell(k, jrank()));
}/*}}}*/

void bigmatpoly_ft::mul2(bigmatpoly_ft & a, bigmatpoly_ft & b)/*{{{*/
{
    ASSERT_ALWAYS(a.n == b.m);
    ASSERT_ALWAYS(a.n1 == b.m1);
    if (check_pre_init()) {
        finish_init(a.ab, a.m, b.n, fti);
    }
    ASSERT_ALWAYS(m);
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    MPI_Datatype mpi_ft;
    MPI_Type_contiguous(tsize, MPI_BYTE, &mpi_ft);
    MPI_Type_commit(&mpi_ft);

    matpoly_ft & lc = my_cell();
    zero();

    for(unsigned int k = 0 ; k < a.n1 ; k++) {
        logline_printf(1, "MUL-round %u/%u", k, a.n1);

        logline_printf(1, "; row");
        /* Cells irank, * receive data from irank, k for a */
        int leader_a = k == (unsigned int) jrank();
        matpoly_ft & xa = a.cell(irank(), k);
        ASSERT_ALWAYS(leader_a == (xa.data != NULL));
        if (!xa.data) xa = matpoly_ft(a.ab, a.m0, a.n0, fti);
        if (leader_a) xa.to_export();
        MPI_Bcast(xa.data, xa.m * xa.n, mpi_ft, k, a.com[1]);
        xa.to_import();

        logline_printf(1, "; col");
        /* Cells *, jrank receive data from k, jrank for b */
        int leader_b = k == (unsigned int) irank();
        matpoly_ft & xb = b.cell(k, jrank());
        ASSERT_ALWAYS(leader_b == (xb.data != NULL));
        if (!xb.data) xb = matpoly_ft(b.ab, b.m0, b.n0, fti);
        if (leader_b) xb.to_export();
        MPI_Bcast(xb.data, xb.m * xb.n, mpi_ft, k, b.com[2]);
        xb.to_import();

        logline_printf(1, "; addmul");
        /* This is one part of the product.
         * XXX This must be parallelized !!!
         */
        addmul(lc, xa, xb);

        logline_printf(1, "; done\n");
    }

    MPI_Type_free(&mpi_ft);
}/*}}}*/
#endif

#if 0
void bigmatpoly_ft::dft(bigmatpoly const & a)
{
    ::dft(my_cell(), a.my_cell());
}

void bigmatpoly_ft::ift(bigmatpoly & a)
{
    a.set_size(a.size); /* this is weird. This probably has side-effects, but those are not too well documented, and relying on them here seems just obnoxious */
    ::ift(a.my_cell(), my_cell());
}

void bigmatpoly_ft::ift_mp(bigmatpoly & a, unsigned int shift)
{
    a.set_size(a.size);
    ::ift_mp(a.my_cell(), my_cell(), shift);
}
#endif

#if 0
void bigmatpoly_mul_caching_adj(bigmatpoly & c, bigmatpoly & a, bigmatpoly & b, unsigned int adj, const struct lingen_substep_schedule * S MAYBE_UNUSED)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);

    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol(fti, a.ab->p, a.size, b.size, a.n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    logline_printf(1, "FT size=%zuk\n", fft_alloc_sizes[0]>>10);

    bigmatpoly_model const & model = (bigmatpoly_model const&) a;
    c = bigmatpoly(a.ab, model, a.m, b.n, csize);
    bigmatpoly_ft ta(c.ab, model, a.m, a.n, fti);
    bigmatpoly_ft tb(c.ab, model, b.m, b.n, fti);
    bigmatpoly_ft tc(c.ab, model, a.m, b.n, fti);
    logline_printf(1, "DFT(A, %u*%u)\n", a.m0, a.n0);
    ta.dft(a);
    logline_printf(1, "DFT(B, %u*%u)\n", b.m0, b.n0);
    tb.dft(b);
    tc.mul2(ta, tb);
    c.size = csize;
    logline_printf(1, "IFT(C, %u*%u)\n", c.m0, c.n0);
    tc.ift(c);
}/*}}}*/

void bigmatpoly_mp_caching_adj(bigmatpoly & c, bigmatpoly & a, bigmatpoly & b, unsigned int adj, const struct lingen_substep_schedule * S MAYBE_UNUSED)/*{{{*/
{
    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol_mp(fti, a.ab->p, MIN(a.size, b.size), MAX(a.size, b.size), a.n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    logline_printf(1, "FT size=%zuk\n", fft_alloc_sizes[0]>>10);

    bigmatpoly_model const & model = (bigmatpoly_model const&) a;
    c = bigmatpoly(a.ab, model, a.m, b.n, MAX(a.size, b.size) - MIN(a.size, b.size) + 1);
    /* The first grave mistake is here. We should not allocate 3*m*n
     * transforms locally. We don't need that much. */
    bigmatpoly_ft ta(c.ab, model, a.m, a.n, fti);
    bigmatpoly_ft tb(c.ab, model, b.m, b.n, fti);
    bigmatpoly_ft tc(c.ab, model, a.m, b.n, fti);
    logline_printf(1, "DFT(A, %u*%u)\n", a.m0, a.n0);
    ta.dft(a);
    logline_printf(1, "DFT(B, %u*%u)\n", b.m0, b.n0);
    tb.dft(b);
    tc.mul2(ta, tb);
    c.size = MAX(a.size, b.size) - MIN(a.size, b.size) + 1;
    logline_printf(1, "IFT-MP(C, %u*%u)\n", c.m0, c.n0);
    tc.ift_mp(c, MIN(a.size, b.size) - 1);
}/*}}}*/
#endif

/* middle product and multiplication are really the same thing, so better
 * avoid code duplication */

struct op_mul {/*{{{*/
    size_t csize;
    static const char * name() { return "MUL"; }
    op_mul(bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, fft_transform_info * fti)
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
struct op_mp {/*{{{*/
    size_t csize;
    static const char * name() { return "MP"; }
    unsigned int shift;
    op_mp(bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, fft_transform_info * fti)
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

template<typename T>
static void mp_or_mul(T& OP, tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, const struct fft_transform_info * fti, const struct lingen_substep_schedule * S)/*{{{*/
{
    bigmatpoly_model const & model(a.get_model());

    if (c.m != a.m || c.n != a.n || c.size != OP.csize)
        c = bigmatpoly(a.ab, model, a.m, b.n, OP.csize);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    MPI_Datatype mpi_ft;
    MPI_Type_contiguous(tsize, MPI_BYTE, &mpi_ft);
    MPI_Type_commit(&mpi_ft);

    ASSERT_ALWAYS(a.get_model().is_square());
    ASSERT_ALWAYS(a.get_model() == b.get_model());

    const unsigned int r = a.m1;

    subdivision mpi_split0(a.m, r);
    subdivision mpi_split1(a.n, r);
    subdivision mpi_split2(b.n, r);
    subdivision shrink_split0(mpi_split0.block_size_upper_bound(), S->shrink0);
    subdivision shrink_split2(mpi_split2.block_size_upper_bound(), S->shrink2);
    /* The order in which we do the transforms is not really our main
     * concern at this point. If sharing makes sense, then probably
     * shrink0 and shrink2 do not. So they're serving opposite purposes.
     */
    /* Declare ta, tb, tc early on so that we don't malloc/free n times.
     */
    matpoly_ft ta,tb,tc;
    const unsigned int nr1 = mpi_split1.block_size_upper_bound();

    /* first, upper bounds on nrs0 and nrs2 */
    unsigned int nrs0 = shrink_split0.block_size_upper_bound();
    unsigned int nrs2 = shrink_split2.block_size_upper_bound();
    tc = matpoly_ft (c.ab, nrs0, nrs2, fti);

    /* We must decide on an ordering beforehand. We cannot do this
     * dynamically because of rounding issues: e.g. for 13=7+6, we
     * will do both 7*6 and 6*7 in the inner loops.
     */
    stats.begin_smallstep("alloc");
    bool inner_is_row_major = nrs0 < nrs2;
    if (inner_is_row_major) {
        ta = matpoly_ft(a.ab, nrs0, r * S->batch, fti);
        tb = matpoly_ft(a.ab, r * S->batch, 1, fti);
    } else {
        ta = matpoly_ft(a.ab, 1, r * S->batch, fti);
        tb = matpoly_ft(a.ab, r * S->batch, nrs2, fti);
    }
    stats.end_smallstep();

    for(unsigned int round0 = 0 ; round0 < S->shrink0 ; round0++) {
        unsigned int i0,i1;
        std::tie(i0, i1) = shrink_split0.nth_block(round0);
        for(unsigned int round2 = 0 ; round2 < S->shrink2 ; round2++) {
            unsigned int j0,j1;
            std::tie(j0, j1) = shrink_split2.nth_block(round2);
            submatrix_range Rc (i0,j0,i1-i0,j1-j0);
            submatrix_range Rct (0, 0,i1-i0,j1-j0);

            /* Now do a subblock */
            tc.zero();
            for(unsigned int k = 0 ; k < nr1 ; k += S->batch) {
                unsigned int ak0mpi,ak1mpi;
                unsigned int bk0mpi,bk1mpi;
                std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(a.jrank());
                std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(b.irank());
                unsigned int ak0 = ak0mpi + k;
                unsigned int ak1 = MIN(ak1mpi, ak0 + S->batch);
                unsigned int bk0 = bk0mpi + k;
                unsigned int bk1 = MIN(bk1mpi, bk0 + S->batch);
                if (inner_is_row_major) {
                    submatrix_range Ra(i0,ak0-ak0mpi,          i1-i0, ak1-ak0);
                    submatrix_range Rat(0,a.jrank()*S->batch,  i1-i0, ak1-ak0);
                    submatrix_range Ratx(0,a.jrank()*S->batch, i1-i0,S->batch);
                    submatrix_range Rbt(b.irank()*S->batch,0,  bk1-bk0, 1);
                    submatrix_range Rbtx(b.irank()*S->batch,0, S->batch, 1);
                    matpoly_ft::view_t ta_loc = ta.view(Rat);
                    matpoly_ft::view_t tb_loc = tb.view(Rbt);

                    stats.begin_smallstep("dft_A", ta_loc.size());
                    ta.zero();  // for safety because of rounding.
                    dft(ta_loc, a.my_cell().view(Ra));
                    stats.end_smallstep();

                    /* XXX not sure about ncalls here */
                    stats.begin_smallstep("dft_A_comm", ta_loc.size());
                    // allgather ta among r nodes.
                    to_export(ta.view(Ratx));
                    /* The data isn't contiguous, so we have to do
                     * several allgather operations.  */
                    for(unsigned int i = i0 ; i < i1 ; i++) {
                        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                ta.part(i-i0,0), S->batch,
                                mpi_ft, a.get_model().com[1]);
                    }
                    to_import(ta);
                    stats.end_smallstep();

                    for(unsigned int j = 0 ; j < j1-j0 ; j++) {
                        submatrix_range Rb(bk0-bk0mpi,j0+j,bk1-bk0,1);
                        stats.begin_smallstep("dft_B", tb_loc.size());
                        tb.zero();
                        dft(tb_loc, b.my_cell().view(Rb));
                        stats.end_smallstep();

                        /* XXX not sure about ncalls here */
                        stats.begin_smallstep("dft_B_comm", tb_loc.size());
                        // allgather tb among r nodes
                        to_export(tb.view(Rbtx));
                        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                tb.data, S->batch,
                                mpi_ft, b.get_model().com[2]);
                        to_import(tb);
                        stats.end_smallstep();

                        stats.begin_smallstep("addmul");
                        // rounding might surprise us.
                        addmul(tc.view(submatrix_range(0,j,i1-i0,1)),
                                ta.view(submatrix_range(0,0,i1-i0,ta.ncols())),
                                tb.view(submatrix_range(0,0,tb.nrows(),1)));
                        stats.end_smallstep();
                    }
                } else {
                    submatrix_range Rb(bk0-bk0mpi,j0,           bk1-bk0, j1-j0);
                    submatrix_range Rbt(b.irank()*S->batch,0,   bk1-bk0, j1-j0);
                    submatrix_range Rbtx(b.irank()*S->batch,0,  S->batch,j1-j0);
                    submatrix_range Rat(0,a.jrank()*S->batch,1, ak1-ak0);
                    submatrix_range Ratx(0,a.jrank()*S->batch,1,S->batch);
                    matpoly_ft::view_t ta_loc = ta.view(Rat);
                    matpoly_ft::view_t tb_loc = tb.view(Rbt);

                    stats.begin_smallstep("dft_B", tb_loc.size());
                    tb.zero();
                    dft(tb_loc, b.my_cell().view(Rb));
                    stats.end_smallstep();

                    /* XXX not sure about ncalls here */
                    stats.begin_smallstep("dft_B_comm", tb_loc.size());
                    // allgather tb among r nodes
                    to_export(tb.view(Rbtx));
                    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                            tb.data, nrs2 * S->batch,
                            mpi_ft, b.get_model().com[2]);
                    to_import(tb);
                    stats.end_smallstep();

                    for(unsigned int i = i0 ; i < i1 ; i++) {
                        stats.begin_smallstep("dft_A", ta_loc.size());
                        ta.zero();
                        submatrix_range Ra(i,ak0-ak0mpi,1,ak1-ak0);
                        dft(ta_loc, a.my_cell().view(Ra));
                        stats.end_smallstep();

                        /* XXX not sure about ncalls here */
                        stats.begin_smallstep("dft_A_comm", ta_loc.size());
                        // allgather ta among r nodes
                        to_export(ta.view(Ratx));
                        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                ta.data, S->batch,
                                mpi_ft, a.get_model().com[1]);
                        to_import(ta);
                        stats.end_smallstep();

                        stats.begin_smallstep("addmul");
                        addmul(tc.view(submatrix_range(i-i0,0,1,j1-j0)),
                                ta.view(submatrix_range(0,0,1,ta.ncols())),
                                tb.view(submatrix_range(0,0,tb.nrows(),j1-j0)));
                        stats.end_smallstep();
                    }
                }
            }
            c.size = OP.csize;
            c.my_cell().size = OP.csize;
            OP.ift(c.my_cell().view(Rc), tc.view(Rct));
        }
    }
    MPI_Type_free(&mpi_ft);
}/*}}}*/

void bigmatpoly_mp_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S)/*{{{*/
{
    struct fft_transform_info fti[1];
    op_mp OP(a, b, adj, fti);
    mp_or_mul(OP, t, c, a, b, fti, S);
} /* }}} */
void bigmatpoly_mul_caching_adj(tree_stats & t, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, const struct lingen_substep_schedule * S)/*{{{*/
{
    struct fft_transform_info fti[1];
    op_mul OP(a, b, adj, fti);
    mp_or_mul(OP, t, c, a, b, fti, S);
} /* }}} */

