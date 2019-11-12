#include "cado.h"
#include <cstdlib>
#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "logline.h"
#include "lingen_matpoly_select.hpp"
#include "lingen_bigmatpoly.hpp"
#include "lingen_mul_substeps.hpp"

int bigmatpoly_model::rank() const
{
    int x;
    MPI_Comm_rank(com[0], &x);
    return x;
}

int bigmatpoly_model::irank() const
{
    int x;
    MPI_Comm_rank(com[2], &x);
    return x;
}

int bigmatpoly_model::jrank() const
{
    int x;
    MPI_Comm_rank(com[1], &x);
    return x;
}


/* {{{  init/zero/clear interface for bigmatpoly */

bigmatpoly_model::bigmatpoly_model(MPI_Comm * comm, unsigned int m, unsigned int n)
{
    m1 = m;
    n1 = n;
    memcpy(com, comm, 3 * sizeof(MPI_Comm));
}

/* This completes the initialization process. This is _not_ a collective
 * operation */
void bigmatpoly::finish_init(abdst_field ab, unsigned int m, unsigned int n, int len)
{
    this->ab = ab;
    ASSERT_ALWAYS(check_pre_init());
    this->m = m;
    this->n = n;
    // m0 = subdivision(m, m1).nth_block_size(irank());;
    // n0 = subdivision(n, n1).nth_block_size(jrank());;
    m0 = subdivision(m, m1).block_size_upper_bound();
    n0 = subdivision(n, n1).block_size_upper_bound();
    size = len;
    my_cell() = matpoly(ab, m0, n0, len);
}

/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigmatpoly_finish_init
 */
bigmatpoly::bigmatpoly(bigmatpoly_model const & model)
    : bigmatpoly_model(model), ab(NULL), m(0), n(0)
{
    m0 = n0 = 0;
    size = 0;
    /* except that finish_init wants this allocated. We could, of course,
     * do the reservation in finish_init, but better have both ctors
     * leave the same post-condition.
     */
    cells.reserve(m1*n1);
    for(unsigned int k = m1*n1; k--;) cells.emplace_back();
}

bigmatpoly::bigmatpoly(abdst_field ab, bigmatpoly_model const & model, unsigned int m, unsigned int n, int len)
    : bigmatpoly_model(model), ab(ab), m(m), n(n)
{
    cells.reserve(m1*n1);
    for(unsigned int k = m1*n1; k--;) cells.emplace_back();
    // m0 = subdivision(m, m1).nth_block_size(irank());;
    // n0 = subdivision(n, n1).nth_block_size(jrank());;
    m0 = subdivision(m, m1).block_size_upper_bound();
    n0 = subdivision(n, n1).block_size_upper_bound();
    /* Either none or all must be non-zero */
    ASSERT_ALWAYS((!m||!n||!len) ^ (m&&n&&len));
    if (!len)
        return;

    /* Allocate only our own cell. If we happen to appear as a left
     * multiplicand someday, we'll need recipient cells for our whole
     * row. If we happen to appear as a right multiplicand, we'll need
     * recipient cells for our whole columns. But there is no reason for
     * both situations (let alone any) to happen. Thus this allocation of
     * recipient cells is dynamic */
    finish_init(ab, m, n, len);
}

bigmatpoly::bigmatpoly(bigmatpoly && a)
    : bigmatpoly_model(a.get_model())
    , ab(a.ab)
    , m(a.m), n(a.n)
    , m0(a.m0), n0(a.n0)
{
    size=a.size;
    a.m0=a.n0=a.m=a.n=a.size=0;
    std::swap(cells, a.cells);
}
bigmatpoly& bigmatpoly::operator=(bigmatpoly&& a)
{
    get_model() = a.get_model();
    ab = a.ab;
    m = a.m;
    n = a.n;
    m0 = a.m0;
    n0 = a.n0;
    size=a.size;
    a.m0=a.n0=a.m=a.n=a.size=0;
    std::swap(cells, a.cells);
    return *this;
}
/* }}} */

#if 0
/* Return a bitmask indicating whether bigmatpoly_provision_{row,col} has
 * been called on this matrix before. bit 0 is for row, bit 1 is for col.
 * If the returned value is zero, then we really have only the local part
 * of this matrix for the moment.
 */
int bigmatpoly::provisioned() const
{
    if (check_pre_init()) return 0;
    unsigned int np_row = 0;
    unsigned int np_col = 0;
    for(unsigned int j = 0 ; j < n1 ; j++)
        np_row += !cell(irank(), j).check_pre_init();
    for(unsigned int i = 0 ; i < m1 ; i++)
        np_col += !cell(i, jrank()).check_pre_init();
    ASSERT_ALWAYS(np_row == 1 || np_row == n1);
    ASSERT_ALWAYS(np_col == 1 || np_col == m1);
    return (np_row>1)+((np_col>1)<<1);
}
#endif

/* We are a left multiplicand. This is a no-op if space for our row has
 * already been allocated */
void bigmatpoly::provision_row()
{
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly & them = cell(irank(), j);
        if (them.check_pre_init())
            them = matpoly(ab, m0, n0, my_cell().capacity());
    }
}

#if 0
/* Rarely useful. We do need it because we resort to a kludgy
 * implementation of scatter_mat, which calls for provisioning on all
 * rows.
 */
void bigmatpoly::unprovision_row()
{
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly & them = cell(irank(), j);
        if (!them.check_pre_init())
            them = matpoly();
    }
}
#endif

/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
void bigmatpoly::provision_col() // bigmatpoly & p(*this)
{
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly & them = cell(i, jrank());
        if (them.check_pre_init())
            them = matpoly(ab, m0, n0, my_cell().capacity());
    }
}

/* Set size to be large enough to receive the given number of
 * coefficients. If space is provisioned for other cells in row or
 * column, allocation is triggered for them as well.
 */
void bigmatpoly::set_size(size_t nsize)
{
    matpoly & me = my_cell();
    ASSERT_ALWAYS(nsize <= me.capacity());
    size = nsize;
    me.set_size(nsize);
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.set_size(nsize);
        ASSERT_ALWAYS(nsize <= them.capacity());
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.set_size(nsize);
        ASSERT_ALWAYS(nsize <= them.capacity());
    }
}
void bigmatpoly::zero_pad(size_t nsize)/*{{{*/
{
    matpoly & me = my_cell();
    size = nsize;
    me.zero_pad(nsize);
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.zero_pad(nsize);
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.zero_pad(nsize);
    }
}

/* If our row or col cells have already been allocated, then reallocate
 * them as well (XXX is it clear or not ?) */
#if 0 /* This function has never been used or needed */
void bigmatpoly_realloc(bigmatpoly & p, int newalloc)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    matpoly_ptr me = bigmatpoly_my_cell(p);
    matpoly_realloc(me, newalloc);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        matpoly_realloc(them, newalloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (!them->x) continue;
        matpoly_realloc(them, newalloc);
    }
}
#endif

/* This zeroes out _our_ cell */
void bigmatpoly::zero()
{
    my_cell().zero();
}

/* okay, it's ugly */
#if 0
void bigmatpoly::swap(bigmatpoly & a)
{
    bigmatpoly x = std::move(*this);
    *this = std::move(a);
    a = std::move(x);
}
#endif
/* }}} */

void bigmatpoly::truncate(bigmatpoly const & src, unsigned int nsize)/*{{{*/
{
    // ASSERT_ALWAYS(src.provisioned());
    // ASSERT_ALWAYS(provisioned());
    if (check_pre_init()) {
        finish_init(src.ab, src.m, src.n, nsize);
    }
    my_cell().truncate(src.my_cell(), nsize);
    size = my_cell().get_size();
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly const & sthem = src.cell(irank(), j);
        matpoly & dthem = cell(irank(), j);
        if (sthem.check_pre_init()) {
            dthem.clear();
        } else if (dthem.check_pre_init()) {
            continue;
        } else {
            dthem.truncate(sthem, nsize);
        }
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly const & sthem = src.cell(i, jrank());
        matpoly & dthem = cell(i, jrank());
        if (sthem.check_pre_init()) {
            dthem.clear();
        } else if (dthem.check_pre_init()) {
            continue;
        } else {
            dthem.truncate(sthem, nsize);
        }
    }
}
/*}}}*/

bool bigmatpoly::high_word_is_clear() const
{
    matpoly const & me = my_cell();
    if (!me.high_word_is_clear()) return false;
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        if (!cell(irank(), j).high_word_is_clear()) return false;
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        if (!cell(i, jrank()).high_word_is_clear()) return false;
    }
    return true;
}

void bigmatpoly::clear_high_word()
{
    matpoly & me = my_cell();
    me.clear_high_word();
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        matpoly & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.clear_high_word();
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        matpoly & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.clear_high_word();
    }
}

void bigmatpoly::coeff_set_zero_loc(unsigned int k)
{
    my_cell().coeff_set_zero(k);
}

int bigmatpoly::coeff_is_zero(unsigned int k) const
{
    int t = my_cell().coeff_is_zero(k);
    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_INT, MPI_MIN, com[0]);
    return t;
}


void bigmatpoly::rshift(bigmatpoly & src, unsigned int k) /*{{{*/
{
    matpoly & me = my_cell();
    if (check_pre_init())
        finish_init(src.ab, src.m, src.n, src.size - k);
    // ASSERT_ALWAYS(provisioned() == 0);
    // ASSERT_ALWAYS(src.provisioned() == 0);
    me.rshift(src.my_cell(), k);
    size = me.get_size();
}
/*}}}*/

#define CHECK_MPI_DATASIZE_FITS(_size0, _size1, _type0, _code) do {	\
    ASSERT_ALWAYS((size_t) _size0 <= (size_t) INT_MAX);			\
    ASSERT_ALWAYS((size_t) _size1 <= (size_t) INT_MAX);			\
    size_t _datasize = (size_t) _size0 * (size_t) _size1;		\
    if (_datasize > (size_t) INT_MAX) {					\
        MPI_Datatype _datatype;						\
        MPI_Type_contiguous(_size1, _type0, &_datatype);		\
        MPI_Type_commit(&_datatype);					\
        int _datasize = _size1;						\
        _code;								\
        MPI_Type_free(&_datatype);					\
    } else {								\
        MPI_Datatype _datatype = _type0;				\
        _code;								\
    }									\
} while (0)

/* {{{ allgather operations */
void bigmatpoly::allgather_row()
{
    provision_row();
    for(unsigned int k = 0 ; k < n1 ; k++) {
        matpoly & data = cell(irank(), k);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long dsize = data.get_size();
        MPI_Bcast(&dsize, 1, MPI_UNSIGNED_LONG, k, com[1]);
        data.set_size(dsize);
        ASSERT_ALWAYS(data.get_size() <= data.capacity());
        ASSERT_ALWAYS((data.m * data.n * data.capacity()) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data.m * data.n, data.data_entry_alloc_size(),
                MPI_BYTE,
                MPI_Bcast(data.x, _datasize, _datatype, k, com[1])
        );
    }
}
void bigmatpoly::allgather_col()
{
    provision_col();
    for(unsigned int k = 0 ; k < m1 ; k++) {
        matpoly & data = cell(k, jrank());
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long dsize = data.get_size();
        MPI_Bcast(&dsize, 1, MPI_UNSIGNED_LONG, k, com[2]);
        data.set_size(dsize);
        ASSERT_ALWAYS(data.get_size() <= data.capacity());
        ASSERT_ALWAYS((data.m * data.n * data.capacity()) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data.m * data.n, data.data_entry_alloc_size(),
                MPI_BYTE,
                MPI_Bcast(data.x, _datasize, _datatype, k, com[2])
        );
    }
}
/* }}} */

bigmatpoly bigmatpoly::mul(bigmatpoly & a, bigmatpoly & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);
    ASSERT_ALWAYS(a.n == b.m);
    ASSERT_ALWAYS(a.n1 == b.m1);
    bigmatpoly c(a.ab, a.get_model(), a.m, b.n, csize);
    a.allgather_row();
    a.allgather_col();
    c.set_size(csize);
    matpoly & lc = c.my_cell();
    lc.zero();
    c.set_size(csize);
    ASSERT_ALWAYS(a.n == b.m);
    unsigned int i = c.irank();
    unsigned int j = c.jrank();
    for(unsigned int k = 0 ; k < a.n1 ; k++) {
        lc.addmul(a.cell(i, k), b.cell(k, j));
    }
    return c;
}/*}}}*/

bigmatpoly bigmatpoly::mp(bigmatpoly & a, bigmatpoly & c) /*{{{*/
{
    unsigned bsize = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    ASSERT_ALWAYS(a.n1 == c.m1);
    bigmatpoly b(a.ab, a.get_model(), a.m, c.n, bsize);

    a.allgather_row();
    c.allgather_col();

    matpoly & lb = b.my_cell();
    lb.zero();
    b.set_size(bsize);
    ASSERT_ALWAYS(lb.get_size() == bsize);
    ASSERT_ALWAYS(lb.capacity() >= lb.get_size());
    unsigned int i = c.irank();
    unsigned int j = c.jrank();
    for(unsigned int k = 0 ; k < a.n1 ; k++) {
        lb.addmp(a.cell(i, k), c.cell(k, j));
    }
    return b;
}/*}}}*/

struct OP_CTX {
    bigmatpoly & c;
    bigmatpoly const & a;
    bigmatpoly const & b;
    MPI_Datatype mpi_entry_a;
    MPI_Datatype mpi_entry_b;
    tree_stats & stats;
    OP_CTX(tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b)
        : c(c), a(a), b(b), stats(stats)
    {
        MPI_Type_contiguous(a.my_cell().data_entry_alloc_size(), MPI_BYTE, &mpi_entry_a);
        MPI_Type_contiguous(b.my_cell().data_entry_alloc_size(), MPI_BYTE, &mpi_entry_b);
        MPI_Type_commit(&mpi_entry_a);
        MPI_Type_commit(&mpi_entry_b);
    }
    ~OP_CTX() {
        MPI_Type_free(&mpi_entry_a);
        MPI_Type_free(&mpi_entry_b);
    }
    inline int a_irank() const { return a.irank(); }
    inline int b_irank() const { return b.irank(); }
    inline int a_jrank() const { return a.jrank(); }
    inline int b_jrank() const { return b.jrank(); }
    inline int mesh_inner_size() const { return a.n1; }
    static const bool uses_mpi = true;
    inline void mesh_checks() const {
        ASSERT_ALWAYS(a.get_model().is_square());
        ASSERT_ALWAYS(a.get_model() == b.get_model());
        ASSERT_ALWAYS(a.irank() == b.irank());
        ASSERT_ALWAYS(a.jrank() == b.jrank());
    }
    void alloc_c_if_needed(size_t size) const {
        if (c.m != a.m || c.n != a.n || c.get_size() != size)
            c = bigmatpoly(a.ab, a.get_model(), a.m, b.n, size);
    }
    inline matpoly const & a_local() { return a.my_cell(); }
    inline matpoly const & b_local() { return b.my_cell(); }
    inline matpoly & c_local() { return c.my_cell(); }
    inline void a_allgather(void * p, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               p, n, mpi_entry_a, a.get_model().com[1]);
    }
    inline void b_allgather(void * p, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               p, n, mpi_entry_b, b.get_model().com[2]);
    }
};
template<typename OP_T> struct mp_or_mul : public OP_CTX { 
    OP_T & OP;
    const lingen_call_companion::mul_or_mp_times * M;
    subdivision mpi_split0;
    subdivision mpi_split1;
    subdivision mpi_split2;
    unsigned int nrs0, nrs2;
    unsigned int b0, b1, b2;
    matpoly a_peers;
    matpoly b_peers;
    /* Declare ta, tb, tc early on so that we don't malloc/free n times.  */
    mp_or_mul(OP_CTX & CTX, OP_T & OP,
            const lingen_call_companion::mul_or_mp_times * M)
        : OP_CTX(CTX)
        , OP(OP)
        , M(M)
        , mpi_split0(a.m, mesh_inner_size())
        , mpi_split1(a.n, mesh_inner_size())
        , mpi_split2(b.n, mesh_inner_size())
        /* first, upper bounds on output block dimensions */
        , nrs0(mpi_split0.block_size_upper_bound())
        , nrs2(mpi_split2.block_size_upper_bound())
        /* By default we use full batching, which costs some memory ! */
        , b0(M ? M->S.batch[0] : nrs0)
        , b1(M ? M->S.batch[1] : a.n)
        , b2(M ? M->S.batch[2] : nrs2)
        , a_peers(a.ab, b0, mesh_inner_size() * b1, a.my_cell().capacity())
        , b_peers(b.ab, mesh_inner_size() * b1, b2, b.my_cell().capacity())
    {
        mesh_checks();

        if (!M) return;

        constexpr const char * opname = OP_T::name;

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
        begin_plan_smallstep_microsteps(opname);
        plan_smallstep("gather_A", M->t_dft_A);
        plan_smallstep("gather_B", M->t_dft_B);
        plan_smallstep("addmul", M->t_conv);
        plan_smallstep("ift_C", M->t_ift_C);
        end_plan_smallstep();
    }
    template<typename... Args>
    inline void begin_smallstep(Args&& ...args) {
        if (M) stats.begin_smallstep(args...);
    }
    template<typename... Args>
    inline void skip_smallstep(Args&& ...args) {
        if (M) stats.skip_smallstep(args...);
    }
    template<typename... Args>
    inline void end_smallstep(Args&& ...args) {
        if (M) stats.end_smallstep(args...);
    }
    template<typename... Args>
    inline void plan_smallstep(Args&& ...args) {
        if (M) stats.plan_smallstep(args...);
    }
    template<typename... Args>
    inline void begin_plan_smallstep_microsteps(Args&& ...args) {
        if (M) stats.begin_plan_smallstep_microsteps(args...);
    }
    template<typename... Args>
    inline void begin_plan_smallstep(Args&& ...args) {
        if (M) stats.begin_plan_smallstep(args...);
    }
    template<typename... Args>
    inline void end_plan_smallstep(Args&& ...args) {
        if (M) stats.end_plan_smallstep(args...);
    }
    inline bool local_smallsteps_done(bool compulsory = false) {
        return M ? stats.local_smallsteps_done(compulsory) : true;
    }


    /* loop0 and loop2 depend on the exact output block. At most we're
     * iterating on, respectively, ceil(ceil(n0/r)), and
     * ceil(ceil(n2/r)).
     */
    subdivision loop0;
    subdivision loop1;
    subdivision loop2;

    /*
     * 
     * Entries in
     * the output block are processed as blocks of size b0*b2. For each
     * these, b1 pairs of input data are used at the same time
     * (that is, b1*(b0+b2)), and it total we collect
     * mesh_inner_size() as many from the MPI peers. The order in
     * which the blocks of size b0*b2 are processed to cover the range of
     * size of nrs0 * nrs2 output values controls the amount of input
     * data we have to fetch in total, namely:
     *  ceil(nrs0/b0)*(b0+b2)*b1 if we process blocks row-major.
     *  ceil(nrs2/b2)*(b0+b2)*b1 if we process blocks col-major.
     *
     * Note though that this is a rather exaggerated notion: we *must*
     * have either b0==nrs0 or n2==nrs2.
     * Therefore, if we heed this adjustment, the number of transforms
     * that are computed to process a block of size b0*b2 in the output
     * is always (b0+b2)*b1, with the processing order being determined
     * by which of b0==nrs0 or b2==nrs2 holds.
     */
    void gather_A(unsigned int i0, unsigned int iloop0, unsigned int iloop1)/*{{{*/
    {
        begin_smallstep("gather_A", b0 * b1);
        unsigned int aj = a_jrank();
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

        /* This is the analogue of the "dft_A" step in the transform
         * case.
         */
        a_peers.zero();  // for safety because of rounding.
        matpoly::copy(a_peers.view(Rat), a_local().view(Ra));

        // allgather ta among r nodes. No serialization needed here.
        /* The data isn't contiguous, so we have to do
         * several allgather operations.  */
        for(unsigned int i = 0 ; i < ii1 - ii0 ; i++) {
            a_allgather(a_peers.part(i, 0), b1);
        }
        end_smallstep();
    }/*}}}*/

    void gather_B(unsigned int j0, unsigned int iloop1, unsigned int iloop2)/*{{{*/
    {
        begin_smallstep("gather_A", b1 * b2);
        unsigned int bi = b_irank();
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

        b_peers.zero();
        matpoly::copy(b_peers.view(Rbt), b_local().view(Rb));

        // allgather tb among r nodes
        b_allgather(b_peers.part(0, 0), b1 * b2);
        end_smallstep();
    }/*}}}*/

    void addmul_for_block(matpoly::view_t & cdst, unsigned int iloop0, unsigned int iloop2)/*{{{*/
    {
        const unsigned int r = mesh_inner_size();
        unsigned int ii0, ii1;
        unsigned int jj0, jj1;
        std::tie(ii0, ii1) = loop0.nth_block(iloop0);
        std::tie(jj0, jj1) = loop2.nth_block(iloop2);

        begin_smallstep("addmul", b0 * b1 * b2 * r);

        // rounding might surprise us.
        submatrix_range Ratxx(0,   0,   ii1 - ii0, r*b1);
        submatrix_range Rbtxx(0,   0,   r*b1,      jj1 - jj0);
        submatrix_range Rct  (cdst.i0 + ii0, cdst.j0 + jj0, ii1 - ii0, jj1 - jj0);

        matpoly::view_t c_loc(cdst.M, Rct);
        matpoly::view_t a_loc = a_peers.view(Ratxx);
        matpoly::view_t b_loc = b_peers.view(Rbtxx);

        OP_T::addcompose(c_loc, a_loc, b_loc);

        end_smallstep();
    }/*}}}*/

    void operator()() {
        constexpr const char * opname = OP_T::name;
        begin_smallstep(opname);

        alloc_c_if_needed(OP.csize);

        const unsigned int r = mesh_inner_size();

        /* The order in which we do the transforms is not really our main
         * concern at this point. If sharing makes sense, then probably
         * shrink0 and shrink2 do not. So they're serving opposite purposes.
         */
        const unsigned int nr1 = mpi_split1.block_size_upper_bound();
        loop1 = subdivision::by_block_size(nr1, b1);

        // unsigned int imax = mpi_split0.nth_block_size(a_irank());
        // unsigned int jmax = mpi_split2.nth_block_size(b_jrank());

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
        unsigned int aj = a_jrank();
        unsigned int bi = b_jrank();
        unsigned int ak0mpi, ak1mpi;
        unsigned int bk0mpi, bk1mpi;
        std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(aj);
        std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(bi);

        /* Prepare the processing of the small blocks of size b0*b2
        */
        unsigned int i0 = 0, i1 = mpi_split0.block_size_upper_bound();
        unsigned int j0 = 0, j1 = mpi_split2.block_size_upper_bound();
        /* Note that i0, i1, j0, j1 are equal for all mpi jobs. Therefore
         * the three loops below are synchronous for all jobs. */
        ASSERT_ALWAYS((i1 - i0) <= nrs0);
        ASSERT_ALWAYS((j1 - j0) <= nrs2);
        loop0 = subdivision::by_block_size(i1 - i0, b0);
        loop2 = subdivision::by_block_size(j1 - j0, b2);

        ASSERT_ALWAYS(loop0.nblocks() == 1 || loop2.nblocks() == 1);
        bool process_blocks_row_major = b0 == nrs0;

        /* Now do a subblock */
        submatrix_range Rc(i0, j0, i1-i0, j1-j0);
        matpoly::view_t cdst = c_local().view(Rc);
        cdst.zero();

        for(unsigned int iloop1 = 0 ; iloop1 < loop1.nblocks() ; iloop1++) {
            if (process_blocks_row_major) {
                ASSERT_ALWAYS(loop0.nblocks() == 1);
                ASSERT_ALWAYS(nrs0 == b0);
                unsigned int iloop0 = 0;
                logline_printf(1, "gather_A (%u*%u)\n", b0, b1);
                gather_A(i0, iloop0, iloop1);
                for(unsigned int iloop2 = 0 ; iloop2 < loop2.nblocks() ; iloop2++) {
                    logline_printf(1, "gather_B (%u*%u)\n", b1, b2);
                    gather_B(j0, iloop1, iloop2);

                    logline_printf(1, "addmul\n");
                    addmul_for_block(cdst, iloop0, iloop2);
                }
                /* adjust counts */
                unsigned int e2 = (iceildiv(nrs2, b2) - loop2.nblocks());
                unsigned int x2 = e2 * b2;
                skip_smallstep("gather_B", b1 * x2);
                skip_smallstep("addmul", b0 * b1 * x2 * r);
            } else {
                ASSERT_ALWAYS(loop2.nblocks() == 1);
                ASSERT_ALWAYS(nrs2 == b2);
                unsigned int iloop2 = 0;
                logline_printf(1, "gather_B (%u*%u)\n", b1, b2);
                gather_B(j0, iloop1, iloop2);
                for(unsigned int iloop0 = 0 ; iloop0 < loop0.nblocks() ; iloop0++) {
                    logline_printf(1, "gather_A (%u*%u)\n", b0, b1);
                    gather_A(i0, iloop0, iloop1);

                    logline_printf(1, "addmul\n");
                    addmul_for_block(cdst, iloop0, iloop2);
                }
                /* adjust counts */
                unsigned int e0 = (iceildiv(nrs0, b0) - loop0.nblocks());
                unsigned int x0 = e0 * b0;
                skip_smallstep("gather_A", x0 * b1);
                skip_smallstep("addmul", x0 * b1 * b2 * r);
            }
        }

        c.set_size(OP.csize);
        /* make it compulsory so that we gain some error reporting */
        ASSERT_ALWAYS(local_smallsteps_done(true));

        end_smallstep();
    }
};

bigmatpoly bigmatpoly::mp(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M)
{
    op_mp<void> op(a, b, UINT_MAX);
    bigmatpoly c(a.get_model());
    OP_CTX CTX(stats, c, a, b);
    mp_or_mul<op_mp<void>>(CTX, op, M)();
    return c;
}

bigmatpoly bigmatpoly::mul(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M)
{
    op_mul<void> op(a, b, UINT_MAX);
    bigmatpoly c(a.get_model());
    OP_CTX CTX(stats, c, a, b);
    mp_or_mul<op_mul<void>>(CTX, op, M)();
    return c;
}

/*
 * gather to node 0, or scatter from node 0, but use "partial" transfer
 * in order to allow some kind of streaming between reading a matrix on
 * node 0 and scaterring it on all the nodes. (and the same for gather /
 * writing to file on node 0).
 */

/* The piece [offset, offset+length[ of the bigmatpoly source is gathered
 * in the matpoly dst on node 0.
 * We assume that all the data structures are already set up properly,
 * and that dst has indeed room for length elements (we set the size, but
 * don't realloc dst).
 */
void bigmatpoly::gather_mat_partial(matpoly & dst,
        size_t dst_k,
        size_t offset, size_t length_raw) const
{
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
#ifdef SELECT_MPFQ_LAYER_u64k1
    static_assert(std::is_same<absrc_vec, const unsigned long *>::value, "uh ?");
#endif

    /* sanity checks, because the code below assumes this. */
    ASSERT_ALWAYS(irank() * (int) n1 + jrank() == rank());
    ASSERT_ALWAYS(length_raw <= (size_t) INT_MAX);

    size_t length = simd * iceildiv(length_raw, simd);

    MPI_Datatype mt;
#ifndef SELECT_MPFQ_LAYER_u64k1
    MPI_Type_contiguous(length * abvec_elt_stride(ab, 1), MPI_BYTE, &mt);
#else
    ASSERT_ALWAYS(length % simd == 0);
    MPI_Type_contiguous(length / simd, MPI_UNSIGNED_LONG, &mt);
#endif
    MPI_Type_commit(&mt);

    matpoly const & me = my_cell();

    subdivision R(m, m1);       /* Row split */
    subdivision C(n, n1);       /* Col split */

    // Node 0 receives data
    if (!rank()) {
        ASSERT_ALWAYS(dst.m == m);
        ASSERT_ALWAYS(dst.n == n);
        ASSERT_ALWAYS(dst.capacity() >= dst_k + length);
        MPI_Request * reqs = new MPI_Request[m * n];
        MPI_Request * req = reqs;
        /* the master receives data from everyone */
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int ii = R.flatten(i1, i0);
                        unsigned int jj = C.flatten(j1, j0);
                        unsigned int peer = i1 * n1 + j1;
                        unsigned int tag = ii * n + jj;
                        ASSERT_ALWAYS(offset + length <= me.capacity());
                        abdst_vec to = abvec_subvec(ab, dst.part(ii, jj), dst_k);
                        absrc_vec from = abvec_subvec_const(ab, me.part(i0, j0), offset);
                        if (peer == 0) {
                            /* talk to ourself */
                            abvec_set(ab, to, from, length);
                        } else {
                            MPI_Irecv(to, 1, mt, peer, tag, com[0], req);
                        }
                        req++;
                    }
                }
            }
        }

        req = reqs;
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int peer = i1 * n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
        delete[] reqs;
    } else {
        // All the other nodes send their data.
        MPI_Request * reqs = new MPI_Request[m0r() * n0r()];
        MPI_Request * req = reqs;
        /* receive. Each job will receive exactly dst.m0 transfers */
        for(unsigned int i0 = 0 ; i0 < m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < n0r() ; j0++) {
                unsigned int ii = R.flatten(irank(), i0);
                unsigned int jj = C.flatten(jrank(), j0);
                unsigned int tag = ii * n + jj;
                absrc_vec from = abvec_subvec_const(ab, me.part(i0, j0), offset);
                /* battle const-deprived MPI prototypes... */
                MPI_Isend((void*)from, 1, mt, 0, tag, com[0], req);
                req++;
            }
        }
        req = reqs;
        for(unsigned int i0 = 0 ; i0 < m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < n0r() ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
        delete[] reqs;
    }
    MPI_Type_free(&mt);
}

/* Exactly the converse of the previous function.
 * Take length element in the src matrix on node 0, and scatter it
 * in dst, with the given offset.
 * We assume that dst has been initialized: all the communicators, mn's,
 * are already set and enough space to accomodate length+offset elements
 * have been already allocated. The only non-data field of dst that is
 * modified is size.
 */
void bigmatpoly::scatter_mat_partial(
        matpoly const & src,
        size_t src_k,
        size_t offset, size_t length_raw)
{
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
#ifdef SELECT_MPFQ_LAYER_u64k1
    static_assert(std::is_same<absrc_vec, const unsigned long *>::value, "uh ?");
#endif
    /* The length is not necessarily aligned on the simd width */
    size_t length = simd * iceildiv(length_raw, simd);

    MPI_Datatype mt;

#ifndef SELECT_MPFQ_LAYER_u64k1
    MPI_Type_contiguous(length * abvec_elt_stride(ab, 1), MPI_BYTE, &mt);
#else
    MPI_Type_contiguous(length / simd, MPI_UNSIGNED_LONG, &mt);
#endif
    MPI_Type_commit(&mt);

    /* sanity check, because the code below assumes this. */
    ASSERT_ALWAYS(irank() * (int) n1 + jrank() == rank());

    matpoly & me = my_cell();

    subdivision R(m, m1);       /* Row split */
    subdivision C(n, n1);       /* Col split */

    if (!rank()) {
        MPI_Request * reqs = new MPI_Request[m * n];
        MPI_Request * req = reqs;
        /* the master sends data to everyone */
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int ii = R.flatten(i1, i0);
                        unsigned int jj = C.flatten(j1, j0);
                        unsigned int peer = i1 * n1 + j1;
                        unsigned int tag = ii * n + jj;
                        absrc_vec from = abvec_subvec_const(ab, src.part(ii, jj), src_k);
                        abdst_vec to = abvec_subvec(ab, me.part(i0, j0), offset);

                        if (peer == 0) {
                            /* talk to ourself */
                            abvec_set(ab, to, from, length);
                        } else {
                            /* battle const-deprived MPI prototypes... */
                            MPI_Isend((void*) from, 1, mt, peer, tag, com[0], req);
                        }
                        req++;
                    }
                }
            }
        }

        req = reqs;
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int peer = i1 * n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
        delete[] reqs;
    } else {
        MPI_Request * reqs = new MPI_Request[m0r() * n0r()];
        MPI_Request * req = reqs;
        for(unsigned int i0 = 0 ; i0 < m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < n0r() ; j0++) {
                unsigned int ii = R.flatten(irank(), i0);
                unsigned int jj = C.flatten(jrank(), j0);
                unsigned int tag = ii * n + jj;
                abdst_vec to = abvec_subvec(ab, me.part(i0, j0), offset);
                MPI_Irecv(to, 1, mt, 0, tag, com[0], req);
                req++;
            }
        }
        req = reqs;
        for(unsigned int i0 = 0 ; i0 < m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < n0r() ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
        delete[] reqs;
    }
    MPI_Type_free(&mt);
}

/* Collect everything into node 0 */
void bigmatpoly::gather_mat(matpoly & dst) const
{
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
    matpoly dst_partial;
#ifndef SELECT_MPFQ_LAYER_u64k1
    size_t length = 100;
#else
    size_t length = 1024;
#endif

    if (!rank()) {
        // Leader should initialize the result matrix
        if (dst.check_pre_init()) {
            dst = matpoly(ab, m, n, size);
        }
        dst.set_size(size);

        // Leader creates a buffer matpoly of size length
        dst_partial = matpoly(ab, m, n, length);
        dst_partial.set_size(length);
    }

    size_t offset = 0;
    while (size > offset) {
        size_t len_raw = MIN(length, (size-offset));
        size_t len = simd * iceildiv(len_raw, simd);

        gather_mat_partial(dst_partial, 0, offset, len);

        // Copy the partial data into dst. This is the place where we
        // could write directly on disk if memory is a concern:
        if (!rank()) {
            for (unsigned int i = 0; i < dst.m; ++i) {
                for (unsigned int j = 0; j < dst.n; ++j) {
                    abdst_vec to = abvec_subvec(ab, dst.part(i, j), offset);
                    absrc_vec from = abvec_subvec_const(ab, dst_partial.part(i, j), 0);
                    abvec_set(ab, to, from, len);
                }
            }
        }
        offset += len;
    }
}

/* Exactly the converse of the previous function. */
void bigmatpoly::scatter_mat(matpoly const & src)
{
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
    matpoly src_partial;
#ifndef SELECT_MPFQ_LAYER_u64k1
    size_t length = 100;
#else
    size_t length = 1024;
#endif

    /* share allocation size. */
    struct {
        unsigned int m;
        unsigned int n;
        size_t size;
        size_t alloc;
    } shell { src.m, src.n, src.get_size(), src.capacity() };

    MPI_Bcast(&shell, sizeof(shell), MPI_BYTE, 0, com[0]);

    /* dst must be in pre-init mode */
    ASSERT_ALWAYS(check_pre_init());

    /* Allocate enough space on each node */
    finish_init(ab, shell.m, shell.n, shell.alloc);
    zero_pad(shell.size);

    if (!rank()) {
        // Leader creates a buffer matpoly of size length
        src_partial = matpoly(src.ab, src.m, src.n, length);
        src_partial.set_size(length);
    }

    size_t offset = 0;
    while (shell.size > offset) {
        size_t len_raw = MIN(length, (shell.size-offset));
        size_t len = simd * iceildiv(len_raw, simd);
        // Copy the partial data into src_partial. This is the place where we
        // could read directly from disk if memory is a concern:
        if (!rank()) {
            for (unsigned int i = 0; i < src.m; ++i) {
                for (unsigned int j = 0; j < src.n; ++j) {
                    abdst_vec to = abvec_subvec(ab, src_partial.part(i, j), 0);
                    absrc_vec from = abvec_subvec_const(ab, src.part(i, j), offset);
                    abvec_set(ab, to, from, len);
                }
            }
        }
        scatter_mat_partial(src_partial, 0, offset, len);
        offset += len;
    }
}

bigmatpoly bigmatpoly::truncate_and_rshift(unsigned int truncated_size, unsigned int shiftcount)
{
    bigmatpoly other(ab, get_model(), m, n, size - shiftcount);
    other.rshift(*this, shiftcount);
    truncate(*this, truncated_size);
    shrink_to_fit();
    std::swap(*this, other);
    return other;
}
