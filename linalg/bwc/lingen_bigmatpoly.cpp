#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <climits>

#include <algorithm>
#include <array>
#include <memory>
#include <tuple>
#include <utility>

#include "lingen_bigmatpoly.hpp"
#include "lingen_matpoly_select.hpp"
#include "lingen_mul_substeps.hpp"
#include "lingen_substep_schedule.hpp"
#include "logline.hpp"
#include "macros.h"
#include "runtime_numeric_cast.hpp"
#include "select_mpi.h"
#include "submatrix_range.hpp"
#include "tree_stats.hpp"

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


/* {{{  init/zero/clear interface for bigmatpoly<is_binary> */

bigmatpoly_model::bigmatpoly_model(MPI_Comm * comm, unsigned int m, unsigned int n)
    : m1(m)
    , n1(n)
{
    memcpy(com, comm, 3 * sizeof(MPI_Comm));
}

/* This completes the initialization process. This is _not_ a collective
 * operation */
template<bool is_binary>
void bigmatpoly<is_binary>::finish_init(typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, int len)
{
    this->ab = ab;
    ASSERT_ALWAYS(check_pre_init());
    this->m = m;
    this->n = n;
    // m0 = subdivision(m, m1).nth_block_size(irank());;
    // n0 = subdivision(n, n1).nth_block_size(jrank());;
    m0 = subdivision(m, m1).block_size_upper_bound();
    n0 = subdivision(n, n1).block_size_upper_bound();
    size = runtime_numeric_cast<size_t>(len);
    my_cell() = matpoly<is_binary>(ab, m0, n0, len);
}

/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigmatpoly_finish_init
 */
template<bool is_binary>
bigmatpoly<is_binary>::bigmatpoly(bigmatpoly_model const & model)
    : bigmatpoly_model(model)
    , cells(m1 * n1)
{
}

template<bool is_binary>
bigmatpoly<is_binary>::bigmatpoly(typename matpoly<is_binary>::arith_hard * ab, bigmatpoly_model const & model, unsigned int m, unsigned int n, int len)
    : bigmatpoly_model(model)
    , ab(ab)
    , m(m)
    , n(n)
    , m0(subdivision(m, m1).block_size_upper_bound())
    , n0(subdivision(n, n1).block_size_upper_bound())
    , cells(m1 * n1)
{
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

template<bool is_binary>
bigmatpoly<is_binary>::bigmatpoly(bigmatpoly && a) noexcept
    : bigmatpoly_model(a.get_model())
    , ab(a.ab)
    , m(a.m), n(a.n)
    , m0(a.m0), n0(a.n0)
    , size(a.size)
    , cells(std::move(a.cells))
{
    a.m0=a.n0=a.m=a.n=a.size=0;
}
template<bool is_binary>
bigmatpoly<is_binary>& bigmatpoly<is_binary>::operator=(bigmatpoly&& a) noexcept
{
    get_model() = a.get_model();
    ab = a.ab;
    m = a.m;
    n = a.n;
    m0 = a.m0;
    n0 = a.n0;
    size=a.size;
    std::swap(cells, a.cells);
    a.m0=a.n0=a.m=a.n=a.size=0;
    return *this;
}
/* }}} */

#if 0
/* Return a bitmask indicating whether bigmatpoly_provision_{row,col} has
 * been called on this matrix before. bit 0 is for row, bit 1 is for col.
 * If the returned value is zero, then we really have only the local part
 * of this matrix for the moment.
 */
template<bool is_binary>
int bigmatpoly<is_binary>::provisioned() const
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
template<bool is_binary>
void bigmatpoly<is_binary>::provision_row()
{
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        auto & them = cell(irank(), j);
        if (them.check_pre_init())
            them = matpoly<is_binary>(ab, m0, n0, my_cell().capacity());
    }
}

#if 0
/* Rarely useful. We do need it because we resort to a kludgy
 * implementation of scatter_mat, which calls for provisioning on all
 * rows.
 */
template<bool is_binary>
void bigmatpoly<is_binary>::unprovision_row()
{
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        auto & them = cell(irank(), j);
        if (!them.check_pre_init())
            them = matpoly();
    }
}
#endif

/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
template<bool is_binary>
void bigmatpoly<is_binary>::provision_col() // bigmatpoly<is_binary> & p(*this)
{
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        auto & them = cell(i, jrank());
        if (them.check_pre_init())
            them = matpoly<is_binary>(ab, m0, n0, my_cell().capacity());
    }
}

/* Set size to be large enough to receive the given number of
 * coefficients. If space is provisioned for other cells in row or
 * column, allocation is triggered for them as well.
 */
template<bool is_binary>
void bigmatpoly<is_binary>::set_size(size_t nsize)
{
    auto & me = my_cell();
    ASSERT_ALWAYS(nsize <= me.capacity());
    size = nsize;
    me.set_size(nsize);
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        auto & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.set_size(nsize);
        ASSERT_ALWAYS(nsize <= them.capacity());
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        auto & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.set_size(nsize);
        ASSERT_ALWAYS(nsize <= them.capacity());
    }
}
template<bool is_binary>
void bigmatpoly<is_binary>::zero_pad(size_t nsize)/*{{{*/
{
    auto & me = my_cell();
    size = nsize;
    me.zero_pad(nsize);
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        auto & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.zero_pad(nsize);
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        auto & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.zero_pad(nsize);
    }
}

/* If our row or col cells have already been allocated, then reallocate
 * them as well (XXX is it clear or not ?) */
#if 0 /* This function has never been used or needed */
void bigmatpoly_realloc(bigmatpoly<is_binary> & p, int newalloc)
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
template<bool is_binary>
void bigmatpoly<is_binary>::zero()
{
    my_cell().zero();
}

/* okay, it's ugly */
#if 0
template<bool is_binary>
void bigmatpoly<is_binary>::swap(bigmatpoly<is_binary> & a)
{
    bigmatpoly x = std::move(*this);
    *this = std::move(a);
    a = std::move(x);
}
#endif
/* }}} */

template<bool is_binary>
void bigmatpoly<is_binary>::truncate(bigmatpoly const & src, unsigned int nsize)/*{{{*/
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
        auto const & sthem = src.cell(irank(), j);
        auto & dthem = cell(irank(), j);
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
        auto const & sthem = src.cell(i, jrank());
        auto & dthem = cell(i, jrank());
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

template<bool is_binary>
bool bigmatpoly<is_binary>::high_word_is_clear() const
{
    auto const & me = my_cell();
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

template<bool is_binary>
void bigmatpoly<is_binary>::clear_high_word()
{
    auto & me = my_cell();
    me.clear_high_word();
    for(unsigned int j = 0 ; j < n1 ; j++) {
        if (j == (unsigned int) jrank()) continue;
        auto & them = cell(irank(), j);
        if (them.check_pre_init()) continue;
        them.clear_high_word();
    }
    for(unsigned int i = 0 ; i < m1 ; i++) {
        if (i == (unsigned int) irank()) continue;
        auto & them = cell(i, jrank());
        if (them.check_pre_init()) continue;
        them.clear_high_word();
    }
}

template<bool is_binary>
void bigmatpoly<is_binary>::coeff_set_zero_loc(unsigned int k)
{
    my_cell().coeff_set_zero(k);
}

template<bool is_binary>
int bigmatpoly<is_binary>::coeff_is_zero(unsigned int k) const
{
    int t = my_cell().coeff_is_zero(k);
    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_INT, MPI_MIN, com[0]);
    return t;
}


template<bool is_binary>
void bigmatpoly<is_binary>::rshift(bigmatpoly & src, unsigned int k) /*{{{*/
{
    auto & me = my_cell();
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
template<bool is_binary>
void bigmatpoly<is_binary>::allgather_row()
{
    provision_row();
    for(unsigned int k = 0 ; k < n1 ; k++) {
        auto & data = cell(irank(), k);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long dsize = data.get_size();
        MPI_Bcast(&dsize, 1, MPI_UNSIGNED_LONG, k, com[1]);
        data.set_size(dsize);
        ASSERT_ALWAYS(data.get_size() <= data.capacity());
        ASSERT_ALWAYS((data.m * data.n * data.capacity()) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data.m * data.n, data.data_entry_alloc_size_in_bytes(),
                MPI_BYTE,
                MPI_Bcast(data.x, _datasize, _datatype, k, com[1])
        );
    }
}
template<bool is_binary>
void bigmatpoly<is_binary>::allgather_col()
{
    provision_col();
    for(unsigned int k = 0 ; k < m1 ; k++) {
        auto & data = cell(k, jrank());
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long dsize = data.get_size();
        MPI_Bcast(&dsize, 1, MPI_UNSIGNED_LONG, k, com[2]);
        data.set_size(dsize);
        ASSERT_ALWAYS(data.get_size() <= data.capacity());
        ASSERT_ALWAYS((data.m * data.n * data.capacity()) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data.m * data.n, data.data_entry_alloc_size_in_bytes(),
                MPI_BYTE,
                MPI_Bcast(data.x, _datasize, _datatype, k, com[2])
        );
    }
}
/* }}} */

template<bool is_binary>
bigmatpoly<is_binary> bigmatpoly<is_binary>::mul(bigmatpoly & a, bigmatpoly & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);
    ASSERT_ALWAYS(a.n == b.m);
    ASSERT_ALWAYS(a.n1 == b.m1);
    bigmatpoly c(a.ab, a.get_model(), a.m, b.n, csize);
    a.allgather_row();
    a.allgather_col();
    c.set_size(csize);
    auto & lc = c.my_cell();
    lc.zero();
    c.set_size(csize);
    ASSERT_ALWAYS(a.n == b.m);
    unsigned int const i = c.irank();
    unsigned int const j = c.jrank();
    for(unsigned int k = 0 ; k < a.n1 ; k++) {
        lc.addmul(a.cell(i, k), b.cell(k, j));
    }
    return c;
}/*}}}*/

template<bool is_binary>
bigmatpoly<is_binary> bigmatpoly<is_binary>::mp(bigmatpoly & a, bigmatpoly & c) /*{{{*/
{
    unsigned const bsize = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    ASSERT_ALWAYS(a.n1 == c.m1);
    bigmatpoly b(a.ab, a.get_model(), a.m, c.n, bsize);

    a.allgather_row();
    c.allgather_col();

    auto & lb = b.my_cell();
    lb.zero();
    b.set_size(bsize);
    ASSERT_ALWAYS(lb.get_size() == bsize);
    ASSERT_ALWAYS(lb.capacity() >= lb.get_size());
    unsigned int const i = c.irank();
    unsigned int const j = c.jrank();
    for(unsigned int k = 0 ; k < a.n1 ; k++) {
        lb.addmp(a.cell(i, k), c.cell(k, j));
    }
    return b;
}/*}}}*/

template<bool is_binary>
struct OP_CTX2 {
    bigmatpoly<is_binary> & c;
    bigmatpoly<is_binary> const & a;
    bigmatpoly<is_binary> const & b;
    MPI_Datatype mpi_entry_a;
    MPI_Datatype mpi_entry_b;
    tree_stats & stats;
    OP_CTX2(tree_stats & stats, bigmatpoly<is_binary> & c, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b)
        : c(c), a(a), b(b), stats(stats)
    {
        MPI_Type_contiguous(a.my_cell().data_entry_alloc_size_in_bytes(), MPI_BYTE, &mpi_entry_a);
        MPI_Type_contiguous(b.my_cell().data_entry_alloc_size_in_bytes(), MPI_BYTE, &mpi_entry_b);
        MPI_Type_commit(&mpi_entry_a);
        MPI_Type_commit(&mpi_entry_b);
    }
    OP_CTX2(OP_CTX2 const&) = delete;
    ~OP_CTX2() {
        MPI_Type_free(&mpi_entry_a);
        MPI_Type_free(&mpi_entry_b);
    }
    int a_irank() const { return a.irank(); }
    int b_irank() const { return b.irank(); }
    int a_jrank() const { return a.jrank(); }
    int b_jrank() const { return b.jrank(); }
    int mesh_inner_size() const { return a.n1; }
    static const bool uses_mpi = true;
    void mesh_checks() const {
        ASSERT_ALWAYS(a.get_model().is_square());
        ASSERT_ALWAYS(a.get_model() == b.get_model());
        ASSERT_ALWAYS(a.irank() == b.irank());
        ASSERT_ALWAYS(a.jrank() == b.jrank());
    }
    void alloc_c_if_needed(size_t size) const {
        if (c.m != a.m || c.n != a.n || runtime_numeric_cast<size_t>(c.get_size()) != size)
            c = bigmatpoly<is_binary>(a.ab, a.get_model(), a.m, b.n, size);
    }
    matpoly<is_binary> const & a_local() { return a.my_cell(); }
    matpoly<is_binary> const & b_local() { return b.my_cell(); }
    matpoly<is_binary> & c_local() { return c.my_cell(); }
    void a_allgather(void * p, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               p, n, mpi_entry_a, a.get_model().com[1]);
    }
    void b_allgather(void * p, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
               p, n, mpi_entry_b, b.get_model().com[2]);
    }
};
template<typename OP_T>
struct mp_or_mul : public OP_CTX2<OP_T::is_binary>
{ 
    static constexpr bool is_binary = OP_T::is_binary;
    OP_T & OP;
    const lingen_call_companion::mul_or_mp_times * M;
    subdivision mpi_split0;
    subdivision mpi_split1;
    subdivision mpi_split2;
    unsigned int nrs0, nrs2;
    unsigned int b0, b1, b2;
    matpoly<is_binary> a_peers;
    matpoly<is_binary> b_peers;
    using OP_CTX2<is_binary>::mesh_inner_size;
    using OP_CTX2<is_binary>::mesh_checks;
    using OP_CTX2<is_binary>::stats;
    using OP_CTX2<is_binary>::a_jrank;
    using OP_CTX2<is_binary>::b_irank;
    using OP_CTX2<is_binary>::b_jrank;
    using OP_CTX2<is_binary>::a;
    using OP_CTX2<is_binary>::b;
    using OP_CTX2<is_binary>::c;
    using OP_CTX2<is_binary>::a_local;
    using OP_CTX2<is_binary>::b_local;
    using OP_CTX2<is_binary>::c_local;
    using OP_CTX2<is_binary>::alloc_c_if_needed;
    using OP_CTX2<is_binary>::a_allgather;
    using OP_CTX2<is_binary>::b_allgather;

    /* Declare ta, tb, tc early on so that we don't malloc/free n times.  */
    mp_or_mul(tree_stats & stats, bigmatpoly<is_binary> & c, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b, OP_T & OP,
            const lingen_call_companion::mul_or_mp_times * M)
        : OP_CTX2<is_binary>(stats, c, a, b)
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
        plan_smallstep("gather_A", M->t_dft_A);
        plan_smallstep("gather_B", M->t_dft_B);
        plan_smallstep("addmul", M->t_conv);
        end_plan_smallstep();
    }
    template<typename... Args>
    void begin_smallstep(Args&& ...args) {
        if (M) stats.begin_smallstep(args...);
    }
    template<typename... Args>
    void skip_smallstep(Args&& ...args) {
        if (M) stats.skip_smallstep(args...);
    }
    template<typename... Args>
    void end_smallstep(Args&& ...args) {
        if (M) stats.end_smallstep(args...);
    }
    template<typename... Args>
    void plan_smallstep(Args&& ...args) {
        if (M) stats.plan_smallstep(args...);
    }
    template<typename... Args>
    void begin_plan_smallstep_microsteps(Args&& ...args) {
        if (M) stats.begin_plan_smallstep_microsteps(args...);
    }
    template<typename... Args>
    void begin_plan_smallstep(Args&& ...args) {
        if (M) stats.begin_plan_smallstep(args...);
    }
    template<typename... Args>
    void end_plan_smallstep(Args&& ...args) {
        if (M) stats.end_plan_smallstep(args...);
    }
    bool local_smallsteps_done(bool compulsory = false) {
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
        unsigned int const aj = a_jrank();
        unsigned int ak0mpi, ak1mpi;
        std::tie(ak0mpi, ak1mpi) = mpi_split1.nth_block(aj);

        unsigned int kk0,kk1;
        std::tie(kk0, kk1) = loop1.nth_block(iloop1);
        ASSERT_ALWAYS((kk1 - kk0) <= b1);
        /* XXX ak0 and co, and esp. ak1-ak0, are *NOT* identical across
         * mpi jobs. All that we have is ak1-ak0 <= b1 and bk1-bk0 <= b1.
         * In the non-mpi case, ak0==bk0 and ak1==bk1 */
        unsigned int const ak0 = ak0mpi + kk0;
        unsigned int const ak1 = std::min(ak1mpi, ak0 + b1);

        unsigned int ii0, ii1;
        std::tie(ii0, ii1) = loop0.nth_block(iloop0);

        submatrix_range const Ra  (i0 + ii0, ak0-ak0mpi, ii1 - ii0, ak1-ak0);
        submatrix_range const Rat (     0,   aj * b1,    ii1 - ii0, ak1-ak0);

        /* This is the analogue of the "dft_A" step in the transform
         * case.
         */
        a_peers.zero_with_size(a.get_size());  // for safety because of rounding.
        matpoly<is_binary>::copy(a_peers.view(Rat), a_local().view(Ra));

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
        begin_smallstep("gather_B", b1 * b2);
        unsigned int const bi = b_irank();
        unsigned int bk0mpi, bk1mpi;
        std::tie(bk0mpi, bk1mpi) = mpi_split1.nth_block(bi);

        unsigned int kk0,kk1;
        std::tie(kk0, kk1) = loop1.nth_block(iloop1);
        ASSERT_ALWAYS((kk1 - kk0) <= b1);
        /* XXX ak0 and co, and esp. ak1-ak0, are *NOT* identical across
         * mpi jobs. All that we have is ak1-ak0 <= b1 and bk1-bk0 <= b1.
         * In the non-mpi case, ak0==bk0 and ak1==bk1 */
        unsigned int const bk0 = bk0mpi + kk0;
        unsigned int const bk1 = std::min(bk1mpi, bk0 + b1);

        unsigned int jj0, jj1;
        std::tie(jj0, jj1) = loop2.nth_block(iloop2);

        submatrix_range const Rb   (bk0-bk0mpi, j0 + jj0, bk1-bk0, jj1 - jj0);
        submatrix_range const Rbt  (bi * b1,      0,      bk1-bk0, jj1 - jj0);

        b_peers.zero_with_size(b.get_size());
        matpoly<is_binary>::copy(b_peers.view(Rbt), b_local().view(Rb));

        // allgather tb among r nodes
        b_allgather(b_peers.part(0, 0), b1 * b2);
        end_smallstep();
    }/*}}}*/

    void addmul_for_block(typename matpoly<is_binary>::view_t & cdst, unsigned int iloop0, unsigned int iloop2)/*{{{*/
    {
        const unsigned int r = mesh_inner_size();
        unsigned int ii0, ii1;
        unsigned int jj0, jj1;
        std::tie(ii0, ii1) = loop0.nth_block(iloop0);
        std::tie(jj0, jj1) = loop2.nth_block(iloop2);

        begin_smallstep("addmul", b0 * b1 * b2 * r);

        // rounding might surprise us.
        submatrix_range const Ratxx(0,   0,   ii1 - ii0, r*b1);
        submatrix_range const Rbtxx(0,   0,   r*b1,      jj1 - jj0);
        submatrix_range const Rct  (cdst.i0 + ii0, cdst.j0 + jj0, ii1 - ii0, jj1 - jj0);

        typename matpoly<is_binary>::view_t const c_loc(cdst.M, Rct);
        typename matpoly<is_binary>::view_t const a_loc = a_peers.view(Ratxx);
        typename matpoly<is_binary>::view_t const b_loc = b_peers.view(Rbtxx);

        OP_T::addcompose(c_loc, a_loc, b_loc);

        end_smallstep();
    }/*}}}*/

    void operator()() {
        begin_smallstep(M->step_name());

        alloc_c_if_needed(OP.csize);
        c.zero_with_size(OP.csize);

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
        unsigned int const aj = a_jrank();
        unsigned int const bi = b_jrank();
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
        bool const process_blocks_row_major = b0 == nrs0;

        /* Now do a subblock */
        submatrix_range const Rc(i0, j0, i1-i0, j1-j0);
        typename matpoly<is_binary>::view_t cdst = c_local().view(Rc);
        cdst.zero();

        for(unsigned int iloop1 = 0 ; iloop1 < loop1.nblocks() ; iloop1++) {
            if (process_blocks_row_major) {
                ASSERT_ALWAYS(loop0.nblocks() == 1);
                ASSERT_ALWAYS(nrs0 == b0);
                unsigned int const iloop0 = 0;
                logline_printf(1, "gather_A (%u*%u)\n", b0, b1);
                gather_A(i0, iloop0, iloop1);
                for(unsigned int iloop2 = 0 ; iloop2 < loop2.nblocks() ; iloop2++) {
                    logline_printf(1, "gather_B (%u*%u)\n", b1, b2);
                    gather_B(j0, iloop1, iloop2);

                    logline_printf(1, "addmul\n");
                    addmul_for_block(cdst, iloop0, iloop2);
                }
                /* adjust counts */
                unsigned int const e2 = (iceildiv(nrs2, b2) - loop2.nblocks());
                unsigned int const x2 = e2 * b2;
                skip_smallstep("gather_B", b1 * x2);
                skip_smallstep("addmul", b0 * b1 * x2 * r);
            } else {
                ASSERT_ALWAYS(loop2.nblocks() == 1);
                ASSERT_ALWAYS(nrs2 == b2);
                unsigned int const iloop2 = 0;
                logline_printf(1, "gather_B (%u*%u)\n", b1, b2);
                gather_B(j0, iloop1, iloop2);
                for(unsigned int iloop0 = 0 ; iloop0 < loop0.nblocks() ; iloop0++) {
                    logline_printf(1, "gather_A (%u*%u)\n", b0, b1);
                    gather_A(i0, iloop0, iloop1);

                    logline_printf(1, "addmul\n");
                    addmul_for_block(cdst, iloop0, iloop2);
                }
                /* adjust counts */
                unsigned int const e0 = (iceildiv(nrs0, b0) - loop0.nblocks());
                unsigned int const x0 = e0 * b0;
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

template<bool is_binary>
bigmatpoly<is_binary> bigmatpoly<is_binary>::mp(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M)
{
    op_mp<is_binary, void> op(a, b, UINT_MAX);
    bigmatpoly c(a.get_model());
    mp_or_mul<op_mp<is_binary, void>>(stats, c, a, b, op, M)();
    return c;
}

template<bool is_binary>
bigmatpoly<is_binary> bigmatpoly<is_binary>::mul(tree_stats & stats, bigmatpoly const & a, bigmatpoly const & b, lingen_call_companion::mul_or_mp_times * M)
{
    op_mul<is_binary, void> op(a, b, UINT_MAX);
    bigmatpoly c(a.get_model());
    mp_or_mul<op_mul<is_binary, void>>(stats, c, a, b, op, M)();
    return c;
}

/*
 * gather to node 0, or scatter from node 0, but use "partial" transfer
 * in order to allow some kind of streaming between reading a matrix on
 * node 0 and scaterring it on all the nodes. (and the same for gather /
 * writing to file on node 0).
 */

#if 0
/* The piece [offset, offset+length[ of the bigmatpoly source is gathered
 * in the matpoly dst on node 0.
 * We assume that all the data structures are already set up properly,
 * and that dst has indeed room for length elements (we set the size, but
 * don't realloc dst).
 */
template<bool is_binary>
void bigmatpoly<is_binary>::gather_mat_partial(matpoly<is_binary> & dst,
        size_t dst_k,
        size_t offset, size_t length_raw) const
{
    constexpr const unsigned int simd = typename matpoly<is_binary>::over_gf2 ? ULONG_BITS : 1;

    choke me;
    /* the compilation must not vary based on LINGEN_BINARY */
#ifdef LINGEN_BINARY
    static_assert(std::is_same<absrc_vec, const unsigned long *>::value, "uh ?");
#endif

    /* sanity checks, because the code below assumes this. */
    ASSERT_ALWAYS(irank() * (int) n1 + jrank() == rank());
    ASSERT_ALWAYS(length_raw <= (size_t) INT_MAX);

    size_t length = simd * iceildiv(length_raw, simd);

    MPI_Datatype mt;
#ifndef LINGEN_BINARY
    MPI_Type_contiguous(length * ab->elt_stride(), MPI_BYTE, &mt);
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


/* Collect everything into node 0 */
template<bool is_binary>
void bigmatpoly<is_binary>::gather_mat(matpoly & dst) const
{
    constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
    matpoly<is_binary> dst_partial;
    size_t length = is_binary ? 1024 : 100;
    length = simd * iceildiv(length, simd);

    if (!rank()) {
        // Leader should initialize the result matrix
        if (dst.check_pre_init()) {
            dst = matpoly<is_binary>(ab, m, n, size);
        }
        dst.set_size(size);
    }

    for(size_t offset = 0 ; offset < size ; ) {
        size_t len = MIN(length, size - offset);
        gather_mat_partial(dst, offset, offset, len);
        offset += len;
    }
}
#endif

template<bool is_binary>
struct scatter_gather_base {/*{{{*/
    protected:
    /* These are just immediately derived from dst.get_model() and src,
     * but we keep copies because they're handy to have */
        typename matpoly<is_binary>::arith_hard * ab;

    size_t batch_length;

    bigmatpoly_model const & _model;

    unsigned int m1;
    unsigned int n1;

    struct shell_t {
        unsigned int m;
        unsigned int n;
        size_t size;
        shell_t(unsigned int m, unsigned int n, size_t size, bigmatpoly_model const & model)
            : m(m), n(n), size(size)
        {
            MPI_Bcast(this, sizeof(*this), MPI_BYTE, 0, model.com[0]);
        }
    } shell;

    subdivision R;       /* Row split */
    subdivision C;       /* Col split */
    unsigned int m0 = R.block_size_upper_bound();
    unsigned int n0 = C.block_size_upper_bound();

    MPI_Datatype mt;

    std::unique_ptr<MPI_Request[]> reqs;

    bigmatpoly_model const & get_model() const { return _model; }
    int rank() const { return get_model().rank(); }

    /* This structure is initialized correctly at root first, and then
     * the ctor bcasts the content to other nodes. This is the reason why
     * we have to put the initialization very early on in the base class.
     * In truth, it is only important for scatter_mat.
     */
    scatter_gather_base(
            typename matpoly<is_binary>::arith_hard * ab,
            bigmatpoly_model const & model,
            unsigned int m, unsigned int n, size_t size)
        /*{{{*/
        : ab(ab)
        , batch_length(roundup_simd(std::min(size_t(is_binary ? 1024 : 100), size)))
        , _model(model)
        , m1(model.m1)
        , n1(model.n1)
        , shell(m, n, size, model)
        , R(shell.m, m1)
        , C(shell.n, n1)
        , m0(R.block_size_upper_bound())
        , n0(C.block_size_upper_bound())
        , reqs(new MPI_Request[rank() ? 1 : (m1 * n1)])
    {
        /* The batch length should normally not be an issue, as we expect
         * that the call thresholds will be such that only smaller sizes will
         * be used.
         */
        MPI_Bcast(&batch_length, sizeof(batch_length), MPI_BYTE, 0, _model.com[0]);

        define_mpi_type();

    }/*}}}*/

    static size_t roundup_simd(size_t x) {
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        return simd * iceildiv(x, simd);
    }

    void define_mpi_type()/*{{{*/
    {
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        if constexpr (is_binary)
            MPI_Type_contiguous(batch_length / simd, MPI_UNSIGNED_LONG, &mt);
        else
            MPI_Type_contiguous(batch_length * ab->elt_stride(), MPI_BYTE, &mt);
        MPI_Type_commit(&mt);
    }/*}}}*/

    void undefine_mpi_type()/*{{{*/
    {
        MPI_Type_free(&mt);
    }/*}}}*/

    ~scatter_gather_base()/*{{{*/
    {
        undefine_mpi_type();
    }/*}}}*/

    /* use_intermediary_tight=true is especially important in the binary
     * case, in order to avoid having very many small sends.
     *
     * use_intermediary_tight=false is the old code that used to work
     * for the prime field case. At least we hope that the two
     * options are equivalent.
     *
     * It might be that use_intermediary_tight=true is The Right
     * Thing anyway. We need to find a way to elide copies in the
     * important special case where src or dst is already tight.
     * Shouldn't be too hard.
     */
    static constexpr const bool use_intermediary_tight = true;

    /* blocking or non-blocking doesn't seem to make a whole lot of
     * difference.
     */
    static constexpr const bool use_nonblocking = false;
};/*}}}*/

template<bool is_binary>
class gather_mat : public scatter_gather_base<is_binary> {/*{{{*/
    using matpoly_type = matpoly<is_binary>;
    using bigmatpoly_type = bigmatpoly<is_binary>;
    bigmatpoly_type const & src;
    matpoly_type & dst;

    /* When use_intermediary_tight = trye, we want to communicate *full*
     * cells from the gathered source matrix to the scattered one.
     * Unfortunately we cannot guarantee that either has tight
     * allocation. Furthermore, the layout in the source matrix does not
     * have contiguous block mathcing the peer sub-blocks. For these two
     * reasons (stride and layout), we need copies at both ends.
     *
     * We have therefore an intermediary source and an intermediary
     * destination, both filled at each loop iteration.
     *
     */
    bigmatpoly_type src_partial;
    matpoly_type dst_partial;

    using scatter_gather_base<is_binary>::get_model;
    using scatter_gather_base<is_binary>::ab;
    using scatter_gather_base<is_binary>::shell;
    using scatter_gather_base<is_binary>::use_intermediary_tight;
    using scatter_gather_base<is_binary>::rank;
    using scatter_gather_base<is_binary>::batch_length;
    using scatter_gather_base<is_binary>::use_nonblocking;
    using scatter_gather_base<is_binary>::mt;
    using scatter_gather_base<is_binary>::reqs;
    using scatter_gather_base<is_binary>::R;
    using scatter_gather_base<is_binary>::C;
    using scatter_gather_base<is_binary>::undefine_mpi_type;
    using scatter_gather_base<is_binary>::define_mpi_type;
    using scatter_gather_base<is_binary>::n0;
    using scatter_gather_base<is_binary>::n1;
    using scatter_gather_base<is_binary>::m1;
    using scatter_gather_base<is_binary>::roundup_simd;

    public:
    gather_mat(matpoly_type & dst, bigmatpoly_type const & src)/*{{{*/
        : scatter_gather_base<is_binary>(src.ab, src.get_model(), src.m, src.n, src.get_size())
       , src(src)
       , dst(dst)
       , src_partial(get_model())
    {
        if (use_intermediary_tight) {
            /* the temporary source and destination */
            if (!rank()) {
                dst_partial = matpoly_type(ab, shell.m, shell.n, batch_length);
                dst_partial.zero_with_size(batch_length);
                ASSERT_ALWAYS(dst_partial.is_tight());
            }
            src_partial = bigmatpoly_type(ab, src, shell.m, shell.n, batch_length);
            src_partial.zero_with_size(batch_length);
            ASSERT_ALWAYS(src_partial.my_cell().is_tight());
        }
    }/*}}}*/

    private:

    /* {{{ building blocks for the old strategy, where we receive one
     * polynomial at a time, avoiding some copies */
    void mini_recv(size_t src_offset, size_t dst_offset, unsigned int i0, unsigned int i1, unsigned int j0, unsigned int j1, MPI_Request * & req)/*{{{*/
    {

        ASSERT_ALWAYS(!use_intermediary_tight);
        unsigned int const ii = R.flatten(i1, i0);
        unsigned int const jj = C.flatten(j1, j0);
        unsigned int const peer = i1 * n1 + j1;
        unsigned int const tag = ii * shell.n + jj;
        auto to = ab->vec_subvec(dst.part(ii, jj), dst_offset);
        auto from = ab->vec_subvec(src.my_cell().part(i0, j0), src_offset);

        /* XXX There's a subtlety here. batch_length and mt are tinkered
         * with by the main_loop code for the last iteration, so that the
         * final write doesn't overflow. Yes it's kludgy.
         */

        if (peer == 0) {
            /* talk to ourself */
            ab->vec_set(to, from, roundup_simd(batch_length));
        } else {
            if (use_nonblocking) {
                MPI_Irecv((void*) to, 1, mt, peer, tag, get_model().com[0], req);
            } else {
                MPI_Recv((void*) to, 1, mt, peer, tag, get_model().com[0], MPI_STATUS_IGNORE);
            }
        }
        req++;
    }/*}}}*/
    void post_mini_recvs(size_t src_offset, size_t dst_offset) {/*{{{*/
        ASSERT_ALWAYS(!use_intermediary_tight);
        MPI_Request * req = reqs.get();
        /* the master sends data to everyone */
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        mini_recv(src_offset, dst_offset, i0, i1, j0, j1, req);
                    }
                }
            }
        }
    }/*}}}*/
    void wait_mini_recvs()/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        if (!use_nonblocking) return;
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int const peer = i1 * n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
    }/*}}}*/
    void post_mini_sends(size_t src_offset)/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        ASSERT_ALWAYS(rank());
        MPI_Request * req = reqs.get();
        for(unsigned int i0 = 0 ; i0 < src.m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < src.n0r() ; j0++) {
                unsigned int const ii = R.flatten(src.irank(), i0);
                unsigned int const jj = C.flatten(src.jrank(), j0);
                unsigned int const tag = ii * shell.n + jj;
                auto from = ab->vec_subvec(src.my_cell().part(i0, j0), src_offset);
                /* battle const-deprived MPI prototypes... */
                if (use_nonblocking) {
                    MPI_Isend((void*) from, 1, mt, 0, tag, get_model().com[0], req);
                } else {
                    MPI_Send((void*) from, 1, mt, 0, tag, get_model().com[0]);
                }
                req++;
            }
        }
    }/*}}}*/
    void wait_mini_sends()/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        if (!use_nonblocking) return;
        MPI_Request * req = reqs.get();
        for(unsigned int i0 = 0 ; i0 < src.m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < src.n0r() ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
    }/*}}}*/
    /* }}} */

    /* {{{ This alternative strategy, with use_intermediary_tight=true, adds
     * extra copies, but uses much larger messages, which is a clear win
     * in some cases (always ?) */
    void copy_tight_to_dst(size_t dst_offset, size_t len, unsigned int i1, unsigned int j1)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        /* dst_partial has the same size as dst, but a different layout,
         * as blocks for all peers are stored contiguously. Our first
         * task is to determine exactly where the block for the current
         * peer starts.
         */
        unsigned int v = 0;
        /* all peers with row index are stored before us */
        v += R.nth_block_start(i1) * shell.n;
        /* and then peers with the same row index, but lower col index */
        v += R.nth_block_size(i1) * C.nth_block_start(j1);
        for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
            for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                unsigned int const ii = R.flatten(i1, i0);
                unsigned int const jj = C.flatten(j1, j0);
                auto to = ab->vec_subvec(dst.part(ii, jj), dst_offset);
                unsigned int const w = v + i0 * C.nth_block_size(j1) + j0;
                unsigned int const wi = w / shell.n;
                unsigned int const wj = w % shell.n;
                auto from = dst_partial.part(wi, wj);
                ab->vec_set(to, from, roundup_simd(len));
            }
        }
    }/*}}}*/

    void post_block_recv(unsigned int i1, unsigned int j1, MPI_Request * & req)
    {
        unsigned int const peer = i1 * n1 + j1;
        typename matpoly_type::ptr rank0_to;

        {
            unsigned int v = 0;
            v += R.nth_block_start(i1) * shell.n;
            v += R.nth_block_size(i1) * C.nth_block_start(j1);
            unsigned int const vi = v / shell.n;
            unsigned int const vj = v % shell.n;
            rank0_to = dst_partial.part(vi, vj);
        }
        unsigned int const blocksize = R.nth_block_size(i1) * C.nth_block_size(j1);

        auto peer_from = src_partial.my_cell().part(0,0);

        unsigned int const tag = peer;

        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        ASSERT_ALWAYS(batch_length % simd == 0);
        if (peer == 0) {
            /* talk to ourself */
            ab->vec_set(rank0_to, peer_from, blocksize*batch_length);
        } else {
            /* battle const-deprived MPI prototypes... */
            if (use_nonblocking) {
                MPI_Irecv((void*) rank0_to, blocksize, mt, peer, tag, get_model().com[0], req);
            } else {
                MPI_Recv((void*) rank0_to, blocksize, mt, peer, tag, get_model().com[0], MPI_STATUS_IGNORE);
            }
        }
        req++;
    }

    void post_recvs()/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < m1 ; i1++) {
            for(unsigned int j1 = 0 ; j1 < n1 ; j1++) {
                post_block_recv(i1, j1, req);
            }
        }
    }/*}}}*/
    void wait_recvs_and_copy_tight_to_dst(size_t dst_offset, size_t len)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < m1 ; i1++) {
            for(unsigned int j1 = 0 ; j1 < n1 ; j1++) {
                unsigned int const peer = i1 * n1 + j1;
                if (use_nonblocking) {
                    if (peer) MPI_Wait(req, MPI_STATUS_IGNORE);
                    req++;
                }
                copy_tight_to_dst(dst_offset, len, i1, j1);
            }
        }
    }/*}}}*/
    void post_sends()/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        ASSERT_ALWAYS (rank());
        auto peer_to = src_partial.my_cell().part(0,0);
        unsigned int const tag = rank();
        unsigned int const blocksize = src.m0r() * src.n0r();
        if (use_nonblocking) {
            MPI_Request req[1];
            MPI_Isend((void*) peer_to, blocksize, mt, 0, tag, get_model().com[0], req);
            MPI_Wait(req, MPI_STATUS_IGNORE);
        } else {
            MPI_Send((void*) peer_to, blocksize, mt, 0, tag, get_model().com[0]);
        }
    }/*}}}*/
    void copy_src_to_tight(size_t src_offset, size_t len)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        /* XXX in dst_partial, the data is stored contiguously ! */
        for (unsigned int i = 0; i < src.m0r(); ++i) {
            for (unsigned int j = 0; j < src.n0r(); ++j) {
                unsigned int const v = i * src.n0r() + j;
                unsigned int const vi = v / n0;
                unsigned int const vj = v % n0;
                auto to = ab->vec_subvec(src_partial.my_cell().part(vi, vj), 0);
                auto from = ab->vec_subvec(src.my_cell().part(i, j), src_offset);
                ab->vec_set(to, from, roundup_simd(len));
            }
        }
    }/*}}}*/
    /* }}} */

    /* {{{ wrap it up */
    void partial(size_t src_offset, size_t dst_offset, size_t len)/*{{{*/
    {
        if (!use_intermediary_tight && len < batch_length) {
            /* XXX Yes this is an ugly hack.  Bear in mind that all
             * transfers, always, have length batch_length. Since the
             * mini- strategy does not have intermediary storage at
             * both ends, it's difficult to avoid overflow. Hence the
             * "len" argument can only be used safely for the other
             * strategy. Here we have to resort to ugly stuff.
             */
            batch_length = roundup_simd(len);
            undefine_mpi_type();
            define_mpi_type();
        }

        ASSERT_ALWAYS(src_offset + len <= src.get_size());
        ASSERT_ALWAYS(dst_offset + len <= dst.get_size() || rank());

        if (!rank()) dst_partial.zero_with_size(batch_length);
        src_partial.zero_with_size(batch_length);

        /* post Isends (everywhere), and Irecvs (at root) */

        if (use_intermediary_tight) {
            ASSERT_ALWAYS(len <= dst_partial.get_size() || rank());
            ASSERT_ALWAYS(len <= src_partial.get_size());
            copy_src_to_tight(src_offset, len);
            if (!rank()) {
                post_recvs();
                wait_recvs_and_copy_tight_to_dst(dst_offset, len);
            } else {
                post_sends();
            }
        } else {
            ASSERT_ALWAYS(len == batch_length);
            if (!rank()) {
                post_mini_recvs(src_offset, dst_offset);
                wait_mini_recvs();
            } else {
                post_mini_sends(src_offset);
                wait_mini_sends();
            }
        }
    }/*}}}*/

    public:

    void main_loop_shifted(size_t src_offset, size_t dst_offset, size_t len)/*{{{*/
    {
        for( ; len >= batch_length ; ) {
            partial(src_offset, dst_offset, batch_length);
            src_offset += batch_length;
            dst_offset += batch_length;
            len -= batch_length;
        }

        if (!len) return;

        partial(src_offset, dst_offset, len);
    }/*}}}*/

    void main_loop()/*{{{*/
    {
        if (!rank()) {
            dst = matpoly_type(ab, shell.m, shell.n, shell.size);
            dst.zero_with_size(shell.size);
        }

        /* This one is quite simple, we have src_offset == dst_offset
         * throughout.
         */
        for(size_t offset = 0 ; offset < shell.size ; ) {
            size_t const len = MIN(batch_length, shell.size - offset);

            partial(offset, offset, len);

            offset += len;
        }
    }/*}}}*/
    /* }}} */
};/*}}}*/

template<bool is_binary>
class scatter_mat : public scatter_gather_base<is_binary> {/*{{{*/
    using matpoly_type = matpoly<is_binary>;
    using bigmatpoly_type = bigmatpoly<is_binary>;
    matpoly_type const & src;
    bigmatpoly_type & dst;

    /* When use_intermediary_tight = trye, we want to communicate *full*
     * cells from the gathered source matrix to the scattered one.
     * Unfortunately we cannot guarantee that either has tight
     * allocation. Furthermore, the layout in the source matrix does not
     * have contiguous block matching the peer sub-blocks. For these two
     * reasons (stride and layout), we need copies at both ends.
     *
     * We have therefore an intermediary source and an intermediary
     * destination, both filled at each loop iteration.
     *
     */
    matpoly_type src_partial;
    bigmatpoly_type dst_partial;

    using scatter_gather_base<is_binary>::get_model;
    using scatter_gather_base<is_binary>::ab;
    using scatter_gather_base<is_binary>::shell;
    using scatter_gather_base<is_binary>::use_intermediary_tight;
    using scatter_gather_base<is_binary>::rank;
    using scatter_gather_base<is_binary>::batch_length;
    using scatter_gather_base<is_binary>::use_nonblocking;
    using scatter_gather_base<is_binary>::mt;
    using scatter_gather_base<is_binary>::reqs;
    using scatter_gather_base<is_binary>::R;
    using scatter_gather_base<is_binary>::C;
    using scatter_gather_base<is_binary>::undefine_mpi_type;
    using scatter_gather_base<is_binary>::define_mpi_type;
    using scatter_gather_base<is_binary>::n0;
    using scatter_gather_base<is_binary>::n1;
    using scatter_gather_base<is_binary>::m1;
    using scatter_gather_base<is_binary>::roundup_simd;

    public:
    scatter_mat(bigmatpoly_type & dst, matpoly_type const & src)/*{{{*/
        /* NOTE that src.ab is not relevant everywhere. We must rely on
         * dst.ab, whence we must be sure that dst.ab != NULL */
        : scatter_gather_base<is_binary>(dst.ab, dst.get_model(), src.m, src.n, src.get_size())
       , src(src)
       , dst(dst)
       , dst_partial(get_model())
    {
        ASSERT_ALWAYS(is_binary || dst.ab != nullptr);
        if (use_intermediary_tight) {
            /* the temporary source and destination */
            if (!rank()) {
                src_partial = matpoly_type(ab, shell.m, shell.n, batch_length);
                src_partial.zero_with_size(batch_length);
                ASSERT_ALWAYS(src_partial.is_tight());
            }
            dst_partial = bigmatpoly_type(ab, dst, shell.m, shell.n, batch_length);
            dst_partial.zero_with_size(batch_length);
            ASSERT_ALWAYS(dst_partial.my_cell().is_tight());
        }

    }/*}}}*/

    private:

    /* {{{ building blocks for the old strategy, where we send one polynomial
     * at a time, avoiding some copies */
    void mini_send(size_t src_offset, size_t dst_offset, unsigned int i0, unsigned int i1, unsigned int j0, unsigned int j1, MPI_Request * & req)/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        unsigned int const ii = R.flatten(i1, i0);
        unsigned int const jj = C.flatten(j1, j0);
        unsigned int const peer = i1 * n1 + j1;
        unsigned int const tag = ii * shell.n + jj;
        auto from = ab->vec_subvec(src.part(ii, jj), src_offset);
        auto to = ab->vec_subvec(dst.my_cell().part(i0, j0), dst_offset);

        /* XXX There's a subtlety here. batch_length and mt are tinkered
         * with by the main_loop code for the last iteration, so that the
         * final write doesn't overflow. Yes it's kludgy.
         */

        if (peer == 0) {
            /* talk to ourself */
            ab->vec_set(to, from, roundup_simd(batch_length));
        } else {
            /* battle const-deprived MPI prototypes... */
            if (use_nonblocking) {
                MPI_Isend((void*) from, 1, mt, peer, tag, get_model().com[0], req);
            } else {
                MPI_Send((void*) from, 1, mt, peer, tag, get_model().com[0]);
            }
        }
        req++;
    }/*}}}*/
    void post_mini_sends(size_t src_offset, size_t dst_offset) {/*{{{*/
        ASSERT_ALWAYS(!use_intermediary_tight);
        MPI_Request * req = reqs.get();
        /* the master sends data to everyone */
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        mini_send(src_offset, dst_offset, i0, i1, j0, j1, req);
                    }
                }
            }
        }
    }/*}}}*/
    void wait_mini_sends()/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        if (!use_nonblocking) return;
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < R.nblocks() ; i1++) {
            for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
                for(unsigned int j1 = 0 ; j1 < C.nblocks() ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                        unsigned int const peer = i1 * n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
    }/*}}}*/
    void post_mini_recvs(size_t dst_offset)/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        ASSERT_ALWAYS(rank());
        MPI_Request * req = reqs.get();
        for(unsigned int i0 = 0 ; i0 < dst.m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < dst.n0r() ; j0++) {
                unsigned int const ii = R.flatten(dst.irank(), i0);
                unsigned int const jj = C.flatten(dst.jrank(), j0);
                unsigned int const tag = ii * shell.n + jj;
                auto to = ab->vec_subvec(dst.my_cell().part(i0, j0), dst_offset);
                if (use_nonblocking) {
                    MPI_Irecv(to, 1, mt, 0, tag, get_model().com[0], req);
                } else {
                    MPI_Recv(to, 1, mt, 0, tag, get_model().com[0], MPI_STATUS_IGNORE);
                }
                req++;
            }
        }
    }/*}}}*/
    void wait_mini_recvs()/*{{{*/
    {
        ASSERT_ALWAYS(!use_intermediary_tight);
        if (!use_nonblocking) return;
        MPI_Request * req = reqs.get();
        for(unsigned int i0 = 0 ; i0 < dst.m0r() ; i0++) {
            for(unsigned int j0 = 0 ; j0 < dst.n0r() ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
    }/*}}}*/
    /* }}} */

    /* {{{ This alternative strategy, with use_intermediary_tight=true, adds
     * extra copies, but uses much larger messages, which is a clear win
     * in some cases (always ?) */
    void copy_src_to_tight(size_t src_offset, size_t len, unsigned int i1, unsigned int j1)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        /* src_partial has the same size as src, but a different layout,
         * as blocks for all peers are stored contiguously. Our first
         * task is to determine exactly where the block for the current
         * peer starts.
         */
        unsigned int v = 0;
        /* all peers with row index are stored before us */
        v += R.nth_block_start(i1) * shell.n;
        /* and then peers with the same row index, but lower col index */
        v += R.nth_block_size(i1) * C.nth_block_start(j1);
        for(unsigned int i0 = 0 ; i0 < R.nth_block_size(i1) ; i0++) {
            for(unsigned int j0 = 0 ; j0 < C.nth_block_size(j1) ; j0++) {
                unsigned int const ii = R.flatten(i1, i0);
                unsigned int const jj = C.flatten(j1, j0);
                auto from = ab->vec_subvec(src.part(ii, jj), src_offset);
                unsigned int const w = v + i0 * C.nth_block_size(j1) + j0;
                unsigned int const wi = w / shell.n;
                unsigned int const wj = w % shell.n;
                auto to = src_partial.part(wi, wj);
                ab->vec_set(to, from, roundup_simd(len));
            }
        }
    }/*}}}*/

    void post_block_send(unsigned int i1, unsigned int j1, MPI_Request * & req)
    {
        unsigned int const peer = i1 * n1 + j1;
        typename matpoly_type::srcptr rank0_from;

        {
            unsigned int v = 0;
            v += R.nth_block_start(i1) * shell.n;
            v += R.nth_block_size(i1) * C.nth_block_start(j1);
            unsigned int const vi = v / shell.n;
            unsigned int const vj = v % shell.n;
            rank0_from = src_partial.part(vi, vj);
        }
        unsigned int const blocksize = R.nth_block_size(i1) * C.nth_block_size(j1);

        auto peer_to = dst_partial.my_cell().part(0,0);

        /* We can send the data now */
        unsigned int const tag = peer;

        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        ASSERT_ALWAYS(batch_length % simd == 0);
        if (peer == 0) {
            /* talk to ourself */
            ab->vec_set(peer_to, rank0_from, blocksize*batch_length);
        } else {
            /* battle const-deprived MPI prototypes... */
            if (use_nonblocking) {
                MPI_Isend((void*) rank0_from, blocksize, mt, peer, tag, get_model().com[0], req);
            } else {
                MPI_Send((void*) rank0_from, blocksize, mt, peer, tag, get_model().com[0]);
            }
        }
        req++;

    }
    void copy_src_to_tight_and_send(size_t src_offset, size_t len)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < m1 ; i1++) {
            for(unsigned int j1 = 0 ; j1 < n1 ; j1++) {
                copy_src_to_tight(src_offset, len, i1, j1);
                post_block_send(i1, j1, req);
            }
        }
    }/*}}}*/
    void wait_sends()/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        if (!use_nonblocking) return;
        MPI_Request * req = reqs.get();
        for(unsigned int i1 = 0 ; i1 < m1 ; i1++) {
            for(unsigned int j1 = 0 ; j1 < n1 ; j1++) {
                unsigned int const peer = i1 * n1 + j1;
                if (peer) MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
    }/*}}}*/
    void post_recvs()/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        ASSERT_ALWAYS (rank());
        auto peer_to = dst_partial.my_cell().part(0,0);
        unsigned int const tag = rank();
        unsigned int const blocksize = dst.m0r() * dst.n0r();
        if (use_nonblocking) {
            MPI_Request req[1];
            MPI_Irecv(peer_to, blocksize, mt, 0, tag, get_model().com[0], req);
            MPI_Wait(req, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(peer_to, blocksize, mt, 0, tag, get_model().com[0], MPI_STATUS_IGNORE);
        }
    }/*}}}*/
    void copy_tight_to_dst(size_t dst_offset, size_t len)/*{{{*/
    {
        ASSERT_ALWAYS(use_intermediary_tight);
        /* XXX in dst_partial, the data is stored contiguously ! */
        for (unsigned int i = 0; i < dst.m0r(); ++i) {
            for (unsigned int j = 0; j < dst.n0r(); ++j) {
                unsigned int const v = i * dst.n0r() + j;
                unsigned int const vi = v / n0;
                unsigned int const vj = v % n0;
                auto from = ab->vec_subvec(dst_partial.my_cell().part(vi, vj), 0);
                auto to = ab->vec_subvec(dst.my_cell().part(i, j), dst_offset);
                ab->vec_set(to, from, roundup_simd(len));
            }
        }
    }/*}}}*/
    /* }}} */

    public:

    /* {{{ wrap it up */
    /* 
     * Take length element in the src matrix on node 0, and scatter it
     * in dst, with the given src and dst offsets.
     * We assume that dst has been initialized: all the communicators, mn's,
     * are already set and enough space to accomodate length+dst_offset elements
     * have been already allocated.
     *
     * dst.size is not modified by this function, we assume that the
     * caller has called zero_with_size() beforehand (for example).
     */
    void partial(size_t src_offset, size_t dst_offset, size_t len)/*{{{*/
    {
        if (!use_intermediary_tight && len < batch_length) {
            /* XXX Yes this is an ugly hack.  Bear in mind that all
             * transfers, always, have length batch_length. Since the
             * mini- strategy does not have intermediary storage at
             * both ends, it's difficult to avoid overflow. Hence the
             * "len" argument can only be used safely for the other
             * strategy. Here we have to resort to ugly stuff.
             */
            batch_length = roundup_simd(len);
            undefine_mpi_type();
            define_mpi_type();
        }

        ASSERT_ALWAYS(src_offset + len <= src.get_size() || rank());
        ASSERT_ALWAYS(dst_offset + len <= dst.get_size());

        dst_partial.zero_with_size(batch_length);
        if (!rank()) src_partial.zero_with_size(batch_length);

        /* post Isends (at root), and Irecvs (elsewhere) */
        if (use_intermediary_tight) {
            ASSERT_ALWAYS(len <= src_partial.get_size() || rank());
            ASSERT_ALWAYS(len <= dst_partial.get_size());
            if (!rank()) {
                copy_src_to_tight_and_send(src_offset, len);
                wait_sends();
            } else {
                post_recvs();
            }
            copy_tight_to_dst(dst_offset, len);
        } else {
            ASSERT_ALWAYS(len == batch_length);
            if (!rank()) {
                post_mini_sends(src_offset, dst_offset);
                wait_mini_sends();
            } else {
                post_mini_recvs(dst_offset);
                wait_mini_recvs();
            }
        }
    }/*}}}*/

    void main_loop_shifted(size_t src_offset, size_t dst_offset, size_t len)/*{{{*/
    {
        for( ; len >= batch_length ; ) {
            partial(src_offset, dst_offset, batch_length);
            src_offset += batch_length;
            dst_offset += batch_length;
            len -= batch_length;
        }

        if (!len) return;

        partial(src_offset, dst_offset, len);
    }/*}}}*/

    void main_loop()/*{{{*/
    {
        /* dst must be in pre-init mode */
        ASSERT_ALWAYS(dst.check_pre_init());

        /* Allocate enough space on each node */
        dst.finish_init(ab, shell.m, shell.n, shell.size);
        dst.zero_with_size(shell.size);

        /* This one is quite simple, we have src_offset == dst_offset
         * throughout.
         */
        for(size_t offset = 0 ; offset < shell.size ; ) {
            size_t const len = MIN(batch_length, shell.size - offset);

            partial(offset, offset, len);

            offset += len;
        }
    }/*}}}*/
    /* }}} */
};/*}}}*/

template<bool is_binary>
void bigmatpoly<is_binary>::scatter_mat_partial(matpoly<is_binary> const & src, size_t src_offset, size_t dst_offset, size_t length)
{
    ::scatter_mat<is_binary>(*this, src).main_loop_shifted(src_offset, dst_offset, length);
}


/* Exactly the converse of the previous function. */
template<bool is_binary>
void bigmatpoly<is_binary>::scatter_mat(matpoly<is_binary> const & src)
{
    ::scatter_mat<is_binary>(*this, src).main_loop();
}

template<bool is_binary>
void bigmatpoly<is_binary>::gather_mat_partial(matpoly<is_binary> & dst, size_t dst_offset, size_t src_offset, size_t length) const
{
    ::gather_mat<is_binary>(dst, *this).main_loop_shifted(src_offset, dst_offset, length);
}


/* Exactly the converse of the previous function. */
template<bool is_binary>
void bigmatpoly<is_binary>::gather_mat(matpoly<is_binary> & dst) const
{
    ::gather_mat<is_binary>(dst, *this).main_loop();
}

template<bool is_binary>
bigmatpoly<is_binary> bigmatpoly<is_binary>::truncate_and_rshift(unsigned int truncated_size, unsigned int shiftcount)
{
    bigmatpoly other(ab, get_model(), m, n, size - shiftcount);
    other.rshift(*this, shiftcount);
    truncate(*this, truncated_size);
    shrink_to_fit();
    std::swap(*this, other);
    return other;
}

#ifdef LINGEN_BINARY
template class bigmatpoly<true>;
#else
template class bigmatpoly<false>;
#endif
