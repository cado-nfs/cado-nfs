#include "cado.h"
#include <cstdlib>
#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "lingen_matpoly_select.hpp"
#include "lingen_bigmatpoly.hpp"

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
        size_t offset, size_t length) const
{
#ifdef SELECT_MPFQ_LAYER_u64k1
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
    static_assert(std::is_same<absrc_vec, const unsigned long *>::value, "uh ?");
#endif

    /* sanity checks, because the code below assumes this. */
    ASSERT_ALWAYS(irank() * (int) n1 + jrank() == rank());
    ASSERT_ALWAYS(length <= (size_t) INT_MAX);

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
#ifndef SELECT_MPFQ_LAYER_u64k1
                            abvec_set(ab, to, from, length);
#else
                            std::copy(from, from + length / simd, to);
#endif
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
        size_t offset, size_t length)
{
#ifdef SELECT_MPFQ_LAYER_u64k1
    constexpr const unsigned int simd = matpoly::over_gf2 ? ULONG_BITS : 1;
    static_assert(std::is_same<absrc_vec, const unsigned long *>::value, "uh ?");
#endif

    MPI_Datatype mt;
#ifndef SELECT_MPFQ_LAYER_u64k1
    MPI_Type_contiguous(length * abvec_elt_stride(ab, 1), MPI_BYTE, &mt);
#else
    ASSERT_ALWAYS(length % simd == 0);
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
#ifndef SELECT_MPFQ_LAYER_u64k1
                            abvec_set(ab, to, from, length);
#else
                            std::copy(from, from + length / simd, to);
#endif
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
        size_t len = MIN(length, (size-offset));
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
    ASSERT_ALWAYS(ab);
    ASSERT_ALWAYS(!src.ab || ab == src.ab);

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
        size_t len = MIN(length, (shell.size-offset));
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
