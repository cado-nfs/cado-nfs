#include "cado.h" // IWYU pragma: keep
#include <limits.h>          // for INT_MAX
#include <stdlib.h>          // for free, malloc, abort
#include <string.h>          // for memcpy, memset
#include "lingen-bigmatpoly.h"
#include "lingen-matpoly.h"  // for matpoly_s, matpoly_part, matpoly_init
#include "macros.h"          // for ASSERT_ALWAYS
#include "mpfq_layer.h"      // for abdst_field, mpfq_pz_dst_field, mpfq_pz_...
#include "timing.h"          // for wct_seconds


/* {{{  init/zero/clear interface for bigmatpoly */
matpoly_ptr bigmatpoly_my_cell(bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    return bigmatpoly_cell(p, irank, jrank);
}

matpoly_srcptr bigmatpoly_my_cell_const(bigmatpoly_srcptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    return bigmatpoly_cell_const(p, irank, jrank);
}

void bigmatpoly_init_model(bigmatpoly_ptr model, MPI_Comm * comm, unsigned int m, unsigned int n)
{
    memset(model, 0, sizeof(bigmatpoly));
    model->m1 = m;
    model->n1 = n;
    memcpy(model->com, comm, 3 * sizeof(MPI_Comm));
}

/* This completes the initialization process. This is _not_ a collective
 * operation */
void bigmatpoly_finish_init(abdst_field ab, bigmatpoly_ptr p, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(bigmatpoly_check_pre_init(p));
    p->m = m;
    p->n = n;
    p->m0 = m / p->m1;
    p->n0 = n / p->n1;
    p->size = len;
    matpoly_ptr me = bigmatpoly_my_cell(p);
    /* For matpoly, init or finish_init are identical */
    matpoly_init(ab, me, p->m0, p->n0, len);
}
/* If m,n,len are zero, then this is exactly equivalent to duplicating
 * the mpi_model. In which case the initialization may be completed later
 * on with bigmatpoly_finish_init
 */
void bigmatpoly_init(abdst_field ab, bigmatpoly_ptr p, bigmatpoly_srcptr model, unsigned int m, unsigned int n, int len)
{
    ASSERT_ALWAYS(model);
    ASSERT_ALWAYS(m % model->m1 == 0);
    ASSERT_ALWAYS(n % model->n1 == 0);
    int irank;
    int jrank;
    MPI_Comm_rank(model->com[2], &irank);
    MPI_Comm_rank(model->com[1], &jrank);
    memset(p, 0, sizeof(bigmatpoly));
    p->m1 = model->m1;
    p->n1 = model->n1;
    memcpy(p->com, model->com, 3 * sizeof(MPI_Comm));

    p->cells = malloc(p->m1*p->n1*sizeof(matpoly));
    memset(p->cells, 0, p->m1*p->n1*sizeof(matpoly));
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
    bigmatpoly_finish_init(ab, p, m, n, len);
}

/* This checks the validity of the m,n,size fields. If the structure
 * corresponds to a lazy allocation state, this function returns 1. If
 * the structure corresponds to something fully initializaed already,
 * this function returns 0. If the fields are filled in an inconsistent
 * manner, abort() is called.
 */
int bigmatpoly_check_pre_init(bigmatpoly_srcptr p)
{
    if (p->m && p->n && p->size)
        return 0;
    if (!p->m && !p->n && !p->size)
        return 1;
    abort();
    return 0;
}


void bigmatpoly_clear_model(bigmatpoly_ptr p)
{
    memset(p, 0, sizeof(bigmatpoly));
}

void bigmatpoly_clear(abdst_field ab, bigmatpoly_ptr p)
{
    for(unsigned int t = 0 ; t < p->m1 * p->n1 ; t++) {
        if (p->cells[t]->x) {
            matpoly_clear(ab, p->cells[t]);
        }
    }
    free(p->cells);
    bigmatpoly_clear_model(p);
}

/* Return a bitmask indicating whether bigmatpoly_provision_{row,col} has
 * been called on this matrix before. bit 0 is for row, bit 1 is for col.
 * If the returned value is zero, then we really have only the local part
 * of this matrix for the moment.
 */
int bigmatpoly_provisioned(bigmatpoly_ptr p)
{
    if (bigmatpoly_check_pre_init(p)) return 0;
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    unsigned int np_row = 0;
    unsigned int np_col = 0;
    for(unsigned int j = 0 ; j < p->n1 ; j++)
        np_row += bigmatpoly_cell(p, irank, j)->x != NULL;
    for(unsigned int i = 0 ; i < p->m1 ; i++)
        np_col += bigmatpoly_cell(p, i, jrank)->x != NULL;
    ASSERT_ALWAYS(np_row == 1 || np_row == p->n1);
    ASSERT_ALWAYS(np_col == 1 || np_col == p->m1);
    return (np_row>1)+((np_col>1)<<1);
}


/* We are a left multiplicand. This is a no-op if space for our row has
 * already been allocated */
void bigmatpoly_provision_row(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    matpoly_srcptr me = bigmatpoly_my_cell(p);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (them->x) continue;
        matpoly_init(ab, them, p->m0, p->n0, me->alloc);
    }
}
/* Rarely useful. We do need it because we resort to a kludgy
 * implementation of scatter_mat, which calls for provisioning on all
 * rows.
 */

void bigmatpoly_unprovision_row(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        matpoly_clear(ab, them);
    }
}

/* We are a right multiplicand. This is a no-op if space for our col has
 * already been allocated */
void bigmatpoly_provision_col(abdst_field ab, bigmatpoly_ptr p)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    matpoly_srcptr me = bigmatpoly_my_cell(p);
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (them->x) continue;
        matpoly_init(ab, them, p->m0, p->n0, me->alloc);
    }
}

/* Set size to be large enough to receive the given number of
 * coefficients. If space is provisioned for other cells in row or
 * column, allocation is triggered for them as well.
 */
void bigmatpoly_set_size(bigmatpoly_ptr p, size_t size)
{
    int irank;
    int jrank;
    MPI_Comm_rank(p->com[2], &irank);
    MPI_Comm_rank(p->com[1], &jrank);
    matpoly_ptr me = bigmatpoly_my_cell(p);
    ASSERT_ALWAYS(size <= me->alloc);
    p->size = size;
    me->size = size;
    for(unsigned int j = 0 ; j < p->n1 ; j++) {
        if (j == (unsigned int) jrank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, irank, j);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
    for(unsigned int i = 0 ; i < p->m1 ; i++) {
        if (i == (unsigned int) irank) continue;
        matpoly_ptr them = bigmatpoly_cell(p, i, jrank);
        if (!them->x) continue;
        them->size = size;
        ASSERT_ALWAYS(size <= them->alloc);
    }
}

/* If our row or col cells have already been allocated, then reallocate
 * them as well (XXX is it clear or not ?) */
#if 0 /* This function has never been used or needed */
void bigmatpoly_realloc(bigmatpoly_ptr p, int newalloc)
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
void bigmatpoly_zero(abdst_field ab, bigmatpoly_ptr p)
{
    matpoly_ptr me = bigmatpoly_my_cell(p);
    matpoly_zero(ab, me);
}

/* okay, it's ugly */
void bigmatpoly_swap(bigmatpoly_ptr a, bigmatpoly_ptr b)
{
    bigmatpoly x;
    memcpy(x, a, sizeof(bigmatpoly));
    memcpy(a, b, sizeof(bigmatpoly));
    memcpy(b, x, sizeof(bigmatpoly));
}

/* }}} */

void bigmatpoly_truncate_loc(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(bigmatpoly_provisioned(src) == 0);
    ASSERT_ALWAYS(bigmatpoly_provisioned(dst) == 0);
    if (bigmatpoly_check_pre_init(dst)) {
        bigmatpoly_finish_init(ab, dst, src->m, src->n, size);
    }
    matpoly_truncate(ab, bigmatpoly_my_cell(dst), bigmatpoly_my_cell(src), size);
    dst->size = bigmatpoly_my_cell(dst)->size;
}
/*}}}*/

void bigmatpoly_coeff_set_zero_loc(abdst_field ab, bigmatpoly_ptr pi, unsigned int k)
{
    matpoly_coeff_set_zero(ab, bigmatpoly_my_cell(pi), k);
}

int bigmatpoly_coeff_is_zero(abdst_field ab, bigmatpoly_srcptr pi, unsigned int k)
{
    int t = matpoly_coeff_is_zero(ab, bigmatpoly_my_cell_const(pi), k);
    MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_INT, MPI_MIN, pi->com[0]);
    return t;
}


void bigmatpoly_rshift(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int k)/*{{{*/
{
    if (bigmatpoly_check_pre_init(dst)) {
        bigmatpoly_finish_init(ab, dst, src->m, src->n, src->size - k);
    }
    ASSERT_ALWAYS(bigmatpoly_provisioned(dst) == 0);
    ASSERT_ALWAYS(bigmatpoly_provisioned(src) == 0);
    matpoly_rshift(ab, bigmatpoly_my_cell(dst), bigmatpoly_my_cell(src), k);
    dst->size = bigmatpoly_my_cell(dst)->size;
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
void bigmatpoly_allgather_row(abdst_field ab, bigmatpoly a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->com[2], &irank);
    MPI_Comm_rank(a->com[1], &jrank);
    bigmatpoly_provision_row(ab, a);
    for(unsigned int k = 0 ; k < a->n1 ; k++) {
        matpoly_ptr data = bigmatpoly_cell(a, irank, k);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->com[1]);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        ASSERT_ALWAYS((data->m * data->n * data->alloc) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data->m * data->n, data->alloc * abvec_elt_stride(ab, 1),
                MPI_BYTE,
                MPI_Bcast(data->x, _datasize, _datatype, k, a->com[1])
        );
    }
}
void bigmatpoly_allgather_col(abdst_field ab, bigmatpoly a)
{
    int irank;
    int jrank;
    MPI_Comm_rank(a->com[2], &irank);
    MPI_Comm_rank(a->com[1], &jrank);
    bigmatpoly_provision_col(ab, a);
    for(unsigned int k = 0 ; k < a->m1 ; k++) {
        matpoly_ptr data = bigmatpoly_cell(a, k, jrank);
        /* XXX: Should we ensure earlier that we agree on the size ? */
        unsigned long size = data->size;
        MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, k, a->com[2]);
        data->size = size;
        ASSERT_ALWAYS(data->size <= data->alloc);
        ASSERT_ALWAYS((data->m * data->n * data->alloc) < (size_t) INT_MAX);
        CHECK_MPI_DATASIZE_FITS(
                data->m * data->n, data->alloc * abvec_elt_stride(ab, 1),
                MPI_BYTE,
                MPI_Bcast(data->x, _datasize, _datatype, k, a->com[2])
        );
    }
}
/* }}} */

double bigmatpoly_mul(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b, int draft)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);

    ASSERT_ALWAYS(a->n == b->m);
    ASSERT_ALWAYS(a->n1 == b->m1);
    if (bigmatpoly_check_pre_init(c)) {
        bigmatpoly_finish_init(ab, c, a->m, b->n, csize);
    }
    ASSERT_ALWAYS(c->m);
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    int irank;
    int jrank;
    MPI_Comm_rank(c->com[2], &irank);
    MPI_Comm_rank(c->com[1], &jrank);
    bigmatpoly_allgather_row(ab, a);
    bigmatpoly_allgather_col(ab, b);
    bigmatpoly_set_size(c, csize);
    // ASSERT_ALWAYS(c->alloc >= c->size);

    matpoly_ptr lc = bigmatpoly_my_cell(c);
    bigmatpoly_zero(ab, c);
    lc->size = c->size;
    ASSERT_ALWAYS(a->n == b->m);

    double tt0=0, tt=0; /* tt=0: placate gcc */
    unsigned long c0 = 0;
    int set = 0;
    double x0 = 0;

    if (draft) { tt=wct_seconds(); tt0 -= tt; }
    for(unsigned int k = 0 ; k < a->n1 && !set ; k++) {
        x0 += matpoly_addmul(ab, lc, 
                bigmatpoly_cell(a, irank, k),
                bigmatpoly_cell(b, k, jrank), draft);
        if (draft) { ++c0; set = (tt0 + (tt = wct_seconds())) >= draft; }
    }

#if 0
    fprintf(stderr, "%s(%zu,%zu) %s after %lu/%u [%f]\n",
            __func__,
            a->size,b->size,
            set ? "gives up" : "does not give up",
            c0, a->n1, tt0+tt);
#endif

    if (!set) return x0;
    tt0 += tt;
    /* The *real* time we've spent is tt0, we should not count it */
    return (tt0 + x0) * a->n1 / c0 - tt0;
}/*}}}*/

double bigmatpoly_mp(abdst_field ab,/*{{{*/
        bigmatpoly b,
        bigmatpoly a,
        bigmatpoly c,
        int draft)
{
    if (bigmatpoly_check_pre_init(b)) {
        bigmatpoly_finish_init(ab, b, a->m, c->n, MAX(a->size, c->size) - MIN(a->size, c->size) + 1);
    }
    ASSERT_ALWAYS(a->n == c->m);
    ASSERT_ALWAYS(a->n1 == c->m1);
    ASSERT_ALWAYS(b->m);
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    int irank;
    int jrank;
    MPI_Comm_rank(b->com[2], &irank);
    MPI_Comm_rank(b->com[1], &jrank);

    bigmatpoly_allgather_row(ab, a);
    bigmatpoly_allgather_col(ab, c);

    matpoly_ptr lb = bigmatpoly_my_cell(b);
    matpoly_zero(ab, lb);
    bigmatpoly_set_size(b, b->size);
    ASSERT_ALWAYS(lb->size == b->size);
    ASSERT_ALWAYS(lb->alloc >= lb->size);

    double tt0=0, tt=0; /* tt=0: placate gcc */
    unsigned long c0 = 0;
    int set = 0;
    double x0 = 0;

    if (draft) { tt=wct_seconds(); tt0 -= tt; }
    for(unsigned int k = 0 ; k < a->n1 && !set ; k++) {
        x0 += matpoly_addmp(ab, lb,
                bigmatpoly_cell(a, irank, k),
                bigmatpoly_cell(c, k, jrank), draft);
        if (draft) { ++c0; set = (tt0 + (tt = wct_seconds())) >= draft; }
    }

#if 0
    fprintf(stderr, "%s(%zu,%zu) %s after %lu/%u [%f]\n",
            __func__,
            a->size,c->size,
            set ? "gives up" : "does not give up",
            c0, a->n1, tt0+tt);
#endif

    if (!set) return x0;
    tt0 += tt;
    /* The *real* time we've spent is tt0, we should not count it */
    return (tt0 + x0) * a->n1 / c0 - tt0;
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
void bigmatpoly_gather_mat_partial(abdst_field ab, matpoly_ptr dst, bigmatpoly_srcptr src,
        size_t offset, size_t length)
{
    int rank;
    int irank;
    int jrank;
    MPI_Comm_rank(src->com[0], &rank);
    MPI_Comm_rank(src->com[2], &irank);
    MPI_Comm_rank(src->com[1], &jrank);

    /* sanity checks, because the code below assumes this. */
    ASSERT_ALWAYS(irank * (int) src->n1 + jrank == rank);
    ASSERT_ALWAYS(length <= (size_t) INT_MAX);

    MPI_Datatype mt;
    MPI_Type_contiguous(length * abvec_elt_stride(ab, 1), MPI_BYTE, &mt);
    MPI_Type_commit(&mt);

    // Node 0 receives data
    if (!rank) {
        ASSERT_ALWAYS(dst->m == src->m);
        ASSERT_ALWAYS(dst->n == src->n);
        ASSERT_ALWAYS(dst->alloc >= length);
        dst->size = length;
        MPI_Request * reqs = malloc(src->m1 * src->n1
                * src->m0 * src->n0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        /* the master receives data from everyone */
        for(unsigned int i1 = 0 ; i1 < src->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < src->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < src->n1 ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < src->n0 ; j0++) {
                        unsigned int ii = i1 * src->m0 + i0;
                        unsigned int jj = j1 * src->n0 + j0;
                        unsigned int peer = i1 * src->n1 + j1;
                        unsigned int tag = ii * src->n + jj;
                        abdst_vec to = matpoly_part(ab, dst, ii, jj, 0);

                        if (peer == 0) {
                            /* talk to ourself */
                            matpoly_srcptr me = bigmatpoly_my_cell_const(src);
                            ASSERT_ALWAYS(offset + length <= me->alloc);
                            absrc_vec from = matpoly_part_const(ab, me, i0, j0, offset);
                            abvec_set(ab, to, from, length);
                        } else {
                            MPI_Irecv(to, 1, mt, peer, tag, src->com[0], req);
                        }
                        req++;
                    }
                }
            }
        }

        req = reqs;
        for(unsigned int i1 = 0 ; i1 < src->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < src->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < src->n1 ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < src->n0 ; j0++) {
                        unsigned int peer = i1 * src->n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
        free(reqs);
    } else {
        // All the other nodes send their data.
        MPI_Request * reqs = malloc(src->m0 * src->n0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        /* receive. Each job will receive exactly dst->m0 transfers */
        matpoly_srcptr me = bigmatpoly_my_cell_const(src);
        for(unsigned int i0 = 0 ; i0 < src->m0 ; i0++) {
            for(unsigned int j0 = 0 ; j0 < src->n0 ; j0++) {
                unsigned int i1 = irank;
                unsigned int j1 = jrank;
                unsigned int ii = i1 * src->m0 + i0;
                unsigned int jj = j1 * src->n0 + j0;
                unsigned int tag = ii * src->n + jj;
                absrc_vec from = matpoly_part_const(ab, me, i0, j0, offset);
                /* battle with const-deprived MPI prototypes... */
                MPI_Isend((void*)from, 1, mt, 0, tag, src->com[0], req);
                req++;
            }
        }
        req = reqs;
        for(unsigned int i0 = 0 ; i0 < src->m0 ; i0++) {
            for(unsigned int j0 = 0 ; j0 < src->n0 ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
        free(reqs);
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
void bigmatpoly_scatter_mat_partial(abdst_field ab,
        bigmatpoly_ptr dst, matpoly_ptr src,
        size_t offset, size_t length)
{
    int rank;
    int irank;
    int jrank;
    MPI_Comm_rank(dst->com[0], &rank);
    MPI_Comm_rank(dst->com[2], &irank);
    MPI_Comm_rank(dst->com[1], &jrank);

    bigmatpoly_set_size(dst, offset+length);

    MPI_Datatype mt;
    MPI_Type_contiguous(length * abvec_elt_stride(ab, 1), MPI_BYTE, &mt);
    MPI_Type_commit(&mt);


    /* sanity check, because the code below assumes this. */
    ASSERT_ALWAYS(irank * (int) dst->n1 + jrank == rank);

    if (!rank) {
        MPI_Request * reqs = malloc(dst->m1 * dst->n1
                * dst->m0 * dst->n0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        /* the master sends data to everyone */
        for(unsigned int i1 = 0 ; i1 < dst->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < dst->n0 ; j0++) {
                        unsigned int ii = i1 * dst->m0 + i0;
                        unsigned int jj = j1 * dst->n0 + j0;
                        unsigned int peer = i1 * dst->n1 + j1;
                        unsigned int tag = ii * dst->n + jj;
                        abdst_vec from = matpoly_part(ab, src, ii, jj, 0);

                        if (peer == 0) {
                            /* talk to ourself */
                            matpoly_ptr me = bigmatpoly_my_cell(dst);
                            abdst_vec to = matpoly_part(ab, me, i0, j0, offset);
                            abvec_set(ab, to, from, length);
                        } else {
                            MPI_Isend(from, 1, mt, peer, tag, dst->com[0], req);
                        }
                        req++;
                    }
                }
            }
        }

        req = reqs;
        for(unsigned int i1 = 0 ; i1 < dst->m1 ; i1++) {
            for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
                for(unsigned int j1 = 0 ; j1 < dst->n1 ; j1++) {
                    for(unsigned int j0 = 0 ; j0 < dst->n0 ; j0++) {
                        unsigned int peer = i1 * dst->n1 + j1;
                        if (peer)
                            MPI_Wait(req, MPI_STATUS_IGNORE);
                        req++;
                    }
                }
            }
        }
        free(reqs);
    } else {
        MPI_Request * reqs = malloc(dst->m0 * dst->n0 * sizeof(MPI_Request));
        MPI_Request * req = reqs;
        matpoly_ptr me = bigmatpoly_my_cell(dst);
        for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
            for(unsigned int j0 = 0 ; j0 < dst->n0 ; j0++) {
                unsigned int i1 = irank;
                unsigned int j1 = jrank;
                unsigned int ii = i1 * dst->m0 + i0;
                unsigned int jj = j1 * dst->n0 + j0;
                unsigned int tag = ii * dst->n + jj;
                abdst_vec to = matpoly_part(ab, me, i0, j0, offset);
                MPI_Irecv(to, 1, mt, 0, tag, dst->com[0], req);
                req++;
            }
        }
        req = reqs;
        for(unsigned int i0 = 0 ; i0 < dst->m0 ; i0++) {
            for(unsigned int j0 = 0 ; j0 < dst->n0 ; j0++) {
                MPI_Wait(req, MPI_STATUS_IGNORE);
                req++;
            }
        }
        free(reqs);
    }
    MPI_Type_free(&mt);
}



/* Collect everything into node 0 */
void bigmatpoly_gather_mat(abdst_field ab, matpoly dst, bigmatpoly src)
{
    matpoly dst_partial;
    size_t length = 100;
    int rank;
    MPI_Comm_rank(src->com[0], &rank);

    if (!rank) {
        // Leader should initialize the result matrix
        if (matpoly_check_pre_init(dst)) {
            matpoly_init(ab, dst, src->m, src->n, src->size);
        }
        dst->size = src->size;

        // Leader creates a buffer matpoly of size length
        matpoly_init(ab, dst_partial, src->m, src->n, length);
        dst_partial->size = length;
    }

    size_t offset = 0;
    while (src->size > offset) {
        size_t len = MIN(length, (src->size-offset));
        bigmatpoly_gather_mat_partial(ab, dst_partial, src, offset, len);

        // Copy the partial data into dst. This is the place where we
        // could write directly on disk if memory is a concern:
        if (!rank) {
            for (unsigned int i = 0; i < dst->m; ++i) {
                for (unsigned int j = 0; j < dst->n; ++j) {
                    abdst_vec to = matpoly_part(ab, dst, i, j, offset);
                    absrc_vec from = matpoly_part(ab, dst_partial, i, j, 0);
                    abvec_set(ab, to, from, len);
                }
            }
        }
        offset += len;
    }

    if (!rank) {
        matpoly_clear(ab, dst_partial);
    }
}


/* Exactly the converse of the previous function. */
void bigmatpoly_scatter_mat(abdst_field ab,
        bigmatpoly_ptr dst, matpoly_ptr src)
{
    matpoly src_partial;
    size_t length = 100;
    int rank;
    MPI_Comm_rank(dst->com[0], &rank);

    /* share allocation size. */
    matpoly shell;
    shell->m = src->m;
    shell->n = src->n;
    shell->size = src->size;
    shell->alloc = src->alloc;

    MPI_Bcast(shell, sizeof(matpoly), MPI_BYTE, 0, dst->com[0]);

    /* dst must be in pre-init mode */
    ASSERT_ALWAYS(bigmatpoly_check_pre_init(dst));

    /* Allocate enough space on each node */
    bigmatpoly_finish_init(ab, dst, shell->m, shell->n, shell->alloc);
    // bigmatpoly_set_size(dst, shell->size);

    if (!rank) {
        // Leader creates a buffer matpoly of size length
        matpoly_init(ab, src_partial, src->m, src->n, length);
        src_partial->size = length;
    }

    size_t offset = 0;
    while (shell->size > offset) {
        size_t len = MIN(length, (shell->size-offset));
        // Copy the partial data into src_partial. This is the place where we
        // could read directly from disk if memory is a concern:
        if (!rank) {
            for (unsigned int i = 0; i < src->m; ++i) {
                for (unsigned int j = 0; j < src->n; ++j) {
                    abdst_vec to = matpoly_part(ab, src_partial, i, j, 0);
                    absrc_vec from = matpoly_part(ab, src, i, j, offset);
                    abvec_set(ab, to, from, len);
                }
            }
        }
        bigmatpoly_scatter_mat_partial(ab, dst, src_partial, offset, len);
        offset += len;
    }

    if (!rank) {
        matpoly_clear(ab, src_partial);
    }
}
