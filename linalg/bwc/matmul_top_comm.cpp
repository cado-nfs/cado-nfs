#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "macros.h"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "matmul_top_vec.hpp"
#include "misc.h"
#include "parallelizing_info.hpp"
#include "select_mpi.h"
#include "timing.h"
#include "verbose.h"

/* Our innermost communication routines are essentially all-gather and
 * reduce-scatter, following the MPI terminology. We provide several
 * choices for doing this. The compile-time definitions here allow to
 * change which is used. The "best" should in theory always be the one
 * provided by the MPI implementation. Unfortunately, it is not always
 * that clear.
 *
 * Note that because we have several threads wanting to do all these
 * operations one after the other, we also have some interest in using
 * non-blockng collectives, or even also some sort of
 * communication-computation overlap.
 */

/* choices for reduce-scatter */
#define RS_CHOICE_STOCK_RS           1
#define RS_CHOICE_STOCK_RSBLOCK      (MPI_VERSION_ATLEAST(2,2) ? 2 : -2)
#define RS_CHOICE_STOCK_IRS          (MPI_VERSION_ATLEAST(3,0) ? 3 : -3)
#define RS_CHOICE_STOCK_IRSBLOCK     (MPI_VERSION_ATLEAST(3,0) ? 4 : -4)
#define RS_CHOICE_MINE               5
#define RS_CHOICE_MINE_DROP_IN       (MPI_VERSION_ATLEAST(3,0) ? 6 : -6)
#define RS_CHOICE_MINE_PARALLEL      7
#define RS_CHOICE_MINE_OVERLAPPING   8  /* TODO */
/* choices for all-gather */
#define AG_CHOICE_STOCK_AG           1
#define AG_CHOICE_STOCK_IAG          (MPI_VERSION_ATLEAST(3,0) ? 2 : -2)

/* The actual performance is affected by the communicator size, the chunk
 * size, and the operation. When we're doing bxor, we have the following
 * measurements for RS.
 *
 * n=2, chunk=174 MB, openmpi-1.8.3, IB FDR.
 * RS_CHOICE_MINE_PARALLEL .25
 * RS_CHOICE_MINE .45
 * RS_CHOICE_MINE_DROP_IN  .48
 * RS_CHOICE_STOCK_RS      .77
 * RS_CHOICE_STOCK_RSBLOCK 1.63    
 * RS_CHOICE_STOCK_IRS        => not measured
 * RS_CHOICE_STOCK_IRSBLOCK        => not measured
 * RS_CHOICE_MINE_OVERLAPPING      => to be implemented
 */

/* _this_ part can be configured */
// #define RS_CHOICE RS_CHOICE_STOCK_IRSBLOCK
// #define AG_CHOICE AG_CHOICE_STOCK_IAG
#define RS_CHOICE RS_CHOICE_MINE_PARALLEL
#define AG_CHOICE (MPI_VERSION_ATLEAST(3,0) ? AG_CHOICE_STOCK_IAG : AG_CHOICE_STOCK_AG)

/* some sanity checking */
#if RS_CHOICE < 0
#error "Choice for reduce-scatter strategy is invalid or not supported"
#endif
#if AG_CHOICE < 0
#error "Choice for reduce-scatter strategy is invalid or not supported"
#endif

/* we no longer use this.
#ifdef  HAVE_MPI3_API
#define MPI_LIBRARY_NONBLOCKING_COLLECTIVES
#endif
*/

/* {{{ mmt_vec_broadcast (generic interface) */
/* mmt_vec_broadcast reads data in mmt.wr[d]->v, and broadcasts it across the
 * communicator mmt.pi->wr[d] ; eventually everybody on the communicator
 * mmt.pi->wr[d] has the data.
 *
 * Note that the combination of mmt_vec_reduce + mmt_vec_broadcast is not the
 * identity (because of the shuffled_product).
 */

/* XXX
 * ``across'' (horizontally) and ``down'' (vertically) are just here for
 * exposition. The binding of this operation to job/thread arrangement
 * is flexible, through argument d.
 * XXX
 */
void
mmt_vec_broadcast(mmt_vec & v)
{
    ASSERT_ALWAYS(v.consistency == 1);

    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
    int err;

    /* communicator wr is in the direction we are broadcasting */
    pi_comm_ptr wr = v.pi->wr[v.d];

    /* communicator xwr is in the other direction */
    pi_comm_ptr xwr = v.pi->wr[!v.d];
    mmt_vec ** xwrpals = v.wrpals[!v.d].get();

    pi_log_op(v.pi->m, "[%s:%d] enter first loop", __func__, __LINE__);
    /* Make sure that no thread on the column is wandering in other
     * places -- when we're leaving reduce, this is important. */
    serialize_threads(wr);

    /* This loop suffers from two-dimensional serializing, so the
     * simplistic macros SEVERAL_THREADS_PLAY_MPI_BEGIN and END do not
     * help here.
     */

    size_t const eblock = mmt_my_own_size_in_items(v);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (!mmt_vec_is_shared(v)) {
        /* not shared: begin by collecting everything on thread 0 */
        mmt_own_vec_set2(v, v.sibling(0), v);
    }
    serialize_threads(v.pi->m);
    if (wr->trank == 0 && xwr->trank == 0) {
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
        MPI_Request * req = (MPI_Request *) malloc(xwr->ncores * sizeof(MPI_Request));
#endif  /* AG_CHOICE == AG_CHOICE_STOCK_IAG */

        for(unsigned int t = 0 ; t < xwr->ncores ; t++) {
            // although the openmpi man page looks funny, I'm assuming that
            // MPI_Allgather wants MPI_IN_PLACE as a sendbuf argument.
            pi_log_op(wr, "[%s:%d] MPI_Allgather (round %u)", __func__, __LINE__, t);
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
            err = MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                    xwrpals[t]->v, v.abase->vec_elt_stride(1) * eblock * wr->ncores, MPI_BYTE, wr->pals, &req[t]);
#elif AG_CHOICE == AG_CHOICE_STOCK_AG
            err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, xwrpals[t]->v, v.abase->vec_elt_stride(1) * eblock * wr->ncores, MPI_BYTE, wr->pals);
#else   /* AG_CHOICE */
#error "Bad AG_CHOICE setting"
#endif  /* AG_CHOICE */
            pi_log_op(wr, "[%s:%d] MPI_Allgather (round %u) done", __func__, __LINE__, t);
            ASSERT_ALWAYS(!err);
        }

#if AG_CHOICE == AG_CHOICE_STOCK_IAG
        for(unsigned int t = 0 ; t < xwr->ncores ; t++) {
            MPI_Wait(&req[t], MPI_STATUS_IGNORE);
        }
        free(req);
#endif  /* AG_CHOICE == AG_CHOICE_STOCK_IAG */
    }
    v.consistency = 2;
    serialize_threads(v.pi->m);
    if (!mmt_vec_is_shared(v)) {
        mmt_full_vec_set(v, v.sibling(0));
        serialize_threads(wr);
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    pi_log_op(v.pi->m, "[%s:%d] trailer", __func__, __LINE__);
}
/* }}} */

/* {{{ mmt_vec_reduce */
/* {{{ various reduce_scatter implementations */
/* {{{ alternative_reduce_scatter */
#if 1 || RS_CHOICE == RS_CHOICE_MINE
// TODO: test!
static void alternative_reduce_scatter [[maybe_unused]] (mmt_vec & v)
{
    pi_comm_ptr wr = v.pi->wr[v.d];
    unsigned int const njobs = wr->njobs;
    unsigned int const rank = wr->jrank;
    MPI_Datatype const t = v.pitype->datatype;

    size_t const eitems = mmt_my_own_size_in_items(v) * wr->ncores;
    if (v.rsbuf_items < eitems) {
        ASSERT_ALWAYS(v.rsbuf_items == 0);
        v.rsbuf[0] = v.abase->alloc_vector(eitems);
        v.rsbuf[1] = v.abase->alloc_vector(eitems);
        v.rsbuf_items = eitems;
    }

    arith_generic::elt *b[2];
    b[0] = v.rsbuf[0].get();
    b[1] = v.rsbuf[1].get();
    v.abase->vec_set_zero(b[0], eitems);

    unsigned int l = (rank + 1) % njobs;
    unsigned int const srank = (rank + 1) % njobs;
    unsigned int const drank = (rank + njobs - 1) % njobs;

    for (unsigned int i = 0; i < njobs; i++) {
        unsigned int const j0 = l * eitems; 
        unsigned int const j1 = j0 +  eitems;
        v.abase->vec_add_and_reduce(b[0],
                v.abase->vec_subvec(v.sibling(0).v, j0),
                j1-j0);
        if (i == njobs - 1)  
            break;
        MPI_Sendrecv(b[0], eitems, t, drank, (i<<16) + rank,
                     b[1], eitems, t, srank, (i<<16) + srank,
                     wr->pals, MPI_STATUS_IGNORE);
        arith_generic::elt * tb = b[0];   b[0] = b[1];  b[1] = tb; 
        l = (l + 1) % njobs;
    }
    /* This is going at the beginning of the memory area, which is weird.
     * But it's the convention used by the MPI call...
     */
    v.abase->vec_set(v.sibling(0).v, b[0], eitems);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE */
/* }}} */
/* {{{ alternative_reduce_scatter_parallel */
#if 1 || RS_CHOICE == RS_CHOICE_MINE_PARALLEL
/* Example data for a factoring matrix (rsa100) of size 135820*135692,
 * split over 2x3 mpi jobs, and 7x5 threads.
 *
 * 2 == mmt.pi->wr[1]->njobs (number of jobs encountered on a vertical axis).
 * 7 == mmt.pi->wr[1]->ncores (number of jobs per core on a vertical axis).
 * 3 == mmt.pi->wr[0]->njobs (number of jobs encountered on an horiz. axis).
 * 5 == mmt.pi->wr[0]->ncores (number of jobs per core on an horiz. axis).
 *
 * matrix is padded to a multiple of 210 = 2*3*5*7, which is * N=135870=210*647
 *
 * we'll call 647 the "small chunk" size.
 *
 * for all jobs/threads, the following relations hold:
 *      mmt.wr[0]->i1 - mmt.wr[0]->i0 == N/14 == 9705 == 15 * 647
 *      mmt.wr[1]->i1 - mmt.wr[1]->i0 == N/15 == 9058 == 14 * 647
 *
 * a mmt_vec_reduce operation, in the context of factoring, is with d==1
 * below. Hence, in fact, we're doing a reduction down a column.
 *
 * the value eitems fed to this function is mmt.pi->wr[d]->ncores (here,
 * 7) times the small chunk size. Here 7*647 == 4529.
 */
/* all threads in mmt.wr[!d], one after another a priori, are going to
 * do alternative_reduce_scatter on their vector v[i]
 */
// TODO: test!
static void alternative_reduce_scatter_parallel [[maybe_unused]] (pi_comm_ptr xr, mmt_vec ** vs)
{
    /* we write all data counts below with comments indicating the typical
     * size in our toy example above */

    /* he have xr->ncores vectors. The vs[] array accesses data from the
     * other peers. Note that the pi structures belong separately to each
     * peer, and so do the embedded communicators. Therefore the proper
     * way to see "our" xr is none other that the given parameter. And for
     * "our" wr, then we'll have to look for our own vector within vs.
     */
    mmt_vec & v = *(vs[xr->trank]);
    arith_generic * ab = v.abase;
    pi_comm_ptr wr = v.pi->wr[v.d];  /* 2 jobs, 7 cores */
    /* what we're going to do will happen completely in parallel over
     * xr->njobs==3 communicators. We're simulating here the split into 5
     * communicators, even though those collide into a unique
     * communicator as far as MPI is concerned.
     */
    unsigned int const njobs = wr->njobs;      /* 2 */
    /* note that we no longer care at this point about what happened at
     * the thread level in our dimension. This is already done, period.
     */
    unsigned int const rank = wr->jrank;
    MPI_Datatype const t = v.pitype->datatype;

    /* If the rsbuf[] buffers have not yet been allocated, it is time to
     * do so now. We also take the opportunity to possibly re-allocate
     * them if because of a larger abase, the corresponding storage has
     * to be expanded.
     */
    size_t const eitems = mmt_my_own_size_in_items(v) * wr->ncores;
    /* notice that we are allocating a temp buffer only for one vector.
     * Of course, since this is a multithreaded routine, each thread in
     * xr is doing so at the same time */
    if (v.rsbuf_items < eitems) {
        /* It's very very ugly, right? We should probably _at least_ let
         * this go through the virtual hierarchy */
        v.rsbuf[0] = v.abase->alloc_vector(eitems);
        v.rsbuf[1] = v.abase->alloc_vector(eitems);
        v.rsbuf_items = eitems;
    }

    ab->vec_set_zero(v.rsbuf[0].get(), eitems);
    serialize_threads(xr);

    unsigned int const srank = (rank + 1) % njobs;
    unsigned int const drank = (rank + njobs - 1) % njobs;

    /* We describe the algorithm for one of the xr->ncores==5 threads.
     * local vector areas [9058] split into [wr->njobs==2] sub-areas of
     * size [eitems (==4529)].
     *
     * on each job, each thread has 9058 == 2*4529 items
     *          more generally: njobs * eitems items.
     *
     * each has a place where it wants to receive data, let's call that
     * "his" zone, even though everyone has data at all places.
     *
     */
    /* scheme for 4 jobs: (input, output. All chunks have size eitmes. AA
     * is the sum a3+a2+a1+a0)
     *
     * j0: a0 b0 c0 d0      j0: AA .. .. ..
     * j1: a1 b1 c1 d1      j1: .. BB .. ..
     * j2: a2 b2 c2 d2      j2: .. .. CC ..
     * j3: a3 b3 c3 d3      j3: .. .. .. DD
     *
     * loop 0: j0: sends       b0 to j3, receives       c1 from j1
     * loop 1: j0: sends    c1+c0 to j3, receives    d2+d1 from j1
     * loop 2: j0: sends d2+d1+d0 to j3, receives a3+a2+a1 from j1
     * loop 3: we only have to adjust.
     */

    for (unsigned int i = 0, s = 0; i < njobs; i++, s=!s) {
        unsigned int const l = (rank + 1 + i) % njobs;
        unsigned int const j0 = l * eitems; 
        unsigned int const j1 = j0 +  eitems;
        ab->vec_add_and_reduce(v.rsbuf[s].get(),
                ab->vec_subvec(v.sibling(0).v, j0),
                j1-j0);
        serialize_threads(xr);

        if (i == njobs - 1)  
            break;

        if (xr->trank == 0) {
            MPI_Request * r = nullptr;
            r = (MPI_Request *) malloc(2 * xr->ncores * sizeof(MPI_Request));
            for(unsigned int w = 0 ; w < xr->ncores ; w++) {
                MPI_Request * rs = r + 2*w;
                MPI_Request * rr = r + 2*w + 1;
                MPI_Isend(vs[w]->rsbuf[s].get(),  eitems, t, drank, 0xb00+w, wr->pals, rs);
                MPI_Irecv(vs[w]->rsbuf[!s].get(), eitems, t, srank, 0xb00+w, wr->pals, rr);
                /*
                MPI_Sendrecv(vs[w]->rsbuf[s].get(), eitems, t, drank, 0xbeef,
                        vs[w]->rsbuf[!s].get(), eitems, t, srank, 0xbeef,
                        wr->pals, MPI_STATUS_IGNORE);
                        */
                // MPI_Waitall(2, r + 2*w, MPI_STATUSES_IGNORE);
            }
            MPI_Waitall(2 * xr->ncores, r, MPI_STATUSES_IGNORE);
            free(r);
        }
        serialize_threads(xr);
    }

    ab->vec_set(v.sibling(0).v, v.rsbuf[(njobs-1)&1].get(), eitems);

    pi_log_op(wr, "[%s:%d] MPI_Reduce_scatter done", __func__, __LINE__);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE_PARALLEL */
/* }}} */
/* {{{ my_MPI_Reduce_scatter_block */
#if 1 || RS_CHOICE == RS_CHOICE_MINE_DROP_IN
// TODO: test!
static int my_MPI_Reduce_scatter_block [[maybe_unused]] (void *sendbuf, void *recvbuf, int recvcount,
                MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, MPI_Comm wr)
{
    ASSERT_ALWAYS(sendbuf == MPI_IN_PLACE);
    int njobs;
    int rank;
    MPI_Comm_size(wr, &njobs);
    MPI_Comm_rank(wr, &rank);

    int tsize;
    MPI_Type_size(datatype, &tsize);

    void * v = recvbuf;
    
    /* This is a deliberate leak. Note that we expect to be serailized
     * here, so there is no concurrency issue with the static data. */
    static size_t rsbuf_items = 0;
    static arith_generic::elt * rsbuf[2] { nullptr, nullptr };

    arith_generic * abase = pi_arith_datatype_get_abase(datatype);

    if (rsbuf_items < (size_t) recvcount) {
        /* FIXME: how do I access abase ? */
        rsbuf[0] = abase->realloc(rsbuf[0], rsbuf_items, recvcount);
        rsbuf[1] = abase->realloc(rsbuf[1], rsbuf_items, recvcount);
        rsbuf_items = recvcount;
    }

    memset(rsbuf[0], 0, recvcount * tsize);

    int const srank = (rank + 1) % njobs;
    int const drank = (rank + njobs - 1) % njobs;

    for (int i = 0, w = 0; i < njobs; i++, w^=1) {
#if MPI_VERSION_ATLEAST(2,2)
        int j0 = ((rank + i + 1) % njobs) * recvcount;
        const void * share = pointer_arith_const(v, j0 * tsize);
        MPI_Reduce_local(share, rsbuf[w], recvcount, datatype, op);
#else
        /* XXX it's miserable, really. We used to have code pretending to
         * do this, but it worked only with unsigned longs and bxors
         * anyway.  In theory, if we ever reach here, it means that we
         * have some real MPI implementation, and not the fakempi.h one.
         * For some obscure reason, fakempi.h cannot advertise an
         * MPI_Reduce_local, so we're really out of luck (at this point,
         * at least. I might get back to it).
         */
        abort();
#endif
        if (i == njobs - 1) {
            memcpy(v, rsbuf[w], recvcount * tsize);
            break;
        }
        MPI_Sendrecv(rsbuf[w],  recvcount, datatype, drank, (i<<16) + rank,
                     rsbuf[!w], recvcount, datatype, srank, (i<<16) + srank,
                     wr, MPI_STATUS_IGNORE);
    }
    return 0;
}
#endif
/* }}} */
/* }}} */

/* mmt_vec_reduce_inner reads data in v (vector for side d), sums it up
 * across the communicator mmt.pi->wr[d], and collects the results in
 * vector v again, except that it's put in thread0's buffer (counting
 * threads in direction d, of course), *AT THE BEGINNING* of the data
 * area (which is surprising).
 *
 * mmt_vec_reduce completes the work by saving the resulting data in vector w
 * (vector in dimension !d).
 *
 * Note that the combination of mmt_vec_reduce + mmt_vec_broadcast is not the
 * identity (because of the shuffled_product).
 */
static void mmt_vec_reduce_inner(mmt_vec & v)
{
    ASSERT_ALWAYS(v.consistency != 2);

    /* reducing across a row is when d == 0 */
    pi_comm_ptr wr = v.pi->wr[v.d];
    pi_comm_ptr xr = v.pi->wr[!v.d];

    pi_log_op(v.pi->m, "[%s:%d] enter first loop", __func__, __LINE__);

    // I don't think that the case of shared vectors has been tested
    // correctly for reduction. Well, to start with, I doubt it really
    // makes a lot of sense anyhow.
    // ASSERT_ALWAYS((v.flags & THREAD_SHARED_VECTOR) == 0);

    if (wr->ncores > 1 && !mmt_vec_is_shared(v)) {
        /* row threads have to sum up their data. Of course it's
         * irrelevant when there is only one such thread...
         *
         * Concerning locking, we have to make sure that everybody on the
         * row has finished its computation task, but besides that,
         * there's no locking until we start mpi stuff.
         */
        pi_log_op(wr, "[%s:%d] serialize_threads", __func__, __LINE__);
        serialize_threads(wr);
        pi_log_op(wr, "[%s:%d] serialize_threads done", __func__, __LINE__);

        /* Our [i0,i1[ range is split into wr->ncores parts. This range
         * represent coordinates which are common to all threads on
         * wr. Corresponding data has to be summed. Amongst the
         * wr->ncores threads, thread k takes the responsibility of
         * summing data in the k-th block (that is, indices
         * i0+k*(i1-i0)/wr->ncores and further). As a convention,
         * the thread which eventually owns the final data is thread 0.
         *
         * Note that one should consider that data in threads other than
         * the destination thread may be clobbered by the operation
         * (although in the present implementation it is not).
         */
        size_t const thread_chunk = wr->njobs * mmt_my_own_size_in_items(v);
        arith_generic::elt * dptr = v.abase->vec_subvec(
                v.sibling(0).v, 
                wr->trank * thread_chunk);
        for(unsigned int w = 1 ; w < wr->ncores ; w++) {
            const arith_generic::elt * sptr = v.abase->vec_subvec(
                    v.sibling(w).v, 
                    wr->trank * thread_chunk);
            v.abase->vec_add_and_reduce(dptr, sptr, thread_chunk);
        }
        pi_log_op(wr, "[%s:%d] thread reduction done", __func__, __LINE__);
    }

    /* Good. Now on each node, thread 0 has the reduced data for the
     * range [i0..i1[ (well, actually, this corresponds to indices
     * distributed in a funny manner in the original matrix, but as far
     * as we care here, it's really the range v.i0..v.i1
     *
     * These data areas must now be reduced or reduce-scattered across
     * nodes. Note that in the case where the MPI library does not
     * support MT operation, we are performing the reduction in place,
     * and copy to the buffer in the other direction is done in a
     * secondary step -- which will be among the reasons which imply a
     * thread serialization before anything interesting can be done after
     * this function.
     */

    /* XXX We have to serialize threads here. At least row-wise, so that
     * the computations above are finished. Easily seen, that's just a
     * couple of lines above.
     *
     * Unfortunately, that's not the whole story. We must also serialize
     * column-wise. Because we're writing on the right vector buffers,
     * which are used as input for the matrix multiplication code --
     * which may be still running at this point because there has been no
     * column serialization in the current procedure yet.
     */

    pi_log_op(v.pi->m, "[%s:%d] secondary loop", __func__, __LINE__);

    pi_log_op(v.pi->m, "[%s:%d] serialize_threads", __func__, __LINE__);
    serialize_threads(wr);
    pi_log_op(v.pi->m, "[%s:%d] serialize_threads done", __func__, __LINE__);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (wr->trank == 0) {
        /* openmpi-1.8.2 does not seem to have a working non-blocking
         * reduce_scatter, at least not a very efficient one. All
         * I've been able to do is to run MPI_Ireduce_scatter with
         * MPI_UNSIGNED_LONG and MPI_BXOR. With custom types it seems
         * to crash. And anyway, it's very inefficient.
         */

        ASSERT((v.i1 - v.i0) % wr->totalsize == 0);
#if RS_CHOICE == RS_CHOICE_STOCK_RS 
        void * dptr = v.sibling(0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            // all recvcounts are equal
            int * rc = malloc(wr->njobs * sizeof(int));
            for(unsigned int k = 0 ; k < wr->njobs ; k++)
                rc[k] = (v.i1 - v.i0) / wr->njobs;
            int err = MPI_Reduce_scatter(dptr, dptr, rc,
                    v.pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            free(rc);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_RSBLOCK
        void * dptr = v.sibling(0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = MPI_Reduce_scatter_block(dptr, dptr,
                    (v.i1 - v.i0) / wr->njobs,
                    v.pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_MINE_DROP_IN
        void * dptr = v.sibling(0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = my_MPI_Reduce_scatter_block(MPI_IN_PLACE, dptr,
                    (v.i1 - v.i0) / wr->njobs,
                    v.pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_IRSBLOCK
        void * dptr = v.sibling(0)->v;
        auto req = pi_shared_array<MPI_Request>(xr, xr->ncores);
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = MPI_Ireduce_scatter_block(dptr, dptr,
                    (v.i1 - v.i0) / wr->njobs,
                    v.pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals, &req[t__]);
            ASSERT_ALWAYS(!err);
            pi_log_op(wr, "[%s:%d] MPI_Reduce_scatter done", __func__, __LINE__);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
        serialize_threads(xr);
        if (xr->trank == 0) {
            for(unsigned int t = 0 ; t < xr->ncores ; t++) {
                MPI_Wait(&req[t], MPI_STATUS_IGNORE);
            }
        }
#elif RS_CHOICE == RS_CHOICE_STOCK_IRS
        void * dptr = v.sibling(0)->v;
        auto req = pi_shared_array<MPI_Request>(xr, xr->ncores);
        int * rc = malloc(wr->njobs * sizeof(int));
        for(unsigned int k = 0 ; k < wr->njobs ; k++)
            rc[k] = (v.i1 - v.i0) / wr->njobs;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = MPI_Ireduce_scatter(dptr, dptr,
                    rc,
                    v.pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals, &req[t__]);
            free(rc);
            ASSERT_ALWAYS(!err);
            pi_log_op(wr, "[%s:%d] MPI_Reduce_scatter done", __func__, __LINE__);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
        serialize_threads(xr);
        if (xr->trank == 0) {
            for(unsigned int t = 0 ; t < xr->ncores ; t++) {
                MPI_Wait(&req[t], MPI_STATUS_IGNORE);
            }
        }
#elif RS_CHOICE == RS_CHOICE_MINE
        /* This strategy exposes code which is really similar to
         * RS_CHOICE_MINE_DROP_IN, with the only exception that we
         * have a slightly different interface. There's no reason for
         * both to stay in the long run.
         */
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            alternative_reduce_scatter(v);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_MINE_PARALLEL
#if 0
        auto vs = pi_shared_array<mmt_vec *>(xr, xr->ncores);
        vs[xr->trank] = &v;
        serialize_threads(xr);
        alternative_reduce_scatter_parallel(xr, vs.get());
#else
        alternative_reduce_scatter_parallel(xr, v.wrpals[!v.d].get());
#endif
#elif RS_CHOICE == RS_CHOICE_MINE_OVERLAPPING
#error "not implemented, but planned"
#endif
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif
    serialize_threads(wr);
}
void
mmt_vec_reduce(mmt_vec & w, mmt_vec & v)
{
    ASSERT_ALWAYS(v.abase == w.abase);
    ASSERT_ALWAYS(v.d != w.d);
    mmt_vec_reduce_inner(v);
    pi_comm_ptr wr = v.pi->wr[v.d];
    // row threads pick what they're interested in in thread0's reduced
    // buffer. Writes are non-overlapping in the mcol buffer here.
    // Different row threads always have different mcol buffers, and
    // sibling col threads write to different locations in their
    // (generally shared) mcol buffer, depending on which row they
    // intersect.

    // Notice that results are being written inside v.all_v[0],
    // just packed at the beginning.

    // row job rj, row thread rt has a col buffer containing
    // picol->totalsize blocks of size eblock.  col job cj, col
    // thread ct has in its row buffer pirow->totalsize blocks, but
    // only virtually. Because of the reduce_scatter operation, only
    // pirow->ncores blocks are here.
    //
    // Thus the k-th block in the row buffer is rather understood as
    // the one of indek rj * pirow->ncores + k in the data the row
    // threads have collectively computed.
    //
    // This block is of interest to the row thread of index k of
    // course, thus we restrict to rt==k
    //
    // Now among the picol->totalsize blocks of the col buffer, this
    // will go to position cj * picol->ncores + ct
 
    size_t const eblock = mmt_my_own_size_in_items(v);
    ASSERT_ALWAYS(mmt_my_own_size_in_items(w) == eblock);

    v.abase->vec_set(
            mmt_my_own_subvec(w),
            /* Note: reduce-scatter packs everything at the beginning in
             * the leader block, which is why we don't have our usual offset
             * here. */
            v.abase->vec_subvec(v.sibling(0).v, wr->trank * eblock),
            eblock);

    // as usual, we do not serialize on exit. Up to the next routine to
    // do so if needed.
    // 
    // what _is_ guaranteed is that in each column, the _leader_ thread
    // has all the necessary data to begin the column broadcast.
    //
    // In most cases, threads must be prevented from starting computation
    // before the leader thread has finished importing data. This means
    // that a column thread serialization is probably needed in most
    // circumstances after this step.
    w.consistency = 1;
}

/* This small variant is used only for the transposition routines. We do
 * not want, in that case, to rely on the data moving back and forth
 * between the left and right vectors, because in full generality, we
 * have no guarantee that they are the same size.
 *
 * Therefore we're merely building upon mmt_vec_reduce_inner, only with
 * the weirdo "pack-at-the-beginning" behaviour removed.
 */
void
mmt_vec_reduce_sameside(mmt_vec & v)
{
    pi_comm_ptr wr = v.pi->wr[v.d];
    mmt_vec_reduce_inner(v);
    size_t const eblock = mmt_my_own_size_in_items(v);
    v.abase->vec_set(
            mmt_my_own_subvec(v),
            v.abase->vec_subvec(v.sibling(0).v, wr->trank * eblock),
            eblock);
    /* This ensures that we've effectively *moved* the data, not copied
     * it */
    if (wr->trank) {
        v.abase->vec_set_zero(
                v.abase->vec_subvec(v.sibling(0).v, wr->trank * eblock),
                eblock);
    }
    serialize_threads(wr);
    v.consistency = 1;
}
/* }}} */

/* {{{ mmt_vec_allreduce */
/* This is only for convenience now. Eventually this will be relevant for
 * block Lanczos.  Note that allreduce is conceptually much simpler.
 * There is no funny permutation to be considered.
 */
void
mmt_vec_allreduce(mmt_vec & v)
{
    ASSERT_ALWAYS(v.consistency != 2);
    /* reducing across a row is when d == 0 */
    pi_comm_ptr wr = v.pi->wr[v.d];

    pi_log_op(v.pi->m, "[%s:%d] enter first loop", __func__, __LINE__);

    serialize_threads(v.pi->m);
    /* sum up row threads, so that only one thread on each row is used
     * for communication */
    size_t const thread_chunk = wr->njobs * mmt_my_own_size_in_items(v);
    if (!mmt_vec_is_shared(v)) {
        arith_generic::elt * dv = v.abase->vec_subvec(v.v,
                wr->trank * thread_chunk);
        for(unsigned int k = 1 ; k < wr->ncores ; k++) {
            arith_generic::elt * sv = v.abase->vec_subvec(
                    v.sibling((wr->trank+k) % wr->ncores).v,
                    wr->trank * thread_chunk);
            v.abase->vec_add_and_reduce(dv, sv, thread_chunk);
        }
    }
    /* Compared to the SEVERAL_THREADS_PLAY_MPI_BEGIN() approach, this
     * one has thread 0 do the work for all other threads, while other
     * threads are waiting.
     */
    SEVERAL_THREADS_PLAY_MPI_BEGIN2(v.pi->m, peer) {
        arith_generic::elt * dv = v.abase->vec_subvec(
                v.mpals.get()[peer]->v,
                v.mpals.get()[peer]->pi->wr[v.d]->trank * thread_chunk);
        MPI_Allreduce(MPI_IN_PLACE,
                dv,
                thread_chunk,
                v.pitype->datatype,
                BWC_PI_SUM->custom,
                wr->pals);
    }
    SEVERAL_THREADS_PLAY_MPI_END2(v.pi->m);
    if (!mmt_vec_is_shared(v)) {
        arith_generic::elt * sv = v.abase->vec_subvec(v.v,
                wr->trank * thread_chunk);
        for(unsigned int k = 1 ; k < wr->ncores ; k++) {
            arith_generic::elt * dv = v.abase->vec_subvec(
                    v.sibling((wr->trank+k) % wr->ncores).v,
                    wr->trank * thread_chunk);
            v.abase->vec_set(dv, sv, thread_chunk);
        }
    }
    v.consistency = 2;
}
/* }}} */

/**********************************************************************/
/* bench code */

static void matmul_top_comm_bench_helper(int * pk, double * pt,
                                  void (*f) (mmt_vec &),
				  mmt_vec & v)
{
    int k;
    double t0, t1;
    int cont;
    t0 = wct_seconds();
    for (k = 0;; k++) {
	t1 = wct_seconds();
	cont = t1 < t0 + 0.25;
	cont = cont && (t1 < t0 + 1 || k < 100);
        pi_allreduce(nullptr, &cont, 1, BWC_PI_INT, BWC_PI_MIN, v.pi->m);
	if (!cont)
	    break;
        /* It's difficult to be faithful to the requirements on
         * consistency here. But apparently 1 pleases both operations
         * tested. */
        v.consistency = 1;
	(*f) (v);
    }
    int target = 10 * k / (t1 - t0);
    ASSERT_ALWAYS(target >= 0);
    target = std::min(target, 100);
    if (target == 0)
        target = 1;
    pi_bcast(&target, 1, BWC_PI_INT, 0, 0, v.pi->m);
    t0 = wct_seconds();
    for (k = 0; k < target; k++) {
        pi_log_op(v.pi->m, "[%s] iter%d/%d", __func__, k, target);
        v.consistency = 1;     /* see above */
        (*f) (v);
    }
    serialize(v.pi->m);
    t1 = wct_seconds();
    *pk = k;
    *pt = t1 - t0;
}


void matmul_top_comm_bench(matmul_top_data & mmt, int d)
{
    /* like matmul_top_mul_comm, we'll call mmt_vec_reduce with !d, and
     * mmt_vec_broadcast with d */
    int k;
    double dt;

    void (*funcs[2])(mmt_vec &) = {
        mmt_vec_broadcast,
        mmt_vec_reduce_sameside,
    };
    const char * text[2] = { "bd", "ra" };

    arith_generic * abase = mmt.abase;

    int const is_shared[2] = {0,0};
    mmt_vec test_vectors[2] {
        { mmt, nullptr, nullptr, 0, is_shared[0], mmt.n[0] },
        { mmt, nullptr, nullptr, 1, is_shared[1], mmt.n[1] },
    };

    size_t datasize[2];
    
    {
        pi_comm_ptr pirow = mmt.pi->wr[!d];
        pi_comm_ptr picol = mmt.pi->wr[d];
        /* within each row, all jobs are concerned with the same range
         * vrow->i0 to vrow->i1. This is split into m=pirow->njobs
         * chunks, and reduce_scatter has m-1 communication rounds where
         * all of the m nodes output and receive one such chunk. All this
         * happens for all column threads, so we multiply by
         * picol->ncores, too.
         *
         * Note that vrow->i1 - vrow->i0 is #rows / picol->totalsize
         *
         * Note also that picol->ncores * #rows / picol->totalsize =
         * #rows / picol->njobs, so that the final thing we compute is
         * really:
         *      #rows / mmt.pi->m->njobs * (pirow->njobs - 1)
         */
        size_t const data_out_ra = abase->vec_elt_stride(
                picol->ncores * (mmt.n[!d] / picol->totalsize) /
                pirow->njobs * (pirow->njobs - 1));

        /* one way to do all-gather is to mimick this, except that each
         * node will output the same chunk at each round. Beyond that,
         * the calculation is similar, and we'll use it as a guide. Note
         * of course that if hardware-level multicast is used, our
         * throughput estimation is way off.
         *
         * as above, this is really:
         *      #cols / mmt.pi->m->njobs * (picol->njobs - 1)
         *
         */
        size_t const data_out_ag = abase->vec_elt_stride(
                pirow->ncores * (mmt.n[d] / pirow->totalsize) /
                picol->njobs * (picol->njobs - 1));

        datasize[0] = data_out_ag;
        datasize[1] = data_out_ra;
    }

    for(int s = 0 ; s < 2 ; s++) {
        /* we have our axis, and the other axis */
        pi_comm_ptr wr = mmt.pi->wr[d ^ s];          /* our axis */
        pi_comm_ptr xr = mmt.pi->wr[d ^ s ^ 1];      /* other axis */
        /* our operation has operated on the axis wr ; hence, we must
         * display data relative to the different indices within the
         * communicator xr.
         */
        matmul_top_comm_bench_helper(&k, &dt, funcs[s], test_vectors[d^s]);
        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_TIMING_GRIDS)) {
            for(unsigned int z = 0 ; z < xr->njobs ; z++) {
                if (xr->jrank == z && wr->jrank == 0) {
                    for(unsigned int w = 0 ; w < xr->ncores ; w++) {
                        if (xr->trank == w && wr->trank == 0) {
                            char buf[16];
                            printf("%s %2d/%d, %s: %2d in %.1fs ; one: %.2fs (xput: %s/s)\n",
                                    wr->th->desc,
                                    xr->jrank * xr->ncores + xr->trank,
                                    xr->totalsize,
                                    text[s],
                                    k, dt, dt/k,
                                    size_disp(datasize[d^s] * k/dt, buf));
                        }
                    }
                }
                serialize(mmt.pi->m);
            }
        }
        serialize(mmt.pi->m);
    }
}


