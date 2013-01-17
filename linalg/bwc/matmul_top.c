#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>             // truncate()
#include <sys/types.h>          // truncate()
#include <errno.h>
#include "bwc_config.h"
#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.h"
#include "debug.h"
#include "filenames.h"
#include "balancing_workhorse.h"
#include "misc.h"

#ifndef CONJUGATED_PERMUTATIONS
#error "Do you really, really want to use arbitrary left and right sigmas ?"
#endif


/* This enables an alternative algorithm instead of the stock
 * reduce_scatter implementation. It's slightly worrying to reach this
 * conclusion, but our simple contender wins with flying colours against
 * mvapich2's.
 */
#define USE_ALTERNATIVE_REDUCE_SCATTER

#define PADD(P, n)      abase_generic_ptr_add(P, n)


///////////////////////////////////////////////////////////////////
/* Start with stuff that does not depend on abase at all -- this
 * provides a half-baked interface */


/* At some point we had this. Not sure it's still useful. */
#define ABASE_UNIVERSAL_READAHEAD_ITEMS 8

#define ALWAYS_ALIGN_LARGE_MALLOCS      16

#ifdef  ALWAYS_ALIGN_LARGE_MALLOCS
#include "misc.h"       /* malloc/free_aligned are in utils */
#define alignable_malloc(s)  malloc_aligned((s), ALWAYS_ALIGN_LARGE_MALLOCS)
#define alignable_free(p,s)  free_aligned((p), (s), ALWAYS_ALIGN_LARGE_MALLOCS)
#else
#define alignable_malloc(s)  malloc((s))
#define alignable_free(p,s)  free((p))
#endif


/* {{{ vector init/clear (generic low-level interface) */
void vec_init_generic(pi_wiring_ptr picol, abase_vbase_ptr abase, mmt_vec_ptr v, int flags, unsigned int n)
{
    v->abase = abase;
    v->stride = v->abase->vec_elt_stride(v->abase, 1);
    v->flags = flags;

    // all_v pointers are stored twice, in a rotating buffer manner, so
    // that it's possible for one thread to address the n other threads
    // without having to compute a modulo
    size_t allocsize = 2 * picol->ncores * sizeof(void *);
    v->all_v = alignable_malloc(allocsize);

    /* Because we sometimes provide extra readahead space, we need to
     * zero out the allocated area. Strictly speaking, it would suffice
     * to zero out the readahead zone only, but it's easier this way.
     */
    if (flags & THREAD_SHARED_VECTOR) {
        void * r;
        if (picol->trank == 0) {
            abase->vec_init(abase, &r, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
            abase->vec_set_zero(abase, r, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
        }
        thread_broadcast(picol, &r, 0);
        v->v = r;
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            v->all_v[t] = v->v;
        }
    } else {
        abase->vec_init(abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
        abase->vec_set_zero(abase, v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
        v->all_v[picol->trank] = v->v;
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            serialize_threads(picol);
            void * r;
            r = v->v;
            thread_broadcast(picol, &r, t);
            v->all_v[t] = r;
        }
    }
    for(unsigned int t = 0 ; t < picol->ncores ; t++) {
        v->all_v[picol->ncores + t] = v->all_v[t];
    }
}

void vec_clear_generic(pi_wiring_ptr picol, mmt_vec_ptr v, unsigned int n MAYBE_UNUSED)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (picol->trank == 0)
            v->abase->vec_clear(v->abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
    } else {
        v->abase->vec_clear(v->abase, &v->v, n + ABASE_UNIVERSAL_READAHEAD_ITEMS);
    }
    size_t allocsize MAYBE_UNUSED;
    allocsize = 2 * picol->ncores * sizeof(void *);
    alignable_free(v->all_v, allocsize);
}
/*}}}*/
/* {{{ vector init/clear (still generic, but for usual-size vectors) */
void matmul_top_vec_init_generic(matmul_top_data_ptr mmt, abase_vbase_ptr abase, mmt_vec_ptr v, int d, int flags)
{
    if (v == NULL) v = mmt->wr[d]->v;
    if (abase == NULL) abase = mmt->abase;
    unsigned int n1 = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n1);
    vec_init_generic(mmt->pi->wr[d], abase, v, flags, n1);
}

void matmul_top_vec_clear_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;
    unsigned int n1 = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    matmul_aux(mmt->mm, MATMUL_AUX_GET_READAHEAD, &n1);
    vec_clear_generic(mmt->pi->wr[d], v, n1);
}
/* }}} */
/* {{{ vector init/clear (top level) */
void matmul_top_vec_init(matmul_top_data_ptr mmt, int d, int flags)
{
    matmul_top_vec_init_generic(mmt, NULL, NULL, d, flags);
}
void matmul_top_vec_clear(matmul_top_data_ptr mmt, int d)
{
    matmul_top_vec_clear_generic(mmt, NULL, d);
}
/* }}} */

/* {{{ broadcast_down (generic interface) */
/* broadcast_down reads data in mmt->wr[d]->v, and broadcasts it across the
 * communicator mmt->pi->wr[d] ; eventually everybody on the communicator
 * mmt->pi->wr[d] has the data.
 *
 * Note that for shuffled product, the combination of reduce_across +
 * broadcast_down is not the identity.
 */
static void
broadcast_down_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;

    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif
    int err;

    mmt_wiring_ptr mcol = mmt->wr[d];
    // mmt_wiring_ptr mrow = mmt->wr[!d];
    
    pi_wiring_ptr picol = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
#endif

    unsigned int ncols = mmt->n[d];
    unsigned int nrows = mmt->n[!d];

    ASSERT(mcol->i0 % picol->totalsize == 0);
    ASSERT(mcol->i1 % picol->totalsize == 0);
    unsigned int z = ncols / mmt->pi->m->totalsize;
    ASSERT(nrows == ncols);
    ASSERT(mcol->i0 % z == 0);
    ASSERT(mcol->i1 % z == 0);
    unsigned int nv = pirow->totalsize;
    unsigned int nh = picol->totalsize;
    // unsigned int cz = ncols / nv; ASSERT(cz == z * nh);
    unsigned int rz = nrows / nh; ASSERT(rz == z * nv);

    /* The intersection of our column-wise input range [i0..i1[ with
     * all the fences existing for the output ranges is described by
     * mcol->x ; note how this intersection is common to a
     * complete column.
     */
    /* We used to support (prior to 2011/01/19) the source vector *not*
     * being shared among threads. This has never been actively checked,
     * so despite the fact that the code probably used to work at some
     * point, it was most probably buggy and has been removed.
     * Furthermore, the corresponding code has not been ported to support
     * shuffled mult. The corresponding code has thus been removed.
     */
    // ASSERT_ALWAYS(v->flags & THREAD_SHARED_VECTOR);

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);
    /* Make sure that no thread on the column is wandering in other
     * places -- when we're leaving reduce, this is important. */
    serialize_threads(picol);

    /* This loop suffers from two-dimensional serializing, so the
     * simplistic macros SEVERAL_THREADS_PLAY_MPI_BEGIN and END do not
     * help here.
     */


#ifndef MPI_LIBRARY_MT_CAPABLE
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
            size_t eblock = (mcol->i1 - mcol->i0) / picol->totalsize;
            if (picol->trank) {
                int pos = picol->jrank * picol->ncores + picol->trank;
                int offset = pos * eblock;
                v->abase->vec_set(v->abase,
                        SUBVEC(v, all_v[0], offset),
                        SUBVEC(v, v, offset),
                        eblock);
            }
            serialize_threads(picol);
        }
        if (picol->trank == 0) {
            for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
                pi_log_op(pirow, "[%s] serialize_threads", __func__);
                serialize_threads(pirow);
                pi_log_op(pirow, "[%s] serialize_threads done", __func__);
                if (t != pirow->trank)
                    continue;   // not our turn.
                // although the openmpi man page looks funny, I'm assuming that
                // MPI_Allgather wants MPI_IN_PLACE as a sendbuf argument.
                pi_log_op(picol, "[%s] MPI_Allgatherv", __func__);
                ASSERT((mcol->i1 - mcol->i0) % picol->totalsize == 0);
                size_t eblock = (mcol->i1 - mcol->i0) /  picol->totalsize;
                err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, v->v, v->stride * eblock * picol->ncores, MPI_BYTE, picol->pals);
                pi_log_op(picol, "[%s] MPI_Allgatherv done", __func__);
                ASSERT_ALWAYS(!err);
            }
        }
        if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
            serialize_threads(picol);
            if (picol->trank) {
                v->abase->vec_set(v->abase,
                        v->v, v->all_v[0], mcol->i1 - mcol->i0);
            }
        }
    } else {
        for(unsigned int t = 0 ; t < pirow->ncores ; t++) {
            pi_log_op(pirow, "[%s] serialize_threads", __func__);
            serialize_threads(pirow);
            pi_log_op(pirow, "[%s] serialize_threads done", __func__);
            if (t != pirow->trank) 
                continue;
            for(unsigned int u = 0 ; u < picol->ncores ; u++) {
                serialize_threads(picol);
                if (u != picol->trank)
                    continue;
                // our turn.
                for(unsigned int ii = mcol->i0 ; ii < mcol->i1 ; ) {
                    unsigned int znum = ii / z;
                    unsigned int k = znum / nv;
                    unsigned int offset_there = (znum % nv) * z;
                    unsigned int count = MIN(mcol->i1 - ii, rz - offset_there);
                    unsigned int offset_me = ii - mcol->i0;
                    /*
                       ASSERT(l < mcol->xlen);
                       ASSERT(mcol->x[l].k == (int) k);
                       ASSERT(mcol->x[l].offset_me == offset_me);
                       ASSERT(mcol->x[l].offset_there == offset_there);
                       ASSERT(mcol->x[l].count == count);
                       */
                    if (k % picol->ncores == picol->trank) {
                        void * ptr = SUBVEC(v, v, offset_me);
                        // XXX This would morally be something like:
                        // v->abase->mpi_bcast(ptr, count, root, picol->pals)
                        // or intermediately something which at leasts
                        // uses the datatype.
                        // size_t siz = count * v->stride;
                        unsigned int root = k / picol->ncores;
                        pi_log_op(picol, "[%s] MPI_Bcast", __func__);
                        // err = MPI_Bcast(ptr, siz, MPI_BYTE, root, picol->pals);
                        err = MPI_Bcast(ptr, count, v->abase->mpi_datatype(v->abase), root, picol->pals);
                        pi_log_op(picol, "[%s] MPI_Bcast done", __func__);
                        ASSERT_ALWAYS(!err);
                    }
                    ii += count;
                }
            }
        }
        if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
            serialize_threads(picol);
            for(unsigned int ii = mcol->i0 ; ii < mcol->i1 ; ) {
                unsigned int znum = ii / z;
                unsigned int k = znum / nv;
                unsigned int offset_there = (znum % nv) * z;
                unsigned int count = MIN(mcol->i1 - ii, rz - offset_there);
                unsigned int offset_me = ii - mcol->i0;
                unsigned int root = k % picol->ncores;
                const void * sptr = SUBVEC(v, all_v[root], offset_me);
                void * dptr = SUBVEC(v, v, offset_me);
                v->abase->vec_set(v->abase, dptr, sptr, count);
                // if (k / picol->ncores == picol->jrank) { }
                ii += count;
            }
        }
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    pi_log_op(mmt->pi->m, "[%s] trailer", __func__);

    ASSERT_ALWAYS(mcol->i1 <= nrows);
}
/* }}} */

/* {{{ backends for load/save */
/* {{{ save, backend */
/* It really something relevant to the pirow communicator. Turn it so */
static void save_vector_toprow_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    if (v == NULL) v = mmt->wr[d]->v;

    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    size_t mysize = v->abase->vec_elt_stride(v->abase, mcol->i1 - mcol->i0);
    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);

    int rc = pi_save_file(pirow, name, iter, v->v, mysize, sizeondisk);

    if (rc == 0) fprintf(stderr, "WARNING: checkpointing failed\n");
}
/* }}} */
// now the exact opposite.

/* {{{ load, backend */
/* The backend is exactly dual to save_vector_toprow_generic.  */
static void load_vector_toprow_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    if (v == NULL) v = mmt->wr[d]->v;

    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    size_t mysize = v->abase->vec_elt_stride(v->abase, mcol->i1 - mcol->i0);
    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);

    int rc = pi_load_file(pirow, name, iter, v->v, mysize, sizeondisk);

    if (rc == 0) fprintf(stderr, "WARNING: checkpointing failed\n");
}
/* }}} */
/* }}} */

/* {{{ generic interfaces for load/save */
/* {{{ load (generic) */
void matmul_top_load_vector_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    int err;
    if (v == NULL) v = mmt->wr[d]->v;

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
#endif
    mmt_wiring_ptr mcol = mmt->wr[d];
    // mmt_wiring_ptr mrow = mmt->wr[!d];

    if (picol->jrank == 0 && picol->trank == 0) {
        load_vector_toprow_generic(mmt, v, name, d, iter, itemsondisk);
    }

    if (picol->jrank == 0) {
        serialize_threads(mmt->pi->m);
    }
    if (picol->trank == 0) {
        SEVERAL_THREADS_PLAY_MPI_BEGIN(pirow) {
            int err = MPI_Barrier(picol->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END;
    }
    // serialize_threads(picol);
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    // after loading, we need to broadcast the data. As in other
    // occcasions, this has to be one one thread at a time unless the MPI
    // library can handle concurrent calls along parallel communicators.

    if (picol->trank == 0) {
        // first, down on columns.
        SEVERAL_THREADS_PLAY_MPI_BEGIN(pirow) {
            err = MPI_Bcast(v->v, mcol->i1 - mcol->i0, v->abase->mpi_datatype(v->abase), 0, picol->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END;
    }

    serialize(mmt->pi->m);

    // then maybe amongst threads if they're not sharing data.
    if ((v->flags & THREAD_SHARED_VECTOR) == 0) {
        // we have picol->ncores threads, which would be delighted to get
        // access to some common data. The data is known to sit in thread
        // zero, so that's relatively easy.

        if (picol->trank) {
            v->abase->vec_set(v->abase, v->v, v->all_v[0], mcol->i1 - mcol->i0);
        }
    }
}
/* }}} */
/* {{{ save (generic) */
void matmul_top_save_vector_generic(matmul_top_data_ptr mmt, mmt_vec_ptr v, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    if (v == NULL) v = mmt->wr[d]->v;

    // we want row 0 to have everything.
    broadcast_down_generic(mmt, v, d);

    // we'll do some funky stuff, so serialize in order to avoid jokes.
    // Penalty is insignificant anyway since I/O are rare.
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    pi_wiring_ptr picol = mmt->pi->wr[d];
    // pi_wiring_ptr pirow = mmt->pi->wr[!d];

    // Now every job row has the complete vector.  We'll do I/O from only
    // one row, that's easier. Pick the topmost one.
    if (picol->jrank == 0 && picol->trank == 0) {
        save_vector_toprow_generic(mmt, v, name, d, iter, itemsondisk);
    }

    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);
}
/* }}} */
/* }}} */

//////////////////////////////////////////////////////////////////////////

#if 0/*{{{*/
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source_generic(matmul_top_data_ptr mmt, size_t stride, mmt_vec_ptr v, int d)
{
    if (v == NULL) v = mmt->wr[d]->v;

    // In conjugation mode, it is possible to fill exactly the data chunk
    // that will eventually be relevant. However, it's easy enough to
    // fill our output vector with garbage, and do broadcast_down
    // afterwards...
    if ((v->flags & THREAD_SHARED_VECTOR) == 0 || mmt->pi->wr[d]->trank == 0)
        abase_generic_random(stride, v->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    broadcast_down_generic(mmt, stride, v, d);
}
#endif/*}}}*/

/* XXX
 * ``across'' (horizontally) and ``down'' (vertically) are just here for
 * exposition. The binding of this operation to job/thread arrangement
 * is flexible, through argument d.
 * XXX
 */
void
broadcast_down(matmul_top_data_ptr mmt, int d)
{
    broadcast_down_generic(mmt, NULL, d);
}

/* {{{ reduce_across */
#ifdef USE_ALTERNATIVE_REDUCE_SCATTER
void alternative_reduce_scatter(matmul_top_data_ptr mmt, mmt_vec_ptr v, int eitems, int d)
{
    pi_wiring_ptr wr = mmt->pi->wr[d];
    int njobs = wr->njobs;
    int rank = wr->jrank;
    MPI_Datatype t = v->abase->mpi_datatype(v->abase);

    /* If the rsbuf[] buffers have not yet been allocated, it is time to
     * do so now. We also take the opportunity to possibly re-allocate
     * them if because of a larger abase, the corresponding storage has
     * to be expanded.
     */
    size_t needed = v->abase->vec_elt_stride(v->abase, 
        mmt->n[d] / mmt->pi->m->totalsize * mmt->pi->wr[d]->ncores);
    if (mmt->wr[d]->rsbuf_size < needed) {
        mmt->wr[d]->rsbuf[0] = malloc(needed);
        mmt->wr[d]->rsbuf[1] = malloc(needed);
        mmt->wr[d]->rsbuf_size = needed;
    }

    void *b[2];
    b[0] = mmt->wr[d]->rsbuf[0];
    b[1] = mmt->wr[d]->rsbuf[1];
    v->abase->vec_set_zero(v->abase, b[0], eitems);

    int l = (rank + 1) % njobs;
    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    for (int i = 0; i < njobs; i++) {
        int j0, j1;
        j0 = l * eitems; 
        j1 = j0 +  eitems;
        v->abase->vec_add(v->abase, b[0], b[0], SUBVEC(v, all_v[0], j0), j1-j0);
        if (i == njobs - 1)  
            break;
        MPI_Sendrecv(b[0], eitems, t, drank, 0xbeef,
                     b[1], eitems, t, srank, 0xbeef,
                     wr->pals, MPI_STATUS_IGNORE);
        void * tb = b[0];   b[0] = b[1];  b[1] = tb; 
        l = (l + 1) % njobs;
    }
    v->abase->vec_set(v->abase, v->all_v[0], b[0], eitems);
}
#endif  /* USE_ALTERNATIVE_REDUCE_SCATTER */

/* reduce_across reads data in mmt->wr[d]->v, sums it up across the
 * communicator mmt->pi->wr[d], and collects the results in the areas pointed
 * to by mmt->wr[!d]->v
 *
 * Note that for shuffled product, the combination of reduce_across +
 * broadcast_down is not the identity.
 */
void
reduce_across(matmul_top_data_ptr mmt, int d)
{
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif

    int err=0;

    /* reducing across a row is when d == 0 */
    pi_wiring_ptr pirow = mmt->pi->wr[d];
#ifndef MPI_LIBRARY_MT_CAPABLE
    pi_wiring_ptr picol = mmt->pi->wr[!d];
#endif

    mmt_wiring_ptr mrow = mmt->wr[d];
    mmt_wiring_ptr mcol = mmt->wr[!d];
    unsigned int nrows = mmt->n[d];
    unsigned int ncols = mmt->n[!d];

    ASSERT(mrow->i0 % pirow->ncores == 0);
    ASSERT(mrow->i1 % pirow->ncores == 0);
    unsigned int z = nrows / mmt->pi->m->totalsize;
    unsigned int z1 = z * pirow->njobs;
    ASSERT(nrows == ncols);
    ASSERT(mrow->i0 % z1 == 0);
    ASSERT(mrow->i1 % z1 == 0);
    unsigned int nv = pirow->totalsize;
    unsigned int nh = picol->totalsize;
    unsigned int cz = ncols / nv; ASSERT(cz == z * nh);
    // unsigned int rz = nrows / nh; ASSERT(rz == z * nv);

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);

    // I don't think that the case of shared vectors has been tested
    // correctly for reduction.
    // ASSERT_ALWAYS((mmt->wr[d]->v->flags & THREAD_SHARED_VECTOR) == 0);

    if (pirow->ncores > 1 && !(mrow->v->flags & THREAD_SHARED_VECTOR)) {
        /* row threads have to sum up their data. Of course it's
         * irrelevant when there is only one such thread...
         *
         * Concerning locking, we have to make sure that everybody on the
         * row has finished its computation task, but besides that,
         * there's no locking until we start mpi stuff.
         */
        pi_log_op(pirow, "[%s] serialize_threads", __func__);
        serialize_threads(pirow);
        pi_log_op(pirow, "[%s] serialize_threads done", __func__);

        /* Our [i0,i1[ range is split into pirow->ncores parts. This range
         * represent coordinates which are common to all threads on
         * pirow. Corresponding data has to be summed. Amongst the
         * pirow->ncores threads, thread k takes the responsibility of
         * summing data in the k-th block (that is, indices
         * i0+k*(i1-i0)/pirow->ncores and further). As a convention, for
         * the shuffled product, the thread which eventually owns the
         * final data is thread 0. In the ``normal case'', the exact
         * thread which eventually owns the final data is not necessarily
         * the same for the whole range, nor does it carry any
         * relationship with the thread doing the computation !
         *
         * Note that one should consider that data in threads other than
         * the destination thread may be clobbered by the operation
         * (although in the present implementation it is not).
         */
        unsigned int ii0 = mrow->i0 + pirow->trank * z1;
        unsigned int ii1 = ii0 + z1;
        if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
            void * dptr = SUBVEC(mrow->v, all_v[0], ii0 - mrow->i0);
            for(uint32_t w = 1 ; w < pirow->ncores ; w++) {
                const void * sptr = SUBVEC(mrow->v, all_v[w], ii0 - mrow->i0);
                mrow->v->abase->vec_add(mrow->v->abase, dptr, dptr, sptr, ii1-ii0);
            }
        } else {
            for(unsigned int ii = ii0 ; ii < ii1 ; ) {
                unsigned int znum = ii / z;
                unsigned int k = znum / nh;
                unsigned int offset_there = (znum % nh) * z;
                unsigned int count = MIN(ii1 - ii, cz - offset_there);
                unsigned int offset_me = ii - mrow->i0;
                /*
                ASSERT(l < mrow->ylen);
                ASSERT(mrow->y[l].k == (int) k);
                ASSERT(mrow->y[l].offset_me == offset_me);
                ASSERT(mrow->y[l].offset_there == offset_there);
                ASSERT(mrow->y[l].count == count);
                */

                unsigned int dst = k % pirow->ncores;
                /* Note that ``offset_there'' does not count here, since
                 * at this point we're not yet flipping to the buffer in
                 * the other direction.
                 */
                void * dptr = SUBVEC(mrow->v, all_v[dst], offset_me);
                // size_t siz = abbytes(mmt->abase, count);
                ASSERT(offset_me + count <= mrow->i1 - mrow->i0);
                /* j indicates a thread number offset -- 0 means
                 * destination, and so on. all_v has been allocated with
                 * wraparound pointers precisely for accomodating the
                 * hack here. */
                for(unsigned int j = 1 ; j < pirow->ncores ; j++) {
                    /* A given data offset is summed by only one of the row
                     * threads. Otherwise, we're clearly creating rubbish
                     * since sptr and dptr lie within mrow->v
                     */
                    const void * sptr = SUBVEC(mrow->v, all_v[dst+j], offset_me);
                    mrow->v->abase->vec_add(mrow->v->abase, dptr, dptr, sptr, count);
                }
                ii += count;
            }
        }
        pi_log_op(pirow, "[%s] thread reduction done", __func__);
    }

    /* Good. Now on each node, there's one thread (thread 0 in the
     * shuffled case) whose buffer contains the reduced data for the
     * range [i0..i1[ (well, actually, this corresponds to indices
     * distributed in a funny manner in the original matrix, but as far
     * as we care here, it's really the range mrow->i0..mrow->i1
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

    pi_log_op(mmt->pi->m, "[%s] secondary loop", __func__);

    pi_log_op(mmt->pi->m, "[%s] serialize_threads", __func__);
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        serialize_threads(pirow);
    } else {
        serialize_threads(mmt->pi->m);
    }
    pi_log_op(mmt->pi->m, "[%s] serialize_threads done", __func__);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (mmt->bal->h->flags & FLAG_SHUFFLED_MUL) {
        if (pirow->trank == 0) {
            for(unsigned int t = 0 ; t < picol->ncores ; t++) {
                pi_log_op(picol, "[%s] serialize_threads", __func__);
                serialize_threads(picol);
                pi_log_op(picol, "[%s] serialize_threads done", __func__);
                if (t != picol->trank)
                    continue;   // not our turn.

                pi_log_op(pirow, "[%s] MPI_Reduce_scatter", __func__);
                ASSERT((mrow->i1 - mrow->i0) % pirow->totalsize == 0);
#ifdef USE_ALTERNATIVE_REDUCE_SCATTER
                int z = (mrow->i1 - mrow->i0) / pirow->njobs;
                alternative_reduce_scatter(mmt, mrow->v, z, d);
#else   /* USE_ALTERNATIVE_REDUCE_SCATTER */
                void * dptr = mrow->v->all_v[0];
                // all recvcounts are equal
                int * rc = malloc(pirow->njobs * sizeof(int));
                for(unsigned int k = 0 ; k < pirow->njobs ; k++) {
                    rc[k] = (mrow->i1 - mrow->i0) / pirow->njobs;
                }
                int z = (mrow->i1 - mrow->i0) / pirow->njobs;
                err = MPI_Reduce_scatter(dptr, dptr, rc,
                        mrow->v->abase->mpi_datatype(mrow->v->abase),
                        mrow->v->abase->mpi_addition_op(mrow->v->abase),
                        pirow->pals);
                ASSERT_ALWAYS(!err);
                free(rc);
#endif  /* USE_ALTERNATIVE_REDUCE_SCATTER */
                pi_log_op(pirow, "[%s] MPI_Reduce_scatter done", __func__);
            }
        }
        serialize_threads(pirow);
        // row threads pick what they're interested in in thread0's reduced
        // buffer. Writes are non-overlapping in the mcol buffer here.
        // Different row threads always have different mcol buffers, and
        // sibling col threads write to different locations in their
        // (generally shared) mcol buffer, depending on which row they
        // intersect.
        
        // Notice that results are being written inside mrow->v->all_v[0],
        // just packed at the beginning.

        size_t eblock = (mrow->i1 - mrow->i0) / pirow->totalsize;

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
        int pos = picol->jrank * picol->ncores + picol->trank;
        int offset = pos * eblock;
        mrow->v->abase->vec_set(mrow->v->abase,
                SUBVEC(mcol->v, v, offset),
                SUBVEC(mrow->v, all_v[0], pirow->trank * eblock),
                eblock);
    } else {
        /* We need two levels of locking :-( */
        for(unsigned int t = 0 ; t < picol->ncores ; t++) {
            pi_log_op(picol, "[%s] serialize_threads", __func__);
            serialize_threads(picol);
            pi_log_op(picol, "[%s] serialize_threads done", __func__);
            if (t != picol->trank)
                continue;   // not our turn.

            // all row threads are likely to engage in a collective operation
            // at one moment or another.
            for(unsigned int ii = mrow->i0 ; ii < mrow->i1 ; ) {
                unsigned int znum = ii / z;
                unsigned int k = znum / nh;
                unsigned int offset_there = (znum % nh) * z;
                unsigned int count = MIN(mrow->i1 - ii, cz - offset_there);
                unsigned int offset_me = ii - mrow->i0;
                /*
                   ASSERT(l < mrow->xlen);
                   ASSERT(mrow->x[l].k == (int) k);
                   ASSERT(mrow->x[l].offset_me == offset_me);
                   ASSERT(mrow->x[l].offset_there == offset_there);
                   ASSERT(mrow->x[l].count == count);
                   */
                unsigned int tdst = k % pirow->ncores;
                serialize_threads(pirow);
                if (pirow->trank == tdst) {
                    unsigned int jdst = k / pirow->ncores;
                    // MPI_Reduce does not know about const-ness !
                    void * sptr = SUBVEC(mrow->v, v, offset_me);
                    void * dptr = SUBVEC(mcol->v, v, offset_there);
                    pi_log_op(pirow, "[%s] MPI_Reduce", __func__);
                    err = MPI_Reduce(sptr, dptr, count,
                            mrow->v->abase->mpi_datatype(mrow->v->abase),
                            mrow->v->abase->mpi_addition_op(mrow->v->abase), 
                            jdst, pirow->pals);
                    pi_log_op(pirow, "[%s] MPI_Reduce done", __func__);
                    ASSERT_ALWAYS(!err);
                }
                ii+=count;
            }
        }
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif

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
}
/* }}} */

/**********************************************************************/
/* Utility stuff for applying standard permutations to vectors.
 *
 * Depending on whether ``shuffled'' product is used or not, the
 * matmul_top_mul_comm routine, which is part of matmul_top_mul, does not
 * only a reduction+broadcast, but induces multiplication by a
 * permutation matrix P
 *
 * The functions below apply the transformations to the input vector.
 *
 * None of these functions is optimized for speed. Their purpose is only
 * normalization w.r.t permutations being used.
 */

/* {{{ allreduce_across */
/* This is only for convenience now. Eventually this will be relevant for
 * block Lanczos.  Note that allreduce is conceptually much simpler.
 * There is no funny permutation to be considered.
 */
void
allreduce_across(matmul_top_data_ptr mmt, int d)
{
    /* reducing across a row is when d == 0 */
    pi_wiring_ptr pirow = mmt->pi->wr[d];

    mmt_wiring_ptr mrow = mmt->wr[d];
    unsigned int nrows = mmt->n[d];
    unsigned int z = nrows / mmt->pi->m->totalsize;
    unsigned int z1 = z * pirow->njobs;

    pi_log_op(mmt->pi->m, "[%s] enter first loop", __func__);

    serialize_threads(mmt->pi->m);
    /* sum up row threads, so that only one thread on each row is used
     * for communication */
    unsigned int ii0 = mrow->i0 + pirow->trank * z1;
    unsigned int ii1 = ii0 + z1;
    if (!(mrow->v->flags & THREAD_SHARED_VECTOR)) {
        for(unsigned int k = 1 ; k < pirow->ncores ; k++) {
            void * dv = SUBVEC(mrow->v, v, ii0 - mrow->i0);
            void * sv = SUBVEC(mrow->v, all_v[pirow->trank+k], ii0 - mrow->i0);
            mrow->v->abase->vec_add(mrow->v->abase, dv, dv, sv, ii1 - ii0);
        }
    }
    SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
        MPI_Allreduce(MPI_IN_PLACE,
                SUBVEC(mrow->v, v, ii0 - mrow->i0),
                z1,
                mrow->v->abase->mpi_datatype(mrow->v->abase),
                mrow->v->abase->mpi_addition_op(mrow->v->abase), 
                pirow->pals);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
    if (!(mrow->v->flags & THREAD_SHARED_VECTOR)) {
        for(unsigned int k = 1 ; k < pirow->ncores ; k++) {
            void * sv = SUBVEC(mrow->v, v, ii0 - mrow->i0);
            void * dv = SUBVEC(mrow->v, all_v[pirow->trank+k], ii0 - mrow->i0);
            mrow->v->abase->vec_set(mrow->v->abase, dv, sv, ii1 - ii0);
        }
    }
}
/* }}} */

/* {{{ utility: matmul_top_zero_vec_area */
static void matmul_top_zero_vec_area(matmul_top_data_ptr mmt, int d)
{
    mmt_wiring_ptr mdst = mmt->wr[d];
    pi_wiring_ptr pidst = mmt->pi->wr[d];
    if (mdst->v->flags & THREAD_SHARED_VECTOR) {
        serialize_threads(pidst);
        if (pidst->trank == 0)
            mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
        serialize_threads(pidst);
    } else {
        mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
    }
}
/* }}} */
/* {{{ utility: matmul_top_restrict_vec_area */
/* On a shared vector which is assumed to be commonly known by all
 * nodes/threads, select only one portion to remain, and zero out the
 * rest (for later use by e.g.  matmul_top_mul_comm or other integrated
 * function).
 */
static void matmul_top_restrict_vec_area(matmul_top_data_ptr mmt, int d)
{
    mmt_wiring_ptr mdst = mmt->wr[d];
    pi_wiring_ptr pidst = mmt->pi->wr[d];
    if (mdst->v->flags & THREAD_SHARED_VECTOR)
        serialize_threads(pidst);
    if (pidst->trank == 0 || !(mdst->v->flags & THREAD_SHARED_VECTOR)) {
        if (pidst->jrank == 0 && pidst->trank == 0) {
            /* ok, we keep the data */
        } else {
            mdst->v->abase->vec_set_zero(mdst->v->abase, mdst->v->v, mdst->i1 - mdst->i0);
        }
    }
}
/* }}} */

/* {{{ apply_permutation -- this applies the permutation S */
/* This is sort of an equivalent to matmul_top_cpu, although since it is
 * meant to work on permutation matrices, it works equally well even if
 * output buffers are shared.  */
void apply_permutation(matmul_top_data_ptr mmt, int d)
{
    mmt_wiring_ptr mrow = mmt->wr[0];
    mmt_wiring_ptr mcol = mmt->wr[1];

    matmul_top_zero_vec_area(mmt, !d);
    serialize_threads(mmt->pi->wr[d]);
    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][0];
        uint32_t v = mmt->mm->twist[i][1];
        ASSERT_ALWAYS(u >= mrow->i0);
        ASSERT_ALWAYS(u < mrow->i1);
        ASSERT_ALWAYS(v >= mcol->i0);
        ASSERT_ALWAYS(v < mcol->i1);
        if (d==0) {
            mcol->v->abase->vec_set(mcol->v->abase,
                    SUBVEC(mcol->v, v, v - mcol->i0),
                    SUBVEC(mrow->v, v, u - mrow->i0),
                    1);
        } else {
            mrow->v->abase->vec_set(mrow->v->abase,
                    SUBVEC(mrow->v, v, u - mrow->i0),
                    SUBVEC(mcol->v, v, v - mcol->i0),
                    1);
        }
    }
}
/* }}} */

/* {{{ This applies the identity matrix (thus copies the vector from [d]
 * to [!d]) */
void apply_identity(matmul_top_data_ptr mmt, int d)
{
    mmt_wiring_ptr msrc = mmt->wr[d];
    mmt_wiring_ptr mdst = mmt->wr[!d];

    matmul_top_zero_vec_area(mmt, !d);
    serialize_threads(mmt->pi->wr[d]);

    for(unsigned int i = msrc->i0 ; i < msrc->i1 ; i++) {
        if (i >= mdst->i0 && i < mdst->i1) {
            mdst->v->abase->vec_set(mdst->v->abase,
                    SUBVEC(mdst->v, v, i - mdst->i0),
                    SUBVEC(msrc->v, v, i - msrc->i0),
                    1);
        }
    }
}
/* }}} */

/* {{{ matmul_top_unapply_S_unapply_P */
// for nullspace=left, bw->dir == 0 : vectors are represented as
// horizontal blocks n*N
// d == 0       v <- v * Sr * Pr
// d == 1       v <- v * Sr^-1 * Pr^-1
//
// for nullspace=right, bw->dir == 1 : vectors are represented as
// vertical blocks N*n
// d == 1       v <- Pr * Sr * v;
// d == 0       v <- Pr^-1 * Sr^-1 * v;
void matmul_top_unapply_S_unapply_P(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    apply_permutation(mmt, d);
    matmul_top_mul_comm(mmt, d);
}
/* }}} */

/* {{{ matmul_top_apply_P_apply_S */
// for nullspace=left, bw->dir == 0 : vectors are represented as
// horizontal blocks n*N
// d == 0       v <- v * Pr * Sr
// d == 1       v <- v * Pr^-1 * Sr^-1
//
// for nullspace=right, bw->dir == 1 : vectors are represented as
// vertical blocks N*n
// d == 1       v <- Sr * Pr * v;
// d == 0       v <- Sr^-1 * Pr^-1 * v;
void matmul_top_apply_P_apply_S(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    pi_wiring_ptr picol = mmt->pi->wr[d];
    matmul_top_restrict_vec_area(mmt, d);
    serialize_threads(picol);
    matmul_top_mul_comm(mmt, !d);
    apply_permutation(mmt, !d);
    allreduce_across(mmt, d);
}
/* }}} */

/* {{{ matmul_top_apply_P */
void matmul_top_apply_P(matmul_top_data_ptr mmt, int d)
{
    // if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) printf("[%s] d == %d\n", __func__, d);
    matmul_top_restrict_vec_area(mmt, d);
    matmul_top_mul_comm(mmt, !d);
    apply_identity(mmt, !d);
    allreduce_across(mmt, d);
}
/* }}} */

/* {{{ matmul_top_unapply_P */
void matmul_top_unapply_P(matmul_top_data_ptr mmt, int d)
{
    apply_identity(mmt, d);
    matmul_top_mul_comm(mmt, d);
}
/* }}} */

/* {{{ matmul_top_unapply_S */
void matmul_top_unapply_S(matmul_top_data_ptr mmt, int d)
{
    apply_permutation(mmt, d);
    allreduce_across(mmt, !d);
    apply_identity(mmt, !d);
    allreduce_across(mmt, d);
}
/* }}} */

/* {{{ matmul_top_apply_S */
void matmul_top_apply_S(matmul_top_data_ptr mmt, int d)
{
    apply_identity(mmt, d);
    allreduce_across(mmt, !d);
    apply_permutation(mmt, !d);
    allreduce_across(mmt, d);
}
/* }}} */

/* {{{ matmul_top_twist_vector */
void matmul_top_twist_vector(matmul_top_data_ptr mmt, int d)
{
    serialize(mmt->pi->m);
    if (d == 1)
        matmul_top_apply_S(mmt, d);
    else
        matmul_top_unapply_S_unapply_P(mmt, d);
    // I know this looks odd, but I have witnessed inconsistencies when
    // matmul_top_mul is called directly after twist_vector, for the case
    // of on-shared vectors. twist_vector() is rare enough to allow an
    // extra serializing call.
    serialize(mmt->pi->m);
}
/* }}} */

/* {{{ matmul_top_untwist_vector */
void matmul_top_untwist_vector(matmul_top_data_ptr mmt, int d)
{
    serialize(mmt->pi->m);
    if (d == 1)
        matmul_top_unapply_S(mmt, d);
    else
        matmul_top_apply_P_apply_S(mmt, d);
    // I know this looks odd, but I have witnessed inconsistencies when
    // matmul_top_mul is called directly after twist_vector, for the case
    // of on-shared vectors. twist_vector() is rare enough to allow an
    // extra serializing call.
    serialize(mmt->pi->m);
}
/* }}} */

/* {{{ Application of permutations to indices */
static void share_u32_table(pi_wiring_ptr wr, uint32_t * g, unsigned int n)
{
    uint32_t ** allg = NULL;
    if (wr->trank == 0)
        allg = malloc(wr->ncores * sizeof(uint32_t *));
    thread_broadcast(wr, (void**) &allg, 0);
    allg[wr->trank] = g;
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                g[k] += allg[t][k];
            }
        }
        /* We expect one given table cell to contain a non-zero value
         * only at one location. Thus it doesn't really matter whether we
         * do +, XOR, OR */
        int err = MPI_Allreduce(MPI_IN_PLACE, g, n*sizeof(uint32_t), MPI_BYTE, MPI_BXOR, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);
    if (wr->trank)
        memcpy(g, allg[0], n * sizeof(uint32_t));
    serialize_threads(wr);
    if (wr->trank == 0)
        free(allg);
}

#if 0
void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d)
{
    mmt_wiring_ptr mrow = mmt->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];
    unsigned int ncols = mmt->n[d];

    /* rebuild the complete permutation S. Slightly expensive, but if we
     * try to do otherwise, we'll go nuts */
    uint32_t * perm = malloc(ncols * sizeof(uint32_t));
    memset(perm, 0, ncols * sizeof(uint32_t));
    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][0];
        uint32_t v = mmt->mm->twist[i][1];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        perm[u] = v;
    }

    share_u32_table(mmt->pi->m, perm, ncols);

    for(unsigned int j = 0 ; j < n ; j++) {
        xs[j] = perm[xs[j]];
    }
    free(perm);
}
#endif

void indices_apply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d MAYBE_UNUSED)
{
    mmt_wiring_ptr mrow = mmt->wr[0];
    mmt_wiring_ptr mcol = mmt->wr[1];

    uint32_t * r = NULL;
    r = malloc((mrow->i1 - mrow->i0) * sizeof(uint32_t));
    memset(r, 0, (mrow->i1 - mrow->i0) * sizeof(uint32_t));

    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][1];
        uint32_t v = mmt->mm->twist[i][0];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        r[v - mrow->i0] = u;
    }
    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t i = xs[j];
        if (i >= mrow->i0 && i < mrow->i1)
            xs[j] = r[i - mrow->i0];
        else
            xs[j] = 0;
    }
    free(r);
    share_u32_table(mmt->pi->m, xs, n);
}

void indices_unapply_S(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d MAYBE_UNUSED)
{
    mmt_wiring_ptr mrow = mmt->wr[0];
    mmt_wiring_ptr mcol = mmt->wr[1];

    uint32_t * c = NULL;
    c = malloc((mcol->i1 - mcol->i0) * sizeof(uint32_t));
    memset(c, 0, (mcol->i1 - mcol->i0) * sizeof(uint32_t));

    for(uint32_t i = 0 ; i < mmt->mm->ntwists ; i++) {
        uint32_t u = mmt->mm->twist[i][1];
        uint32_t v = mmt->mm->twist[i][0];
        ASSERT_ALWAYS(v >= mrow->i0);
        ASSERT_ALWAYS(v < mrow->i1);
        ASSERT_ALWAYS(u >= mcol->i0);
        ASSERT_ALWAYS(u < mcol->i1);
        c[u - mcol->i0] = v;
    }
    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t i = xs[j];
        if (i >= mcol->i0 && i < mcol->i1)
            xs[j] = c[i - mcol->i0];
        else
            xs[j] = 0;
    }
    free(c);
    share_u32_table(mmt->pi->m, xs, n);
}

void indices_apply_P(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d)
{
    if (! (mmt->bal->h->flags & FLAG_SHUFFLED_MUL)) {
        return;
    }
    pi_wiring_ptr picol = mmt->pi->wr[d];
    pi_wiring_ptr pirow = mmt->pi->wr[!d];
    unsigned int ncols = mmt->n[d];
    unsigned int nv = pirow->totalsize;
    unsigned int nh = picol->totalsize;

    for(unsigned int j = 0 ; j < n ; j++) {
        uint32_t v = xs[j];
        /* Now apply the permutation P to v */
        uint32_t z = ncols / (nh*nv);
        uint32_t q = v / z;
        uint32_t r = v % z;
        uint32_t q1 = q / nv;
        uint32_t q0 = q % nv;
        xs[j] = (q0 * nh + q1) * z + r;
    }
}
/* }}} */

/**********************************************************************/
#if 0
/* no longer used -- was only used by prep.
 * It's not buggy, but making this work in a context where we have
 * multiple threads is tricky.
 */
void matmul_top_fill_random_source(matmul_top_data_ptr mmt, int d)
{
    matmul_top_fill_random_source_generic(mmt, mmt->vr->stride, NULL, d);
}
#endif

#if 0
static void mmt_debug_writeout(matmul_top_data_ptr mmt, int d, const char * name)
{
    // serialize(mmt->pi->m);
    debug_write(mmt->wr[d]->v->v,
            abbytes(mmt->abase, mmt->wr[d]->i1 - mmt->wr[d]->i0),
            "%s.j%u.t%u", name, mmt->pi->m->jrank, mmt->pi->m->trank);
    // serialize(mmt->pi->m);
}
#endif

/* Takes data in mmt->wd[d]->v, and compute the corresponding partial result in
 * mmt->wr[!d]->v.
 */
void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int d)
{
#ifndef NDEBUG
    unsigned int di_in = mmt->wr[d]->i1 - mmt->wr[d]->i0;
    unsigned int di_out = mmt->wr[!d]->i1 - mmt->wr[!d]->i0;
    ASSERT(mmt->mm->dim[1] == (d ? di_in : di_out));
    ASSERT(mmt->mm->dim[0] == (d ? di_out : di_in));
#endif
    mmt_wiring_ptr mrow = mmt->wr[!d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    ASSERT_ALWAYS((mrow->v->flags & THREAD_SHARED_VECTOR) == 0);

    pi_log_op(mmt->pi->m, "[%s] enter matmul_mul", __func__);

    /* Note that matmul_init copies the calling abase argument to the
     * lower-level mm structure. It can quite probably be qualified as a
     * flaw.
     */
    matmul_mul(mmt->mm, mrow->v->v, mcol->v->v, d);
}

/* This takes partial results in the areas mmt->wr[!d]->v, and puts the
 * collected and re-broadcasted results in the areas mmt->wd[d]->v
 *
 * Note that for the shuffled product, this is not equivalent to a trivial
 * operation.
 */
void matmul_top_mul_comm(matmul_top_data_ptr mmt, int d)
{
    pi_wiring_ptr picol = mmt->pi->wr[d];
    mmt_wiring_ptr mcol = mmt->wr[d];

    pi_log_op(mmt->pi->m, "[%s] enter reduce_across", __func__);
    reduce_across(mmt, !d);
    pi_log_op(mmt->pi->m, "[%s] enter broadcast_down", __func__);
    broadcast_down(mmt, d);

    /* If we have shared input data for the column threads, then we'd
     * better make sure it has arrived completely, because while all
     * threads will need the data, only one is actually importing it.
     */
    if (mcol->v->flags & THREAD_SHARED_VECTOR) {
        pi_log_op(picol, "[%s] serialize threads", __func__);
        serialize_threads(picol);
    }
}


/* Vector I/O is done by only one job, one thread. It incurs a
 * significant amount of memory allocation, but this is done relatively
 * rarely.
 *
 * Vectors, as handled within the core routines, are permuted. He have
 * coordinates (v_{\sigma(0)}...v_{\sigma(n-1)}), where \sigma is the
 * permutation given by the *_perm files. Notice that we assume that we
 * have conjugated permutations on the two sides. This means that no
 * data manipulation is required within the critical loops (this would
 * imply a fuzzy communication pattern, boiling down to essentially
 * using allreduce instead of reduce).
 */

/* As always, comments are written with one given point of view in mind.
 * The vector wich gets saved is always the source vector. Here, our
 * arbitrary description choice is that we iterate M^i times vector,
 * hence the source vector is on the right: d == 1.
 */

// comments, variable names and so on are names which for simplicity
// reflect the situation d == 1 (save right vector). The stable state we
// start with is after a matrix-times-vector product (or just before
// one): The left vector has been computed, and the data for indices
// [i0..i1[ is on the nodes seeing those indices vertically as well. Data
// is therefore found in the right vector area, as after the reduce step.
//
// Doing a broadcast_down columns will ensure that each row contains
// the complete data set for our vector.

/* The first function is the back-end, run only on the top row.
 */
void matmul_top_save_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    matmul_top_save_vector_generic(mmt, NULL, name, d, iter, itemsondisk);
}

void matmul_top_load_vector(matmul_top_data_ptr mmt, const char * name, int d, unsigned int iter, unsigned int itemsondisk)
{
    matmul_top_load_vector_generic(mmt, NULL, name, d, iter, itemsondisk);
}



/**********************************************************************/
static void mmt_finish_init(matmul_top_data_ptr mmt, int const *, param_list pl, int optimized_direction);
static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list pl, int optimized_direction);

static void mmt_fill_fields_from_balancing(matmul_top_data_ptr mmt, param_list pl)
{
    mmt->n[0] = mmt->bal->trows;
    mmt->n[1] = mmt->bal->tcols;
    mmt->n0[0] = mmt->bal->h->nrows;
    mmt->n0[1] = mmt->bal->h->ncols;
   
    // nslices[0] is the number of slices in a row -- this is the number of
    // ``vertical slices'', in balance.c parlance.
    unsigned int nslices[2];
    nslices[1] = mmt->bal->h->nh;
    nslices[0] = mmt->bal->h->nv;

    int ok = 1;

    for(int d = 0 ; d < 2 ; d++)
        ok = ok && mmt->pi->wr[d]->totalsize == nslices[d];

    if (!ok) {
        fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                "on-disk data %ux%u\n",
                mmt->pi->wr[1]->njobs,
                mmt->pi->wr[0]->njobs,
                mmt->pi->wr[1]->ncores,
                mmt->pi->wr[0]->ncores,
                nslices[1],
                nslices[0]);
        serialize(mmt->pi->m);
        exit(1);
    }

    // mmt->ncoeffs_total = 0;
    unsigned int ix[2];
    // ix[1] is something between 0 and nslices[1]-1, i.e. between 0 and nh-1
    // ix[0] is something between 0 and nslices[0]-1, i.e. between 0 and nv-1

    // fences[0] has 1+nslices[1] stop, which is nh+1
    // fences[1] has 1+nslices[0] stop, which is nv+1
    //
    // ix[1] and ix[1]+1 are thus indices within fences[0]
    for(int d = 0 ; d < 2 ; d++)  {
        unsigned int * f = malloc((1 + nslices[d])*sizeof(unsigned int));
        // we use UINT_MAX markers for sanity checking.
        memset(f, -1, (1+nslices[d]) * sizeof(unsigned int));
        f[nslices[d]] = mmt->n[d];
        uint32_t quo = mmt->n[d] / nslices[d];
        unsigned int rem = mmt->n[d] % nslices[d];
        for(unsigned int k = 0 ; k < nslices[d] ; k++) {
            f[k] = k * quo + MIN(rem, k);
        }
        mmt->fences[d^1] = f;
    }

    for(int d = 0 ; d < 2 ; d++)  {
        ix[d^1] = mmt->pi->wr[1^d]->jrank * mmt->pi->wr[1^d]->ncores;
        ix[d^1]+= mmt->pi->wr[1^d]->trank;
        mmt->wr[d]->i0 = mmt->fences[d][ix[d^1]];
        mmt->wr[d]->i1 = mmt->fences[d][ix[d^1] + 1];
        ASSERT_ALWAYS(mmt->wr[d]->i0 < mmt->wr[d]->i1);
        ASSERT_ALWAYS(mmt->wr[d]->i1 <= mmt->n[d]);
        ASSERT_ALWAYS((mmt->wr[d]->i1 - mmt->wr[d]->i0) % mmt->pi->wr[d]->totalsize == 0);

        // the mm layer has an ncoeffs field as well, which is used,
        // while to the contrary the one in the mmt structure is totally
        // unused (and is now gone)
        // mmt->ncoeffs = ncoeffs;
        // mmt->ncoeffs_total += ncoeffs;
    }

    char * base;
    {
        const char * mname;
        // get basename of matrix file. Remove potential suffix, .txt or .bin
        mname = param_list_lookup_string(pl, "matrix");
        ASSERT_ALWAYS(mname);
        char * last_slash = strrchr(mname, '/');
        if (last_slash) {
            mname = last_slash + 1;
        }
        const char * suffixes[] = { ".txt", ".bin" };
        unsigned int l = strlen(mname);
        for(unsigned int k = 0 ; k < sizeof(suffixes) / sizeof(suffixes[0]) ; k++) {
            unsigned int ls = strlen(suffixes[k]);
            if (ls > l) continue;
            if (strcmp(mname+l-ls,suffixes[k]) == 0)
                l -= ls;
        }
        base = strndup(mname, l);
    }

    /* The name choice is only indicative. Depending on the value of
     * the nullspace argument, the matrices in h?.v? may in fact be
     * transposed because this is what the cache building code prefers to
     * receive.
     */
    int rc = asprintf(&mmt->locfile, "%s.%dx%d.%08"PRIx32".h%d.v%d", base, nslices[1], nslices[0], mmt->bal->h->checksum, ix[1], ix[0]);
    free(base);
    ASSERT_ALWAYS(rc >= 0);
}

void matmul_top_init(matmul_top_data_ptr mmt,
        abase_vbase_ptr abase,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        int const * flags,
        param_list pl,
        int optimized_direction)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    mmt->pi = pi;
    mmt->mm = NULL;

    if (mmt->pi->m->trank == 0) abase->mpi_ops_init(abase);

    // n[]
    // ncoeffs_total, ncoeffs,
    // locfile,
    // fences[]
    // wr[]->*
    // are filled in later within read_info_file
    //

    const char * tmp = param_list_lookup_string(pl, "balancing");
    if (!tmp) {
        fprintf(stderr, "Missing parameter balancing=\n");
        exit(1);
    }

    if (pi->m->jrank == 0 && pi->m->trank == 0)
        balancing_read_header(mmt->bal, tmp);
    global_broadcast(pi->m, mmt->bal, sizeof(mmt->bal), 0, 0);
    serialize(mmt->pi->m);

    // after that, we need to check for every node if the cache file can
    // be found. If none is found, rebuild the cache files with an
    // equivalent of mf_dobal.
    mmt_finish_init(mmt, flags, pl, optimized_direction);
}

/* Some work has to be done in order to fill the remaining fields in the
 * matmul_top structure.
 */
static void mmt_finish_init(matmul_top_data_ptr mmt, int const * flags, param_list pl, int optimized_direction)
{
    mmt_fill_fields_from_balancing(mmt, pl);
#ifndef  CONJUGATED_PERMUTATIONS
    choke me;
#endif
    // TODO: reuse the ../bw-matmul/matmul/matrix_base.cpp things for
    // displaying communication info.

    if (!mmt->pi->interleaved) {
        matmul_top_read_submatrix(mmt, pl, optimized_direction);
    } else {
        /* Interleaved threads will share their matrix data. The first
         * thread to arrive will do the initialization, and the second
         * one will just grab the pointer. The trick is to be able to
         * pick the pointer in the right location ! We add a generic
         * pointer dictionary feature in the parallelizing_info
         * interface for this purpose.
         */

#define MMT_MM_MAGIC_KEY        0xaa000000UL

        if (mmt->pi->interleaved->idx == 0) {
            matmul_top_read_submatrix(mmt, pl, optimized_direction);
            pi_store_generic(mmt->pi, MMT_MM_MAGIC_KEY, mmt->pi->m->trank, mmt->mm);
        } else {
            mmt->mm = pi_load_generic(mmt->pi, MMT_MM_MAGIC_KEY, mmt->pi->m->trank);
        }
    }

    /* NOTE: flags default to zero */
    matmul_top_vec_init(mmt, 0, flags ? flags[0] : 0);
    matmul_top_vec_init(mmt, 1, flags ? flags[1] : 0);


}

static int export_cache_list_if_requested(matmul_top_data_ptr mmt, param_list pl)
{
    const char * cachelist = param_list_lookup_string(pl, "export_cachelist");
    if (!cachelist) return 0;

    char * myline = NULL;
    int rc;
    rc = asprintf(&myline, "%s %s", mmt->pi->nodename, mmt->mm->cachefile_name);
    ASSERT_ALWAYS(rc >= 0);
    ASSERT_ALWAYS(myline != NULL);
    char ** tlines = NULL;
    if (mmt->pi->m->trank == 0) {
        tlines = malloc(mmt->pi->m->ncores * sizeof(const char *));
    }
    thread_broadcast(mmt->pi->m, (void**) &tlines, 0);
    tlines[mmt->pi->m->trank] = myline;
    serialize_threads(mmt->pi->m);
    int len = 0;
    if (mmt->pi->m->trank == 0) {
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
            int s = strlen(tlines[j]);
            if (s >= len) len = s;
        }
        MPI_Allreduce(MPI_IN_PLACE, &len, 1, MPI_INT, MPI_MAX, mmt->pi->m->pals);
        char * info = malloc(mmt->pi->m->totalsize * (len + 1));
        char * mybuf = malloc(mmt->pi->m->ncores * (len+1));
        memset(mybuf, 0, mmt->pi->m->ncores * (len+1));
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
            memcpy(mybuf + j * (len+1), tlines[j], strlen(tlines[j])+1);
        }
        MPI_Allgather(mybuf, mmt->pi->m->ncores * (len+1), MPI_BYTE,
                info, mmt->pi->m->ncores * (len+1), MPI_BYTE,
                mmt->pi->m->pals);
        if (mmt->pi->m->jrank == 0) {
            FILE * f = fopen(cachelist, "w");
            DIE_ERRNO_DIAG(f == NULL, "fopen", cachelist);
            for(unsigned int j = 0 ; j < mmt->pi->m->njobs ; j++) {
                unsigned int j0 = j * mmt->pi->m->ncores;
                fprintf(f, "get-cache %s", info + j0 * (len + 1));
                for(unsigned int k = 1 ; k < mmt->pi->m->ncores ; k++) {
                    char * t = strchr(info + (j0 + k) * (len + 1), ' ');
                    ASSERT_ALWAYS(t);
                    fprintf(f, "%s", t);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        free(info);
        free(mybuf);
        free(tlines);
    }
    serialize_threads(mmt->pi->m);
    free(myline);
    serialize(mmt->pi->m);

    return 1;
}


static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, param_list pl, int optimized_direction)
{
    const char * impl = param_list_lookup_string(pl, "mm_impl");
    int rebuild = 0;
    param_list_parse_int(pl, "rebuild_cache", &rebuild);

    mmt->mm = matmul_init(mmt->abase,
            mmt->wr[0]->i1 - mmt->wr[0]->i0,
            mmt->wr[1]->i1 - mmt->wr[1]->i0,
            mmt->locfile, impl, pl, optimized_direction); 

    mmt->mm->nslices[1] = mmt->bal->h->nh;
    mmt->mm->nslices[0] = mmt->bal->h->nv;

    // *IF* we need to do a collective read of the matrix, we need to
    // provide the pointer *now*.
    unsigned int sqread = 0;
    param_list_parse_uint(pl, "sequential_cache_read", &sqread);

    int cache_loaded = 0;

    if (export_cache_list_if_requested(mmt, pl)) {
        /* This has the effect of terminating the program if we are being
         * called from dispatch */
        return;
    }

    if (!rebuild) {
        if (sqread) {
            for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
                serialize_threads(mmt->pi->m);
                if (j == mmt->pi->m->trank) {
                    cache_loaded = matmul_reload_cache(mmt->mm);
                }
            }
        } else {
            cache_loaded = matmul_reload_cache(mmt->mm);
        }
    }

    if (!global_data_eq(mmt->pi, &cache_loaded, sizeof(int))) {
        if (mmt->pi->m->trank == 0) {
        if (mmt->pi->m->jrank == 0) {
            fprintf(stderr, "Fatal error: cache files not present at expected locations\n");
        }
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
            fprintf(stderr, "[%s] J%uT%u: cache %s: %s\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank,
                    mmt->pi->m->trank,
                    mmt->mm->cachefile_name,
                    cache_loaded ? "ok" : "not ok");
        }
        SEVERAL_THREADS_PLAY_MPI_END;
        serialize(mmt->pi->m);
        abort();
    }

    unsigned int sqb = 0;
    param_list_parse_uint(pl, "sequential_cache_build", &sqb);

    matrix_u32 m;
    memset(m, 0, sizeof(matrix_u32));
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (!cache_loaded) {
        if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
            fprintf(stderr, "Matrix dispatching starts\n");
        }
        m->mfile = param_list_lookup_string(pl, "matrix");
        m->bfile = param_list_lookup_string(pl, "balancing");
        // the mm layer is informed of the higher priority computations
        // that will take place. Depending on the implementation, this
        // may cause the direct or transposed ordering to be preferred.
        // Thus we have to read this back from the mm structure.
        m->transpose = mmt->mm->store_transposed;
        m->withcoeffs = param_list_lookup_string(pl, "prime") != NULL;
        balancing_get_matrix_u32(mmt->pi, pl, m);

        int ssm = 0;
        param_list_parse_int(pl, "save_submatrices", &ssm);
        if (ssm) {
            fprintf(stderr, "DEBUG: creating %s\n", mmt->locfile);
            FILE * f = fopen(mmt->locfile, "w");
            fwrite(m->p, sizeof(uint32_t), m->size, f);
            fclose(f);
        }
    }


    if (!sqb) {
        if (!cache_loaded) {
            // everybody does it in parallel
            fprintf(stderr,"[%s] J%uT%u building cache for %s\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank,
                    mmt->pi->m->trank,
                    mmt->locfile);
            matmul_build_cache(mmt->mm, m);
            matmul_save_cache(mmt->mm);
        }
    } else {
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores + 1 ; j++) {
            serialize_threads(mmt->pi->m);
            if (cache_loaded) continue;
            if (j == mmt->pi->m->trank) {
                fprintf(stderr,"[%s] J%uT%u building cache for %s\n",
                        mmt->pi->nodenumber_s,
                        mmt->pi->m->jrank,
                        mmt->pi->m->trank,
                        mmt->locfile);
                matmul_build_cache(mmt->mm, m);
            } else if (j == mmt->pi->m->trank + 1) {
                matmul_save_cache(mmt->mm);
            }
        }
    }
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    my_pthread_mutex_lock(mmt->pi->m->th->m);
    fprintf(stderr, "[%s] J%uT%u uses cache file %s\n",
            mmt->pi->nodenumber_s,
            mmt->pi->m->jrank, mmt->pi->m->trank,
            /* cache for mmt->locfile, */
            mmt->mm->cachefile_name);
    my_pthread_mutex_unlock(mmt->pi->m->th->m);
}


void matmul_top_clear(matmul_top_data_ptr mmt)
{
    serialize(mmt->pi->m);

    if (mmt->pi->m->trank == 0) mmt->abase->mpi_ops_clear(mmt->abase);

    matmul_top_vec_clear(mmt,0);
    matmul_top_vec_clear(mmt,1);
#ifdef USE_ALTERNATIVE_REDUCE_SCATTER
    for(int d = 0 ; d < 2 ; d++)  {
        free(mmt->wr[d]->rsbuf[0]); mmt->wr[d]->rsbuf[0] = NULL;
        free(mmt->wr[d]->rsbuf[1]); mmt->wr[d]->rsbuf[1] = NULL;
        mmt->wr[d]->rsbuf_size = 0;
        // both are expected to hold storage for:
        // (mmt->n[d] / mmt->pi->m->totalsize * mmt->pi->wr[d]->ncores))
        // elements, corresponding to the largest abase encountered.
    }
#endif  /* USE_ALTERNATIVE_REDUCE_SCATTER */
    serialize(mmt->pi->m);
    if (!mmt->pi->interleaved) {
        matmul_clear(mmt->mm);
    } else if (mmt->pi->interleaved->idx == 1) {
        matmul_clear(mmt->mm);
        /* group 0 is the first to leave, thus it doesn't to freeing.
         */
    }
    free(mmt->fences[0]);
    free(mmt->fences[1]);
    free(mmt->locfile);
}

