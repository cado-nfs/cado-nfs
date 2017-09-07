#include "cado.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>             // truncate()
#include <sys/types.h>          // truncate()
#include <errno.h>
#include <sys/stat.h>
#include "bwc_config.h"
#include "matmul.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "intersections.h"
#include "balancing_workhorse.h"
#include "portability.h"
#include "misc.h"
#include "random_matrix.h"
#include "cheating_vec_init.h"
#include "mf_bal.h"

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
#define RS_CHOICE_STOCK_IRSBLOCK     (MPI_VERSION_ATLEAST(3,0) ? 3 : -3)
#define RS_CHOICE_MINE               4
#define RS_CHOICE_MINE_DROP_IN       (MPI_VERSION_ATLEAST(3,0) ? 5 : -5)
#define RS_CHOICE_MINE_PARALLEL      6
#define RS_CHOICE_MINE_OVERLAPPING   7  /* TODO */
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
 * RS_CHOICE_STOCK_IRSBLOCK        => seg fault
 * RS_CHOICE_MINE_OVERLAPPING      => to be implemented
 */

/* _this_ part can be configured */
// #define RS_CHOICE RS_CHOICE_STOCK_RS
// #define AG_CHOICE AG_CHOICE_STOCK_AG
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

///////////////////////////////////////////////////////////////////
/* Start with stuff that does not depend on abase at all -- this
 * provides a half-baked interface */

/* At some point we had this. Not sure it's still useful. */
#define ABASE_UNIVERSAL_READAHEAD_ITEMS 8

static permutation_data_ptr permutation_data_alloc();
static void permutation_data_free(permutation_data_ptr a);
static void permutation_data_push(permutation_data_ptr a, unsigned int u[2]);
#ifdef  MVAPICH2_NUMVERSION
static void permutation_data_ensure(permutation_data_ptr a, size_t n);
#endif


void matmul_top_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "matrix",
            "the matrix file (binary)");
    param_list_decl_usage(pl, "balancing",
            "the matrix balancing file, as computed by mf_bal");
    param_list_decl_usage(pl, "random_matrix",
            "characteristics of a random matrix to be used for staged runs.");

    param_list_decl_usage(pl, "rebuild_cache",
            "force rebuilding matrix caches");
    param_list_decl_usage(pl, "export_cachelist",
            "print the per-node needed cache files and exit");
    param_list_decl_usage(pl, "save_submatrices",
            "after dispatching, save a copy of the local uncompressed matrix before creating the cache file");
    param_list_decl_usage(pl, "sequential_cache_build",
            "build the cache files sequentially on each node");
    param_list_decl_usage(pl, "sequential_cache_read",
            "read the cache files sequentially on each node");
    param_list_decl_usage(pl, "balancing_options",
            "options to pass to the balancing subprogram (see mf_bal_adjust_from_option_string)");
    balancing_decl_usage(pl);
    matmul_decl_usage(pl);
}

void matmul_top_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "matrix");
    param_list_lookup_string(pl, "balancing");
    param_list_lookup_string(pl, "random_matrix");
    param_list_lookup_string(pl, "rebuild_cache");
    param_list_lookup_string(pl, "export_cachelist");
    param_list_lookup_string(pl, "save_submatrices");
    param_list_lookup_string(pl, "sequential_cache_build");
    param_list_lookup_string(pl, "sequential_cache_read");
    param_list_lookup_string(pl, "balancing_options");
    balancing_lookup_parameters(pl);
    matmul_lookup_parameters(pl);
}

static int is_char2(mpfq_vbase_ptr abase)
{
    mpz_t p;
    mpz_init(p);
    abase->field_characteristic(abase, p);
    int char2 = mpz_cmp_ui(p, 2) == 0;
    mpz_clear(p);
    return char2;
}


/* Some info about distributed vectors.
 *
 * See linalg/bwc/GUIDED_TOUR_OF_SOURCES for documentation about what's
 * in the mmt_vec type.
 *
 * Here are some pre- and post- conditions about consistency for the
 * common routines:
 *
 * matmul_top_mul_cpu:
 *      input: fully consistent
 *      output: inconsistent
 * reduce (all variants)
 *      input: inconsistent
 *      output: partially consistent
 * allreduce
 *      input: inconsistent
 *      output: fully consistent
 * broadcast
 *      input: partially consistent
 *      output: fully consistent
 * matmul_top_mul_comm:
 *      input: inconsistent
 *      output: fully consistent
 *
 * The following are non-critical. So for easiness, we require and
 * provide full consistency.
 *
 * apply_identity:
 *      input: fully consistent
 *      output: fully consistent
 *
 */


/* {{{ vector init/clear */

/* this is for a vector which will be of interest to a group of threads
 * and jobs in direction d */
void mmt_vec_init(matmul_top_data_ptr mmt, mpfq_vbase_ptr abase, pi_datatype_ptr pitype, mmt_vec_ptr v, int d, int flags, unsigned int n)
{
    ASSERT_ALWAYS(v != NULL);
    if (abase == NULL) abase = mmt->abase;
    if (pitype == NULL) pitype = mmt->pitype;
    memset(v, 0, sizeof(mmt_vec));
    v->pi = mmt->pi;
    v->d = d;
    v->abase = abase;
    v->pitype = pitype;
    v->n = n;

    ASSERT_ALWAYS(n % mmt->pi->m->totalsize == 0);

    pi_comm_ptr wr = mmt->pi->wr[d];
    pi_comm_ptr xwr = mmt->pi->wr[!d];

    /* now what is the size which we are going to allocate locally */
    n /= xwr->totalsize;
    v->i0 = n * (xwr->jrank * xwr->ncores + xwr->trank);
    v->i1 = v->i0 + n;

    /* Look for readahead settings for all submatrices */
    n += ABASE_UNIVERSAL_READAHEAD_ITEMS;
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        matmul_aux(mmt->matrices[i]->mm, MATMUL_AUX_GET_READAHEAD, &n);
    }

    if (flags & THREAD_SHARED_VECTOR) {
        if (wr->trank == 0) {
            cheating_vec_init(abase, &v->v, n);
            abase->vec_set_zero(abase, v->v, n);
        }
        pi_thread_bcast(&v->v, sizeof(void*), BWC_PI_BYTE, 0, wr);
        v->siblings = NULL;
    } else {
        cheating_vec_init(abase, &v->v, n);
        abase->vec_set_zero(abase, v->v, n);
        v->siblings = shared_malloc(wr, wr->ncores * sizeof(void *));
        v->siblings[wr->trank] = v;
    }
    /* Vectors begin initialized to zero, so we have full consistency */
    v->consistency = 2;
    serialize_threads(v->pi->m);

    // pi_log_op(v->pi->m, "Hello, world");
    /* fill wrpals and mpals */
    v->wrpals[0] = shared_malloc(v->pi->wr[0], v->pi->wr[0]->ncores * sizeof(void *));
    v->wrpals[0][v->pi->wr[0]->trank] = v;
    serialize_threads(v->pi->m);
    v->wrpals[1] = shared_malloc(v->pi->wr[1], v->pi->wr[1]->ncores * sizeof(void *));
    v->wrpals[1][v->pi->wr[1]->trank] = v;
    serialize_threads(v->pi->m);
    v->mpals = shared_malloc(v->pi->m, v->pi->m->ncores * sizeof(void *));
    v->mpals[v->pi->m->trank] = v;
    serialize_threads(v->pi->m);

}

void mmt_vec_clear(matmul_top_data_ptr mmt, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    pi_comm_ptr wr = mmt->pi->wr[v->d];
    serialize_threads(wr);
    if (v->rsbuf[0]) free(v->rsbuf[0]);
    if (v->rsbuf[1]) free(v->rsbuf[1]);
    unsigned int n = v->i1 - v->i0;
    /* see above */
    n += ABASE_UNIVERSAL_READAHEAD_ITEMS;
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        matmul_aux(mmt->matrices[i]->mm, MATMUL_AUX_GET_READAHEAD, &n);
    }
    if (v->siblings) {
        cheating_vec_clear(v->abase, &v->v, n);
        shared_free(wr, v->siblings);
    } else {
        if (wr->trank == 0)
            cheating_vec_clear(v->abase, &v->v, n);
    }
    shared_free(v->pi->wr[0], v->wrpals[0]);
    shared_free(v->pi->wr[1], v->wrpals[1]);
    shared_free(v->pi->m, v->mpals);
    memset(v, 0, sizeof(mmt_vec));
}
/* }}} */

/* my "own" offset is the added offset within my locally stored data area
 * which represents the data range I am the owner of. This data range
 * correspond to the index range v->i0 + offset to v->i0 + offset + size
 */
size_t mmt_my_own_offset_in_items(mmt_vec_ptr v)
{
    pi_comm_ptr wr = v->pi->wr[v->d];
    size_t eblock = (v->i1 - v->i0) /  wr->totalsize;
    int pos = wr->jrank * wr->ncores + wr->trank;
    return pos * eblock;
}

size_t mmt_my_own_offset_in_bytes(mmt_vec_ptr v)
{
    return v->abase->vec_elt_stride(v->abase, mmt_my_own_offset_in_items(v));
}

void * mmt_my_own_subvec(mmt_vec_ptr v)
{
    return v->abase->vec_subvec(v->abase, v->v, mmt_my_own_offset_in_items(v));
}

size_t mmt_my_own_size_in_items(mmt_vec_ptr v)
{
    pi_comm_ptr wr = v->pi->wr[v->d];
    size_t eblock = (v->i1 - v->i0) /  wr->totalsize;
    return eblock;
}

size_t mmt_my_own_size_in_bytes(mmt_vec_ptr v)
{
    return v->abase->vec_elt_stride(v->abase, mmt_my_own_size_in_items(v));
}

/* This copies **ONLY** the data we are supposed to own from v to w.
 * mmt_own_vec_set2 is slightly special, since in this case we might in
 * this way be picking data from vector areas owned by other blocks (or
 * writing there). Therefore we convey the info about which vector piece
 * we care about with the first argument z.
 *
 */
void mmt_own_vec_set2(mmt_vec_ptr z, mmt_vec_ptr w, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(z != NULL);
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(w != NULL);
    if (v == w) return;
    ASSERT_ALWAYS(z->d == v->d);
    ASSERT_ALWAYS(z->d == w->d);
    size_t off = mmt_my_own_offset_in_items(z);
    size_t sz = mmt_my_own_size_in_items(z);
    ASSERT_ALWAYS(sz == mmt_my_own_size_in_items(v));
    ASSERT_ALWAYS(sz == mmt_my_own_size_in_items(w));
    v->abase->vec_set(v->abase,
            w->abase->vec_subvec(w->abase, w->v, off),
            v->abase->vec_subvec(v->abase, v->v, off),
            sz);
}
void mmt_own_vec_set(mmt_vec_ptr w, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v->abase == w->abase);
    mmt_own_vec_set2(v, w, v);
    w->consistency = 1;
}
void mmt_full_vec_set(mmt_vec_ptr w, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(w != NULL);
    /* DO **NOT** early-quit when v==w, because we might be calling this
     * with v and w being siblings, maybe equal for one of the threads.
     */
    // same remark as above
    // ASSERT_ALWAYS(v->abase == w->abase);
    ASSERT_ALWAYS(v->d == w->d);
    ASSERT_ALWAYS(mmt_my_own_size_in_items(v) == mmt_my_own_size_in_items(w));
    if (w->siblings) {
        if (w->v != v->v) {
            v->abase->vec_set(v->abase, w->v, v->v, v->i1 - v->i0);
        }
    } else {
        ASSERT_ALWAYS(v->siblings == NULL);
        if (w->v != v->v) {
            if (w->pi->wr[w->d]->trank == 0) {
                v->abase->vec_set(v->abase, w->v, v->v, v->i1 - v->i0);
            }
        }
        serialize_threads(w->pi->wr[w->d]);
    }
    if (v != w)
        w->consistency = v->consistency;
}

void mmt_full_vec_set_zero(mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    if (v->siblings) {
        v->abase->vec_set_zero(v->abase, v->v, v->i1 - v->i0);
    } else {
        serialize_threads(v->pi->wr[v->d]);
        if (v->pi->wr[v->d]->trank == 0)
            v->abase->vec_set_zero(v->abase, v->v, v->i1 - v->i0);
    }
    v->consistency = 2;
    serialize_threads(v->pi->wr[v->d]);
}

void mmt_vec_downgrade_consistency(mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v->consistency == 2);
    size_t erase[2][2];
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
        serialize_threads(v->pi->wr[v->d]);
    if (v->siblings) {
        erase[0][0] = 0;
        erase[0][1] = off;
        erase[1][0] = off + sz;
        erase[1][1] = v->i1 - v->i0;
    } else {
        /* There are no siblings, which means that this vector is shared
         * across all threads in this direction. Let only one thread do
         * the job.
         */
        if (v->pi->wr[v->d]->trank == 0) {
            erase[0][0] = 0;
            /* because we are rank 0, this is the minimal offset for this set
             * of threads */
            erase[0][1] = off;
            erase[1][0] = off + sz * v->pi->wr[v->d]->ncores;
            erase[1][1] = v->i1 - v->i0;
        } else {
            erase[0][0] = 0;
            erase[0][1] = 0;
            erase[1][0] = 0;
            erase[1][1] = 0;
        }
    }
    for(int i = 0 ; i < 2 ; i++) {
        if (erase[i][1] != erase[i][0]) {
            v->abase->vec_set_zero(v->abase,
                    v->abase->vec_subvec(v->abase, v->v, erase[i][0]),
                    erase[i][1] - erase[i][0]);
        }
    }
    v->consistency = 1;
    serialize_threads(v->pi->wr[v->d]);
}


#if 0
// looks utterly bogus, so...
/* On a shared vector which is assumed to be commonly known by all
 * nodes/threads, select only one portion to remain, and zero out the
 * rest (for later use by e.g.  matmul_top_mul_comm or other integrated
 * function).
 */
static void mmt_own_vec_clear_complement(matmul_top_data_ptr mmt, int d)
{
    mmt_comm_ptr mdst = mmt->wr[d];
    pi_comm_ptr pidst = mmt->pi->wr[d];
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
#endif
void mmt_vec_clear_padding(mmt_vec_ptr v, size_t unpadded, size_t padded)
{
    /* This can be applied no matter what the consistency argument says
     * */
    serialize(v->pi->m);
    if (unpadded >= padded) return;

    size_t s0 = unpadded >= v->i0 ? (unpadded - v->i0) : 0;
    size_t s1 = padded >= v->i0 ? (padded - v->i0) : 0;
    s0 = MIN(s0, v->i1 - v->i0);
    s1 = MIN(s1, v->i1 - v->i0);

    if (s1 - s0)
        v->abase->vec_set_zero(v->abase,
                v->abase->vec_subvec(v->abase, v->v, s0), s1-s0);

    serialize(v->pi->m);
}

mmt_vec_ptr mmt_vec_sibling(mmt_vec_ptr v, unsigned int i)
{
    if (v->siblings) {
        return v->siblings[i];
    } else {
        return v;
    }
}

/* {{{ mmt_vec_broadcast (generic interface) */
/* mmt_vec_broadcast reads data in mmt->wr[d]->v, and broadcasts it across the
 * communicator mmt->pi->wr[d] ; eventually everybody on the communicator
 * mmt->pi->wr[d] has the data.
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
mmt_vec_broadcast(mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(v->consistency == 1);

    // broadcasting down columns is for common agreement on a right
    // source vector, once a left output has been computed.
    //
    // broadcasting down a column is when d == 1.
    int err;

    /* communicator wr is in the direction we are broadcasting */
    pi_comm_ptr wr = v->pi->wr[v->d];

    /* communicator xwr is in the other direction */
    pi_comm_ptr xwr = v->pi->wr[!v->d];
    mmt_vec_ptr * xwrpals = v->wrpals[!v->d];

    pi_log_op(v->pi->m, "[%s:%d] enter first loop", __func__, __LINE__);
    /* Make sure that no thread on the column is wandering in other
     * places -- when we're leaving reduce, this is important. */
    serialize_threads(wr);

    /* This loop suffers from two-dimensional serializing, so the
     * simplistic macros SEVERAL_THREADS_PLAY_MPI_BEGIN and END do not
     * help here.
     */

    size_t eblock = mmt_my_own_size_in_items(v);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (v->siblings) {
        /* not shared: begin by collecting everything on thread 0 */
        mmt_own_vec_set2(v, v->siblings[0], v);
    }
    serialize_threads(v->pi->m);
    if (wr->trank == 0 && xwr->trank == 0) {
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
        MPI_Request * req = malloc(xwr->ncores * sizeof(MPI_Request));
#endif  /* AG_CHOICE == AG_CHOICE_STOCK_IAG */

        for(unsigned int t = 0 ; t < xwr->ncores ; t++) {
            // although the openmpi man page looks funny, I'm assuming that
            // MPI_Allgather wants MPI_IN_PLACE as a sendbuf argument.
            pi_log_op(wr, "[%s:%d] MPI_Allgather (round %u)", __func__, __LINE__, t);
#if AG_CHOICE == AG_CHOICE_STOCK_IAG
            err = MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                    xwrpals[t]->v, v->abase->vec_elt_stride(v->abase, 1) * eblock * wr->ncores, MPI_BYTE, wr->pals, &req[t]);
#elif AG_CHOICE == AG_CHOICE_STOCK_AG
            err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, xwrpals[t]->v, v->abase->vec_elt_stride(v->abase, 1) * eblock * wr->ncores, MPI_BYTE, wr->pals);
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
    v->consistency = 2;
    serialize_threads(v->pi->m);
    if (v->siblings) {
        mmt_full_vec_set(v, v->siblings[0]);
        serialize_threads(wr);
    }
#else   /* MPI_LIBRARY_MT_CAPABLE */
    /* Code deleted 20110119, as I've never been able to have enough
     * trust in an MPI implementation to check this */
    ASSERT_ALWAYS(0);
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    pi_log_op(v->pi->m, "[%s:%d] trailer", __func__, __LINE__);
}
/* }}} */

/* {{{ generic interfaces for load/save */
/* {{{ load */
int mmt_vec_load_stream(pi_file_handle f, mmt_vec_ptr v, unsigned int itemsondisk)
{
    ASSERT_ALWAYS(v != NULL);
    serialize(v->pi->m);
    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);
    void * mychunk = mmt_my_own_subvec(v);
    size_t mysize = mmt_my_own_size_in_bytes(v);
    ssize_t s = pi_file_read(f, mychunk, mysize, sizeondisk);
    int ok =  s >= 0 && (size_t) s == sizeondisk;
    v->consistency = ok;
    /* not clear it's useful, but well. */
    if (ok) mmt_vec_broadcast(v);
    serialize_threads(v->pi->m);
    return ok;
}
int mmt_vec_load(mmt_vec_ptr v, const char * filename, unsigned int itemsondisk)
{
    ASSERT_ALWAYS(v != NULL);
    serialize(v->pi->m);

    pi_file_handle f;
    int ok = pi_file_open(f, v->pi, v->d, filename, "rb");
    if (ok) {
        ok = mmt_vec_load_stream(f, v, itemsondisk);
        pi_file_close(f);
    }
    if (!ok) {
        if (v->pi->m->trank == 0 && v->pi->m->jrank == 0) {
            fprintf(stderr, "ERROR: failed to load %s\n", filename);
        }
        if (v->pi->m->trank == 0)
            MPI_Abort(v->pi->m->pals, EXIT_FAILURE);
    }
    serialize_threads(v->pi->m);
    return ok;
}
/* }}} */
/* {{{ save */
int mmt_vec_save_stream(pi_file_handle f, mmt_vec_ptr v, unsigned int itemsondisk)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(v->consistency == 2);
    serialize_threads(v->pi->m);
    size_t sizeondisk = v->abase->vec_elt_stride(v->abase, itemsondisk);
    void * mychunk = mmt_my_own_subvec(v);
    size_t mysize = mmt_my_own_size_in_bytes(v);
    ssize_t s = pi_file_write(f, mychunk, mysize, sizeondisk);
    serialize_threads(v->pi->m);
    return s >= 0 && (size_t) s == sizeondisk;
}
int mmt_vec_save(mmt_vec_ptr v, const char * filename, unsigned int itemsondisk)
{
    serialize_threads(v->pi->m);

    pi_file_handle f;
    int ok = pi_file_open(f, v->pi, v->d, filename, "wb");
    if (ok) {
        ok = mmt_vec_save_stream(f, v, itemsondisk);
        pi_file_close(f);
    }
    if (!ok) {
        if (v->pi->m->trank == 0 && v->pi->m->jrank == 0) {
            fprintf(stderr, "WARNING: failed to save %s\n", filename);
            unlink(filename);
        }
    }
    serialize_threads(v->pi->m);
    return ok;
}
/* }}} */
/* }}} */

void mmt_vec_reduce_mod_p(mmt_vec_ptr v)
{
    /* if it so happens that there's no reduce function defined, it may
     * correspond to a case where we have nothing to do, like over GF(2)
     * -- where unreduced elements are just the same, period */
    if (!v->abase->reduce)
        return;
    void * ptr = mmt_my_own_subvec(v);
    void * tmp;
    v->abase->vec_ur_init(v->abase, &tmp, 1);
    for(size_t i = 0 ; i < mmt_my_own_size_in_items(v) ; i++) {
        v->abase->elt_ur_set_elt(v->abase,
                tmp,
                v->abase->vec_coeff_ptr_const(v->abase, ptr, i));
        v->abase->reduce(v->abase,
                v->abase->vec_coeff_ptr(v->abase, ptr, i),
                tmp);
    }
    v->abase->vec_ur_clear(v->abase, &tmp, 1);
}



//////////////////////////////////////////////////////////////////////////

void matmul_top_mul(matmul_top_data_ptr mmt, mmt_vec * v, struct timing_data * tt)/*{{{*/
{
    /* Do all matrices in turn.
     *
     * We represent M as * M0*M1*..*M_{n-1}.
     *
     * The input vector is v[0], and the result is put in v[0] again.
     *
     * The direction in which we apply the product is given by v[0]->d.
     * For v[0]->d == 0, we do v[0]*M. For v[0]->d==1, we do M*v[0]
     *
     * We use temporaries as follows.
     * For v[0]->d == 0:
     *  v[1] <- v[0] * M0 ; v[1]->d == 1
     *  v[2] <- v[1] * M1 ; v[2]->d == 0
     *  if n is odd:
     *  v[n] <- v[n-1] * M_{n-1} ; v[n]->d == 1
     *  if n is even:
     *  v[0] <- v[n-1] * M_{n-1} ; v[0]->d == 0
     *
     * For v[0]->d == 1:
     *  v[1] <- M0 * v[0] ; v[1]->d == 0
     *  v[2] <- M1 * v[1] ; v[2]->d == 1
     *  if n is odd:
     *  v[n] <- M_{n-1} * v[n-1] ; v[n]->d == 0
     *  if n is even:
     *  v[0] <- M_{n-1} * v[n-1] ; v[0]->d == 1
     *
     * This has the consequence that v must hold exactly n+(n&1) vectors.
     *
     * Appropriate communication operations are run after each step.
     *
     * If tt is not NULL, it should be a timing_data structure holding
     * exactly 4*n timers (or only 4, conceivably). Timers are switched
     * exactly that many times.
     *
     * If the mmt->pi->interleaving setting is on, we interleave
     * computations and communications. We do 2n flips. Communications
     * are forbidden both before and after calls to this function (in the
     * directly adjacent code fragments before the closest flip() call,
     * that is).
     */

    int d = v[0]->d;
    int nmats_odd = mmt->nmatrices & 1;
    int midx = (d ? (mmt->nmatrices - 1) : 0);
    for(int l = 0 ; l < mmt->nmatrices ; l++) {
        mmt_vec_ptr src = v[l];
        int last = l == (mmt->nmatrices - 1);
        int lnext = last && !nmats_odd ? 0 : (l+1);
        mmt_vec_ptr dst = v[lnext];

        ASSERT_ALWAYS(src->consistency == 2);
        matmul_top_mul_cpu(mmt, midx, d, dst, src);
        ASSERT_ALWAYS(dst->consistency == 0);

        timing_next_timer(tt);
        /* now measuring jitter */
        pi_interleaving_flip(mmt->pi);
        serialize(mmt->pi->m);
        timing_next_timer(tt);

        /* Now we can resume MPI communications. */
        if (last && nmats_odd) {
            ASSERT_ALWAYS(lnext == mmt->nmatrices);
            matmul_top_mul_comm(v[0], dst);
        } else {
            mmt_vec_allreduce(dst);
        }
        timing_next_timer(tt);
        /* now measuring jitter */
        pi_interleaving_flip(mmt->pi);
        serialize(mmt->pi->m);

        timing_next_timer(tt);
        midx += d ? -1 : 1;
    }
    ASSERT_ALWAYS(v[0]->consistency == 2);
}
/*}}}*/

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
    // fill our output vector with garbage, and do mmt_vec_broadcast
    // afterwards...
    if ((v->flags & THREAD_SHARED_VECTOR) == 0 || mmt->pi->wr[d]->trank == 0)
        mpfq_generic_random(stride, v->v, mmt->wr[d]->i1 - mmt->wr[d]->i0);

    // reconcile all cells which correspond to the same vertical block.
    mmt_vec_broadcast(mmt, v, d);
}
#endif/*}}}*/

/* {{{ mmt_vec_reduce */
/* {{{ various reduce_scatter implementations */
/* {{{ alternative_reduce_scatter */
#if 1 || RS_CHOICE == RS_CHOICE_MINE
void alternative_reduce_scatter(mmt_vec_ptr v)
{
    pi_comm_ptr wr = v->pi->wr[v->d];
    int njobs = wr->njobs;
    int rank = wr->jrank;
    MPI_Datatype t = v->pitype->datatype;

    size_t eitems = mmt_my_own_size_in_items(v) * wr->ncores;
    size_t needed = v->abase->vec_elt_stride(v->abase, eitems);

    if (v->rsbuf_size < needed) {
        ASSERT_ALWAYS(v->rsbuf_size == 0);
        v->rsbuf[0] = realloc(v->rsbuf[0], needed);
        v->rsbuf[1] = realloc(v->rsbuf[1], needed);
        v->rsbuf_size = needed;
    }

    void *b[2];
    b[0] = v->rsbuf[0];
    b[1] = v->rsbuf[1];
    v->abase->vec_set_zero(v->abase, b[0], eitems);

    int l = (rank + 1) % njobs;
    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    for (int i = 0; i < njobs; i++) {
        int j0, j1;
        j0 = l * eitems; 
        j1 = j0 +  eitems;
        v->abase->vec_add(v->abase, b[0], b[0],
                v->abase->vec_subvec(v->abase, mmt_vec_sibling(v, 0)->v, j0),
                j1-j0);
        if (i == njobs - 1)  
            break;
        MPI_Sendrecv(b[0], eitems, t, drank, (i<<16) + rank,
                     b[1], eitems, t, srank, (i<<16) + srank,
                     wr->pals, MPI_STATUS_IGNORE);
        void * tb = b[0];   b[0] = b[1];  b[1] = tb; 
        l = (l + 1) % njobs;
    }
    /* This is going at the beginning of the memory area, which is weird.
     * But it's the convention used by the MPI call...
     */
    v->abase->vec_set(v->abase, mmt_vec_sibling(v, 0)->v, b[0], eitems);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE */
/* }}} */
/* {{{ alternative_reduce_scatter_parallel */
#if 1 || RS_CHOICE == RS_CHOICE_MINE_PARALLEL
/* Example data for a factoring matrix (rsa100) of size 135820*135692,
 * split over 2x3 mpi jobs, and 7x5 threads.
 *
 * 2 == mmt->pi->wr[1]->njobs (number of jobs encountered on a vertical axis).
 * 7 == mmt->pi->wr[1]->ncores (number of jobs per core on a vertical axis).
 * 3 == mmt->pi->wr[0]->njobs (number of jobs encountered on an horiz. axis).
 * 5 == mmt->pi->wr[0]->ncores (number of jobs per core on an horiz. axis).
 *
 * matrix is padded to a multiple of 210 = 2*3*5*7, which is * N=135870=210*647
 *
 * we'll call 647 the "small chunk" size.
 *
 * for all jobs/threads, the following relations hold:
 *      mmt->wr[0]->i1 - mmt->wr[0]->i0 == N/14 == 9705 == 15 * 647
 *      mmt->wr[1]->i1 - mmt->wr[1]->i0 == N/15 == 9058 == 14 * 647
 *
 * a mmt_vec_reduce operation, in the context of factoring, is with d==1
 * below. Hence, in fact, we're doing a reduction down a column.
 *
 * the value eitems fed to this function is mmt->pi->wr[d]->ncores (here,
 * 7) times the small chunk size. Here 7*647 == 4529.
 */
/* all threads in mmt->wr[!d], one after another a priori, are going to
 * do alternative_reduce_scatter on their vector v[i]
 */
void alternative_reduce_scatter_parallel(pi_comm_ptr xr, mmt_vec_ptr * vs)
{
    /* we write all data counts below with comments indicating the typical
     * size in our toy example above */

    /* he have xr->ncores vectors. The vs[] array accesses data from the
     * other peers. Note that the pi structures belong separately to each
     * peer, and so do the embedded communicators. Therefore the proper
     * way to see "our" xr is none other that the given parameter. And for
     * "our" wr, then we'll have to look for our own vector within vs.
     */
    mmt_vec_ptr v = vs[xr->trank];
    mpfq_vbase_ptr ab = v->abase;
    pi_comm_ptr wr = v->pi->wr[v->d];  /* 2 jobs, 7 cores */
    /* what we're going to do will happen completely in parallel over
     * xr->njobs==3 communicators. We're simulating here the split into 5
     * communicators, even though those collide into a unique
     * communicator as far as MPI is concerned.
     */
    int njobs = wr->njobs;      /* 2 */
    /* note that we no longer care at this point about what happened at
     * the thread level in our dimension. This is already done, period.
     */
    int rank = wr->jrank;
    MPI_Datatype t = v->pitype->datatype;

    /* If the rsbuf[] buffers have not yet been allocated, it is time to
     * do so now. We also take the opportunity to possibly re-allocate
     * them if because of a larger abase, the corresponding storage has
     * to be expanded.
     */
    size_t eitems = mmt_my_own_size_in_items(v) * wr->ncores;
    size_t needed = v->abase->vec_elt_stride(v->abase, eitems);

    /* notice that we are allocating a temp buffer only for one vector.
     * Of course, since this is a multithreaded routine, each thread in
     * xr is doing so at the same time */
    if (v->rsbuf_size < needed) {
        v->rsbuf[0] = realloc(v->rsbuf[0], needed);
        v->rsbuf[1] = realloc(v->rsbuf[1], needed);
        v->rsbuf_size = needed;
    }

    ab->vec_set_zero(ab, v->rsbuf[0], eitems);
    serialize_threads(xr);

    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

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

    for (int i = 0, s = 0; i < njobs; i++, s=!s) {
        int l = (rank + 1 + i) % njobs;
        int j0 = l * eitems; 
        int j1 = j0 +  eitems;
        ab->vec_add(ab, v->rsbuf[s], v->rsbuf[s],
                ab->vec_subvec(ab, mmt_vec_sibling(v, 0)->v, j0),
                j1-j0);
        serialize_threads(xr);

        if (i == njobs - 1)  
            break;

        if (xr->trank == 0) {
            MPI_Request * r = NULL;
            r = malloc(2 * xr->ncores * sizeof(MPI_Request));
            for(unsigned int w = 0 ; w < xr->ncores ; w++) {
                MPI_Request * rs = r + 2*w;
                MPI_Request * rr = r + 2*w + 1;
                MPI_Isend(vs[w]->rsbuf[s],  eitems, t, drank, 0xb00+w, wr->pals, rs);
                MPI_Irecv(vs[w]->rsbuf[!s], eitems, t, srank, 0xb00+w, wr->pals, rr);
                /*
                MPI_Sendrecv(vs[w]->rsbuf[s], eitems, t, drank, 0xbeef,
                        vs[w]->rsbuf[!s], eitems, t, srank, 0xbeef,
                        wr->pals, MPI_STATUS_IGNORE);
                        */
                // MPI_Waitall(2, r + 2*w, MPI_STATUSES_IGNORE);
            }
            MPI_Waitall(2 * xr->ncores, r, MPI_STATUSES_IGNORE);
            free(r);
        }
        serialize_threads(xr);
    }

    ab->vec_set(ab, mmt_vec_sibling(v, 0)->v, v->rsbuf[(njobs-1)&1], eitems);

    pi_log_op(wr, "[%s:%d] MPI_Reduce_scatter done", __func__, __LINE__);
}
#endif  /* RS_CHOICE == RS_CHOICE_MINE_PARALLEL */
/* }}} */
/* {{{ my_MPI_Reduce_scatter_block */
#if 1 || RS_CHOICE == RS_CHOICE_MINE_DROP_IN
int my_MPI_Reduce_scatter_block(void *sendbuf, void *recvbuf, int recvcount,
                MPI_Datatype datatype, MPI_Op op, MPI_Comm wr)
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
    static size_t rsbuf_size = 0;
    static unsigned long * rsbuf[2];

    size_t needed = recvcount * tsize;

    if (rsbuf_size < needed) {
        rsbuf[0] = realloc(rsbuf[0], needed);
        rsbuf[1] = realloc(rsbuf[1], needed);
        rsbuf_size = needed;
    }

    memset(rsbuf[0], 0, recvcount * tsize);

    int srank = (rank + 1) % njobs;
    int drank = (rank + njobs - 1) % njobs;

    for (int i = 0, w = 0; i < njobs; i++, w^=1) {
        int j0 = ((rank + i + 1) % njobs) * recvcount;
        const void * share = pointer_arith_const(v, j0 * tsize);
#if MPI_VERSION_ATLEAST(2,2)
        MPI_Reduce_local(share, rsbuf[w], recvcount, datatype, op);
#else
        {
            ASSERT_ALWAYS(datatype == MPI_UNSIGNED_LONG);
            ASSERT_ALWAYS(op == MPI_BXOR);
            const unsigned long * a = share;
            unsigned long * b = rsbuf[w];
            for(int k = 0 ; k < recvcount ; k++) {
                b[k] ^= a[k];
            }
        }
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
 * across the communicator mmt->pi->wr[d], and collects the results in
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
void
mmt_vec_reduce_inner(mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(v->consistency != 2);

    /* reducing across a row is when d == 0 */
    pi_comm_ptr wr = v->pi->wr[v->d];
    pi_comm_ptr xr = v->pi->wr[!v->d];

    pi_log_op(v->pi->m, "[%s:%d] enter first loop", __func__, __LINE__);

    // I don't think that the case of shared vectors has been tested
    // correctly for reduction. Well, to start with, I doubt it really
    // makes a lot of sense anyhow.
    // ASSERT_ALWAYS((v->flags & THREAD_SHARED_VECTOR) == 0);

    if (wr->ncores > 1 && v->siblings) {
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
        size_t thread_chunk = wr->njobs * mmt_my_own_size_in_items(v);
        void * dptr = v->abase->vec_subvec(v->abase,
                mmt_vec_sibling(v, 0)->v, 
                wr->trank * thread_chunk);
        for(unsigned int w = 1 ; w < wr->ncores ; w++) {
            const void * sptr = v->abase->vec_subvec(v->abase,
                    mmt_vec_sibling(v, w)->v, 
                    wr->trank * thread_chunk);
            v->abase->vec_add(v->abase, dptr, dptr, sptr, thread_chunk);
        }
        pi_log_op(wr, "[%s:%d] thread reduction done", __func__, __LINE__);
    }

    /* Good. Now on each node, thread 0 has the reduced data for the
     * range [i0..i1[ (well, actually, this corresponds to indices
     * distributed in a funny manner in the original matrix, but as far
     * as we care here, it's really the range v->i0..v->i1
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

    pi_log_op(v->pi->m, "[%s:%d] secondary loop", __func__, __LINE__);

    pi_log_op(v->pi->m, "[%s:%d] serialize_threads", __func__, __LINE__);
    serialize_threads(wr);
    pi_log_op(v->pi->m, "[%s:%d] serialize_threads done", __func__, __LINE__);

#ifndef MPI_LIBRARY_MT_CAPABLE
    if (wr->trank == 0) {
        /* openmpi-1.8.2 does not seem to have a working non-blocking
         * reduce_scatter, at least not a very efficient one. All
         * I've been able to do is to run MPI_Ireduce_scatter with
         * MPI_UNSIGNED_LONG and MPI_BXOR. With custom types it seems
         * to crash. And anyway, it's very inefficient.
         */

        ASSERT((v->i1 - v->i0) % wr->totalsize == 0);
#if RS_CHOICE == RS_CHOICE_STOCK_RS 
        void * dptr = mmt_vec_sibling(v, 0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            // all recvcounts are equal
            int * rc = malloc(wr->njobs * sizeof(int));
            for(unsigned int k = 0 ; k < wr->njobs ; k++)
                rc[k] = (v->i1 - v->i0) / wr->njobs;
            int err = MPI_Reduce_scatter(dptr, dptr, rc,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            free(rc);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_RSBLOCK
        void * dptr = mmt_vec_sibling(v, 0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = MPI_Reduce_scatter_block(dptr, dptr,
                    (v->i1 - v->i0) / wr->njobs,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_MINE_DROP_IN
        void * dptr = mmt_vec_sibling(v, 0)->v;
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = my_MPI_Reduce_scatter_block(MPI_IN_PLACE, dptr,
                    (v->i1 - v->i0) / wr->njobs,
                    v->pitype->datatype,
                    BWC_PI_SUM->custom,
                    wr->pals);
            ASSERT_ALWAYS(!err);
        }
        SEVERAL_THREADS_PLAY_MPI_END();
#elif RS_CHOICE == RS_CHOICE_STOCK_IRSBLOCK
        void * dptr = mmt_vec_sibling(v, 0)->v;
        MPI_Request * req = shared_malloc(xr, xr->ncores * sizeof(MPI_Request));
        SEVERAL_THREADS_PLAY_MPI_BEGIN(xr) {
            int err = MPI_Ireduce_scatter(dptr, dptr,
                    (v->i1 - v->i0) / wr->njobs,
                    v->pitype->datatype,
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
        shared_free(xr, req);
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
        mmt_vec_ptr * vs = shared_malloc(xr, xr->ncores * sizeof(mmt_vec_ptr));
        vs[xr->trank] = v;
        serialize_threads(xr);
        alternative_reduce_scatter_parallel(xr, vs);
        shared_free(xr, vs);
#elif RS_CHOICE == RS_CHOICE_MINE_OVERLAPPING
#error "not implemented, but planned"
#endif
    }
    serialize_threads(wr);
}
void
mmt_vec_reduce(mmt_vec_ptr w, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(w != NULL);
    ASSERT_ALWAYS(v->abase == w->abase);
    ASSERT_ALWAYS(v->d != w->d);
    mmt_vec_reduce_inner(v);
    pi_comm_ptr wr = v->pi->wr[v->d];
    // row threads pick what they're interested in in thread0's reduced
    // buffer. Writes are non-overlapping in the mcol buffer here.
    // Different row threads always have different mcol buffers, and
    // sibling col threads write to different locations in their
    // (generally shared) mcol buffer, depending on which row they
    // intersect.

    // Notice that results are being written inside v->all_v[0],
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
 
    size_t eblock = mmt_my_own_size_in_items(v);
    ASSERT_ALWAYS(mmt_my_own_size_in_items(w) == eblock);

    v->abase->vec_set(v->abase,
            mmt_my_own_subvec(w),
            /* Note: reduce-scatter packs everything at the beginning in
             * the leader block, which is why we don't have our usual offset
             * here. */
            v->abase->vec_subvec(v->abase,
                mmt_vec_sibling(v, 0)->v, wr->trank * eblock),
            eblock);
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
    w->consistency = 1;
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
mmt_vec_reduce_sameside(mmt_vec_ptr v)
{
    pi_comm_ptr wr = v->pi->wr[v->d];
    mmt_vec_reduce_inner(v);
    size_t eblock = mmt_my_own_size_in_items(v);
    v->abase->vec_set(v->abase,
            mmt_my_own_subvec(v),
            v->abase->vec_subvec(v->abase,
                                mmt_vec_sibling(v, 0)->v, wr->trank * eblock),
            eblock);
    /* This ensures that we've effectively *moved* the data, not copied
     * it */
    if (wr->trank) {
        v->abase->vec_set_zero(v->abase,
                v->abase->vec_subvec(v->abase,
                    mmt_vec_sibling(v, 0)->v, wr->trank * eblock),
                eblock);
    }
    serialize_threads(wr);
    v->consistency = 1;
}
/* }}} */

/* {{{ mmt_vec_allreduce */
/* This is only for convenience now. Eventually this will be relevant for
 * block Lanczos.  Note that allreduce is conceptually much simpler.
 * There is no funny permutation to be considered.
 */
void
mmt_vec_allreduce(mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    ASSERT_ALWAYS(v->consistency != 2);
    /* reducing across a row is when d == 0 */
    pi_comm_ptr wr = v->pi->wr[v->d];

    pi_log_op(v->pi->m, "[%s:%d] enter first loop", __func__, __LINE__);

    serialize_threads(v->pi->m);
    /* sum up row threads, so that only one thread on each row is used
     * for communication */
    size_t thread_chunk = wr->njobs * mmt_my_own_size_in_items(v);
    if (v->siblings) {
        void * dv = v->abase->vec_subvec(v->abase, v->v,
                wr->trank * thread_chunk);
        for(unsigned int k = 1 ; k < wr->ncores ; k++) {
            void * sv = v->abase->vec_subvec(v->abase,
                    mmt_vec_sibling(v, (wr->trank+k) % wr->ncores)->v,
                    wr->trank * thread_chunk);
            v->abase->vec_add(v->abase, dv, dv, sv, thread_chunk);
        }
    }
    /* Compared to the SEVERAL_THREADS_PLAY_MPI_BEGIN() approach, this
     * one has thread 0 do the work for all other threads, while other
     * threads are waiting.
     */
    SEVERAL_THREADS_PLAY_MPI_BEGIN2(v->pi->m, peer) {
        void * dv = v->abase->vec_subvec(v->abase, v->mpals[peer]->v,
                v->mpals[peer]->pi->wr[v->d]->trank * thread_chunk);
        MPI_Allreduce(MPI_IN_PLACE,
                dv,
                thread_chunk,
                v->pitype->datatype,
                BWC_PI_SUM->custom,
                wr->pals);
    }
    SEVERAL_THREADS_PLAY_MPI_END2(v->pi->m);
    if (v->siblings) {
        void * sv = v->abase->vec_subvec(v->abase, v->v,
                wr->trank * thread_chunk);
        for(unsigned int k = 1 ; k < wr->ncores ; k++) {
            void * dv = v->abase->vec_subvec(v->abase,
                    mmt_vec_sibling(v, (wr->trank+k) % wr->ncores)->v,
                    wr->trank * thread_chunk);
            v->abase->vec_set(v->abase, dv, sv, thread_chunk);
        }
    }
    v->consistency = 2;
}
/* }}} */

/**********************************************************************/
/* bench code */

void matmul_top_comm_bench_helper(int * pk, double * pt,
                                  void (*f) (mmt_vec_ptr),
				  mmt_vec_ptr v)
{
    int k;
    double t0, t1;
    int cont;
    t0 = wct_seconds();
    for (k = 0;; k++) {
	t1 = wct_seconds();
	cont = t1 < t0 + 0.25;
	cont = cont && (t1 < t0 + 1 || k < 100);
        pi_allreduce(NULL, &cont, 1, BWC_PI_INT, BWC_PI_MIN, v->pi->m);
	if (!cont)
	    break;
        /* It's difficult to be faithful to the requirements on
         * consistency here. But apparently 1 pleases both operations
         * tested. */
        v->consistency = 1;
	(*f) (v);
    }
    int target = 10 * k / (t1 - t0);
    ASSERT_ALWAYS(target >= 0);
    if (target > 100)
	target = 100;
    if (target == 0)
        target = 1;
    pi_bcast(&target, 1, BWC_PI_INT, 0, 0, v->pi->m);
    t0 = wct_seconds();
    for (k = 0; k < target; k++) {
        pi_log_op(v->pi->m, "[%s] iter%d/%d", __func__, k, target);
        v->consistency = 1;     /* see above */
        (*f) (v);
    }
    serialize(v->pi->m);
    t1 = wct_seconds();
    *pk = k;
    *pt = t1 - t0;
}


void matmul_top_comm_bench(matmul_top_data_ptr mmt, int d)
{
    /* like matmul_top_mul_comm, we'll call mmt_vec_reduce with !d, and
     * mmt_vec_broadcast with d */
    int k;
    double dt;

    void (*funcs[2])(mmt_vec_ptr) = {
        mmt_vec_broadcast,
        mmt_vec_reduce_sameside,
    };
    const char * text[2] = { "bd", "ra" };

    mpfq_vbase_ptr abase = mmt->abase;

    mmt_vec test_vectors[2];
    int is_shared[2] = {0,0};
    mmt_vec_init(mmt, NULL, NULL, test_vectors[0], 0, is_shared[0], mmt->n[0]);
    mmt_vec_init(mmt, NULL, NULL, test_vectors[1], 1, is_shared[1], mmt->n[1]);

    size_t datasize[2];
    
    {
        pi_comm_ptr pirow = mmt->pi->wr[!d];
        pi_comm_ptr picol = mmt->pi->wr[d];
        /* within each row, all jobs are concerned with the same range
         * vrow->i0 to vrow->i1. This is split into m=pirow->njobs
         * chunks, and reduce_scatter has m-1 communication rounds where
         * all of the m nodes output and receive one such chunk. All this
         * happens for all column threads, so we multiply by
         * picol->ncores, too.
         *
         * Note that vrow->i1 - vrow->i0 is #rows / picol->totalsize
         */
        size_t data_out_ra = abase->vec_elt_stride(abase,
                picol->ncores * (mmt->n[!d] / picol->totalsize) /
                pirow->njobs * (pirow->njobs - 1));

        /* one way to do all-gather is to mimick this, except that each
         * node will output the same chunk at each round. Beyond that,
         * the calculation is similar, and we'll use it as a guide. Note
         * of course that if hardware-level multicast is used, our
         * throughput estimation is way off.
         */
        size_t data_out_ag = abase->vec_elt_stride(abase,
                pirow->ncores * (mmt->n[d] / pirow->totalsize) /
                picol->njobs * (picol->njobs - 1));

        datasize[0] = data_out_ag;
        datasize[1] = data_out_ra;
    }

    for(int s = 0 ; s < 2 ; s++) {
        /* we have our axis, and the other axis */
        pi_comm_ptr wr = mmt->pi->wr[d ^ s];          /* our axis */
        pi_comm_ptr xr = mmt->pi->wr[d ^ s ^ 1];      /* other axis */
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
                serialize(mmt->pi->m);
            }
        }
        serialize(mmt->pi->m);
    }
    mmt_vec_clear(mmt, test_vectors[0]);
    mmt_vec_clear(mmt, test_vectors[1]);
}


/**********************************************************************/
/* Utility stuff for applying standard permutations to vectors. */

/* We have four permutations defined in the context of bwc.
 *
 * Sr -- row permutation
 * Sc -- column permutation
 * P  -- shuffled product permutation
 * T  -- decorrelating permutation
 *
 * As a convention, in the code as well as in the documentation below, we
 * freely handle the following equivalent representations of a
 * permutation:
 *  - a function f,
 *  - a list of images [f(i), * 0<=i<n], 
 *  - a list of pairs [<i,f(i)>, * 0<=i<n], 
 *  - a matrix F of size n*n, with only entries at row i and column f(i) being 1
 *
 * We let Mpad be the "padded" matrix, where #rows and #columns are both
 * multiples of nh*nv
 *
 * T is a constant permutation on the columns of the matrix M. It is
 * mostly working as a preconditioner, meant to eliminate some nasty
 * correlation effects between row and column weights. In every respect,
 * the matrix we work with consistently is the matrix Mpad*T. T is computed
 * by balancing_pre_shuffle and balancing_pre_unshuffle.
 *
 * Sr and Sc are defined in the balancing file (as list of images). They
 * are such that the matrix Mtwisted = Sr*Mpad*T*Sc^-1 is well balanced
 * across the different jobs and threads. For square matrices, depending
 * on how the balancing was computed, one may be implicitly defined from
 * the other (but not equal, for the reason below).
 *
 * P is defined only for square matrices. It reflects the fact that
 * although the matrix is "geographically" stored as Mtwisted =
 * Sr*Mpad*T*Sc^-1 in the jobs and threads, the matmul_top_mul code
 * multiplies by a matrix which is:
 *      - for matrix times vector: Sc*Mpad*T*Sc^-1
 *              (that is, v<-v*Transpose(Sc*Mpad*T*Sc^-1) )
 *      - for vector times matrix: P*Sc*Mpad*T*Sc^-1*P^-1
 *
 * Implicit definition of Sr for square matrices
 * =============================================
 *
 * Sc, read as the colperm[] array in the balancing file, is such that column
 * i of the twisted matrix is in fact column Sc(i) in the original matrix,
 * which means that we build Mtwisted as [...]*M*Sc^-1
 *
 * for square matrices, with FLAG_REPLICATE, the Sr permutation is chosen
 * implicitly.
 *
 * here's how we compute the (forward) row permutation in the code (here, xr
 * is the colperm array).
 *              ix = (i * nv + j) * elem;
 *              iy = (j * nh + i) * elem;
 *              for(unsigned int k = 0 ; k < elem ; k++) {
 *                  m->fw_rowperm[xr[iy+k]] = ix+k;
 *
 * We denote by Sr the permutation, implicitly computed here, such that row i
 * of the twisted matrix is row Sr(i) in the original matrix -- so that
 * Mtwisted is in fact Sr*M*Sc^-1. Here, fw_rowperm is the inverse: row i in
 * the original matrix goes to row fw_rowperm[i] in the twisted matrix. So
 * that fw_rowperm == Sr^-1.
 *
 * Our code does fw_rowperm(Sc(P(x)))=x, for P the permutation which sends
 * sub-block nv*i+j to sub-block nh*j+i. Writing as operations on row vectors
 * (as magma does), this gives:
 *      P * Sc * (Sr^-1) = id
 *      Sr = P * Sc
 *
 * So in this case the twisted matrix is in fact Sr*Mpad*Sc^-1, and Sr = P*Sc;
 *
 * Implicit definition of Sc for square matrices
 * =============================================
 *
 * There is no code doing this at the moment, but one could imagine
 * defining this as well. In such a situation, we write down what Sc
 * would be.
 *
 * Matrix is geographically stored as Sr*Mpad*T*Sc^-1 in the jobs and
 * threads, and the matmul_top_mul code multiplies by a matrix which is:
 *      - for matrix times vector: P^-1*Sr*Mpad*T*Sc^-1
 *      - for vector times matrix: Sr*Mpad*T*Sc^-1*P^-1
 * Therefore we want Sc^-1*P^-1 == Sr^-1, whence Sc == P^-1 * Sr
 *
 * Action of P with matmul_top_mul_comm
 * ====================================
 *
 * After v*Mtwisted, reduce_sameside of the result (in direction 1) followed
 * by broadcast (in direction 0) produces a twisted vector. Since v*Mtwisted
 * is in direction 1, then blocks of the vector are to be read in column-major
 * order. If we reduce then broadcast, then sub-block nh*j+i will go to
 * sub-block nv*i+j. This means that v*Mtwisted will be transformed into
 * v*Mtwisted*P^-1
 *
 * After v*Transpose(Mtwisted), reduce_sameside on the result which is in
 * direction 0 transforms it to v*transpose(Mtwisted)*P (here P is
 * transpose of P^-1: we transpose indices blocks in the other
 * direction), which is thus v*transpose(P^-1*Mtwisted).
 *
 * Conclusion: with Mtwisted = Sr*Mpad*Sc^-1, and Sr = P*Sc, when we do
 * matmul_top_mul_cpu followed by matmul_top_mul_comm, we multiply by:
 *      - for matrix times vector: Sc*Mpad*Sc^-1
 *              (that is, v<-v*Transpose(Sc*Mpad*T*Sc^-1) )
 *      - for vector times matrix: P*Sc*Mpad*Sc^-1*P^-1
 *
 * Applying P to vectors
 * =====================
 *
 * We can emulate the multiplications by P and P^-1 with appropriate
 * combinations of apply_identity and matmul_top_mul_comm. So we give a few
 * details.
 *
 * Let v be a distributed vector. With mmt_apply_identity, the output is
 * distributed in the other direction, but not consistent. If we do
 * mmt_vec_reduce_sameside (or mmt_vec_allreduce, which does more), the
 * resulting vector (looking at locally-owned pieces only) is the same as
 * v, in the other direction.
 *
 * And then, matmul_top_mul_comm on this resulting vector does the P
 * action as above. Therefore:
 *
 * For v->d == 0, applying in sequence the functions mmt_apply_identity
 * matmul_top_mul_comm, then we get v*P^-1
 *
 * For v->d == 1, the same sequence produces v*Transpose(P^-1)==v*P
 *
 * Doing the converse is feasible. For v->d==0, if we do instead
 * matmul_top_mul_comm then mmt_apply_identity, we get v*P
 *
 * For v->d==1, if we do matmul_top_mul_comm then mmt_apply_identity, we
 * get v*P^-1
 *
 */

void mmt_apply_identity(mmt_vec_ptr w, mmt_vec_ptr v)
{
    /* input: fully consistent */
    /* output: inconsistent ! 
     * Need mmt_vec_allreduce or mmt_vec_reduce_sameside, or
     * matmul_top_mul_comm, depending on what we want to do. */
    ASSERT_ALWAYS(v->consistency = 2);
    ASSERT_ALWAYS(w->abase == v->abase);
    ASSERT_ALWAYS(v->d != w->d);
    ASSERT_ALWAYS(v->n == w->n);

    mpfq_vbase_ptr A = v->abase;

    serialize_threads(w->pi->m);
    mmt_full_vec_set_zero(w);
    serialize_threads(w->pi->m);

    unsigned int v_off, w_off;
    unsigned int how_many = intersect_two_intervals(&v_off, &w_off,
            v->i0, v->i1, w->i0, w->i1);

    A->vec_set(A,
            A->vec_coeff_ptr(A,  w->v, w_off),
            A->vec_coeff_ptr(A,  v->v, v_off),
            how_many);
    w->consistency = 1;
}

void mmt_vec_apply_or_unapply_P_inner(matmul_top_data_ptr mmt, mmt_vec_ptr y, int apply)
{
    ASSERT_ALWAYS(y->consistency == 2);
    mmt_vec yt;
    mmt_vec_init(mmt, y->abase, y->pitype, yt, !y->d, 0, y->n);
    if ((apply ^ y->d) == 0) {
        // y->d == 0: get v*P^-1
        // y->d == 1: get v*P
        mmt_apply_identity(yt, y);
        matmul_top_mul_comm(y, yt);
    } else {
        // y->d == 0: get v*P
        // y->d == 1: get v*P^-1
        mmt_vec_downgrade_consistency(y);
        matmul_top_mul_comm(yt, y);
        mmt_apply_identity(y, yt);
        mmt_vec_allreduce(y);
    }
    mmt_vec_clear(mmt, yt);
    ASSERT_ALWAYS(y->consistency == 2);
}

void mmt_vec_unapply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y)
{
    mmt_vec_apply_or_unapply_P_inner(mmt, y, 0);
}

void mmt_vec_apply_P(matmul_top_data_ptr mmt, mmt_vec_ptr y)
{
    mmt_vec_apply_or_unapply_P_inner(mmt, y, 1);
}

/* apply == 1 for apply, apply == 0 for unapply *
 * apply == 1 d == 0 Sr defined:   v <- v * Sr
 * apply == 1 d == 0 Sr implicit:  v <- v * Sc
 * apply == 1 d == 1 Sc defined:   v <- v * Sc
 * apply == 1 d == 1 Sc implicit:  v <- v * Sr
 * apply == 0 d == 0 Sr defined:   v <- v * Sr^-1
 * apply == 0 d == 0 Sr implicit:  v <- v * Sc^-1
 * apply == 0 d == 1 Sc defined:   v <- v * Sc^-1
 * apply == 0 d == 1 Sc implicit:  v <- v * Sr^-1
 *
 * See that when, say, Sr is implicitly defined (to P*Sc), this function
 * only applies Sc, not P !
 */
void mmt_vec_apply_or_unapply_S_inner(matmul_top_data_ptr mmt, int midx, mmt_vec_ptr y, int apply)
{
    ASSERT_ALWAYS(y->consistency == 2);
    /* input: fully consistent */
    /* output: fully consistent */
    int d = y->d;
    mpfq_vbase_ptr A = y->abase;

    serialize_threads(y->pi->m);

    /* We'll have two vectors of size n[d], one named y in direction d,
     * and one named yt in direction !d.
     * 
     *
     * In the permutation mmt->perm[d], the pairs (i,j) are such that,
     * given two vectors of size n[d], one named y in direction d,
     * and one named yt in direction !d, we have:
     *  i in [y->i0..y->i1[
     *  j in [yt->i0..yt->i1[
     */
    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
    permutation_data_ptr s = Mloc->perm[d];

    /* For square matrices, we'll use the other permutation transparently
     * with this piece of code. Note though that when we do so, applying
     * the permutation actually goes in the opposite direction. */
    int xd = d;
    if (!s && (Mloc->bal->h->flags & FLAG_REPLICATE)) {
        ASSERT_ALWAYS(Mloc->n[0] == Mloc->n[1]);
        s = Mloc->perm[!d];
        xd = !d;
    }
    if (!s) {
        /* could well be that we have nothing to do */
        return;
    }

    if ((apply^d^xd) == 0) {
        /*
         * apply == 0 d == 0 Sr defined:   v <- v * Sr^-1
         * apply == 0 d == 1 Sc defined:   v <- v * Sc^-1
         * apply == 1 d == 0 Sr implicit:  v <- v * Sc
         * apply == 1 d == 1 Sc implicit:  v <- v * Sr
         */
        mmt_vec yt;
        mmt_vec_init(mmt, A, y->pitype, yt, !d, 0, y->n);

        mmt_apply_identity(yt, y);
        mmt_vec_allreduce(yt);
        mmt_full_vec_set_zero(y);
        serialize_threads(y->pi->m);
        for(unsigned int k = 0 ; k < s->n ; k++) {
            if (s->x[k][d^xd] < y->i0 || s->x[k][d^xd] >= y->i1)
                continue;
            if (s->x[k][d^xd^1] < yt->i0 || s->x[k][d^xd^1] >= yt->i1)
                continue;
            A->vec_set(A,
                    A->vec_coeff_ptr(A,  y->v, s->x[k][d^xd]  - y->i0),
                    A->vec_coeff_ptr(A, yt->v, s->x[k][d^xd^1] - yt->i0),
                    1);
        }
        y->consistency = 1;
        serialize_threads(y->pi->m);
        mmt_vec_allreduce(y);
        mmt_vec_clear(mmt, yt);
    } else {
        /*
         * apply == 1 d == 0 Sr defined:   v <- v * Sr
         * apply == 1 d == 1 Sc defined:   v <- v * Sc
         * apply == 0 d == 0 Sr implicit:  v <- v * Sc^-1
         * apply == 0 d == 1 Sc implicit:  v <- v * Sr^-1
         */
        mmt_vec yt;
        mmt_vec_init(mmt, A, y->pitype, yt, !d, 0, y->n);
        for(unsigned int k = 0 ; k < s->n ; k++) {
            if (s->x[k][d^xd] < y->i0 || s->x[k][d^xd] >= y->i1)
                continue;
            if (s->x[k][d^xd^1] < yt->i0 || s->x[k][d^xd^1] >= yt->i1)
                continue;
            A->vec_set(A,
                    A->vec_coeff_ptr(A, yt->v, s->x[k][d^xd^1] - yt->i0),
                    A->vec_coeff_ptr(A,  y->v, s->x[k][d^xd] -  y->i0),
                    1);
        }
        yt->consistency = 1;
        mmt_vec_allreduce(yt);
        mmt_apply_identity(y, yt);
        mmt_vec_allreduce(y);
        mmt_vec_clear(mmt, yt);
    }
    serialize_threads(y->pi->m);
    ASSERT_ALWAYS(y->consistency == 2);
}

/* multiply v by Sr^-1 if v->d == 0, by Sc^-1 if v->d == 1 */
void mmt_vec_unapply_S(matmul_top_data_ptr mmt, int midx, mmt_vec_ptr y)
{
    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
    mmt_vec_apply_or_unapply_S_inner(mmt, midx, y, 0);
    if ((Mloc->bal->h->flags & FLAG_REPLICATE) && !Mloc->perm[y->d]) {
        if (y->d == 0) {
            /* implicit Sr^-1 is Sc^-1*P^-1 */
            mmt_vec_unapply_P(mmt, y);
        } else {
            /* implicit Sc^-1 is Sr^-1*P */
            mmt_vec_apply_P(mmt, y);
        }
    }
}

/* multiply v by Sr if v->d == 0, by Sc if v->d == 1 */
void mmt_vec_apply_S(matmul_top_data_ptr mmt, int midx, mmt_vec_ptr y)
{
    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
    if ((Mloc->bal->h->flags & FLAG_REPLICATE) && !Mloc->perm[y->d]) {

        if (y->d == 0) {
            /* implicit Sr is P * Sc */
            mmt_vec_apply_P(mmt, y);
        } else {
            /* implicit Sc is P^-1 * Sr */
            mmt_vec_unapply_P(mmt, y);
        }
    }
    mmt_vec_apply_or_unapply_S_inner(mmt, midx, y, 1);
}

/* for square matrix products, the inner loops do
 *   v <- v * (Sr=P*Sc)*Mpad*Sc^-1*P^-1              (for vector times matrix)
 *   v <- v * Transpose(P^-1*(Sr=P*Sc)*Mpad*Sc^-1)   (for matrix times vector)
 * while for non-square (no FLAG_REPLICATE), it's simply by
 *      Sr*Mpad*Sc^-1
 *
 * therefore the rules for twisting and untwisting are as follows.
 *
 * twisting v->d == 0
 *      we assume we want here to change v to become a good *input* for a
 *      vector times matrix operation. Therefore we do:
 *              v <- v * Sr^-1
 * twisting v->d == 1
 *      v <- v * Sc^-1
 *
 */

void mmt_vec_twist(matmul_top_data_ptr mmt, mmt_vec_ptr y)
{
    mmt_vec_unapply_S(mmt, y->d == 0 ? 0 : (mmt->nmatrices-1), y);
}

void mmt_vec_untwist(matmul_top_data_ptr mmt, mmt_vec_ptr y)
{
    mmt_vec_apply_S(mmt, y->d == 0 ? 0 : (mmt->nmatrices-1), y);
}

/* {{{ mmt_vec_{un,}appy_T -- this applies the fixed column
 * permutation which we use unconditionally in bwc to avoid correlation
 * of row and column weights.
 */
    // pshuf indicates two integers a,b such that the COLUMN i of the input
    // matrix is in fact mapped to column a*i+b mod n in the matrix we work
    // with. pshuf_inv indicates the inverse permutation. a and b do
void mmt_vec_apply_or_unapply_T_inner(matmul_top_data_ptr mmt, mmt_vec_ptr y, int apply)
{
    if (y->d == 0) return;
    matmul_top_matrix_ptr Mloc = mmt->matrices[mmt->nmatrices - 1];
    ASSERT_ALWAYS(y->consistency == 2);
    serialize_threads(y->pi->m);
    mmt_vec yt;
    mmt_vec_init(mmt, y->abase, y->pitype, yt, !y->d, 0, y->n);
    for(unsigned int i = y->i0 ; i < y->i1 ; i++) {
        unsigned int j;
        if (apply) {
            j = balancing_pre_shuffle(Mloc->bal, i);
        } else {
            j = balancing_pre_unshuffle(Mloc->bal, i);
        }
        if (j >= yt->i0 && j < yt->i1) {
            y->abase->vec_set(y->abase,
                    y->abase->vec_coeff_ptr(y->abase, yt->v, j - yt->i0),
                    y->abase->vec_coeff_ptr_const(y->abase, y->v, i - y->i0),
                    1);
        }
    }
    yt->consistency = 0;
    mmt_vec_allreduce(yt);
    mmt_apply_identity(y, yt);
    mmt_vec_allreduce(y);
}
void mmt_vec_unapply_T(matmul_top_data_ptr mmt, mmt_vec_ptr v)
{
    mmt_vec_apply_or_unapply_T_inner(mmt, v, 0);
}

void mmt_vec_apply_T(matmul_top_data_ptr mmt, mmt_vec_ptr v)
{
    mmt_vec_apply_or_unapply_T_inner(mmt, v, 1);
}
/* }}} */

/* {{{ Application of permutations to indices */

/* Given indices which relate to a vector in direction d, modify them in
 * place so that they correspond to coefficients of the matching twisted vector.
 *
 * So we want: v_i == twist(v)_{twist(i)}
 *
 * for d == 0, twist(v) = v * Sr^-1, so that index i in v goes to
 * position S(i) in in the twisted vector. For implicit Sr, we have to
 * take into account the fact that Sr=P*Sc, hence S(i) is in fact
 * Sc[P[i]].
 *
 * The algorithm is simple: only one thread knows about each single
 * (i,Sr(i)) pair. So this one changes the relevant xs[] values, leaving
 * the other to zero. And then we do a global allreduce() on the xs[]
 * values.
 *
 * (here, "relevant" may mean something which includes the P
 * permutation).
 */

/* returns the two intervals such that for all pairs (i,j) in
 * mmt->perm->x, we have ii[0] <= i < ii[1], and jj[0] <= j < jj[1]
 */
static void get_local_permutations_ranges(matmul_top_data_ptr mmt, int d, unsigned int ii[2], unsigned int jj[2])
{
    int pos[2];

    for(int dir = 0 ; dir < 2 ; dir++)  {
        pi_comm_ptr piwr = mmt->pi->wr[dir];
        pos[dir] = piwr->jrank * piwr->ncores + piwr->trank;
    }

    size_t e = mmt->n[d] / mmt->pi->m->totalsize;
    ii[0] = e *  pos[!d]    * mmt->pi->wr[d]->totalsize;
    ii[1] = e * (pos[!d]+1) * mmt->pi->wr[d]->totalsize;
    jj[0] = e *  pos[d]     * mmt->pi->wr[!d]->totalsize;
    jj[1] = e * (pos[d]+1)  * mmt->pi->wr[!d]->totalsize;
}

void indices_twist(matmul_top_data_ptr mmt, uint32_t * xs, unsigned int n, int d)
{
    int midx = d == 0 ? 0 : (mmt->nmatrices - 1);
    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
    /* d == 1: twist(v) = v*Sc^-1
     * d == 0: twist(v) = v*Sr^-1
     */
    unsigned int ii[2], jj[2];

    if (Mloc->perm[d]) {
        /* explicit S */
        /* coordinate S[i] in the original vector becomes i in the
         * twisted vector.
         */
        get_local_permutations_ranges(mmt, d, ii, jj);
        unsigned int * r = malloc((jj[1] - jj[0]) * sizeof(unsigned int));
        memset(r, 0, (jj[1] - jj[0]) * sizeof(unsigned int));
        for(size_t i = 0 ; i < Mloc->perm[d]->n ; i++) {
            ASSERT_ALWAYS(Mloc->perm[d]->x[i][0] >= ii[0]);
            ASSERT_ALWAYS(Mloc->perm[d]->x[i][0] <  ii[1]);
            ASSERT_ALWAYS(Mloc->perm[d]->x[i][1] >= jj[0]);
            ASSERT_ALWAYS(Mloc->perm[d]->x[i][1] <  jj[1]);
            // r[Mloc->perm[d]->x[i][0] - ii[0]] = Mloc->perm[d]->x[i][1];
            r[Mloc->perm[d]->x[i][1] - jj[0]] = Mloc->perm[d]->x[i][0];
        }
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = xs[k];
            if (j >= jj[0] && j < jj[1])
                xs[k] = r[j - jj[0]];
            else
                xs[k] = 0;
        }
        free(r);
        pi_allreduce(NULL, xs, n * sizeof(uint32_t), BWC_PI_BYTE, BWC_PI_BXOR, mmt->pi->m);
    } else if (Mloc->perm[!d] && (Mloc->bal->h->flags & FLAG_REPLICATE)) {
        ASSERT_ALWAYS(Mloc->n[0] == Mloc->n[1]);
        /* implicit S -- first we get the bits about the S in the other
         * direction, because the pieces we have are for the other
         * ranges, which is a bit disturbing...
         */
        get_local_permutations_ranges(mmt, !d, ii, jj);
        unsigned int * r = malloc((jj[1] - jj[0]) * sizeof(unsigned int));
        memset(r, 0, (jj[1] - jj[0]) * sizeof(unsigned int));
        for(size_t i = 0 ; i < Mloc->perm[!d]->n ; i++) {
            ASSERT_ALWAYS(Mloc->perm[!d]->x[i][0] >= ii[0]);
            ASSERT_ALWAYS(Mloc->perm[!d]->x[i][0] <  ii[1]);
            ASSERT_ALWAYS(Mloc->perm[!d]->x[i][1] >= jj[0]);
            ASSERT_ALWAYS(Mloc->perm[!d]->x[i][1] <  jj[1]);
            r[Mloc->perm[!d]->x[i][1] - jj[0]] = Mloc->perm[!d]->x[i][0];
        }
        /* nh and nv are the same for all submatrices, really */
        unsigned int nn[2] = { Mloc->bal->h->nh, Mloc->bal->h->nv };
        /* for d == 0, we have implicit Sr = P * Sc.
         * for d == 1, we have implicit Sc = P^-1 * Sc
         */
        size_t z = Mloc->bal->trows / (nn[0]*nn[1]);
        for(unsigned int k = 0 ; k < n ; k++) {
            unsigned int j = xs[k];
            /* for d == 0, index i goes to P^-1[Sc^-1[i]] */
            if (j >= jj[0] && j < jj[1]) {
                unsigned int i = r[j - jj[0]];
                unsigned int qz = i / z;
                unsigned int rz = i % z;
                /* P    sends sub-block nv*i+j to sub-block nh*j+i */
                /* P^-1 sends sub-block nh*j+i to sub-block nv*i+j */
                unsigned int qh = qz / nn[d];
                unsigned int rh = qz % nn[d];
                ASSERT_ALWAYS(qz == qh * nn[d] + rh);
                ASSERT_ALWAYS(i == (qh * nn[d] + rh)*z + rz);
                xs[k] = (rh * nn[!d] + qh)*z + rz;
            } else {
                xs[k] = 0;
            }
        }
        free(r);
        pi_allreduce(NULL, xs, n * sizeof(uint32_t), BWC_PI_BYTE, BWC_PI_BXOR, mmt->pi->m);
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

/* For d == 0: do w = v * M
 * For d == 1: do w = M * v
 *
 * We do not necessarily have v->d == d, although this will admittedly be
 * the case most often:
 * - in block Wiedemann, we have a vector split in some direction (say
 *   d==1 when we wanna solve Mw=0), we compute w=Mv, and then there's
 *   the matmul_top_mul_comm step which moves stuff to v again.
 * - in block Lanczos, say we start from v->d == 1 again. We do w=M*v,
 *   so that w->d==0. But then we want to compute M^T * w, which is w^T *
 *   M. So again, w->d == 0 is appropriate with the second product being
 *   in the direction foo*M.
 * - the only case where this does not necessarily happen so is when we
 *   have several matrices.
 */
void matmul_top_mul_cpu(matmul_top_data_ptr mmt, int midx, int d, mmt_vec_ptr w, mmt_vec_ptr v)
{
    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
    ASSERT_ALWAYS(v->consistency == 2);
    ASSERT_ALWAYS(w->abase == v->abase);
    unsigned int di_in  = v->i1 - v->i0;
    unsigned int di_out = w->i1 - w->i0;
    ASSERT_ALWAYS(Mloc->mm->dim[!d] == di_out);
    ASSERT_ALWAYS(Mloc->mm->dim[d] == di_in);

    ASSERT_ALWAYS(w->siblings); /* w must not be shared */

    pi_log_op(mmt->pi->m, "[%s:%d] enter matmul_mul", __func__, __LINE__);

    /* Note that matmul_init copies the calling abase argument to the
     * lower-level mm structure. It can quite probably be qualified as a
     * flaw.
     */
    matmul_mul(Mloc->mm, w->v, v->v, d);
    w->consistency = 0;
}

/* This takes partial results in w, and puts the
 * collected and re-broadcasted results in the areas mmt->wd[d]->v
 *
 * Note that for the shuffled product, this is not equivalent to a trivial
 * operation.
 */
void matmul_top_mul_comm(mmt_vec_ptr v, mmt_vec_ptr w)
{
    /* this takes inconsistent input.
     * XXX if we have fully consistent input, then a reduce() is much
     * undesired !
     */
    ASSERT_ALWAYS(w->consistency != 2);
    pi_log_op(v->pi->m, "[%s:%d] enter mmt_vec_reduce", __func__, __LINE__);
    mmt_vec_reduce(v, w);
    ASSERT_ALWAYS(v->consistency == 1);
    pi_log_op(v->pi->m, "[%s:%d] enter mmt_vec_broadcast", __func__, __LINE__);
    mmt_vec_broadcast(v);
    ASSERT_ALWAYS(v->consistency == 2);

    /* If we have shared input data for the column threads, then we'd
     * better make sure it has arrived completely, because while all
     * threads will need the data, only one is actually importing it.
     */
    if (!v->siblings) {
        pi_log_op(v->pi->wr[v->d], "[%s:%d] serialize threads", __func__, __LINE__);
        serialize_threads(v->pi->wr[v->d]);
    }
}


/* Vector I/O is done by only one job, one thread. It incurs a
 * significant amount of memory allocation, but this is done relatively
 * rarely.
 *
 * Vectors, as handled within the core routines, are permuted. He have
 * coordinates (v_{\sigma(0)}...v_{\sigma(n-1)}), where \sigma is the
 * permutation given by the *_perm files. For block Wiedemann, we assume
 * that we have conjugated permutations on the two sides. This means that
 * no data manipulation is required within the critical loops (this would
 * imply a fuzzy communication pattern, boiling down to essentially using
 * allreduce instead of reduce).
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
// Doing a mmt_vec_broadcast columns will ensure that each row contains
// the complete data set for our vector.

void mmt_vec_set_random_through_file(mmt_vec_ptr v, const char * filename, unsigned int itemsondisk, gmp_randstate_t rstate)
{
    /* FIXME: this generates the complete vector on rank 0, saves it, and
     * loads it again. But I'm a bit puzzled by the choice of saving a
     * number of items which is n0[d]. Seems to me that this is in fact
     * incorrect, we want n0[!d] here.
     */
    mpfq_vbase_ptr A = v->abase;
    pi_datatype_ptr A_pi = v->pitype;
    parallelizing_info_ptr pi = v->pi;

    if (pi->m->trank == 0) {
        void * y;
        cheating_vec_init(A, &y, v->n);
        A->vec_set_zero(A, y, v->n);
        /* we generate garbage for the padding coordinates too, but
         * that's not really an issue since that does not get written,
         * thanks to itemsondisk conveying the info about the correct
         * size
         */
        A->vec_random(A, y, v->n, rstate);
        int err = MPI_Bcast(y, v->n, A_pi->datatype, 0, pi->m->pals);
        ASSERT_ALWAYS(!err);
        FILE * f = fopen(filename, "wb");
        ASSERT_ALWAYS(f);
        int rc = fwrite(y, A->vec_elt_stride(A,1), itemsondisk, f);
        ASSERT_ALWAYS(rc == (int) itemsondisk);
        fclose(f);
        cheating_vec_clear(A, &y, v->n);
    }
    mmt_vec_load(v, filename, itemsondisk);
}

void mmt_vec_set_random_inconsistent(mmt_vec_ptr v, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(v != NULL);
    v->abase->vec_random(v->abase, v->v, v->i1 - v->i0, rstate);
    v->consistency=1;
}

void mmt_vec_truncate(matmul_top_data_ptr mmt, mmt_vec_ptr v)
{
    ASSERT_ALWAYS(v != NULL);
    if (mmt->n0[v->d] >= v->i0 && mmt->n0[v->d] < v->i1) {
        v->abase->vec_set_zero(v->abase, 
                v->abase->vec_subvec(v->abase, 
                    v->v, mmt->n0[v->d] - v->i0),
                v->i1 - mmt->n0[v->d]);
    }
}

void mmt_vec_set_x_indices(mmt_vec_ptr y, uint32_t * gxvecs, int m, unsigned int nx)
{
    int shared = !y->siblings;
    mpfq_vbase_ptr A = y->abase;
    mmt_full_vec_set_zero(y);
    void * dummy;
    cheating_vec_init(A, &dummy, 1);
    if (!shared || y->pi->wr[y->d]->trank == 0) {
        for(int j = 0 ; j < m ; j++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                uint32_t i = gxvecs[j*nx+k];
                // set bit j of entry i to 1.
                if (i < y->i0 || i >= y->i1)
                    continue;
                {
                    /* We're not doing to do set_ui_at unconditionally,
                     * because then we would be doing rubbish in presence
                     * of duplicate coordinates -- which can happen.
                     *
                     * So set_ui_at on something which we know is zero
                     * first, then add.
                     */
                    A->set_zero(A, dummy);
                    A->set_ui_at(A, dummy, j, 1);
                    A->add(A,
                            A->vec_coeff_ptr(A, y->v, i - y->i0), 
                            A->vec_coeff_ptr(A, y->v, i - y->i0),
                            dummy);
                }
            }
        }
    }
    cheating_vec_clear(A, &dummy, 1);
    y->consistency=2;
    if (shared)
        serialize_threads(y->pi->wr[y->d]);
}



/**********************************************************************/
static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, int midx, param_list_ptr pl, int optimized_direction);

/* returns an allocated string holding the name of the midx-th submatrix */
static char * matrix_list_get_item(param_list_ptr pl, const char * key, int midx)
{
    char * res = NULL;
    char ** mnames;
    int nmatrices;
    int rc = param_list_parse_string_list_alloc(pl, key, &mnames, &nmatrices, ",");
    if (rc == 0)
        return NULL;
    ASSERT_ALWAYS(midx < nmatrices);
    for(int i = 0 ; i < nmatrices ; i++) {
        if (i == midx) {
            res = mnames[i];
        } else {
            free(mnames[i]);
        }
    }
    free(mnames);
    return res;
}

/* return an allocated string with the name of a balancing file for this
 * matrix and this mpi/thr split.
 */
static char* matrix_get_derived_balancing_filename(const char * matrixname, parallelizing_info_ptr pi)
{
    /* input is NULL in the case of random matrices */
    if (!matrixname) return NULL;
    unsigned int nh = pi->wr[1]->totalsize;
    unsigned int nv = pi->wr[0]->totalsize;
    char * copy = strdup(matrixname);
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int rc = asprintf(&t, "%s.%ux%u.bin", t, nh, nv);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    return t;
}

static char* matrix_get_derived_cache_filename_stem(const char * matrixname, parallelizing_info_ptr pi, uint32_t checksum)
{
    /* input is NULL in the case of random matrices */
    if (!matrixname) return NULL;
    unsigned int nh = pi->wr[1]->totalsize;
    unsigned int nv = pi->wr[0]->totalsize;
    char * copy = strdup(matrixname);
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int pos[2];
    for(int d = 0 ; d < 2 ; d++)  {
        pi_comm_ptr wr = pi->wr[d];
        pos[d] = wr->jrank * wr->ncores + wr->trank;
    }
    int rc = asprintf(&t, "%s.%ux%u.%08" PRIx32 ".h%d.v%d", t, nh, nv, checksum, pos[1], pos[0]);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    return t;
}

static char* matrix_get_derived_submatrix_filename(const char * matrixname, parallelizing_info_ptr pi)
{
    /* input is NULL in the case of random matrices */
    if (!matrixname) return NULL;
    unsigned int nh = pi->wr[1]->totalsize;
    unsigned int nv = pi->wr[0]->totalsize;
    char * copy = strdup(matrixname);
    char * t;
    if (strlen(copy) > 4 && strcmp((t = copy + strlen(copy) - 4), ".bin") == 0) {
        *t = '\0';
    }
    if ((t = strrchr(copy, '/')) == NULL) { /* basename */
        t = copy;
    } else {
        t++;
    }
    int pos[2];
    for(int d = 0 ; d < 2 ; d++)  {
        pi_comm_ptr wr = pi->wr[d];
        pos[d] = wr->jrank * wr->ncores + wr->trank;
    }
    int rc = asprintf(&t, "%s.%ux%u.h%d.v%d.bin", t, nh, nv, pos[1], pos[0]);
    ASSERT_ALWAYS(rc >=0);
    free(copy);
    return t;
}

static void matmul_top_init_fill_balancing_header(matmul_top_data_ptr mmt, int i, param_list_ptr pl)
{
    parallelizing_info_ptr pi = mmt->pi;
    matmul_top_matrix_ptr Mloc = mmt->matrices[i];

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        if (!Mloc->mname) {
            random_matrix_fill_fake_balancing_header(Mloc->bal, pi, param_list_lookup_string(pl, "random_matrix"));
        } else {
            if (access(Mloc->bname, R_OK) != 0) {
                if (errno == ENOENT) {
                    printf("Creating balancing file %s\n", Mloc->bname);
                    struct mf_bal_args mba = {
                        .mfile = Mloc->mname,
                        .bfile = Mloc->bname,
                        .nh = pi->wr[1]->totalsize,
                        .nv = pi->wr[0]->totalsize,
                        .do_perm = { MF_BAL_PERM_AUTO, MF_BAL_PERM_AUTO },
                    };
                    mf_bal_adjust_from_option_string(&mba, param_list_lookup_string(pl, "balancing_options"));
                    /* withcoeffs being a switch for param_list, it is
                     * clobbered by the configure_switch mechanism */
                    mba.withcoeffs = !is_char2(mmt->abase);
                    mf_bal(&mba);
                } else {
                    fprintf(stderr, "Cannot access balancing file %s: %s\n", Mloc->bname, strerror(errno));
                }
            }
            balancing_read_header(Mloc->bal, Mloc->bname);
        }
    }
    pi_bcast(Mloc->bal, sizeof(balancing), BWC_PI_BYTE, 0, 0, mmt->pi->m);

    /* check that balancing dimensions are compatible with our run */
    int ok = 1;
    ok = ok && mmt->pi->wr[0]->totalsize == Mloc->bal->h->nv;
    ok = ok && mmt->pi->wr[1]->totalsize == Mloc->bal->h->nh;
    if (ok) return;

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        fprintf(stderr, "Matrix %d, %s: balancing file %s"
                " has dimensions %ux%u,"
                " this conflicts with the current run,"
                " which expects dimensions (%ux%u)x(%ux%u).\n",
                i, Mloc->mname, Mloc->bname,
                Mloc->bal->h->nh, Mloc->bal->h->nv,
                mmt->pi->wr[1]->njobs,
                mmt->pi->wr[1]->ncores,
                mmt->pi->wr[0]->njobs,
                mmt->pi->wr[0]->ncores);
    }
    serialize(mmt->pi->m);
    exit(1);
}


static void matmul_top_init_prepare_local_permutations(matmul_top_data_ptr mmt, int i)
{
    matmul_top_matrix_ptr Mloc = mmt->matrices[i];
    /* Here, we get a copy of the rowperm and colperm.
     *
     * For each (job,thread), two pairs of intervals are defined.
     *
     * for row indices: [i0[0], i1[0][ = [mrow->i0, mrow->i1[ and [Xi0[0], Xi1[0][
     * for col indices: [i0[1], i1[1][ = [mcol->i0, mcol->i1[ and [Xi0[1], Xi1[1][
     *
     * The parts we get are:
     *      [i, rowperm[i]] for    i in [mrow->i0, mrow->i1[
     *                         and rowperm[i] in [X0, X1[
     *      [j, colperm[j]] for    j in [mcol->i0, mcol->i1[
     *                         and colperm[i] in [Y0, Y1[
     * where X0 is what would be mcol->i0 if the matrix as many columns as
     * rows, and ditto for X1,Y0,Y1.
     */

    unsigned int rowperm_items=0;
    unsigned int colperm_items=0;

    /* Define a complete structure for the balancing which is shared
     * among threads, so that we'll be able to access it from all threads
     * simultaneously. We will put things in bal_tmp->rowperm and
     * bal_tmp->colperm, but beyond that, the header part will be wrong
     * at non-root nodes.
     */
    balancing_ptr bal_tmp = shared_malloc_set_zero(mmt->pi->m, sizeof(balancing));

    if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
        if (Mloc->bname)
            balancing_read(bal_tmp, Mloc->bname);
        /* It's fine if we have nothing. This just means that we'll have
         * no balancing to deal with (this occurs only for matrices
         * which are generated at random on the fly). */
        rowperm_items = bal_tmp->rowperm != NULL ? bal_tmp->trows : 0;
        colperm_items = bal_tmp->colperm != NULL ? bal_tmp->tcols : 0;
    }
    pi_bcast(&rowperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt->pi->m);
    pi_bcast(&colperm_items, 1, BWC_PI_UNSIGNED, 0, 0, mmt->pi->m);

    if (mmt->pi->m->trank == 0) {
        if (rowperm_items) {
            ASSERT_ALWAYS(rowperm_items == Mloc->bal->trows);
            if (mmt->pi->m->jrank != 0)
                bal_tmp->rowperm = malloc(Mloc->bal->trows * sizeof(uint32_t));
            MPI_Bcast(bal_tmp->rowperm, Mloc->bal->trows * sizeof(uint32_t), MPI_BYTE, 0, mmt->pi->m->pals);
        }
        if (colperm_items) {
            ASSERT_ALWAYS(colperm_items == Mloc->bal->tcols);
            if (mmt->pi->m->jrank != 0)
                bal_tmp->colperm = malloc(Mloc->bal->tcols * sizeof(uint32_t));
            MPI_Bcast(bal_tmp->colperm, Mloc->bal->tcols * sizeof(uint32_t), MPI_BYTE, 0, mmt->pi->m->pals);
        }
    }
    serialize_threads(mmt->pi->m);      /* important ! */

    uint32_t * balperm[2] = { bal_tmp->rowperm, bal_tmp->colperm };
    for(int d = 0 ; d < 2 ; d++)  {
        unsigned int ii[2];
        unsigned int jj[2];
        get_local_permutations_ranges(mmt, d, ii, jj);

        if (!balperm[d]) continue;

        Mloc->perm[d] = permutation_data_alloc();
#ifdef  MVAPICH2_NUMVERSION
        /* apparently mvapich2 frowns on realloc() */
        permutation_data_ensure(Mloc->perm[d],ii[1] - ii[0]);
#endif

        /* now create the really local permutation */
        for(unsigned int i = ii[0] ; i < ii[1] ; i++) {
            unsigned int j = balperm[d][i];
            if (j >= jj[0] && j < jj[1]) {
                unsigned int ij[2] = { i, j };
                permutation_data_push(Mloc->perm[d], ij);
            }
        }
#if 0
        const char * text[2] = { "left", "right", };
        printf("[%s] J%uT%u does %zu/%u permutation pairs for %s vectors\n",
                mmt->pi->nodenumber_s,
                mmt->pi->m->jrank, mmt->pi->m->trank,
                Mloc->perm[d]->n, d ? Mloc->bal->tcols : Mloc->bal->trows,
                text[d]);
#endif
    }

    serialize_threads(mmt->pi->m);      /* important ! */

    if (mmt->pi->m->trank == 0) {
        if (bal_tmp->colperm) free(bal_tmp->colperm);
        if (bal_tmp->rowperm) free(bal_tmp->rowperm);
    }

    shared_free(mmt->pi->m, bal_tmp);
}

void matmul_top_init(matmul_top_data_ptr mmt,
        mpfq_vbase_ptr abase,
        /* matmul_ptr mm, */
        parallelizing_info_ptr pi,
        param_list_ptr pl,
        int optimized_direction)
{
    memset(mmt, 0, sizeof(*mmt));

    mmt->abase = abase;
    mmt->pitype = pi_alloc_mpfq_datatype(pi, abase);
    mmt->pi = pi;
    mmt->matrices = NULL;

    int nbals = param_list_get_list_count(pl, "balancing");
    mmt->nmatrices = param_list_get_list_count(pl, "matrix");
    const char * random_description = param_list_lookup_string(pl, "random_matrix");


    if (random_description) {
        if (nbals || mmt->nmatrices) {
            fprintf(stderr, "random_matrix is incompatible with balancing= and matrix=\n");
            exit(EXIT_FAILURE);
        }
    } else if (nbals && !mmt->nmatrices) {
        fprintf(stderr, "missing parameter matrix=\n");
        exit(EXIT_FAILURE);
    } else if (!nbals && mmt->nmatrices) {
        /* nbals == 0 is a hint towards taking the default balancing file
         * names, that's it */
    } else if (nbals != mmt->nmatrices) {
        fprintf(stderr, "balancing= and matrix= have inconsistent number of items\n");
        exit(EXIT_FAILURE);
    }

    if (random_description)
        mmt->nmatrices = 1;

    mmt->matrices = malloc(mmt->nmatrices * sizeof(matmul_top_matrix));
    memset(mmt->matrices, 0, mmt->nmatrices * sizeof(matmul_top_matrix));

    serialize_threads(mmt->pi->m);

    /* The initialization goes through several passes */
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        matmul_top_matrix_ptr Mloc = mmt->matrices[i];
        Mloc->mname = matrix_list_get_item(pl, "matrix", i);
        Mloc->bname = matrix_list_get_item(pl, "balancing", i);
        if (!Mloc->bname) {
            Mloc->bname = matrix_get_derived_balancing_filename(Mloc->mname, mmt->pi);
        }
        ASSERT_ALWAYS((Mloc->bname != NULL) == !random_description);

        matmul_top_init_fill_balancing_header(mmt, i, pl);

        Mloc->n[0] = Mloc->bal->trows;
        Mloc->n[1] = Mloc->bal->tcols;
        Mloc->n0[0] = Mloc->bal->h->nrows;
        Mloc->n0[1] = Mloc->bal->h->ncols;
        Mloc->locfile = matrix_get_derived_cache_filename_stem(Mloc->mname, mmt->pi, Mloc->bal->h->checksum);

    }

    mmt->n[0] = mmt->matrices[0]->n[0];
    mmt->n0[0] = mmt->matrices[0]->n0[0];
    mmt->n[1] = mmt->matrices[mmt->nmatrices-1]->n[1];
    mmt->n0[1] = mmt->matrices[mmt->nmatrices-1]->n0[1];

    /* in the second loop below, get_local_permutations_ranges uses
     * mmt->n[], so we do it in a second pass.
     * Now, given that double matrices only barely work at the moment,
     * I'm not absolutely sure that it's really needed.
     */
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        matmul_top_matrix_ptr Mloc = mmt->matrices[i];

        matmul_top_init_prepare_local_permutations(mmt, i);

        if (!mmt->pi->interleaved) {
            matmul_top_read_submatrix(mmt, i, pl, optimized_direction );
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
                matmul_top_read_submatrix(mmt, i, pl, optimized_direction);
                pi_store_generic(mmt->pi, MMT_MM_MAGIC_KEY + i, mmt->pi->m->trank, Mloc->mm);
            } else {
                Mloc->mm = pi_load_generic(mmt->pi, MMT_MM_MAGIC_KEY + i, mmt->pi->m->trank);
            }
        }
    }
}

unsigned int matmul_top_rank_upper_bound(matmul_top_data_ptr mmt)
{
    unsigned int r = MAX(mmt->n0[0], mmt->n0[1]);
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        matmul_top_matrix_ptr Mloc = mmt->matrices[i];
        r = MAX(r, Mloc->bal->h->nrows - Mloc->bal->h->nzrows);
        r = MAX(r, Mloc->bal->h->ncols - Mloc->bal->h->nzcols);
    }
    return r;
}


static int export_cache_list_if_requested(matmul_top_matrix_ptr Mloc, parallelizing_info_ptr pi, param_list_ptr pl)
{
    const char * cachelist = param_list_lookup_string(pl, "export_cachelist");
    if (!cachelist) return 0;

    char * myline = NULL;
    int rc;
    rc = asprintf(&myline, "%s %s", pi->nodename, Mloc->mm->cachefile_name);
    ASSERT_ALWAYS(rc >= 0);
    ASSERT_ALWAYS(myline != NULL);
    char ** tlines = NULL;
    tlines = shared_malloc_set_zero(pi->m, pi->m->ncores * sizeof(const char *));
    tlines[pi->m->trank] = myline;
    serialize_threads(pi->m);

    /* Also, just out of curiosity, try to see what we have currently */
    struct stat st[1];
    int * has_cache = shared_malloc_set_zero(pi->m, pi->m->totalsize * sizeof(int));
    rc = stat(Mloc->mm->cachefile_name, st);
    unsigned int mynode = pi->m->ncores * pi->m->jrank;
    has_cache[mynode + pi->m->trank] = rc == 0;
    serialize_threads(pi->m);

    int len = 0;
    if (pi->m->trank == 0) {
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                has_cache, pi->m->ncores, MPI_INT, pi->m->pals);

        for(unsigned int j = 0 ; j < pi->m->ncores ; j++) {
            int s = strlen(tlines[j]);
            if (s >= len) len = s;
        }
        MPI_Allreduce(MPI_IN_PLACE, &len, 1, MPI_INT, MPI_MAX, pi->m->pals);
        char * info = malloc(pi->m->totalsize * (len + 1));
        char * mybuf = malloc(pi->m->ncores * (len+1));
        memset(mybuf, 0, pi->m->ncores * (len+1));
        for(unsigned int j = 0 ; j < pi->m->ncores ; j++) {
            memcpy(mybuf + j * (len+1), tlines[j], strlen(tlines[j])+1);
        }
        MPI_Allgather(mybuf, pi->m->ncores * (len+1), MPI_BYTE,
                info, pi->m->ncores * (len+1), MPI_BYTE,
                pi->m->pals);
        if (pi->m->jrank == 0) {
            FILE * f = fopen(cachelist, "wb");
            DIE_ERRNO_DIAG(f == NULL, "fopen", cachelist);
            for(unsigned int j = 0 ; j < pi->m->njobs ; j++) {
                unsigned int j0 = j * pi->m->ncores;
                fprintf(f, "get-cache ");
                for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            for(unsigned int j = 0 ; j < pi->m->njobs ; j++) {
                unsigned int j0 = j * pi->m->ncores;
                fprintf(f, "has-cache ");
                for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
                    char * t = info + (j0 + k) * (len + 1);
                    char * q = strchr(t, ' ');
                    ASSERT_ALWAYS(q);
                    if (!k) {
                        *q='\0';
                        fprintf(f, "%s", t);
                        *q=' ';
                    }
                    if (!has_cache[pi->m->ncores * j + k]) continue;
                    fprintf(f, "%s", q);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        free(info);
        free(mybuf);
    }
    serialize_threads(pi->m);
    shared_free(pi->m, tlines);
    shared_free(pi->m, has_cache);
    serialize_threads(pi->m);
    free(myline);
    serialize(pi->m);

    return 1;
}

static unsigned int local_fraction(unsigned int padded, unsigned int normal, pi_comm_ptr wr)
{
    ASSERT_ALWAYS(padded % wr->totalsize == 0);
    unsigned int i = wr->jrank * wr->ncores + wr->trank;
    unsigned int quo = padded / wr->totalsize;
    return MIN(normal - i * quo, quo);
}


static void matmul_top_read_submatrix(matmul_top_data_ptr mmt, int midx, param_list_ptr pl, int optimized_direction)
{
    int rebuild = 0;
    param_list_parse_int(pl, "rebuild_cache", &rebuild);
    int can_print = (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0);

    matmul_top_matrix_ptr Mloc = mmt->matrices[midx];

    Mloc->mm = matmul_init(mmt->abase,
            Mloc->n[0] / mmt->pi->wr[1]->totalsize,
            Mloc->n[1] / mmt->pi->wr[0]->totalsize,
            Mloc->locfile, NULL  /* means: choose mm_impl from pl */,
            pl, optimized_direction);

    // *IF* we need to do a collective read of the matrix, we need to
    // provide the pointer *now*.
    unsigned int sqread = 0;
    param_list_parse_uint(pl, "sequential_cache_read", &sqread);

    int cache_loaded = 0;

    if (export_cache_list_if_requested(Mloc, mmt->pi, pl)) {
        /* If we are being called from dispatch, once all submatrices
         * have had their list of required files printed, the program
         * will exit. */
        return;
    }

    if (!rebuild) {
        if (can_print) {
            printf("Now trying to load matrix cache files\n");
        }
        if (sqread) {
            for(unsigned int j = 0 ; j < mmt->pi->m->ncores ; j++) {
                serialize_threads(mmt->pi->m);
                if (j == mmt->pi->m->trank) {
                    cache_loaded = matmul_reload_cache(Mloc->mm);
                }
            }
        } else {
            cache_loaded = matmul_reload_cache(Mloc->mm);
        }
    }

    if (!pi_data_eq(&cache_loaded, 1, BWC_PI_INT, mmt->pi->m)) {
        if (can_print) {
            fprintf(stderr, "Fatal error: cache files not present at expected locations\n");
        }
        SEVERAL_THREADS_PLAY_MPI_BEGIN(mmt->pi->m) {
            fprintf(stderr, "[%s] J%uT%u: cache %s: %s\n",
                    mmt->pi->nodenumber_s,
                    mmt->pi->m->jrank,
                    mmt->pi->m->trank,
                    Mloc->mm->cachefile_name,
                    cache_loaded ? "ok" : "not ok");
        }
        SEVERAL_THREADS_PLAY_MPI_END();
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
        // the mm layer is informed of the higher priority computations
        // that will take place. Depending on the implementation, this
        // may cause the direct or transposed ordering to be preferred.
        // Thus we have to read this back from the mm structure.
        m->bfile = Mloc->bname;
        m->mfile = Mloc->mname;
        m->transpose = Mloc->mm->store_transposed;
        m->withcoeffs = !is_char2(mmt->abase);
        if (!(Mloc->mname)) {
            if (can_print) {
                printf("Begin creation of fake matrix data in parallel\n");
            }
            /* Mloc->mm->dim[0,1] contains the dimensions of the padded
             * matrix. This is absolutely fine in the normal case. But in
             * the case of staged matrices, it's a bit different. We must
             * make sure that we generate matrices which have zeroes in
             * the padding area.
             */
            random_matrix_get_u32(mmt->pi, pl, m,
                    local_fraction(Mloc->n[0], Mloc->n0[0], mmt->pi->wr[1]),
                    local_fraction(Mloc->n[1], Mloc->n0[1], mmt->pi->wr[0]));
        } else {
            if (can_print) {
                printf("Matrix dispatching starts\n");
            }
            balancing_get_matrix_u32(mmt->pi, pl, m);

            int ssm = 0;
            param_list_parse_int(pl, "save_submatrices", &ssm);
            if (ssm) {
                char * submat = matrix_get_derived_submatrix_filename(Mloc->mname, mmt->pi);
                fprintf(stderr, "DEBUG: creating %s\n", submat);
                FILE * f = fopen(submat, "wb");
                fwrite(m->p, sizeof(uint32_t), m->size, f);
                fclose(f);
                free(submat);
            }
        }
    }


    if (!sqb) {
        if (!cache_loaded) {
            // everybody does it in parallel
            if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                printf("[%s] J%uT%u building cache for %s\n",
                        mmt->pi->nodenumber_s,
                        mmt->pi->m->jrank,
                        mmt->pi->m->trank,
                        Mloc->locfile);
            matmul_build_cache(Mloc->mm, m);
            matmul_save_cache(Mloc->mm);
        }
    } else {
        for(unsigned int j = 0 ; j < mmt->pi->m->ncores + 1 ; j++) {
            serialize_threads(mmt->pi->m);
            if (cache_loaded) continue;
            if (j == mmt->pi->m->trank) {
                if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
                    printf("[%s] J%uT%u building cache for %s\n",
                            mmt->pi->nodenumber_s,
                            mmt->pi->m->jrank,
                            mmt->pi->m->trank,
                            Mloc->locfile);
                matmul_build_cache(Mloc->mm, m);
            } else if (j == mmt->pi->m->trank + 1) {
                matmul_save_cache(Mloc->mm);
            }
        }
    }
    /* see remark in raw_matrix_u32.h about data ownership for type
     * matrix_u32 */

    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO)) {
        my_pthread_mutex_lock(mmt->pi->m->th->m);
        printf("[%s] J%uT%u uses cache file %s\n",
                mmt->pi->nodenumber_s,
                mmt->pi->m->jrank, mmt->pi->m->trank,
                /* cache for mmt->locfile, */
                Mloc->mm->cachefile_name);
        my_pthread_mutex_unlock(mmt->pi->m->th->m);
    }
}

void matmul_top_report(matmul_top_data_ptr mmt, double scale)
{
    for(int midx = 0 ; midx < mmt->nmatrices ; midx++) {
        matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
        matmul_report(Mloc->mm, scale);
    }
}

void matmul_top_clear(matmul_top_data_ptr mmt)
{
    pi_free_mpfq_datatype(mmt->pi, mmt->pitype);
    serialize_threads(mmt->pi->m);
    serialize(mmt->pi->m);
    serialize_threads(mmt->pi->m);

    for(int midx = 0 ; midx < mmt->nmatrices ; midx++) {
        matmul_top_matrix_ptr Mloc = mmt->matrices[midx];
        for(int d = 0 ; d < 2 ; d++)  {
            permutation_data_free(Mloc->perm[d]);
            // both are expected to hold storage for:
            // (mmt->n[d] / mmt->pi->m->totalsize * mmt->pi->wr[d]->ncores))
            // elements, corresponding to the largest abase encountered.
        }
        serialize(mmt->pi->m);
        if (!mmt->pi->interleaved) {
            matmul_clear(Mloc->mm);
        } else if (mmt->pi->interleaved->idx == 1) {
            matmul_clear(Mloc->mm);
            /* group 0 is the first to leave, thus it doesn't to freeing.
            */
        }
        free(Mloc->locfile);
        free(Mloc->mname);
        free(Mloc->bname);
    }
    serialize(mmt->pi->m);
    free(mmt->matrices);
}

static permutation_data_ptr permutation_data_alloc()
{
    permutation_data_ptr a = malloc(sizeof(permutation_data));
    a->n = 0;
    a->alloc = 0;
    a->x = NULL;
    return a;
}

static void permutation_data_free(permutation_data_ptr a)
{
    if (!a) return;
    a->n = 0;
    a->alloc = 0;
    if (a->x) free(a->x);
    free(a);
}

static pthread_mutex_t pp = PTHREAD_MUTEX_INITIALIZER;

static void permutation_data_push(permutation_data_ptr a, unsigned int u[2])
{
    if (a->n >= a->alloc) {
        pthread_mutex_lock(&pp);
        a->alloc = a->n + 16 + a->alloc / 4;
        a->x = realloc(a->x, a->alloc * sizeof(unsigned int[2]));
        pthread_mutex_unlock(&pp);
        ASSERT_ALWAYS(a->x != NULL);
    }
    a->x[a->n][0] = u[0];
    a->x[a->n][1] = u[1];
    a->n++;
}

#ifdef  MVAPICH2_NUMVERSION
static void permutation_data_ensure(permutation_data_ptr a, size_t n)
{
    if (n >= a->alloc) {
        pthread_mutex_lock(&pp);
        a->alloc = n;
        a->x = realloc(a->x, n * sizeof(unsigned int[2]));
        pthread_mutex_unlock(&pp);
        ASSERT_ALWAYS(a->x != NULL);
    }
}
#endif
