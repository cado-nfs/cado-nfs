#include "cado.h"
#include "bwc_config.h"
#include "xymats.h"

/* This is used as a very simple reduction routine, mostly for the small
 * m*n matrices. It has nothing to do with the more involved reduce
 * operations relevant to vectors. Here, the data we're dealing with is
 * tiny.
 *
 * This operation serializes threads.
 */
#if 0
void allreduce_generic_threadlevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (v->flags & THREAD_SHARED_VECTOR)
        return; /* If the vector is shared, reduction is pointless */

    /* Do it star-like. Easy enough. */
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                abadd(abase, v->v + aboffset(abase,k),
                        v->all_v[t] + aboffset(abase,k));
            }
        }
    }
    serialize_threads(wr);
    if (wr->trank > 0) {
        abcopy(abase, v->v, v->all_v[0], n);
    }
}
#endif

static void allreduce_generic_mpilevel(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (wr->trank == 0) {
        size_t siz = abbytes(abase, n);
        int err = MPI_Allreduce(MPI_IN_PLACE, v->v, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);

    if (wr->trank > 0 && (v->flags & THREAD_SHARED_VECTOR) == 0) {
        /* recover the reduced data from thread zero */
        abcopy(abase, v->v, v->all_v[0], n);
    }
}

/* This combines both, with fewer thread locks */
void allreduce_generic(abobj_t abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        allreduce_generic_mpilevel(abase, v, wr, n);
    }

    /* Do it star-like. Easy enough. */
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            for(unsigned int k = 0 ; k < n ; k++) {
                abadd(abase, v->v + aboffset(abase,k),
                        v->all_v[t] + aboffset(abase,k));
            }
        }
        size_t siz = abbytes(abase, n);
        int err = MPI_Allreduce(MPI_IN_PLACE, v->v, siz, MPI_BYTE, MPI_BXOR, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);
    if (wr->trank > 0) {
        abcopy(abase, v->v, v->all_v[0], n);
    }
}

void broadcast_generic_threadlevel(mmt_generic_vec_ptr v, pi_wiring_ptr wr, size_t siz, unsigned int t0)
{
    if (v->flags & THREAD_SHARED_VECTOR)
        return; /* If the vector is shared, reduction is pointless */

    serialize_threads(wr);
    /* Do it star-like. Easy enough. */
    if (wr->trank != t0) // valgrind is picky sometimes.
        memcpy(v->v, v->all_v[t0], siz);
}

void broadcast_generic_mpilevel(mmt_generic_vec_ptr v, pi_wiring_ptr wr, size_t siz, unsigned int j0)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (wr->trank == 0) {
            int err = MPI_Bcast(v->v, siz, MPI_BYTE, j0, wr->pals);
            ASSERT_ALWAYS(!err);
        }
        return;
    }

    SEVERAL_THREADS_PLAY_MPI_BEGIN(wr) {
        int err = MPI_Bcast(v->v, siz, MPI_BYTE, j0, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
}

void broadcast_generic(mmt_generic_vec_ptr v, pi_wiring_ptr wr, size_t siz, unsigned int j0, unsigned int t0)
{
    broadcast_generic_threadlevel(v, wr, siz, t0);
    broadcast_generic_mpilevel(v, wr, siz, j0);
}
