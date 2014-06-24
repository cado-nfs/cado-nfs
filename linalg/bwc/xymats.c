#include "cado.h"
#include "bwc_config.h"
#include "xymats.h"

#if 0
void allreduce_generic_threadlevel(abobj abase, mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
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

static void allreduce_generic_mpilevel(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (wr->trank == 0) {
        int err = MPI_Allreduce(MPI_IN_PLACE,
                v->v, n,
                v->abase->mpi_datatype(v->abase),
                v->abase->mpi_addition_op(v->abase),
                wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);

    if (wr->trank > 0 && (v->flags & THREAD_SHARED_VECTOR) == 0) {
        /* recover the reduced data from thread zero */
        v->abase->vec_set(v->abase, v->v, v->all_v[0], n);
    }
}

/* This combines both, with fewer thread locks */
void allreduce_generic(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        allreduce_generic_mpilevel(v, wr, n);
        abort(); // weird. I'm adding a return here, because I'm a
        // bit puzzled by its absence. Probably this code path is unused.
        return;
    }

    /* Do it star-like. Easy enough. */
    serialize_threads(wr);
    if (wr->trank == 0) {
        for(unsigned int t = 1 ; t < wr->ncores ; t++) {
            v->abase->vec_add(v->abase, v->v, v->v, v->all_v[t], n);
        }
        int err = MPI_Allreduce(MPI_IN_PLACE,
                v->v, n,
                v->abase->mpi_datatype(v->abase),
                v->abase->mpi_addition_op(v->abase),
                wr->pals);
        ASSERT_ALWAYS(!err);
    }
    serialize_threads(wr);
    if (wr->trank > 0) {
        v->abase->vec_set(v->abase, v->v, v->all_v[0], n);
    }
    serialize_threads(wr);
}

void broadcast_generic_threadlevel(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n, unsigned int t0)
{
    if (v->flags & THREAD_SHARED_VECTOR)
        return; /* If the vector is shared, reduction is pointless */

    serialize_threads(wr);
    /* Do it star-like. Easy enough. */
    if (wr->trank != t0) // valgrind is picky sometimes.
        v->abase->vec_set(v->abase, v->v, v->all_v[t0], n);
}

void broadcast_generic_mpilevel(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n, unsigned int j0)
{
    if (v->flags & THREAD_SHARED_VECTOR) {
        if (wr->trank == 0) {
            int err = MPI_Bcast(v->v,
                    n,
                    v->abase->mpi_datatype(v->abase),
                    j0, wr->pals);
            ASSERT_ALWAYS(!err);
        }
        return;
    }

    SEVERAL_THREADS_PLAY_MPI_BEGIN(wr) {
        int err = MPI_Bcast(v->v, n,
                v->abase->mpi_datatype(v->abase),
                j0, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
}

void broadcast_generic(mmt_vec_ptr v, pi_wiring_ptr wr, unsigned int n, unsigned int j0, unsigned int t0)
{
    broadcast_generic_threadlevel(v, wr, n, t0);
    broadcast_generic_mpilevel(v, wr, n, j0);
}
