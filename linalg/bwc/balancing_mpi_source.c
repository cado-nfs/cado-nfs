#include "cado.h"
#include <string.h>

#include "select_mpi.h"
#include "balancing_mpi_source.h"

size_t mpi_source_get(mpi_source_ptr s, uint32_t ** p, size_t avail)
{
    // fprintf(stderr, "mpi_source_get (over == %d/%d) !\n", s->over, s->nparallel);
    if (*p == NULL || avail == 0) {
        if (s->buf == NULL) {
            s->buf = malloc(s->size * sizeof(uint32_t));
        }
        *p = s->buf;
    } else {
        ASSERT_ALWAYS(avail >= s->size);
    }
    if (s->over >= s->nparallel) {
        // fprintf(stderr, "mpi_source_get called on closed source\n");
        return 0;
    }
    MPI_Status x;
    memset(&x, 0, sizeof(x));
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // fprintf(stderr, "[%d] MPI_Recv(%zu,%d)\n", rank, s->size, s->peer);
    // fprintf(stderr, "MPI_Recv(%zu,%d)\n", s->size, s->peer);
    struct {
        size_t sz;
        size_t over;
    } sz_info[1];
    ASSERT_ALWAYS(sizeof(sz_info) == s->tailsize * sizeof(uint32_t));
    MPI_Recv(*p, s->size * sizeof(uint32_t), MPI_BYTE,
            s->peer, MPI_ANY_TAG, MPI_COMM_WORLD, &x);
    memcpy(sz_info, *p + s->size - s->tailsize, sizeof(sz_info));
    s->over += sz_info->over != 0;

    s->b->pos += sz_info->sz;
#ifdef HAVE_MPI
    s->tag = x.MPI_TAG;
#else
    abort(); // I believe that by now, we should never arrive here.
    s->tag = 0;
#endif
    // fprintf(stderr, "receiving tag %d, length %zu (over: %d)\n",
            // s->tag, sz_info->sz, sz_info->over != 0);
    return sz_info->sz;
}

data_source_ptr mpi_source_alloc(int peer, size_t queue_size)
{
    mpi_source_ptr s = malloc(sizeof(mpi_source));
    memset(s, 0, sizeof(mpi_source));
    s->b->get = (size_t(*)(void*,uint32_t**,size_t)) mpi_source_get;
    s->tag = 0;
    s->peer = peer;
    s->size = queue_size;
    s->nparallel = 1;   // at least we have a default value...
    ASSERT_ALWAYS(sizeof(size_t) % sizeof(uint32_t) == 0);
    s->tailsize = 2*iceildiv(sizeof(size_t),sizeof(uint32_t));
    return (data_source_ptr) s;
}

void mpi_source_free(data_source_ptr q)
{
    mpi_source_ptr p = (mpi_source_ptr) q;
    if (p->buf)
        free(p->buf);
}

