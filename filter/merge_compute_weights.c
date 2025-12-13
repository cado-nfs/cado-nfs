#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <pthread.h>

#include "macros.h"
#include "merge_compute_weights.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "misc.h"
#include "omp_proxy.h"
#include "typedefs.h"

/* 
 * This does a pass on the matrix data (all rows), and collects the
 * following info
 *
 * mat->wt[]
 * mat->rem_ncols
 * mat->tot_weight
 *
 * In the general case, column weights in mat->wt need only be computed
 * up to cwmax + 1 since we only need to know whether the weights are <=
 * cwmax or not).
 *
 * However, this is not true for the shrink case, where a full count is
 * needed in order to accurately compute the density.
 *
 */

/* each thread has a bucket array that consists of one bucket of
 * uint16_t's for each region of 2^16 indices, and this bucket can hold
 * up to B updates
 *
 * this means nthreads * ncols / 2^16 * B * 2 bytes
 *
 * say if B=2^16, ncols=2G, ntreads=128, this means 512G in total,
 * which seems acceptable for such a large problem (rsa-250 is in
 * this ballpark). And it is always possible to decrease B. Note also
 * that we don't have to pre-allocate all arrays.
 *
 */

#define bucket_B        (1 << 16)
struct bucket_s {
    uint16_t * data;
    size_t size;
    size_t alloc;
    size_t last_flush;
    enum { ASYNC_FLUSH, PRIVATE_CACHE } mode;
};
typedef struct bucket_s bucket[1];
typedef struct bucket_s * bucket_ptr;
typedef const struct bucket_s * bucket_srcptr;

void bucket_init(bucket_ptr bu)
{
    bu->data = NULL;
    bu->size = bu->alloc = 0;
    bu->last_flush = 0;
    bu->mode = ASYNC_FLUSH;
}

void bucket_clear(bucket_ptr bu)
{
    free(bu->data);
    bu->data = NULL;
    bu->size = bu->alloc = 0;
    bu->last_flush = 0;
    bu->mode = ASYNC_FLUSH;
}

void bucket_flush_raw(bucket_ptr bu, col_weight_t * dst, unsigned int max)
{
    max = MIN(max, bu->size);
    if (bu->mode == PRIVATE_CACHE) {
        for(unsigned int i = 0 ; i < max ; i++) {
            col_weight_t c = dst[i] + bu->data[i];
            if (c < dst[i])
                dst[i] = UMAX(col_weight_t);
            else
                dst[i] = c;
        }
    } else {
        for(size_t j = 0 ; j < bu->size ; j++) {
            uint16_t i = bu->data[j];
            col_weight_t c = dst[i] + 1;
            if (c < dst[i])
                dst[i] = UMAX(col_weight_t);
            else
                dst[i] = c;
        }
    }
    /* note that we don't even call clear(), nor zero out the
     * contents of cache[] !
     */
}
void bucket_flush(bucket_ptr bu, col_weight_t * dst, size_t when, int T, pthread_mutex_t * mm)
{
    {
        pthread_mutex_lock(mm);
        bucket_flush_raw(bu, dst, UINT_MAX);
        pthread_mutex_unlock(mm);
    }
    if (bu->mode == PRIVATE_CACHE) {
        ASSERT(bu->size == (1 << 16));
        memset(bu->data, 0, (1 << 16) * sizeof(uint16_t));
    } else {
        bu->size = 0;
    }
    /* if this holds, it means that every 16 rows on average, one of
     * the threads is scheduling a flush for this zone, which we
     * consider is a bit too much.
     */
    if (bu->mode == ASYNC_FLUSH && when < bu->last_flush + 16 * T) {
        /* convert this bucket to a cache */
        fprintf(stderr,
                "# thread %d converts writes to %p from async-flush to private-cache\n",
                omp_get_thread_num(), dst);
        bu->mode = PRIVATE_CACHE;
        if (bu->alloc < (1 << 16)) {
            bu->alloc = 1 << 16;
            CHECKED_REALLOC(bu->data, bu->alloc, uint16_t);
        }
        bu->size = 1 << 16;
        memset(bu->data, 0, (1 << 16) * sizeof(uint16_t));
    }
    bu->last_flush = when; 
}

static inline void bucket_push(bucket_ptr bu, uint16_t j, col_weight_t * dst, size_t when, int T, pthread_mutex_t * mm) {
    if (bu->mode == PRIVATE_CACHE) {
        /* still flush every 65536 lines */
        if ((when >> 16) > (bu->last_flush >> 16))
            bucket_flush(bu, dst, when, T, mm);
        bu->data[j]++;
    } else {
        if (bu->size >= bu->alloc) {
            bu->alloc = MAX(16, 2 * bu->alloc);
            CHECKED_REALLOC(bu->data, bu->alloc, uint16_t);
        }
        bu->data[bu->size++] = j;
        if (bu->size == bucket_B)
            bucket_flush(bu, dst, when, T, mm);
    }
}

void compute_weights_backend (filter_matrix_t *mat, index_t j0)
{
    ASSERT_ALWAYS(j0 == 0);

    uint64_t empty_cols = 0;
    uint64_t tot_weight = 0;
    
    size_t number_of_slots = iceildiv(mat->ncols, 1 << 16);

    pthread_mutex_t * locks = malloc(number_of_slots * sizeof(pthread_mutex_t));
    for(size_t i = 0 ; i < number_of_slots ; i++)
        pthread_mutex_init(&locks[i], NULL);

    bucket ** all_buckets = malloc(omp_get_max_threads() * sizeof(bucket *));

    memset (mat->wt, 0, mat->ncols * sizeof (col_weight_t));

#pragma omp parallel reduction(+: empty_cols, tot_weight)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();

        all_buckets[tid] = malloc(number_of_slots * sizeof(bucket));

        for(size_t i = 0 ; i < number_of_slots ; i++)
            bucket_init(all_buckets[tid][i]);

        bucket * BA = all_buckets[tid];

#pragma omp for schedule(guided)
        for (index_t i = 0; i < mat->nrows; i++) {
            if (mat->rows[i] == NULL) /* row was discarded */
                continue;
            tot_weight += matLengthRow(mat, i);
            for (index_t l = matLengthRow (mat, i); l >= 1; l--) {
                index_t j = matCell (mat, i, l);
                index_t jh = j >> 16;
                uint16_t jl = j;
                /* This flushes evey so often */
                bucket_push(BA[jh], jl, mat->wt + (jh << 16), i, T, &locks[jh]);
            }
        }

#pragma omp for schedule(static) /* slightly better than guided */
        for (index_t jh = 0; jh < number_of_slots ; jh++) {
            col_weight_t * where = mat->wt + (jh << 16);
            unsigned int max = MIN(1 << 16, mat->ncols - (jh << 16));
            for(int t = 0 ; t < T ; t++)
                bucket_flush_raw(all_buckets[t][jh], where, max);
            for(index_t jl = 0 ; jl < max ; jl++)
                empty_cols += where[jl] == 0;
        }

        for(size_t i = 0 ; i < number_of_slots ; i++)
            bucket_clear(all_buckets[tid][i]);

        free(all_buckets[tid]);
    }

    mat->rem_ncols = mat->ncols - empty_cols;
    mat->tot_weight = tot_weight;

    free(all_buckets);
    for(size_t i = 0 ; i < number_of_slots ; i++)
        pthread_mutex_destroy(&locks[i]);
    free(locks);
}
