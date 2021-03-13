#include "cado.h"
#include <mutex>
#include <vector>
#include <utility>
#include "macros.h"
#include "omp_proxy.h"
#include "sparse.h"
#include "merge_compute_weights.h"

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

#if 0
/* reference code */
void compute_weights_backend (filter_matrix_t *mat, index_t j0)
{
    uint64_t empty_cols = 0;
    uint64_t tot_weight = 0;
    col_weight_t *Wt[omp_get_max_threads()];
#pragma omp parallel reduction(+: empty_cols, tot_weight)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();

        /* we allocate an array of size mat->ncols, but the first j0 entries are unused */
        if (tid == 0)
            Wt[0] = mat->wt; /* trick: we use wt for Wt[0] */
        else
            Wt[tid] = malloc (mat->ncols * sizeof (col_weight_t));
        memset (Wt[tid] + j0, 0, (mat->ncols - j0) * sizeof (col_weight_t));

        /* Thread k accumulates weights in Wt[k].
           We only consider ideals of index >= j0, and put the weight of ideal j,
           j >= j0, in Wt[k][j]. */

        /* using a dynamic or guided schedule here is crucial, since during
           merge, the distribution of row lengths is no longer uniform
           (including discarded rows) */
        col_weight_t *Wtk = Wt[tid];
#pragma omp for schedule(guided)
        for (index_t i = 0; i < mat->nrows; i++) {
            if (mat->rows[i] == NULL) /* row was discarded */
                continue;
            for (index_t l = matLengthRow (mat, i); l >= 1; l--) {
                index_t j = matCell (mat, i, l);
                if (j < j0) /* assume ideals are sorted by increasing order */
                    break;
                else
                    Wtk[j]++;
            }
        }

        /* Thread k accumulates in Wt[0] the weights for the k-th block of columns,
           Wt[0][j] = Wt[0][j] + Wt[1][j] + ... + Wt[nthreads-1][j] */
        col_weight_t *Wt0 = Wt[0];
#pragma omp for schedule(static) /* slightly better than guided */
        for (index_t i = j0; i < mat->ncols; i++) {
            col_weight_t val = Wt0[i];
            for (int t = 1; t < T; t++)
                val += Wt[t][i];

            Wt0[i] = val;
            empty_cols += val == 0;
            tot_weight += val;
        }

        if (tid > 0)     /* start from 1 since Wt[0] = mat->wt + j0 should be kept */
            free (Wt[tid]);
    }

    mat->rem_ncols = mat->ncols - empty_cols;
    mat->tot_weight = tot_weight;
}
#endif

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

struct bucket : public std::vector<uint16_t> {
    private:
    size_t last_flush = 0;
    enum { ASYNC_FLUSH, PRIVATE_CACHE } mode = ASYNC_FLUSH;
    static constexpr const size_t B = 65536;
    public:
    template<typename... Args> bucket(Args&& ...args) : std::vector<uint16_t>(std::forward<Args>(args)...) {}
    void flush_raw(col_weight_t * dst, unsigned int max = UINT_MAX)
    {
        if (mode == PRIVATE_CACHE) {
            for(unsigned int i = 0 ; i < max && i < size() ; i++) {
                col_weight_t c = dst[i] + (*this)[i];
                if (c < dst[i])
                    dst[i] = std::numeric_limits<col_weight_t>::max();
                else
                    dst[i] = c;
            }
        } else {
            for(auto i : *this) {
                col_weight_t c = dst[i] + 1;
                if (c < dst[i])
                    dst[i] = std::numeric_limits<col_weight_t>::max();
                else
                    dst[i] = c;
            }
        }
        /* note that we don't even call clear(), nor zero out the
         * contents of cache[] !
         */
    }
    void flush(col_weight_t * dst, size_t when, int T, std::mutex & mm) {
        {
            std::lock_guard<std::mutex> dummy(mm);
            flush_raw(dst);
        }
        if (mode == PRIVATE_CACHE) {
            assign(1 << 16, 0);
        } else {
            clear();
        }
        /* if this holds, it means that every 16 rows on average, one of
         * the threads is scheduling a flush for this zone, which we
         * consider is a bit too much.
         */
        if (mode == ASYNC_FLUSH && when < last_flush + 16 * T) {
            /* convert this bucket to a cache */
            fprintf(stderr,
                    "# thread %d converts writes to %p from async-flush to private-cache\n",
                    omp_get_thread_num(), dst);
            mode = PRIVATE_CACHE;
            assign(1<<16, 0);
        }
        last_flush = when; 
    }
    inline void push(uint16_t j, col_weight_t * dst, size_t when, int T, std::mutex & mm) {
        if (mode == PRIVATE_CACHE) {
            /* still flush every 65536 lines */
            if ((when >> 16) > (last_flush >> 16))
                flush(dst, when, T, mm);
            (*this)[j]++;
        } else {
            push_back(j);
            if (size() == B)
                flush(dst, when, T, mm);
        }
    }
};

constexpr const size_t bucket::B;

typedef std::vector<bucket> bucket_array;

void compute_weights_backend (filter_matrix_t *mat, index_t j0)
{
    ASSERT_ALWAYS(j0 == 0);

    uint64_t empty_cols = 0;
    uint64_t tot_weight = 0;
    
    size_t number_of_slots = iceildiv(mat->ncols, 1 << 16);

    std::vector<std::mutex> locks(number_of_slots);

    std::vector<bucket_array> all_buckets(omp_get_max_threads());

    memset (mat->wt, 0, mat->ncols * sizeof (col_weight_t));

#pragma omp parallel reduction(+: empty_cols, tot_weight)
    {
        int T = omp_get_num_threads();
        int tid = omp_get_thread_num();

        all_buckets[tid] = bucket_array(number_of_slots);
        bucket_array & BA = all_buckets[tid];

#pragma omp for schedule(guided)
        for (index_t i = 0; i < mat->nrows; i++) {
            if (mat->rows[i] == NULL) /* row was discarded */
                continue;
            tot_weight += matLengthRow (mat, i);
            for (index_t l = matLengthRow (mat, i); l >= 1; l--) {
                index_t j = matCell (mat, i, l);
                index_t jh = j >> 16;
                uint16_t jl = j;
                /* This flushes evey so often */
                BA[jh].push(jl, mat->wt + (jh << 16), i, T, locks[jh]);
            }
        }

#pragma omp for schedule(static) /* slightly better than guided */
        for (index_t jh = 0; jh < number_of_slots ; jh++) {
            col_weight_t * where = mat->wt + (jh << 16);
            unsigned int max = MIN(1 << 16, mat->ncols - (jh << 16));
            for(int t = 0 ; t < T ; t++)
                all_buckets[t][jh].flush_raw(where, max);
            for(index_t jl = 0 ; jl < max ; jl++)
                empty_cols += where[jl] == 0;
        }
    }

    mat->rem_ncols = mat->ncols - empty_cols;
    mat->tot_weight = tot_weight;
}
