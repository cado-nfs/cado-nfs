/* sparse matrix "transposition"

Copyright 2019 Charles Bouillaguet.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado.h" // IWYU pragma: keep
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "omp_proxy.h"
/* the following should come after cado.h, which sets -Werror=all */
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
#include "typedefs.h"
#include "transpose.h"
#include "memory.h" // free_aligned

// #define BIG_BROTHER

/* OK, so we can do this the easy way, or the hard way. */

#ifdef TRANSPOSE_EASY_WAY
/* simple, cheap and dirty. 
   This is similar to algorithm 2 in "Parallel Transposition of Sparse Data Structures", 
   by Hao Wang, Weifeng Liu, Kaixi Hou, and Wu-chun Feng.
    
   It is in fact a direct parallelization of "distribution sorting" using atomic memory accesses.
   This works decently well, because there are very few conflicts.

   The row pointers MUST point to the END of each row */
void transpose(uint64_t nnz, index_t *Ai, index_t *Aj, index_t Rn, index_t *Rp, index_t *Ri)
{
        (void) Rn;
        /* dispatch entries */
        #pragma omp parallel for schedule(static)
        for (uint64_t k = 0; k < nnz; k++) {
                index_t i = Ai[k];
                index_t j = Aj[k];
                index_t s;
                #pragma omp atomic capture
                s = --Rp[j];
                Ri[s] = i;
        }
}
#else
/* The hard way. 
   Uses a parallel radix sort with a software write-combining buffer.
   Relies on aligned_alloc (OK in C11) and OpenMP. */

#define MAX_RADIX_BITS 10   /* was experimentally found to be OK */
#define MAX_PASSES 4 

/* cache-resident buffer for (i, j) pairs. One such entry per output bucket.
   Invariants:  row[CACHELINE_SIZE - 1] contains COUNT[...] for this bucket,
        col[CACHELINE_SIZE - 1] contains offset of the first entry in this buffer. 
   */

/* cache-resident buffer for (i, j) pairs. One such entry per output bucket.
   Invariants:  row[CACHELINE_SIZE - 1] contains COUNT[...] for this bucket,
        col[CACHELINE_SIZE - 1] contains offset of the first entry in this buffer. 
   */
struct cacheline_t {
    index_t row[CACHELINE_SIZE];
    index_t col[CACHELINE_SIZE];
};


/* Allocate the buffer */
static struct cacheline_t * wc_alloc()
{
  return malloc_aligned (sizeof(struct cacheline_t) *  (1 << MAX_RADIX_BITS), 64);
}

/* Setup the buffer for a new pass */
static inline void wc_prime(struct cacheline_t *buffer, const index_t *COUNT, int n_buckets)
{
    for (int i = 0; i < n_buckets; i++ ) {
        buffer[i].row[CACHELINE_SIZE - 1] = COUNT[i];
        buffer[i].col[CACHELINE_SIZE - 1] = COUNT[i] & (CACHELINE_SIZE - 1);
    }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this bucket is filled to the end */
static inline void wc_flush(struct cacheline_t * self, index_t count, index_t start, index_t *OUTi, index_t *OUTj)
{
    index_t target = count & ~(CACHELINE_SIZE - 1);
    if (start != 0) {   /* incomplete flush */
        for (int i = start; i < CACHELINE_SIZE; i++) {
            OUTi[target + i] = self->row[i];
            OUTj[target + i] = self->col[i];
        }
    } else {            /* complete cache line flush */
        store_nontemp_64B(OUTi + target, self->row);
        store_nontemp_64B(OUTj + target, self->col);
    }
    self->col[CACHELINE_SIZE - 1] = 0;
}

/* push an (i,j) pair into the buffer */
static inline void wc_push(index_t i, index_t j, struct cacheline_t * buffer, index_t bucket_idx, index_t *OUTi, index_t *OUTj)
{
    struct cacheline_t *self = buffer + bucket_idx;
    index_t count = self->row[CACHELINE_SIZE - 1];
    index_t start = self->col[CACHELINE_SIZE - 1];
    index_t slot = count & (CACHELINE_SIZE - 1);
    self->row[slot] = i;
    self->col[slot] = j;
    if (slot == CACHELINE_SIZE - 1)
        wc_flush(self, count, start, OUTi, OUTj);
    self->row[CACHELINE_SIZE - 1] = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_purge(struct cacheline_t * buffer, int n_buckets, index_t *OUTi, index_t *OUTj)
{
    for (int i = 0; i < n_buckets; i++) {
        index_t count = buffer[i].row[CACHELINE_SIZE - 1];
        index_t target = count & ~(CACHELINE_SIZE - 1);
        index_t start = buffer[i].col[CACHELINE_SIZE - 1];
        for (index_t j = target + start; j < count; j++) {
            OUTi[j] = buffer[i].row[j - target];
            OUTj[j] = buffer[i].col[j - target];
        }
    }
}


struct half_cacheline_t {
    index_t row[CACHELINE_SIZE];
};

/* Setup the buffer for a new pass */
static inline void wc_half_prime(struct half_cacheline_t *buffer, char *start, const index_t *COUNT, int n_buckets)
{
    for (int i = 0; i < n_buckets; i++ ) {
        buffer[i].row[CACHELINE_SIZE - 1] = COUNT[i];
        start[i] = COUNT[i] & (CACHELINE_SIZE - 1);
    }
}

/* Transfer data stored in this buffer to the output arrays.
   Assumption: this buckget is filled to the end */
static inline void wc_half_flush(struct half_cacheline_t * self, index_t count, char start, index_t *OUTi)
{
    index_t target = count & ~(CACHELINE_SIZE - 1);
    if (start == 0) {   /* complete cache line flush */
        store_nontemp_64B(OUTi + target, self->row);
    } else {            /* incomplete flush */
        for (int i = start; i < CACHELINE_SIZE; i++)
            OUTi[target + i] = self->row[i];
    }
}

/* push an (i,j) pair into the buffer */
static inline void wc_half_push(index_t i, struct half_cacheline_t * buffer, char *start, index_t bucket_idx, index_t *OUTi)
{
    struct half_cacheline_t *self = buffer + bucket_idx;
    index_t count = self->row[CACHELINE_SIZE - 1];
    index_t slot = count & (CACHELINE_SIZE - 1);
    self->row[slot] = i;
    if (slot == CACHELINE_SIZE - 1) {
        wc_half_flush(self, count, start[bucket_idx], OUTi);
        start[bucket_idx] = 0;
    }
    self->row[CACHELINE_SIZE - 1] = count + 1;
}

/* flush all buffer entries to the OUT arrays */
static inline void wc_half_purge(struct half_cacheline_t * buffer, char *start, int n_buckets, index_t *OUTi)
{
    for (int i = 0; i < n_buckets; i++) {
        index_t count = buffer[i].row[CACHELINE_SIZE - 1];
        index_t target = count & ~(CACHELINE_SIZE - 1);
        index_t start_ptr = start[i];
        for (index_t j = target + start_ptr; j < count; j++)
            OUTi[j] = buffer[i].row[j - target];
    }
}

struct ctx_t {
    int bits;
    int n_passes;   /* last pass is done... first, in parallel. */
    int radix[MAX_PASSES];
    index_t shift[MAX_PASSES];
    int n_buckets[MAX_PASSES];
    int pCOUNT[MAX_PASSES];
    index_t mask[MAX_PASSES];
    index_t *OUTi[MAX_PASSES];
    index_t *OUTj[MAX_PASSES];
    int seq_count_size;
    int par_count_size;
};



/* prepare the ctx object with information needed for all passes */
static void planification(struct ctx_t *ctx, index_t Rn, index_t nnz, index_t *scratch, index_t *Ri)
{
    int bits = 0;
    int tmp = Rn;
    while (tmp > 0) {
        bits++;
        tmp >>= 1;
    }
    int n = ceil((double) bits / MAX_RADIX_BITS);
    ctx->bits = bits;
    ctx->n_passes = n;

    /* the first pass is done on the most significant bits */
    ctx->radix[0] = ceil((double) bits / n);
    bits -= ctx->radix[0];
    ctx->shift[0] = bits;

    /* other passes are done on least significant bits first */
    int s_shift = 0;
    for (int p = 1; p < n; p++) {
        ctx->shift[p] = s_shift;
        int r = ceil((double) bits / (n - p)); /* bits in p-th pass */
        ctx->radix[p] = r;
        bits -= r;
        s_shift += r;
    }
    ASSERT (bits == 0);

    /* common to all passes */
    int s_count = 0;
    for (int p = 0; p < n; p++) {
        if (p == 1)
            s_count = 0;
        int size = 1 << ctx->radix[p];
        ctx->n_buckets[p] = size;
        ctx->pCOUNT[p] = s_count;
        s_count += size;
        ctx->mask[p] = size - 1;
        if ((n - p) & 1) {
            ctx->OUTi[p] = Ri;
            ctx->OUTj[p] = scratch;
        } else {
            ctx->OUTi[p] = scratch + ((nnz | 63) + 1);
            ctx->OUTj[p] = scratch + 2*((nnz | 63) + 1);
        }

        /* check alignment: the hardcoded value of 63 corresponds to the
           L1 cache size on modern cpus */
        unsigned long check MAYBE_UNUSED = (unsigned long) ctx->OUTi[p];
        ASSERT ((check & 63) == 0);
        check = (unsigned long) ctx->OUTj[p];
        ASSERT ((check & 63) == 0);
    }
    ctx->par_count_size = ctx->n_buckets[0];
    ctx->seq_count_size = s_count;
}


/* returns k such that buckets [0:k] are non-empty. */
static int partitioning(
            const struct ctx_t *ctx, 
            const index_t *Ai, 
            const index_t *Aj, 
            index_t nnz,
            struct cacheline_t *buffer,
            index_t *tCOUNT,
            index_t *gCOUNT)
{
    index_t mask = ctx->mask[0];
    index_t size = ctx->n_buckets[0];
    index_t shift = ctx->shift[0];
    index_t *OUTi = ctx->OUTi[0];
    index_t *OUTj = ctx->OUTj[0];
    ASSERT (size == mask + 1);

    /* parallel partitioning */
    int t = omp_get_thread_num();
    int T = omp_get_num_threads();

    index_t* COUNT = tCOUNT + t * size;
    memset(COUNT, 0, size * sizeof(*COUNT));

    #pragma omp for schedule(static)
    for (index_t k = 0; k < nnz; k++) {
        index_t j = Aj[k];
        index_t b = (j >> shift) & mask;
        COUNT[b]++;
    }
 
    index_t last = 0;
    #pragma omp master
    {
        index_t s = 0;
        for (index_t i = 0; i < size; i++) {
            gCOUNT[i] = s;
            for (int t = 0; t < T; t++) {
                index_t w = tCOUNT[t * size + i];
                tCOUNT[t * size + i] = s;
                s += w;
            }
            if (s > gCOUNT[i])
                last = i;
        }
        gCOUNT[size] = s;
        ASSERT (s == nnz);
    }

    #pragma omp barrier

    wc_prime(buffer, COUNT, size);
    #pragma omp for schedule(static) nowait
    for (index_t k = 0; k < nnz; k++) {
        index_t i = Ai[k];
        index_t j = Aj[k];
        index_t b = (j >> shift) & mask;
        wc_push(i, j, buffer, b, OUTi, OUTj);
    }
    wc_purge(buffer, size, OUTi, OUTj);

    /* correctness check
    #pragma omp barrier
    for (index_t b = 0; b < size; b++)
        for (index_t k = gCOUNT[b]; k < gCOUNT[b + 1]; k++) {
            index_t i = OUTi[k];
            index_t j = OUTj[k];
            index_t c = (j >> shift) & mask;
            if (c != b) {
                printf("discrepancy in bucket %d for index %d: found (%d, %d) which should be in bucket %d\n", 
                    b, k, i, j, c);
                ASSERT(0);
            }
        }
    */
    return last + 1;
}


static void histogram(struct ctx_t *ctx, const index_t *Aj, index_t lo, 
                        index_t hi, int n, index_t **W)
{
    index_t *shift = ctx->shift;
    index_t *mask = ctx->mask;

    switch (n) {
    case 2:
        for (index_t k = lo; k < hi; k++) {
            index_t j = Aj[k];
            int q1 = j & mask[1];
            W[1][q1]++;
        }
        break;
    case 3:
        for (index_t k = lo; k < hi; k++) {
            index_t j = Aj[k];
            int q1 = j & mask[1];
            int q2 = (j >> shift[2]) & mask[2];
            W[1][q1]++;
            W[2][q2]++;
        }
        break;
    case 4:
        for (index_t k = lo; k < hi; k++) {
            index_t j = Aj[k];
            int q1 = j & mask[1];
            int q2 = (j >> shift[2]) & mask[2];
            int q3 = (j >> shift[3]) & mask[3];
            W[1][q1]++;
            W[2][q2]++;
            W[3][q3]++;
        }
        break;
#if __SIZEOF_INDEX__ == 8
    case 5:
        #pragma omp for schedule(static)
        for (index_t k = lo; k < hi; k++) {
            index_t j = Aj[k];
            int q1 = j & mask[1];
            int q2 = (j >> shift[2]) & mask[2];
            int q3 = (j >> shift[3]) & mask[3];
            int q4 = (j >> shift[4]) & mask[4];
            W[1][q1]++;
            W[2][q2]++;
            W[3][q3]++;
            W[4][q4]++;
        }
        break;
#endif
    default:
        printf("Ask the programmer to hardcode more passes in (radix) transpose...\n");
        ASSERT_ALWAYS(0);
    }
}


/* sequentially transpose a single bucket */
void transpose_bucket(struct ctx_t *ctx, struct cacheline_t *buffer, index_t lo, index_t hi)
{
    int n = ctx->n_passes;
    index_t csize = ctx->seq_count_size;
    index_t COUNT[csize];
    index_t *W[n];
    
    for (int p = 1; p < n; p++)
        W[p] = COUNT + ctx->pCOUNT[p];
    
    index_t *INi = ctx->OUTi[0];
    index_t *INj = ctx->OUTj[0];

    // printf("histograming [%d:%d] with %d passes\n", lo, hi, n);
    memset(COUNT, 0, csize * sizeof(*COUNT));
    histogram(ctx, INj, lo, hi, n, W);

    for (int p = 1; p < n; p++) {
        index_t shift = ctx->shift[p];
        index_t mask = ctx->mask[p];
        index_t *OUTi = ctx->OUTi[p];
        index_t *OUTj = ctx->OUTj[p];

        /* prefix-sum */
        index_t s = lo;
        index_t size = ctx->n_buckets[p];
        for (index_t i = 0; i < size; i++) {
            index_t w = W[p][i];
            W[p][i] = s;
            s += w;
        }

        if (p < n - 1) {   
            /* full pass */
            wc_prime(buffer, W[p], size);
            for (index_t k = lo; k < hi; k++) {
                index_t i = INi[k];
                index_t j = INj[k];
                int b = (j >> shift) & mask;
                wc_push(i, j, buffer, b, OUTi, OUTj);
            }
            wc_purge(buffer, size, OUTi, OUTj);
        } else {        
            /* last pass: we don't need to write the j values */
            struct half_cacheline_t *half_buffer = (struct half_cacheline_t *) buffer;
            char start[size];
            wc_half_prime(half_buffer, start, W[p], size);
            for (index_t k = lo; k < hi; k++) {
                index_t i = INi[k];
                index_t j = INj[k];
                int b = (j >> shift) & mask;
                wc_half_push(i, half_buffer, start, b, OUTi);
            }
            wc_half_purge(half_buffer, start, size, OUTi);
        }

        /* check 
        for (index_t b = 0; b < size - 1; b++)
            for (index_t k = W[p][b]; k < W[p][b + 1]; k++) {
                index_t i = OUTi[k];
                index_t j = OUTj[k];
                index_t c = (j >> shift) & mask;
                if (c != b) {
                    printf("pass %d, discrepancy in bucket %d for index %d: found (%d, %d) which should be in bucket %d\n", 
                        p, b, k, i, j, c);
                    ASSERT(0);
                }
            }
        */
        INi = OUTi;
        INj = OUTj;
    }
}


void transpose(uint64_t nnz, index_t *Ai, index_t *Aj, index_t Rn, index_t *Rp, index_t *Ri)
{
    (void) Rp;
    
    struct ctx_t ctx;
    index_t *scratch = malloc_aligned (3 * ((nnz | 63) + 1) * sizeof(index_t), 64);
    planification(&ctx, Rn, nnz, scratch, Ri);


#ifdef BIG_BROTHER
    printf("$$$     transpose:\n");
    printf("$$$        bits: %d\n", ctx.bits);
    printf("$$$        passes: \n");
    for (int p = 0; p < ctx.n_passes; p++)
        printf("$$$         - %d \n", ctx.radix[p]);
#endif

    int size = ctx.par_count_size;
    int T = omp_get_max_threads();
    index_t tCOUNT[T * size];
    index_t gCOUNT[size + 1];

    int non_empty;

    #pragma omp parallel
    {
        struct cacheline_t * buffer = wc_alloc();
        #ifdef BIG_BROTHER
                double start = wct_seconds();
        #endif
        int tmp = partitioning(&ctx, Ai, Aj, nnz, buffer, tCOUNT, gCOUNT);

        #pragma omp master
        {
                non_empty = tmp;
                #ifdef BIG_BROTHER
                        printf("$$$        buckets: %d\n", non_empty);
                        printf("$$$        partitioning-wct: %.2f\n", wct_seconds() - start);
                        index_t smallest = nnz;
                        index_t biggest = 0;
                        for (int i = 0; i < non_empty; i++) {
                                index_t size = gCOUNT[i + 1] - gCOUNT[i];
                                smallest = MIN(smallest, size);
                                biggest = MAX(biggest, size);
                        }
                        printf("$$$        smallest-bucket: %d\n", smallest);
                        printf("$$$        avg-bucket: %" PRId64 "\n", nnz / non_empty);
                        printf("$$$        biggest-bucket: %d\n", biggest);
                #endif
        }

        #pragma omp barrier
        #ifdef BIG_BROTHER
        double sub_start = wct_seconds(); 
        #endif

        if (ctx.n_passes > 1)
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < non_empty; i++)
                transpose_bucket(&ctx, buffer, gCOUNT[i], gCOUNT[i + 1]);
        free_aligned (buffer);

        #pragma omp master
        {
                #ifdef BIG_BROTHER
                        printf("$$$        buckets-wct: %.2f\n", wct_seconds() - sub_start);
                #endif
        }

    }
    free_aligned (scratch);
}      
#endif
