#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <pthread.h>
#include "bwc_config.h"
#include "matmul.h"
#include "matmul-threaded.h"
#include "matmul-common.h"
#include "abase.h"
#include "worker-threads.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-threaded"
#define MM_MAGIC_FAMILY        0xa002UL
#define MM_MAGIC_VERSION       0x1005UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* These are just reference values. Everything except MM_DGRP_SIZE is
 * tunable (but for some of them, it forces to rebuild the cache file:
 * DENSIFY_TOLERANCE SGRP_SIZE BLOCKSIZE) */
#define MM_DENSIFY_TOLERANCE    0.25     /* Set to zero to skip dense blocks */
#define MM_DGRP_SIZE            32
#define MM_NTHREADS             4
#define MM_SGRP_SIZE            64
#define MM_STORE_TRANSPOSED     1
#define MM_PREFETCH_OFFSET_1    32
#define MM_PREFETCH_OFFSET_2    16
#define MM_PREFETCH_OFFSET_3    16
#define MM_BLOCKSIZE            8
#define MM_CACHE_LINE_SIZE      64

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

struct matmul_threaded_data_s {
    struct matmul_public_s public_[1];

    /* Now our private fields */
    abobj_t xab;

    struct {
        uint32_t n;
        uint32_t * q;
        uint32_t weight;
        /* This is protected by the mutex, and used for serializing the
         * writing of the data area */
        uint32_t processed;
        pthread_cond_t recheck_please;
    } dense[1];
    struct {
        uint32_t n;
        uint32_t * len;
        uint32_t tlen;
        uint32_t * q;
    } sparse[1];

    unsigned int sgrp_size;
    unsigned int nthreads;
    unsigned int blocksize;
    unsigned int store_transposed;
    unsigned int off1;
    unsigned int off2;
    unsigned int off3;
    double densify_tolerance;

    struct worker_threads_group * tg;

    void * dst;
    const void * src;
    int d;
};

extern void matmul_threaded_mul_sub_dense(struct matmul_threaded_data_s * mm, abt * dst, abt const * src, int d, int i);
extern void matmul_threaded_mul_sub_sparse(struct matmul_threaded_data_s * mm, abt * dst, abt const * src, int d, int i);

void matmul_threaded_clear(struct matmul_threaded_data_s * mm)
{
    matmul_common_clear(mm->public_);

    worker_threads_clear(mm->tg);
    pthread_cond_destroy(&mm->dense->recheck_please);

    free(mm->dense->q);
    free(mm->sparse->len);
    free(mm->sparse->q);
    free(mm);
}

struct matmul_threaded_data_s * matmul_threaded_init(abobj_ptr xx MAYBE_UNUSED, param_list pl, int optimized_direction)
{
    struct matmul_threaded_data_s * mm;
    mm = malloc(sizeof(struct matmul_threaded_data_s));
    memset(mm, 0, sizeof(struct matmul_threaded_data_s));
    abobj_init_set(mm->xab, xx);

    int suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    mm->public_->store_transposed = suggest;
    if (pl) {
        param_list_parse_uint(pl, "mm_store_transposed", 
                &mm->public_->store_transposed);
        if (mm->public_->store_transposed != (unsigned int) suggest) {
            fprintf(stderr, "Warning, mm_store_transposed"
                    " overrides suggested matrix storage ordering\n");
        }           
    }   


    mm->nthreads = MM_NTHREADS;
    mm->densify_tolerance = MM_DENSIFY_TOLERANCE;
    mm->sgrp_size = MM_SGRP_SIZE;
    mm->blocksize = MM_CACHE_LINE_SIZE;
    mm->off1 = MM_PREFETCH_OFFSET_1;
    mm->off2 = MM_PREFETCH_OFFSET_2;
    mm->off3 = MM_PREFETCH_OFFSET_3;

    if (pl) {
        param_list_parse_uint(pl, "cache_line_size", &mm->blocksize);
        param_list_parse_uint(pl, "mm_threaded_nthreads", &mm->nthreads);
        param_list_parse_double(pl, "mm_threaded_densify_tolerance", &mm->densify_tolerance);
        param_list_parse_uint(pl, "mm_threaded_sgroup_size", &mm->sgrp_size);
        param_list_parse_uint(pl, "mm_threaded_offset1", &mm->off1);
        param_list_parse_uint(pl, "mm_threaded_offset2", &mm->off2);
        param_list_parse_uint(pl, "mm_threaded_offset3", &mm->off3);
    }

    mm->blocksize /= abbytes(mm->xab, 1);
    
    return mm;
}

void matmul_threaded_blocks_info(struct matmul_threaded_data_s * mm)
{
    printf("// %u dense blocks of %d %ss\n",
            mm->dense->n, MM_DGRP_SIZE, rowcol[mm->store_transposed]);
    unsigned long sparse_coeffs = mm->public_->ncoeffs - mm->dense->weight;
    unsigned long padding = mm->sparse->tlen - sparse_coeffs;
    printf("// %u sparse blocks of %d %ss ; %lu coeffs + %lu padding (+%.2f%%)\n",
            mm->sparse->n, mm->sgrp_size,
            rowcol[mm->store_transposed],
            sparse_coeffs, padding,
            (100.0 * padding / sparse_coeffs));
}

void matmul_threaded_build_cache(struct matmul_threaded_data_s * mm, uint32_t * data)
{
    if (!data)
        data = matmul_common_read_stupid_data(mm->public_);

    unsigned int nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    unsigned int ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    /* We give some notations and example values.  */

    // unsigned int n = mm->nthreads;      /* e.g. 4 */
    unsigned int g = mm->sgrp_size;    /* e.g. 32 */
    // unsigned int b = mm->blocksize;     /* e.g. 16 */
    ASSERT_ALWAYS(g % (mm->nthreads*mm->blocksize) == 0);
    size_t sparse_readahead = mm->nthreads * MAX(mm->off1, mm->off2) * mm->blocksize;

    /* count dense and sparse blocks. Since one entry in the sparse part
     * occupies a 32-bit location, we consider that it's preferrable to
     * use a dense mode when it's no heavier. For 32-bit consecutive
     * rows, we then store a sequence of N 32-bit integers, whose bit
     * indicates entries in the N columns of the block.
     *
     * 32 is somehow a constant here. But we may allow some waste
     * storage, say up to a factor of 2. That's called
     * mm->densify_tolerance
     */
    uint32_t * ptr = data;
    unsigned int i = 0;
    /* Count dense blocks */
    for( ; i < nrows_t ; mm->dense->n++) {
        unsigned int di;
        unsigned int weight = 0;
        uint32_t * q = ptr;
        for(di = 0 ; di < MM_DGRP_SIZE && i+di < nrows_t ; di++) {
            weight += *q;
            q += 1 + *q;
        }
        int dense = weight * mm->densify_tolerance >= ncols_t;
        printf("%ss %u..%u ; total %u coeffs = n%ss(%u) * %.2f ;%s dense\n",
                rowcol[mm->public_->store_transposed],
                i, i+di-1, weight,
                rowcol[!mm->public_->store_transposed],
                mm->public_->dim[!mm->public_->store_transposed],
                (double) weight / ncols_t,
                dense ? "" : " not");
        if (!dense)
            break;
        mm->public_->ncoeffs += weight;
        mm->dense->weight += weight;
        ptr = q;
        i += di;
    }
    /* Count sparse blocks */
    for( ; i < nrows_t ; mm->sparse->n++) {
        unsigned long slice_coeffs = 0;
        uint32_t mlen = 0;
        unsigned int di;
        for(di = 0 ; di < g && i + di < nrows_t ; di++) {
            if (*ptr >= mlen)
                mlen = *ptr;
            slice_coeffs += *ptr;
            ptr += 1 + *ptr;
        }
        mm->public_->ncoeffs += slice_coeffs;
        mm->sparse->tlen += g * mlen;
        i += di;
    }

    matmul_threaded_blocks_info(mm);

    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, mm->blocksize);

    mm->dense->q = malloc(mm->dense->n * dblock_stride * sizeof(uint32_t));
    mm->sparse->q = malloc((mm->sparse->tlen + sparse_readahead) * sizeof(uint32_t));
    mm->sparse->len = malloc(mm->sparse->n * sizeof(uint32_t));

    /* reset pointer to beginning, start to fill dense blocks */
    ptr = data;
    i = 0;
    memset(mm->dense->q, 0, mm->dense->n * dblock_stride * sizeof(uint32_t));
    for(unsigned int blocknum = 0 ; blocknum < mm->dense->n ; blocknum++) {
        uint32_t * q = mm->dense->q + blocknum * dblock_stride;
        unsigned int di;
        for(di = 0 ; di < MM_DGRP_SIZE && i + di < nrows_t ; di++) {
            uint32_t rlen = *ptr++;
            for(unsigned int j = 0 ; j < rlen ; j++) {
                q[*ptr++] ^= ((uint32_t) 1) << di;
            }
        }
        i += di;
    }

    /* Now proceed with sparse blocks */
    uint32_t * sq = mm->sparse->q;
    /* First fill padding data */
    for(unsigned int j = 0 ; j < mm->sparse->tlen ; j++)
        sq[j] = ncols_t;
    for(unsigned int blocknum = 0 ; blocknum < mm->sparse->n ; blocknum++) {
        unsigned int di;
        /* We'll discover the max size in the strip once again */
        uint32_t mlen = 0;

        /* striding for the entries in one row is g==sgrp_size.
         * Entry g*l+k corresponds to index l in row k.  A
         * complete block of g indices is covered by the n threads in
         * g/(nb) blocks of b rows (g/nb blocks from each thread).
         */

        for(di = 0 ; di < g ; di++) {
            uint32_t rlen;
            
            if (i + di < nrows_t) {
                rlen = *ptr++;
            } else {
                rlen = 0;
            }

            ASSERT(sq-mm->sparse->q+rlen*g <= mm->sparse->tlen);

            for(unsigned int j = 0 ; j < rlen ; j++) sq[di + j * g] = *ptr++;

            if (rlen >= mlen) mlen = rlen;
        }

        mm->sparse->len[blocknum] = mlen * g;
        i += di;
        sq += mlen * g;
    }

    memset(sq, 0, sparse_readahead * sizeof(uint32_t));

    free(data);

    mm->tg = worker_threads_init(mm->nthreads);
    pthread_cond_init(&mm->dense->recheck_please, NULL);
}

int matmul_threaded_reload_cache(struct matmul_threaded_data_s * mm)
{
    FILE * f;

    f = matmul_common_reload_cache_fopen(abbytes(mm->xab, 1), mm->public_, MM_MAGIC);
    if (f == NULL) return 0;

    MATMUL_COMMON_READ_ONE32(mm->dense->weight, f);
    MATMUL_COMMON_READ_ONE32(mm->sgrp_size, f);
    MATMUL_COMMON_READ_ONE32(mm->dense->n, f);
    MATMUL_COMMON_READ_ONE32(mm->sparse->n, f);
    MATMUL_COMMON_READ_ONE32(mm->sparse->tlen, f);

    unsigned int ncols_t = mm->public_->dim[!mm->public_->store_transposed];
    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, mm->blocksize);

    /* sparse_readahead depends on data which is initialized by us, not
     * recovered from the data file */
    size_t sparse_readahead;
    sparse_readahead = mm->nthreads * MAX(mm->off1, mm->off2) * mm->blocksize;
    mm->dense->q = malloc(mm->dense->n * dblock_stride * sizeof(uint32_t));
    memset(mm->dense->q, 0, mm->dense->n * dblock_stride * sizeof(uint32_t));
    mm->sparse->q = malloc((mm->sparse->tlen + sparse_readahead) * sizeof(uint32_t));
    mm->sparse->len = malloc(mm->sparse->n * sizeof(uint32_t));

    for(unsigned int i = 0 ; i < mm->dense->n ; i++) {
        MATMUL_COMMON_READ_MANY32(mm->dense->q + i * dblock_stride, ncols_t, f);
    }
    MATMUL_COMMON_READ_MANY32(mm->sparse->len, mm->sparse->n, f);
    MATMUL_COMMON_READ_MANY32(mm->sparse->q, mm->sparse->tlen, f);

    /* fill up the readahead zone properly */
    memset(mm->sparse->q + mm->sparse->tlen, 0, sparse_readahead * sizeof(uint32_t));

    fclose(f);

    matmul_threaded_blocks_info(mm);
    mm->tg = worker_threads_init(mm->nthreads);
    pthread_cond_init(&mm->dense->recheck_please, NULL);

    return 1;
}

void matmul_threaded_save_cache(struct matmul_threaded_data_s * mm)
{
    FILE * f;

    f = matmul_common_save_cache_fopen(abbytes(mm->xab, 1), mm->public_, MM_MAGIC);

    MATMUL_COMMON_WRITE_ONE32(mm->dense->weight, f);
    MATMUL_COMMON_WRITE_ONE32(mm->sgrp_size, f);
    MATMUL_COMMON_WRITE_ONE32(mm->dense->n, f);
    MATMUL_COMMON_WRITE_ONE32(mm->sparse->n, f);
    MATMUL_COMMON_WRITE_ONE32(mm->sparse->tlen, f);

    unsigned int ncols_t = mm->public_->dim[!mm->public_->store_transposed];
    size_t dblock_stride = next_multiple_of_powerof2(ncols_t, mm->blocksize);

    for(unsigned int i = 0 ; i < mm->dense->n ; i++) {
        MATMUL_COMMON_WRITE_MANY32(mm->dense->q + i * dblock_stride, ncols_t, f);
    }
    MATMUL_COMMON_WRITE_MANY32(mm->sparse->len, mm->sparse->n, f);
    MATMUL_COMMON_WRITE_MANY32(mm->sparse->q, mm->sparse->tlen, f);
    fclose(f);
}

void matmul_threaded_mul_sub_dense(struct matmul_threaded_data_s * mm, abt * dst, abt const * src, int d, int throff)
{
    if (mm->dense->n == 0) {
        /* We have nothing to do. In this special case, the task of
         * clearing the destination area is handed over to the sparse
         * block case */
        return;
    }

    ASM_COMMENT("dense multiplication code");
    abobj_ptr x = mm->xab;
    unsigned int nthr = mm->nthreads;      /* e.g. 4 */
    unsigned int b = mm->blocksize;     /* e.g. 16 */
    unsigned int off3 = mm->off3;

    uint32_t * A = mm->dense->q + throff * b;
    uint32_t i0 = throff * b;

    abt * rx = abinitf(x, 2);
    abzero(x, rx, 1);

    abt * dst0 = dst;
    abt const * src0 = src;

    if (d == !mm->public_->store_transposed) {
        /* t == 1 is well suited for d == 0. In this case, dst is
         * progressively filled. The same holds for t == 0 and d == 1 */

        /* We expect to have very few dense blocks, so it's acceptable to
         * have temporary storage for all of them on the stack */
        abt * tdst = abinitf(x, MM_DGRP_SIZE * mm->dense->n);
        abt * rb = tdst;
        for(unsigned int blocknum = 0 ; blocknum < mm->dense->n ; blocknum++) {
            src = src0 + aboffset(x, i0);
            abzero(x, rb, MM_DGRP_SIZE);
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < mm->public_->dim[d] ; i += nthr * b) {
                for(uint32_t j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off3, 0);
                    __builtin_prefetch(src + aboffset(x, nthr * b * off3), 0);
                    uint32_t rA = *A++;
                    abcopy(x, rx + aboffset(x, 1), src, 1);
                    src += aboffset(x, 1);
                    for(unsigned int k = 0 ; k < MM_DGRP_SIZE ; k++) {
                        abadd(x, rb + aboffset(x, k), rx + aboffset(x, rA&1));
                        rA >>= 1;
                    }
                }
                A += (nthr-1) * b;
                src += aboffset(x, (nthr-1) * b);
            }
            rb += aboffset(x, MM_DGRP_SIZE);
            ASM_COMMENT("end of critical loop");
        }

        /* Now the painful section. Experimentally, it does not seem
         * expensive to adopt this crude approach. */
        pthread_mutex_lock(&mm->tg->mu); 
        if (mm->dense->processed++ == 0) {
            abcopy(x, dst, tdst, MM_DGRP_SIZE * mm->dense->n);
        } else {
            for(unsigned int k = 0 ; k < MM_DGRP_SIZE * mm->dense->n ; k++) {
                abadd(x, dst + k, tdst + k);
            }
        }
        pthread_mutex_unlock(&mm->tg->mu);

#if 0
        /* There could be a way to seemingly avoid requiring a mutex
         * here, although in the end it probably wouldn't change much.
         * However, in any case we have to zero out the dst area in a
         * locked manner beforehand, and this implies serialization as
         * well.
         *
         * For reference, plan B uses the atomic xor operation
         * __sync_xor_and_fetch, which is a gcc builtin (corresponds to
         * using the lock prefix on x86*).  We're cheating with the abase
         * interface here. A lot.  Because depending on the active abase,
         * we may not be allowed to use __sync_xor_and_fetch. Thus we
         * have to cast everything to unsigned longs */

        unsigned int dsize_ulongs;
        dsize_ulongs = aboffset(x, MM_DGRP_SIZE);
        dsize_ulongs *= mm->dense->n;
        dsize_ulongs *= sizeof(abt);
        ASSERT(dsize_ulongs % sizeof(unsigned long) == 0);
        dsize_ulongs /= sizeof(unsigned long);

        for(unsigned int k = 0 ; k < dsize_ulongs ; k++) {
            __sync_xor_and_fetch(
                    ((unsigned long*)dst) + k,
                    ((unsigned long*)tdst) [k]);
        }
#endif
        abclearf(x, tdst, MM_DGRP_SIZE * mm->dense->n);
    } else {
        for(uint32_t i = i0 ; i < mm->public_->dim[!d] ; i += nthr * b) {
            abzero(x, dst0 + aboffset(x, i), b);
        }
        for(unsigned int blocknum = 0 ; blocknum < mm->dense->n ; blocknum++) {
            dst = dst0 + aboffset(x, i0);
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < mm->public_->dim[!d] ; i += nthr * b) {
                for(uint32_t j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + nthr * b * off3, 0);
                    __builtin_prefetch(dst + aboffset(x, nthr * b * off3), 1);
                    uint32_t rA = *A++;
                    abzero(x, rx + aboffset(x, 1), 1);
                    for(unsigned int k = 0 ; k < MM_DGRP_SIZE ; k++) {
                        abadd(x, rx + aboffset(x, rA&1), src + aboffset(x, k));
                        rA >>= 1;
                    }
                    abadd(x, dst, rx + aboffset(x, 1));
                    dst += aboffset(x, 1);
                }
                A += (nthr-1) * b;
                dst += aboffset(x, (nthr-1) * b);
            }
            src += aboffset(x, MM_DGRP_SIZE);
            ASM_COMMENT("end of critical loop");
        }
        pthread_mutex_lock(&mm->tg->mu); 
        /* It's a pity, but if some other thread is already doing the
         * sparse blocks at this point (it can't be anything else than
         * thread 0), we must be sure that everybody here is done with
         * the dense blocks.
         */
        if (++mm->dense->processed == mm->nthreads) {
            pthread_cond_signal(&mm->dense->recheck_please);
        }
        pthread_mutex_unlock(&mm->tg->mu);
    }
    abclearf(x, rx, 2);
    ASM_COMMENT("end of dense multiplication code");
}

void matmul_threaded_mul_sub_sparse(struct matmul_threaded_data_s * mm, abt * dst, abt const * src, int d, int throff)
{
    ASM_COMMENT("sparse multiplication code");
    abobj_ptr x = mm->xab;
    unsigned int nthr MAYBE_UNUSED = mm->nthreads;
    unsigned int g = mm->sgrp_size;
    unsigned int b = mm->blocksize;
    unsigned int off1 = mm->off1;
    unsigned int off2 = mm->off2;

    uint32_t * A = mm->sparse->q + throff * b;
    uint32_t i0 = throff * b;


    /* If d == 0, output has ncols == dim[1] entries
     *             input has nrows == dim[0] entries
     */

    /*
    abt * dst0 = dst;
    abt const * src0 = src;
    unsigned int ra[2];
    ra[0] = mm->public_->dim[0];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[0]));
    ra[1] = mm->public_->dim[1];
    matmul_threaded_aux(mm, MATMUL_AUX_GET_READAHEAD, &(ra[1]));
    */

    size_t sparse_readahead;
    sparse_readahead = mm->nthreads * MAX(mm->off1, mm->off2) * mm->blocksize;

    if (d == !mm->public_->store_transposed) {
        dst += aboffset(x, mm->dense->n * MM_DGRP_SIZE);
#define xxxHARDCODE
#ifdef HARDCODE
        ASSERT_ALWAYS(b == MM_BLOCKSIZE);
        ASSERT_ALWAYS(g == MM_SGRP_SIZE);
        ASSERT_ALWAYS(nthr == MM_NTHREADS);
        ASSERT_ALWAYS(off1 == MM_PREFETCH_OFFSET_1);
        ASSERT_ALWAYS(off2 == MM_PREFETCH_OFFSET_2);
#define b       MM_BLOCKSIZE
#define g       MM_SGRP_SIZE
#define nthr    MM_NTHREADS
#define off1    MM_PREFETCH_OFFSET_1
#define off2    MM_PREFETCH_OFFSET_2
#endif
        for(unsigned int blocknum = 0 ; blocknum < mm->sparse->n ; blocknum++) {
            // number of blocks
            uint32_t rlen= mm->sparse->len[blocknum];
            ASM_COMMENT("critical loop");
            for(uint32_t i = i0 ; i < g ; i += nthr * b) {
                abzero(x, dst + aboffset(x, i), b);
            }
            for(uint32_t i = i0 ; i < rlen ; i+= nthr * b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    // ASSERT((A - mm->sparse->q + nthr * b * off1) < mm->sparse->tlen + sparse_readahead);
                    // ASSERT((A - mm->sparse->q + nthr * b * off2) < mm->sparse->tlen + sparse_readahead);
                    __builtin_prefetch(A + nthr * b * off1, 0);
                    __builtin_prefetch(src + aboffset(x, A[nthr * b * off2]), 0);
                    uint32_t z = (i + j) & (g-1);
                    uint32_t k = *A++;
                    /*
                    ASSERT((dst - dst0) < aboffset(x, ra[!d]-z));
                    ASSERT((src - src0) < aboffset(x, 1 + mm->public_->dim[d]-k));
                    */
                    abadd(x, dst + aboffset(x, z), src + aboffset(x, k));
                }
                A += (nthr-1) * b;
            }
            ASM_COMMENT("end of critical loop");

            dst += aboffset(x, g);
        }
#ifdef  HARDCODE
#undef b
#undef g
#undef nthr
#undef off1
#undef off2
#endif
    } else {
        /* Then it's not so fun. In fact it's even terribly slow, and
         * bug-prone (concurrent writes all over the place). Therefore we
         * prefer to stay single-threaded here.
         *
         * Note that there _could_ be a way around the problem, using
         * atomic xors for abadd. But it isn't worth the trouble.
         */
        if (throff != 0)
            return;

        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        src += aboffset(x, mm->dense->n * MM_DGRP_SIZE);

        if (mm->dense->n == 0) {
            abzero(x, dst, mm->public_->dim[!d]);
        } else {
            /* We must test whether all threads are done with the dense
             * blocks */
            if (mm->dense->processed != mm->nthreads) {
                pthread_mutex_lock(&mm->tg->mu);
                for( ; mm->dense->processed != mm->nthreads ; ) {
                    pthread_cond_wait(&mm->dense->recheck_please, &mm->tg->mu);
                }
                pthread_mutex_unlock(&mm->tg->mu);
            }
        }
        for(unsigned int blocknum = 0 ; blocknum < mm->sparse->n ; blocknum++) {
            ASM_COMMENT("critical loop (transposed mult, SLOW)");
            uint32_t rlen= mm->sparse->len[blocknum];
            for(uint32_t i = i0 ; i < rlen ; i+= b) {
                for(unsigned int j = 0 ; j < b ; j++) {
                    __builtin_prefetch(A + b * off1, 0);
                    __builtin_prefetch(dst + aboffset(x, A[b * off2]), 1);
                    uint32_t z = (i + j) % g;
                    uint32_t k = *A++;
                    /*
                    ASSERT((dst - dst0) < aboffset(x, 1 + mm->public_->dim[!d]-k));
                    ASSERT((src - src0) < aboffset(x, ra[d]-z));
                    */
                    abadd(x, dst + aboffset(x, k), src + aboffset(x, z));
                }
            }
            ASM_COMMENT("end of critical loop (transposed mult, SLOW)");
            src += aboffset(x, g);
        }
        /* The extra readahead location which is used for non-transposed
         * mult is here an extra scrap location. It's innocuous to write
         * there of course, but in the context of checking, it's better
         * if we keep a zero there because it may become readahead again
         * for the next turn.
         */
        abzero(x, dst + aboffset(x, mm->public_->dim[!d]), 1);
    }
    ASM_COMMENT("end of sparse multiplication code");
}

void matmul_thread_worker_routine(struct worker_threads_group * tg MAYBE_UNUSED, int i, struct matmul_threaded_data_s * mm)
{
    /* Do our job with mutexes released */
    matmul_threaded_mul_sub_dense(mm, mm->dst, mm->src, mm->d, i);
    matmul_threaded_mul_sub_sparse(mm, mm->dst, mm->src, mm->d, i);
}

void matmul_threaded_mul(struct matmul_threaded_data_s * mm, abt * dst, abt const * src, int d)
{
    mm->dense->processed = 0;
    mm->dst = dst;
    mm->src = src;
    mm->d = d;

    /* unleash worker threads */
    worker_threads_do(mm->tg, (worker_func_t) &matmul_thread_worker_routine, mm);

    mm->public_->iteration[d]++;
}

void matmul_threaded_report(struct matmul_threaded_data_s * mm MAYBE_UNUSED, double scale MAYBE_UNUSED) {
}

void matmul_threaded_auxv(struct matmul_threaded_data_s * mm, int op, va_list ap)
{
    if (op == MATMUL_AUX_GET_READAHEAD) {
        unsigned int * res = va_arg(ap, unsigned int *);
        // First, we allow at least one extra.
        // But more subtly, we need to make sure that irrespective of the
        // number of dense and sparse blocks, the last sparse block will
        // have plenty of space to end. For example if nrows=40 and there
        // is one dense block of 32 rows, with sgrp_size=64, we need to
        // reach 96 entries, not merely 64. So if s and d are the sparse
        // and dense group sizes, we need to reach at least (x&~(d-1)+s)
        unsigned int start = *res;
        *res = (start & ~(MM_DGRP_SIZE-1)) + mm->sgrp_size;
        if (*res > start) {
            return;
        } else {
            *res = next_multiple_of_powerof2(1+*res, mm->sgrp_size);
        }
    }
}

void matmul_threaded_aux(struct matmul_threaded_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_threaded_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
