/* Manage the in-memory data for the matrix */

/* It's in C++ because the STL is handy, but that's really all there is
 * to it ; a conversion to C would not be extremely difficult */

#include "cado.h"
#include <cstdio>
#include <cstdlib>
#include <cerrno>

#include <climits>
#include <cmath>

#include <algorithm>

#include <stdint.h>
#include <inttypes.h>
#include <string.h>

// C++ headers.
// #include <string>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>    // sort
#include <iostream>     // cout
#include "utils.h"
#include "bwc_config.h"
using namespace std;

#include "matmul-bucket.h"
#include "abase.h"

/* Make sure that the assembly function is only called if it matches
 * correctly the abase header !! */
#if defined(HAVE_GCC_STYLE_AMD64_ASM) && defined(ABASE_U64_H_) && !defined(DISABLE_ASM)
// disabling one particular assembly code is done by simply disabling the
// header file (and optionally removing the assembly file from the link
// list is the CMakeLists.txt file, but in reality it's not needed. This
// file being C++, the names are mangled if no header file marks them as
// having C linkage. So commenting out is enough).
#include "matmul-sub-small1.h"
#include "matmul-sub-small2.h"
#include "matmul-sub-large-fbd.h"
#include "matmul-sub-large-fbi.h"
#include "matmul-sub-vsc-dispatch.h"
#include "matmul-sub-vsc-combine.h"
// #define ENABLE_ASM
#endif

#include "readmat-easy.h"
#include "matmul-common.h"

// take only a portion of the L1 cache.
#define L1_CACHE_SIZE   28000
// for 8-bytes abt values, this gives 3500 items.

// make sure that it's something each core can use (i.e. divide by two
// the L2 cache size for a dual-core w/ shared L2).
#define L2_CACHE_SIZE   1600000

/* To what should we compare the average dj value for determining when we
 * should switch from dense (small) to sparse (large) blocks */
#define DJ_CUTOFF1   1.5

/* Same, but for stopping large blocks and go to the next type
 * (currently, deferred blocks)
 */
#define DJ_CUTOFF2   3.0

/* How many buckets in large slices, and in large slices which are
 * contained in small slices. */
#define LSL_NBUCKETS_MAX      256

/* How many large slices MINIMUM must go in a huge slice ? Below this
 * value, we keep going producing large slices. */
#define HUGE_MPLEX_MIN        2

/* Here, I've seen 97 slices be handled best with 2 mplexed rounds. So
 * let's fix 64. */
#define HUGE_MPLEX_MAX        64

#define xxxDEBUG_BUCKETS

/* This implementation builds upon the ``sliced'' variant. 
 * The format used with the -sliced implementation is still used, but
 * there is also another format used for very dense strips. We refer to
 * matmul-sliced.cpp for the first format, which is called ``small1''
 * here (and there also).
 *
 * The ``small2'' format is specialised for denser strips.  Within an
 * horizontal strip of less than 4096 rows, when two non-zero columns are
 * never more than 16 positions apart, then the information (dj,i) is
 * stored in 16 bits instead of 32. This reduces the memory footprint
 * quite considerably, since many coefficients are in dense strips.
 * Speed-wise, small2 strips provide a nice speed-up sometimes (not
 * always).
 *
 * When we reach sparser areas, we use another method. Dense (small2 and
 * small1) slices having been discarded, the rest is split in ``large''
 * slices of LSL_NBUCKETS_MAX*256 rows at most (in fact we arrange for the
 * slices to have equal size, so it's a bit less).  For each such large
 * slice, we use scratch space equal to the L2 cache size for performing
 * matrix multiplication. While reading column vector data, we copy
 * coefficient to (up to) LSL_NBUCKETS_MAX ``buckets'', each representing
 * 256 possible row indices. The dispatching information, telling in
 * which bucket a given coefficient has to be copied, is an 8-bit integer
 * stored in a table called ``main''. Each entry in table ``main'' also
 * has a second 8-bit integer indicating the offset to the next useful
 * column entry (because we're considering very sparse blocks, it is not
 * at all uncommon to see average values ~ 10 here, even when
 * LSL_NBUCKETS_MAX==256, so that 65536 rows are considered
 * simultaneously).
 *
 * The number of columns considered while filling the scratch space is
 * limited by the scratch space size. Once it's filled, we ``apply the
 * buckets'', so as to be able to reuse our scratch space buffer for
 * another batch of columns. Thus a ``large slice'' is typically split
 * into several ``vertical blocks''. Extremely sparse large slices will
 * have only one vertical block.
 *
 * ``Applying the buckets'' means taking a list of column coefficients
 * which have been directed to the bucket because it is known that they
 * affect at least one of the 256 corresponding row indices on output. So
 * in addition to the N coefficients stored in the bucket, we have a
 * list of N 8-bit integers indicating which row index is to be modified.
 *
 * [the section below is superseded by the next one]
 * That's it for the basic idea. Now on top of that, we use an extra
 * indirection level for very sparse blocks. Several ``large blocks'' are
 * considered together as part of a ``huge block''. Call this number
 * ``nlarge''. Source coefficients are first dispatched in nlarge ``big
 * buckets'' of appropriate size. The number of columns considered is
 * guided by the fact that we want none of these big buckets to exceed
 * the L2 cache size.  Afterwards, we do a second bispatching from each
 * of these big buckets (in turn) to (up to) LSL_NBUCKETS_MAX small
 * buckets in the same manner as above (except that we no longer need the
 * information on the link to the next useful column.
 *
 * [documentation of the new format]
 * After ``large'' blocks, we consider ``deferred'' blocks, also referred
 * to as (vertical) staircase-shaped blocks. The idea is the following.
 * The portion of the matrix under study is split into several vertical
 * blocks of at most 2^16 entries. In the other direction, we make splits
 * according to the average row density. A first batch of row has some
 * indicative density X, the next one X/2, and so on (the ratio 2 being a
 * tunable constant) until the sparsest block at the bottom.  Such
 * horizontal blocks need not have the same size. 
 *
 * Vertical strips are processed in order. Buffers corresponding to the
 * sub-matrices in the different sub-matrices created by the horizontal
 * splitting are filled. We thus populate scratch space with coefficients
 * from the source vector. After each such processing of a vertical
 * strip, we consider whether we have to perform flushing w.r.t. the
 * horizontal strips. The buffers corresponding to the densest horizontal
 * strip are flushed after each vstrip processing. Those which are twice
 * sparser are flushed twice less often. Thus the next vstrip occupies
 * some more buffer shape. This pattern goes until the sparsest block,
 * which is probably flushed only once. The total amount of buffer space
 * needed is proportional to the number of coefficients in an area having
 * a staircase shape. Assume we have N vstrips of width W. Assume
 * horizontal blocks have a number of rows R0, R1, R2, ... corresponding
 * to densities X, X/2, X/4, ... The formula
 * for the tbuf space is:
 *   max weight of the N matrices of size R0*W in h. strip #0
 * + max weight of the N/2 matrices of size R1*2W in h. strip #1
 * + max weight of the N/4 matrices of size R2*4W in h. strip #2
 * ...
 * It should be noted that the numbers added up in this formula are
 * expected to be of the same magnitude asymptotically. For real-life
 * examples though, it's pretty common to observe somewhat degenerate
 * cases.
 */

/* This extension is used to distinguish between several possible
 * implementations of the matrix-times-vector product.
 */
#define MM_EXTENSION   "-bucket"

/* MM_MAGIC is stored in files, so as to check compatibility. The upper
 * word (MM_MAGIC_FAMILY) correspond to the implementation itself, the
 * lower one (MM_MAGIC_VERSION) to the n-th binary incompatible change
 * (make sure to bump it) */
#define MM_MAGIC_FAMILY        0xa003UL
#define MM_MAGIC_VERSION       0x1015UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

template<typename T> inline T * ptrbegin(vector<T>& v) { return &v.front(); }
template<typename T> inline T const * ptrbegin(vector<T> const & v) { return &v.front(); }
template<typename T> inline T * ptrend(vector<T>& v) { return v.size() + &v.front(); }
template<typename T> inline T const * ptrend(vector<T> const & v) { return v.size() + &v.front(); }

#if 0
static unsigned int idiotic_sum(void * p, unsigned int nbytes)
{
    uint8_t * q = (uint8_t *) p;
    unsigned int res = 0x12345678;
    for( ; nbytes-- ; q++) {
        unsigned int z = *q;
        z ^= z << 24;
        res = (((res << 8) ^ (res * z)) - (z ^ 0xabababab)) ^ res >> 13;
    }
    return res;
}
#endif

#define SLICE_TYPE_NONE           0
#define SLICE_TYPE_SMALL1 1
#define SLICE_TYPE_SMALL1_VBLOCK         2
#define SLICE_TYPE_SMALL2         3

#define SLICE_TYPE_LARGE_ENVELOPE 4

/* When coefficients go through several passes, they are considered as
 * one slice for each pass. Note that the ordering of some slices is
 * mandated by the design, e.g. slices of type LARGE_ASB must follow
 * immediately their parent LARGE_FBI
 */

#if 0
#define SLICE_TYPE_LARGE_FBI      5
#define SLICE_TYPE_LARGE_ASB      6
#endif

#define SLICE_TYPE_HUGE_ENVELOPE  7

#if 0
#define SLICE_TYPE_HUGE_FBI       8
#define SLICE_TYPE_HUGE_FBD       9
#define SLICE_TYPE_HUGE_ASB       10
#endif

#define SLICE_TYPE_DEFER_ENVELOPE       11
#define SLICE_TYPE_DEFER_COLUMN         12
#define SLICE_TYPE_DEFER_DIS            13
#define SLICE_TYPE_DEFER_ROW            14
#define SLICE_TYPE_DEFER_CMB        15

#define SLICE_TYPE_MAX  16


static inline const char* slice_name(int s) { 
    switch(s) {
        case SLICE_TYPE_SMALL1_VBLOCK: return "small1v";
        case SLICE_TYPE_SMALL1: return "small1";
        case SLICE_TYPE_SMALL2: return "small2";
        case SLICE_TYPE_LARGE_ENVELOPE: return "large";
        case SLICE_TYPE_HUGE_ENVELOPE: return "huge";
        case SLICE_TYPE_DEFER_ENVELOPE: return "defer";
        case SLICE_TYPE_DEFER_COLUMN: return "d-col";
        case SLICE_TYPE_DEFER_DIS: return "d-dis";
        case SLICE_TYPE_DEFER_ROW: return "d-row";
        case SLICE_TYPE_DEFER_CMB: return "d-cmb";
        default: return "other";
    }
}

struct slice_header_t {
    uint16_t t;
    uint8_t nchildren;  // w.r.t the tree-like shape of the block structure.
    uint8_t pad[5];     // make the structure 256-bit wide.
    uint32_t i0, i1;
    uint32_t j0, j1;
    uint64_t ncoeffs;
};

/* Ok, there's only one field for the moment. But there might be more
 * eventually */
struct slice_runtime_stats {
    double t;
    slice_runtime_stats() : t(0) {}
};

struct matmul_bucket_methods {
    int small1;
    int small2;
    int large;
    int huge;
    int vsc;
    inline matmul_bucket_methods(const char * desc = NULL) {
        small1=small2=large=vsc=huge=0;
        if (!desc) {
            /* default configuration */
            small1=small2=large=vsc=1;
            return;
        } 
        const char * p = desc;
        for( ; ; ) {
            const char * q = strchr(p, ',');
            int n = q ? q-p : strlen(p);
            if (strncmp(p, "small1", n) == 0) {
                ASSERT_ALWAYS(!small1);
                small1=1;
            } else if (strncmp(p, "small2", n) == 0) {
                ASSERT_ALWAYS(!small2);
                small2=1;
            } else if (strncmp(p, "large", n) == 0) {
                ASSERT_ALWAYS(!large);
                large=1;
            } else if (strncmp(p, "huge", n) == 0) {
                ASSERT_ALWAYS(!huge && !vsc);
                huge=1;
            } else if (strncmp(p, "vsc", n) == 0) {
                ASSERT_ALWAYS(!huge && !vsc);
                vsc=1;
            } else {
                fprintf(stderr, "Parse error: %s\n", p);
            }
            if (!q) break;
            p = q + 1;
        }
        ASSERT_ALWAYS(small1||!small2);
        ASSERT_ALWAYS(small1||small2||large||vsc||huge);
    }
};

struct matmul_bucket_data_s {
    /* repeat the fields from the public_ interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    abobj_t xab;
    size_t npack;
    size_t scratch1size;
    size_t scratch2size;
    size_t scratch3size;
    abt * scratch1;
    abt * scratch2;
    abt * scratch3;
    vector<uint16_t> t16;       /* For small (dense) slices */
    vector<uint8_t> t8;         /* For large (sparse) slices */
    vector<unsigned int> aux;   /* Various descriptors -- fairly small */
    /* headers are the first thing found in memory */
    vector<slice_header_t> headers;
    vector<slice_runtime_stats> slice_timings;
    matmul_bucket_methods methods;
};

void matmul_bucket_clear(struct matmul_bucket_data_s * mm)
{
    if (mm->scratch1) abclear(mm->xab, mm->scratch1, mm->scratch1size);
    if (mm->scratch2) abclear(mm->xab, mm->scratch2, mm->scratch2size);
    if (mm->scratch3) abclear(mm->xab, mm->scratch3, mm->scratch3size);
    matmul_common_clear(mm->public_);
    // delete properly calls the destructor for members as well.
    delete mm;
}

static void mm_finish_init(struct matmul_bucket_data_s * mm);

struct matmul_bucket_data_s * matmul_bucket_init(abobj_ptr xx MAYBE_UNUSED, param_list pl, int optimized_direction)
{
    struct matmul_bucket_data_s * mm;
    mm = new matmul_bucket_data_s;

    /* We first make sure that everything gets to zero ; except that it
     * may be more complicated than just calling memset, because we have
     * some vector() types in there...
     */
    memset(mm->public_, 0, sizeof(struct matmul_public_s));
    // memset(mm, 0, sizeof(struct matmul_bucket_data_s));
    mm->scratch1size = 0; mm->scratch1 = NULL;
    mm->scratch2size = 0; mm->scratch2 = NULL;
    mm->scratch3size = 0; mm->scratch3 = NULL;

    abobj_init_set(mm->xab, xx);

    unsigned int npack = L1_CACHE_SIZE;
    if (pl) param_list_parse_uint(pl, "l1_cache_size", &npack);
    npack /= abbytes(mm->xab,1);
    mm->npack = npack;

    unsigned int scratch1size = L2_CACHE_SIZE/2;
    if (pl) param_list_parse_uint(pl, "l2_cache_size", &scratch1size);
    scratch1size /= abbytes(mm->xab,1);
    mm->scratch1size = scratch1size;

    unsigned int scratch2size;
    scratch2size = scratch1size * HUGE_MPLEX_MAX; // (1 << 24)/abbytes(mm->xab,1));
    mm->scratch2size = scratch2size;

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

    const char * tmp = NULL;

    if (pl)
        tmp = param_list_lookup_string(pl, "matmul_bucket_methods");

    mm->methods = matmul_bucket_methods(tmp);

    return mm;
}

/* This moves an element at the tail of a list with no copy, transferring
 * the ownership to the container argument */
template<typename T>
void transfer(list<T> * ctr, T * elem)
{
    ctr->push_back(T());
    swap(ctr->back(), *elem);
}

struct builder {
    uint32_t * data[2];
    // rowhead is a pointer which is naturally excpected to *move* within
    // the data[0] array.
    uint32_t * rowhead;
    uint32_t nrows_t;
    uint32_t ncols_t;
    const char * rowname;
    const char * colname;
    /* Statistics on the prime parameter affecting performance of
     * subroutines. All are computed as (total sum, #samples) */
    struct matmul_bucket_data_s * mm;
};

void builder_init(builder * mb, struct matmul_bucket_data_s * mm, uint32_t * data)
{
    memset(mb, 0, sizeof(struct builder));
    if (!data) data = matmul_common_read_stupid_data(mm->public_);
    mb->data[0] = data;
    mb->data[1] = NULL;

    mb->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    mb->ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    mb->rowname = rowcol[ mm->public_->store_transposed];
    mb->colname = rowcol[!mm->public_->store_transposed];
    mb->rowhead = mb->data[0];
    mb->mm = mm;
}

void builder_clear(builder * mb)
{
    free(mb->data[0]);
    memset(mb, 0, sizeof(struct builder));
}

/************************************/
/* Elementary (!) buliding routines */

/* there's room for a slice type where data is stored row-major. It would
 * make it possible to avoid most stores, for a nice performance gain. We
 * might even consider accumulating four rows at a time to four different
 * registers, so that the loads are compensated. The in-memory format
 * must account for that.
 */

#define SMALL_SLICES_I_BITS     12
#define SMALL_SLICES_DJ_BITS    (16 - SMALL_SLICES_I_BITS)

/* {{{ small slices */

struct small_slice_t {
    uint32_t i0;
    uint32_t i1;
    unsigned int ncoeffs;
    double dj_avg;
    uint32_t dj_max;
    int is_small2;

    typedef std::vector<std::pair<uint16_t, uint16_t> > Lv_t;
    typedef Lv_t::const_iterator Lvci_t;
    typedef Lv_t::iterator Lvi_t;

    typedef std::deque<Lv_t> Lu_t;
    typedef Lu_t::const_iterator Luci_t;
    typedef Lu_t::iterator Lui_t;

    Lu_t Lu;
};

int builder_do_small_slice(builder * mb, struct small_slice_t * S, uint32_t i0, uint32_t i1)
{
    S->i0 = i0;
    S->i1 = i1;
    S->ncoeffs = 0;
    ASSERT_ALWAYS(i1-i0 <= (1 << SMALL_SLICES_I_BITS) );

    /* Make enough vstrips */
    S->Lu.assign(iceildiv(mb->ncols_t, 1UL << 16), small_slice_t::Lv_t());
    uint32_t * ptr0 = mb->rowhead;
    /* We're doing a new slice */
    for(uint32_t i = i0 ; i < i1 ; i++) {
        for(unsigned int j = 0 ; j < *mb->rowhead ; j++) {
            uint32_t jj = mb->rowhead[1+j];
            S->Lu[jj>>16].push_back(std::make_pair(jj%(1<<16),i-i0));
            S->ncoeffs++;
        }
        mb->rowhead += 1 + *mb->rowhead;
    }
    /* L is the list of (column index, row index) of all
     * coefficients in the current horizontal slice */

    /* Convert all j indices to differences */
    S->dj_max = 0;
    S->dj_avg = mb->ncols_t / (double) S->ncoeffs;

    typedef small_slice_t::Lui_t Lui_t;
    typedef small_slice_t::Lvci_t Lvci_t;
    uint32_t lu_offset = 0;
    uint32_t j = 0;
    for(Lui_t lup = S->Lu.begin() ; lup != S->Lu.end() ; ++lup, lu_offset+=1<<16) {
        std::sort(lup->begin(), lup->end());
        for(Lvci_t lvp = lup->begin() ; lvp != lup->end() ; ++lvp) {
            uint32_t dj = (lu_offset + lvp->first) - j;
            j = lu_offset + lvp->first;
            if (dj > S->dj_max) { S->dj_max = dj; }
        }
    }

    /* There are two possible reasons for this slice to be discarded.
     * Either the max dj is too large -- larger than 16 bits -- or the
     * average dj is too large -- larger than some guessed value.
     */
    int keep1 = S->dj_max < (1UL << 16) && S->dj_avg < DJ_CUTOFF1;
    int keep2 = S->dj_max < (1 << SMALL_SLICES_DJ_BITS) && S->dj_avg < DJ_CUTOFF1;

    if (!mb->mm->methods.small2) keep2=0;

    S->is_small2 = keep2;

    if (!keep1 && !keep2) {
        mb->rowhead = ptr0;
    }
    return keep1 || keep2;
}

/* }}} */

/* {{{ large slices */

struct large_slice_vblock_t {
    uint32_t j0;
    uint32_t j1;
    unsigned int n;             // number of coeffs.
    unsigned int pad;           // (half the) number of padding coeffs.
    vector<uint8_t> t8c;        // read: contribution to t8.
    vector<unsigned int> auxc;  // read: contribution to aux.
};

struct large_slice_t {
    slice_header_t hdr[1];
    double dj_avg;
    uint32_t dj_max;
    list<large_slice_vblock_t> vbl;
};

struct large_slice_raw_t {
    vector<uint8_t> ind[LSL_NBUCKETS_MAX];
    vector<uint8_t> main;
    vector<uint32_t> col_sizes;
    vector<uint32_t> pad_sizes;
};

void split_large_slice_in_vblocks(builder * mb, large_slice_t * L, large_slice_raw_t * R, unsigned int scratch1size)
{
    /* Now split into vslices */
    uint8_t * mp = ptrbegin(R->main);
    uint8_t * ip[LSL_NBUCKETS_MAX];
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        ip[k] = ptrbegin(R->ind[k]);
    }
    unsigned int vblocknum = 0;
    double vbl_ncols_variance = 0;
    for(uint32_t j = 0 ; j < mb->ncols_t ; vblocknum++) {
        uint32_t j0 = j;
        size_t n = 0;
        size_t np = 0;
        for( ; j < mb->ncols_t ; j++) {
            size_t dn = R->col_sizes[j];
            size_t dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > scratch1size) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t j1 = j;

        large_slice_vblock_t V;
        V.j0 = j0;
        V.j1 = j1;
        /* First the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * main, and nf bytes in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.auxc.reserve(258);
        V.auxc.push_back(V.j1 - V.j0);
        V.t8c.resize(3 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        memcpy(q, mp, 2 * nf * sizeof(uint8_t));
        if (vblocknum) {
            q[0] = 0;
        }
        unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ind_sizes[mp[2 * k + 1]]++;
        }
        mp += 2 * nf;
        V.auxc.push_back(nf);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ASSERT(ip[k]+ind_sizes[k] <= ptrend(R->ind[k]));
            memcpy(q, ip[k], ind_sizes[k] * sizeof(uint8_t));
            q += ind_sizes[k];
            ip[k] += ind_sizes[k];
            V.auxc.push_back(ind_sizes[k]);

            // mb->asb_avg[0] += ind_sizes[k];
            // mb->asb_avg[1] += MIN(L->hdr->i1 - L->hdr->i0 - k * 256, 256);
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));
#if 0
        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");
#endif
        transfer(&(L->vbl), &V);
    }
    double vbl_ncols_mean = mb->ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    printf(" %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);

}

/* Do the partial transpose. Since we're handling sparse areas, it's
 * quite a bit faster than doing the full-blown transpose (even though
 * transposing may be needed anyway because the input matrix is in the
 * wrong order, it saves us hassles to only view the matrix in one given
 * way.
 *
 * Note though that there _is_ some loss ; when testing, we no longer
 * have access to the pre-computed transposed file, since we recompute
 * the transposed parts each time.
 */

static uint32_t * do_partial_transpose(builder * mb, vector<uint32_t> & cs, uint32_t i0, uint32_t i1) /*{{{*/
{
    uint32_t * cols;
    uint32_t * qptr;

    uint32_t * ptr = mb->rowhead;
    for(uint32_t i = i0 ; i < i1 ; i++) {
        uint32_t w = *ptr++;
        for( ; w-- ; ) {
            uint32_t j = *ptr++;
            cs[j]++;
        }
    }
    cols = new uint32_t[ptr - mb->rowhead - (i1 - i0)];
    qptr = cols;

    vector<uint32_t *> colptrs;
    for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
        colptrs.push_back(qptr);
        qptr += cs[j];
    }
    ASSERT(qptr-cols == ptr - mb->rowhead - (ptrdiff_t) (i1 - i0));
    ptr = mb->rowhead;
    for(uint32_t i = i0 ; i < i1 ; i++) {
        uint32_t w = *ptr++;
        for( ; w-- ; ) {
            uint32_t j = *ptr++;
            *colptrs[j]++ = i - i0;
        }
    }
    colptrs.clear();
    mb->rowhead = ptr;
    return cols;
}/*}}}*/

int builder_do_large_slice(builder * mb, struct large_slice_t * L, uint32_t i0, uint32_t i1, unsigned int scratch1size)
{
    memset(L->hdr, 0, sizeof(slice_header_t));
    L->hdr->t = SLICE_TYPE_LARGE_ENVELOPE;
    L->hdr->j0 = 0;
    L->hdr->j1 = mb->ncols_t;
    L->hdr->i0 = i0;
    L->hdr->i1 = i1;
    L->hdr->ncoeffs = 0;

    L->dj_max = 0;

    large_slice_raw_t R[1];

    R->col_sizes.assign(mb->ncols_t, 0);
    R->pad_sizes.assign(mb->ncols_t, 0);

    uint32_t * ptr0 = mb->rowhead;
    uint32_t * cols = do_partial_transpose(mb, R->col_sizes, i0, i1);
    uint32_t * qptr = cols;

    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
        uint32_t len  = R->col_sizes[j];

        if (!len) continue;

        uint32_t diff = j - last_j;
        if (diff > L->dj_max) {
            L->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->main.push_back(255);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            R->main.push_back(0);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->main.push_back(diff);
        last_j = j;
        for( ; len-- ; ) {
            uint32_t i = *qptr++;
            uint8_t w = i / 256;
            R->main.push_back(w);
            R->main.push_back(0);
            R->ind[w].push_back((uint8_t) i);
            L->hdr->ncoeffs ++;
        }
        R->main.pop_back();
    }
    ASSERT(qptr-cols == mb->rowhead - ptr0 - (ptrdiff_t) (i1 - i0));
    qptr = NULL;
    delete[] cols;
    cols = NULL;

    /* The average dj is quite easy */
    L->dj_avg = mb->ncols_t / (double) L->hdr->ncoeffs;

    printf(" w=%"PRIu64", avg dj=%.1f, max dj=%u, bucket hit=1/%.1f",
            L->hdr->ncoeffs, L->dj_avg, L->dj_max, LSL_NBUCKETS_MAX * L->dj_avg);

    if (L->dj_avg > DJ_CUTOFF2) {
        printf("-> too sparse");
        if (mb->nrows_t - i0 < HUGE_MPLEX_MIN * LSL_NBUCKETS_MAX * 256) {
            printf("; kept because of short tail\n");
        } else {
            /* We won't keep this slice */
            printf("\n");
            mb->rowhead = ptr0;
            return 0;
        }
    }

    printf("\n");

    split_large_slice_in_vblocks(mb, L, R, scratch1size);

    L->hdr->nchildren = L->vbl.size();

    return 1;
}

/* }}} */

/* {{{ huge slices */

struct huge_slice_t {
    slice_header_t hdr[1];
    double dj_avg;
    uint32_t dj_max;
    unsigned int nlarge;
    list<large_slice_vblock_t> vbl;
};

struct huge_subslice_raw_t {
    vector<uint8_t> ind[LSL_NBUCKETS_MAX];
    vector<uint8_t> main;
};

struct huge_subslice_ptrblock_t {
    uint8_t * ip[LSL_NBUCKETS_MAX];
    uint8_t * mp;
};

struct huge_slice_raw_t {
    vector<uint8_t> super;
    vector<huge_subslice_raw_t> subs;
    vector<uint32_t> col_sizes;
    vector<uint32_t> pad_sizes;
};

void split_huge_slice_in_vblocks(builder * mb, huge_slice_t * H, huge_slice_raw_t * R, unsigned int scratch2size)/*{{{*/
{
    /* Now split into vslices */
    // unsigned int lsize = iceildiv(H->hdr->i1 - H->hdr->i0, H->nlarge);
    uint8_t * sp = ptrbegin(R->super);
    vector<huge_subslice_ptrblock_t> ptrs(R->subs.size());
    for(unsigned int i = 0 ; i < R->subs.size() ; i++) {
        ptrs[i].mp = ptrbegin(R->subs[i].main);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ptrs[i].ip[k] = ptrbegin(R->subs[i].ind[k]);
        }
    }
    double vbl_ncols_variance = 0;
    unsigned int vblocknum = 0;
    for(uint32_t j = 0 ; j < mb->ncols_t ; vblocknum++) {
        uint32_t j0 = j;
#if 0
        /* First way to advance to the next j -- compare the total number
         * of coeffs to the scratch2size argument */
        size_t n = 0;
        size_t np = 0;
        for( ; j < mb->ncols_t ; j++) {
            size_t dn = R->col_sizes[j];
            size_t dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > scratch2size) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t j1 = j;
#else
        /* However, there's another way. We might as well arrange so that
         * no mega bucket exceeds the L2 cache. */
        uint32_t j1 = j;
        uint8_t * sp0 = sp;
        uint8_t * spc = sp0;
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        if (vblocknum) {
            sp0[0]=0;
        }
        {
            vector<unsigned int> Lsizes(H->nlarge, 0);
            for( ; j < mb->ncols_t && sp != ptrend(R->super) ; ) {
                unsigned int dj = *sp;
                j += dj;
                if (dj) {
                    /* beginning of a column. */
                    spc = sp;
                    j1 = j;
                }
                sp++;
                unsigned int k = *sp++;
                Lsizes[k]++;
                if (Lsizes[k] > scratch2size / HUGE_MPLEX_MAX) {
                    Lsizes[k]--;
                    sp = spc;
                    break;
                }
            }
            if (sp == ptrend(R->super)) {
                spc = sp;
                j1 = mb->ncols_t;
            }
        }
        /* TODO: abandon the n/np distinction. It's useful for debugging
         * at best, and it would not hurt to have it displayed only at the
         * upper level.
         */
        size_t n = 0;
        size_t np = 0;
        for(j = j0 ; j < j1 ; j++) {
            n += R->col_sizes[j];
            np += R->pad_sizes[j];
        }
        ASSERT_ALWAYS((n+2*np)*2 == (size_t) (spc - sp0));
        sp = sp0;
#endif

        /* That's the same type, although it contains more stuff. */
        large_slice_vblock_t V;
        V.j0 = j0;
        V.j1 = j1;
        V.auxc.reserve(1 + H->nlarge * (1 + LSL_NBUCKETS_MAX));
        V.auxc.push_back(V.j1 - V.j0);
        /* First the super block, and then for each contained large
         * slice, the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * super, nf in mains, and nf in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.t8c.resize(4 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        memcpy(q, sp, 2 * nf * sizeof(uint8_t));
        /* TODO: If the solution above is to be kept, then of course we'd
         * rather reuse this thing above, since it has already been
         * computed.
         */
        vector<unsigned int> L_sizes(H->nlarge, 0);
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ASSERT(sp[2 * k + 1] < H->nlarge);
            L_sizes[sp[2 * k + 1]]++;
        }
        sp += 2 * nf;
        V.auxc.push_back(nf);
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            V.auxc.push_back(L_sizes[l]);
        }
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
            memcpy(q, ptrs[l].mp, L_sizes[l] * sizeof(uint8_t));
            q += L_sizes[l];
            for(unsigned int k = 0 ; k < L_sizes[l] ; k++) {
                ind_sizes[ptrs[l].mp[k]]++;
            }
            ptrs[l].mp += L_sizes[l];

            // unsigned int i_size = MIN(H->hdr->i1 - H->hdr->i0 - l * lsize, lsize);

            for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
                ASSERT(ptrs[l].ip[k]-(ptrbegin(R->subs[l].ind[k]))+(ptrdiff_t) ind_sizes[k] <= (ptrdiff_t) R->subs[l].ind[k].size());
                memcpy(q, ptrs[l].ip[k], ind_sizes[k] * sizeof(uint8_t));
                q += ind_sizes[k];
                ptrs[l].ip[k] += ind_sizes[k];
                V.auxc.push_back(ind_sizes[k]);
                // mb->asb_avg[0] += ind_sizes[k];
                // mb->asb_avg[1] += MIN(i_size - k * 256, 256);
            }
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));

        /* TODO: have some sort of verbosity level */
#if 0
        printf(" vbl%u", vblocknum);
        printf(": %ss %u..%u ", mb->colname, j0, j1);
        printf("; w %u", V.n);
        printf("; avg dj %.1f", (j1 - j0) / (double) V.n);
        if (V.pad) printf("; pad 6*%u", V.pad);
        printf("\n");
#endif

        transfer(&(H->vbl), &V);
    }
    double vbl_ncols_mean = mb->ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    printf(" %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);
}/*}}}*/

int builder_do_huge_slice(builder * mb, struct huge_slice_t * H, uint32_t i0, uint32_t i1, unsigned int scratch2size)
{
    memset(H->hdr, 0, sizeof(slice_header_t));
    H->hdr->t = SLICE_TYPE_HUGE_ENVELOPE;
    H->hdr->j0 = 0;
    H->hdr->j1 = mb->ncols_t;
    H->hdr->i0 = i0;
    H->hdr->i1 = i1;
    H->hdr->ncoeffs = 0;

    H->dj_max = 0;

    huge_slice_raw_t R[1];

    R->col_sizes.assign(mb->ncols_t, 0);
    R->pad_sizes.assign(mb->ncols_t, 0);

    uint32_t * ptr0 = mb->rowhead;
    uint32_t * cols = do_partial_transpose(mb, R->col_sizes, i0, i1);
    uint32_t * qptr = cols;

    /* How many large slices in this huge slice ? */

    H->nlarge = iceildiv(i1 - i0, LSL_NBUCKETS_MAX * 256);
    ASSERT(H->nlarge >= HUGE_MPLEX_MIN);
    ASSERT(H->nlarge <= HUGE_MPLEX_MAX);
    unsigned int lsize = iceildiv(i1 - i0, H->nlarge);
    ASSERT(lsize <= LSL_NBUCKETS_MAX * 256);
    ASSERT(lsize * H->nlarge >= i1 - i0);
    R->subs.assign(H->nlarge, huge_subslice_raw_t());

    printf(" (%u*%u) ", H->nlarge, lsize);
    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    uint32_t next_dot = mb->ncols_t / H->nlarge;
    for(uint32_t j = 0 ; j < mb->ncols_t ; j++) {
        uint32_t len  = R->col_sizes[j];

        if (j > next_dot) {
            printf(".");
            fflush(stdout);
            next_dot += mb->ncols_t / H->nlarge;
        }

        if (!len)
            continue;

        uint32_t diff = j - last_j;
        if (diff > H->dj_max) {
            H->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->super.push_back(255);    // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            R->super.push_back(0);      // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->super.push_back(diff);
        last_j = j;
        for( ; len-- ; ) {
            uint32_t i = *qptr++;
            uint8_t w = i / lsize;
            uint8_t wq = (i % lsize) / 256;
            uint8_t wr = i % lsize;
            R->super.push_back(w);
            R->super.push_back(0);
            R->subs[w].main.push_back(wq);
            R->subs[w].ind[wq].push_back(wr);
            H->hdr->ncoeffs ++;
        }
        R->super.pop_back();
    }
    printf("\n");
    ASSERT(qptr-cols == mb->rowhead - ptr0 - (ptrdiff_t) (i1 - i0));
    qptr = NULL;
    delete[] cols;
    cols = NULL;

    /* The average dj is quite easy */
    H->dj_avg = mb->ncols_t / (double) H->hdr->ncoeffs;

    printf(" w=%"PRIu64", avg dj=%.1f, max dj=%u, bucket block hit=1/%.1f\n",
            H->hdr->ncoeffs, H->dj_avg, H->dj_max,
            H->nlarge * H->dj_avg);

    if (0) {
        /* Don't keep -- well, here this never occurs */
        mb->rowhead = ptr0;
        return 0;
    }

    split_huge_slice_in_vblocks(mb, H, R, scratch2size);

    H->hdr->nchildren = H->vbl.size();

    return 1;
}

/* }}} */

/* {{{ vertical staircase slices */

struct vsc_step_t {
    unsigned int defer;
    uint32_t nrows;
    unsigned int density_upper_bound;
    unsigned long tbuf_space;
};

static inline unsigned int when_flush(unsigned int k, unsigned int nv, unsigned int defer)
{
    k -= k % defer;
    k += defer - 1;
    k = min(k, nv-1);
    return k;
}

static inline int flush_here(unsigned int k, unsigned int nv, unsigned int defer) {
    return k == when_flush(k,nv,defer);
}

struct vsc_sub_slice_t {
    slice_header_t hdr[1];
    vector<uint16_t> x;
    vector<uint8_t> c;
};

struct vsc_middle_slice_t {
    slice_header_t hdr[1];
    vector<vsc_sub_slice_t> sub;
};

struct vsc_slice_t {
    /* This header is mostly a placeholder, and is not here to reflect an
     * additional pass on the data. */
    slice_header_t hdr[1];
    vector<vsc_step_t> steps;
    vector<vsc_middle_slice_t> dispatch;
    vector<uint32_t> transpose_help;
    unsigned long tbuf_space;
};

/* read batches of ROW_BATCH rows, in order to have something which is less
 * sensible to variations incurred by splitting the matrix */
#define ROW_BATCH 256
#define FLUSHFREQ_RATIO 1.7

/* {{{ decide on our flushing periods (this step is data independent) */
static vector<unsigned int> flush_periods(unsigned int nvstrips)
{
    vector<unsigned int> fp;
    for(unsigned int nf = 1 ;  ; ) {
        unsigned int quo = min(255u, iceildiv(nvstrips, nf));
        if (fp.empty() || quo < fp.back()) {
            fp.push_back(quo);
        }
        if (quo == 1)
            break;
        // next ?
        nf = max(nf + 1, (unsigned int) (nf * FLUSHFREQ_RATIO));
    }
    std::reverse(fp.begin(), fp.end());
    printf("Considering cells to be flushed every");
    for(unsigned int s = 0 ; s < fp.size() ; s++) {
        printf(" %u", fp[s]);
    }
    printf(" vstrips\n");
    ASSERT_ALWAYS(fp.back() == min(255u, nvstrips));
    ASSERT_ALWAYS(fp.front() == 1);

    return fp;
}
/* }}} */

/* {{{ read the row weights */
static vector<unsigned long> rowblock_weights(builder * mb, struct vsc_slice_t * V)
{
    vector<unsigned long> blockweight;
    uint32_t * ptr = mb->rowhead;
    /* It's a bit unfortunate, but since here we do speedy parsing of the
     * row weights, we're taking a shortcut which is valid only if our
     * data set spans the entire column range */
    ASSERT_ALWAYS(V->hdr->j0 == 0 && V->hdr->j1 == mb->ncols_t);
    for(uint32_t i = V->hdr->i0 ; i < V->hdr->i1 ; ) {
        unsigned long bw = 0;
        for(uint32_t di = 0 ; di < ROW_BATCH && i < V->hdr->i1 ; i++, di++) {
            unsigned int w = *ptr++;
            ptr += w;
            bw += w;
        }
        blockweight.push_back(bw);
    }
    return blockweight;
}
/* }}} */

/*{{{*/
static void
compute_staircase(builder * mb, struct vsc_slice_t * V)
{
    unsigned int nvstrips = V->dispatch.size();

    vector<unsigned long> blockweight = rowblock_weights(mb, V);
    vector<unsigned int> flushperiod = flush_periods(nvstrips);

    /* {{{ first dispatching pass */
    uint32_t ii = V->hdr->i0;
    int s = 0;
    for(unsigned int k = 0 ; ; k++) {
        if (k == blockweight.size() ||
                (s < (int) flushperiod.size() - 1 &&
                 blockweight[k] < ROW_BATCH * nvstrips / flushperiod[s]))
        {
            uint32_t old_ii = ii;
            ii = min(V->hdr->i0 + k * ROW_BATCH, V->hdr->i1);
            uint32_t di = ii - old_ii;
            if (di) {
                vsc_step_t step;
                step.defer = flushperiod[s];
                step.nrows = di;
                step.tbuf_space = 0;    // to be set later.
                step.density_upper_bound = (s == 0) ?
                    UINT_MAX : iceildiv(nvstrips, flushperiod[s-1]);
                V->steps.push_back(step);
            }
            s++;
        }
        if (k == blockweight.size())
            break;
    }
    /* }}} */
    /* {{{ second pass: do some merges for smallish areas */
    ii = V->hdr->i0;
    for(unsigned int k = 0 ; k < V->steps.size() ; ) {
        unsigned long w = V->steps[k].nrows;
        uint32_t old_ii = ii;
        int merge_with_next = 0;
        int merge_with_previous = 0;
        if ((k+1) < V->steps.size())
            merge_with_next = w <= 8192; // || w <= V->steps[k+1].nrows / 10;
        if (k > 0)
            merge_with_previous = w <= 8192; // || w <= V->steps[k-1].nrows / 10;
        if (!merge_with_previous && !merge_with_next) {
            ii += w;
            k++;
            continue;
        }
        /* Otherwise it's a bit ridiculous to have a separate block for
         * such a small number of rows. So by default, we choose to merge
         * it with the next block, or the previous one if we're reaching
         * the end. */
        V->steps.erase(V->steps.begin() + k);
        if (merge_with_next) {
            printf("strip %"PRIu32"+%lu merged with next\n",
                    old_ii, w);
            ASSERT_ALWAYS(k < V->steps.size());
            V->steps[k].nrows += w;
            // don't increment k here.
            // don't increment ii either.
        } else {
            ASSERT_ALWAYS(merge_with_previous);
            printf("strip %"PRIu32"+%lu merged with previous\n",
                    old_ii, w);
            V->steps[k-1].nrows += w;
            // don't increment k here.
            ii += w;
        }
    }
    /* }}} */
}
/*}}}*/

/*{{{*/
static
void vsc_fill_buffers(builder * mb, struct vsc_slice_t * V)
{
    unsigned int nvstrips = V->dispatch.size();
    uint32_t width = iceildiv(V->hdr->j1 - V->hdr->j0, nvstrips);
    uint32_t * ptr = mb->rowhead;
    uint32_t i = V->hdr->i0;
    V->tbuf_space = 0;

    for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
        /* This strip is scheduled to be flushed every V->steps[s].defer
         * vstrips. */
        unsigned int defer = V->steps[s].defer;
        uint32_t old_ii = i;
        uint32_t ii = i + V->steps[s].nrows;
        ASSERT_ALWAYS(ii <= V->hdr->i1);
        for( ; i < ii ; i++) {
            unsigned int w = *ptr++;
            for( ; w-- ; ) {
                uint32_t j = *ptr++;
                unsigned int d = j / width;
                V->dispatch[d].sub[s].hdr->ncoeffs++;
                V->dispatch[d].sub[s].x.push_back(j % width);
                unsigned int fidx = when_flush(d,nvstrips,defer);
                V->dispatch[fidx].sub[s].c.push_back(1 + (d % defer));
            }
            for(unsigned int d = 0 ; d < nvstrips ; d+= defer) {
                /* row change */
                unsigned int fidx = when_flush(d,nvstrips,defer);
                V->dispatch[fidx].sub[s].c.push_back(0);
            }
        }
        /*
        printf("Post counts\n");
        for(unsigned int d = 0 ; d < nvstrips ; d++) {
            printf(" (nrows=%u) [%u].sub[%u] : %lu coeffs ; xsize=%zu, csize=%zu\n", V->steps[s].nrows, d, s, V->dispatch[d].sub[s].hdr->ncoeffs, V->dispatch[d].sub[s].x.size(), V->dispatch[d].sub[s].c.size());
        }
        */
        unsigned int acc = 0;
        for(unsigned int d = 0 ; d < nvstrips ; d++) {
            acc += V->dispatch[d].sub[s].hdr->ncoeffs;
            if (!flush_here(d,nvstrips,defer))
                continue;
            ASSERT(V->dispatch[d].sub[s].c.size() == acc + V->steps[s].nrows);
            acc = 0;
        }
        unsigned long m = 0;
        unsigned long cm = 0;
        for(unsigned int d = 0 ; d < nvstrips ; d++) {
            if (d % defer == 0) cm = 0;
            cm += V->dispatch[d].sub[s].hdr->ncoeffs;
            if (cm >= m) m = cm;
        }
        printf("Rows %"PRIu32"+%"PRIu32, old_ii, V->steps[s].nrows);
        if (V->steps[s].density_upper_bound != UINT_MAX) {
            printf(": d < %u", V->steps[s].density_upper_bound);
        }
        printf("; %u flushes (every %u), tbuf space %lu\n", iceildiv(nvstrips, defer), defer, m);
        m += 1; // add space for one dummy pointer.
        V->steps[s].tbuf_space = m;
        V->tbuf_space += m;
    }
    for(unsigned int d = 0 ; d < nvstrips ; d++) {
        for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
            V->dispatch[d].hdr->ncoeffs+=V->dispatch[d].sub[s].hdr->ncoeffs;
        }
        V->hdr->ncoeffs+=V->dispatch[d].hdr->ncoeffs;
    }

    printf("Total tbuf space %lu (%lu MB)\n",
            V->tbuf_space, (V->tbuf_space * sizeof(abt)) >> 20);
}
/*}}}*/

int builder_prepare_vsc_slices(builder * mb, struct vsc_slice_t * V, uint32_t i0)
{
    memset(V->hdr, 0, sizeof(slice_header_t));
    V->hdr->t = SLICE_TYPE_DEFER_ENVELOPE;
    V->hdr->i0 = i0;
    V->hdr->i1 = mb->nrows_t;
    V->hdr->j0 = 0;
    V->hdr->j1 = mb->ncols_t;
    V->hdr->ncoeffs = 0;

    unsigned int nvstrips = iceildiv(V->hdr->j1 - V->hdr->j0, 1 << 16);
    V->hdr->nchildren = nvstrips;

    uint32_t width = iceildiv(V->hdr->j1 - V->hdr->j0, nvstrips);
    printf("%u vstrips of width %"PRIu32"\n", nvstrips, width);

    // prepare the dispatch slice headers
    for(uint32_t j = V->hdr->j0 ; j < V->hdr->j1 ; j += width) {
        vsc_middle_slice_t v;
        memset(v.hdr, 0, sizeof(slice_header_t));
        v.hdr->t = SLICE_TYPE_DEFER_COLUMN;
        v.hdr->i0 = V->hdr->i0;
        v.hdr->i1 = V->hdr->i1;
        v.hdr->j0 = j;
        v.hdr->j1 = min(j + width, V->hdr->j1);
        /* ncoeffs is set later, as is the data array x */
        v.hdr->ncoeffs = 0;
        V->dispatch.push_back(v);
    }
    ASSERT(V->dispatch.size() == nvstrips);

    compute_staircase(mb, V);
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        V->dispatch[k].hdr->nchildren = V->steps.size();
        uint32_t i0 = V->hdr->i0;
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            vsc_sub_slice_t w;
            memset(w.hdr, 0, sizeof(slice_header_t));
            w.hdr->t = SLICE_TYPE_DEFER_DIS;
            w.hdr->i0 = i0;
            w.hdr->i1 = (i0 += V->steps[l].nrows);
            w.hdr->j0 = V->dispatch[k].hdr->j0;
            w.hdr->j1 = V->dispatch[k].hdr->j1;
            w.hdr->ncoeffs = 0;
            w.hdr->nchildren = 0;
            V->dispatch[k].sub.push_back(w);
        }
    }
    vsc_fill_buffers(mb, V);

    return 0;
}

/* }}} */

/*************************************************************************/
/* Pushing slices to mm ; all these routines clear the given slice stack */

/* {{{ small slices */
void builder_push_small_slice(struct matmul_bucket_data_s * mm, small_slice_t * S)
{
    unsigned int ncols_t = mm->public_->dim[!mm->public_->store_transposed];
    if (S->is_small2) {
        slice_header_t hdr[1];
        memset(hdr, 0, sizeof(slice_header_t));
        hdr->t = SLICE_TYPE_SMALL2;
        hdr->i0 = S->i0;
        hdr->i1 = S->i1;
        hdr->j0 = 0;
        hdr->j1 = ncols_t;
        hdr->ncoeffs = S->ncoeffs;
        hdr->nchildren = 0;

        /* small2 slices are smaller, and faster -- but more
         * restrictive.  */

        /* Because a small2 slice has data in 16bits chunks, we'll
         * have to ensure proper alignment in the end. */
        ASSERT_ALWAYS((mm->t16.size() & (2-1)) == 0);

        unsigned int j = 0;
        uint32_t lu_offset = 0;
        for( ; ! S->Lu.empty() ; S->Lu.pop_front(), lu_offset+=1<<16) {
            small_slice_t::Lv_t lv;
            swap(lv, S->Lu.front());
            typedef small_slice_t::Lvci_t Lvci_t;
            for(Lvci_t lvp = lv.begin() ; lvp != lv.end() ; ++lvp) {
                uint32_t dj = (lu_offset + lvp->first) - j;
                ASSERT(dj < (1 << SMALL_SLICES_DJ_BITS) );
                ASSERT(lvp->second < (1 << SMALL_SLICES_I_BITS) );
                mm->t16.push_back((dj << SMALL_SLICES_I_BITS) | lvp->second);
                j = lu_offset + lvp->first;
            }
        }

        // align.
        for( ; (mm->t16.size() & (2-1)) ; )
            mm->t16.push_back(0);

        mm->headers.push_back(*hdr);
    } else {
        /* How many vertical parts -- this merely has to do with index
         * wraparound, since we're processing data column-major anyway.
         */
        slice_header_t ehdr[1];
        vector<slice_header_t> hdrs;
        memset(ehdr, 0, sizeof(slice_header_t));
        ehdr->t = SLICE_TYPE_SMALL1;
        ehdr->i0 = S->i0;
        ehdr->i1 = S->i1;
        ehdr->j0 = 0;
        ehdr->j1 = ncols_t;
        ehdr->ncoeffs = S->ncoeffs;

        unsigned int lu_offset = 0;
        for( ; ! S->Lu.empty() ; S->Lu.pop_front()) {
            small_slice_t::Lv_t lv;
            swap(lv, S->Lu.front());
            typedef small_slice_t::Lvci_t Lvci_t;

            slice_header_t hdr[1];
            memset(hdr, 0, sizeof(slice_header_t));
            hdr->t = SLICE_TYPE_SMALL1_VBLOCK;
            hdr->i0 = S->i0;
            hdr->i1 = S->i1;
            hdr->j0 = lu_offset;
            hdr->nchildren = 0;
            hdr->ncoeffs = lv.size();
            for(Lvci_t lvp = lv.begin() ; lvp != lv.end() ; ++lvp) {
                mm->t16.push_back(lvp->first);
                mm->t16.push_back(lvp->second);
            }
            lu_offset += (1UL << 16);
            hdr->j1 = min(lu_offset, ncols_t);
            hdrs.push_back(*hdr);
        }
        ehdr->nchildren = hdrs.size();
        mm->headers.push_back(*ehdr);
        mm->headers.insert(mm->headers.end(), hdrs.begin(), hdrs.end());
    }
    mm->public_->ncoeffs += S->ncoeffs;
}

/* }}} */

/* {{{ large slices */
void builder_push_large_slice(struct matmul_bucket_data_s * mm, large_slice_t * L)
{
    mm->headers.push_back(*L->hdr);
    for( ; ! L->vbl.empty() ; L->vbl.pop_front()) {
        large_slice_vblock_t & V(L->vbl.front());
        mm->aux.insert(mm->aux.end(), V.auxc.begin(), V.auxc.end());
        unsigned t8_size =  V.t8c.size();
        mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
        // large slices, but not huge slices, may put an odd number
        // of coefficients in t8
        if (t8_size & 1) { mm->t8.push_back(0); }
    }
    mm->public_->ncoeffs += L->hdr->ncoeffs;
}

/* }}} */

/* {{{ huge slices */

void builder_push_huge_slice(struct matmul_bucket_data_s * mm, huge_slice_t * H)
{
    mm->headers.push_back(*H->hdr);
    mm->aux.push_back(H->nlarge);
    for( ; ! H->vbl.empty() ; H->vbl.pop_front()) {
        large_slice_vblock_t & V(H->vbl.front());
        mm->aux.insert(mm->aux.end(), V.auxc.begin(), V.auxc.end());
        unsigned t8_size =  V.t8c.size();
        ASSERT_ALWAYS((t8_size & 1) == 0);
        mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
    }
    mm->public_->ncoeffs += H->hdr->ncoeffs;
}
/* }}} */

/******************************************/
/* Iteratively call the building routines */

/* {{{ small slices */
void builder_do_all_small_slices(builder * mb, uint32_t * p_i0, unsigned int npack)
{
    /* npack is a guess for the expected size of small slices ; they are
     * arranged later to all have approximately equal size.
     */
    unsigned int nslices = iceildiv(mb->nrows_t, npack);
    unsigned int s;
    for(s = 0 ; s < nslices ; s++) {
        uint32_t i0 =  s    * mb->nrows_t / nslices;
        uint32_t i1 = (s+1) * mb->nrows_t / nslices;

        small_slice_t S[1];

        printf("Ssl%u: %ss %u+%u...", s, mb->rowname, i0, i1-i0);
        fflush(stdout);

        int keep = builder_do_small_slice(mb, S, i0, i1);

        printf(" w %" PRIu32 " ; avg dj %.1f ; max dj %u%s\n",
                S->ncoeffs, S->dj_avg, S->dj_max,
                S->is_small2 ? " [packed]" : "");
        fflush(stdout);

        if (!keep) {
            printf("Switching to large slices. Ssl%u to be redone\n", s);
            break;
        }
        builder_push_small_slice(mb->mm, S);
        // transfer(Sq, S);
    }
    *p_i0 = s * mb->nrows_t / nslices;
}
/* }}} */

/* {{{ large slices */
void builder_do_all_large_slices(builder * mb, uint32_t * p_i0, unsigned int scratch1size)
{
    unsigned int rem_nrows = mb->nrows_t - *p_i0;
    unsigned int nlarge_slices = iceildiv(rem_nrows, LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nlarge_slices ; s++) {
        large_slice_t L[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nlarge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nlarge_slices;

        printf("Lsl%u %ss %u+%u", s, mb->rowname, i0, i1-i0);
        fflush(stdout);

        int keep = builder_do_large_slice(mb, L, i0, i1, scratch1size);

        if (!keep) {
            printf("Switching to huge slices. Lsl%u to be redone\n", s);
            break;
        }

        // transfer(Lq, L);
        builder_push_large_slice(mb->mm, L);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

/* {{{ huge slices */
void builder_do_all_huge_slices(builder * mb, uint32_t * p_i0, unsigned int scratch2size)
{
    unsigned int rem_nrows = mb->nrows_t - *p_i0;
    unsigned int nhuge_slices = iceildiv(rem_nrows, HUGE_MPLEX_MAX * LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nhuge_slices ; s++) {
        huge_slice_t H[1];
        uint32_t i0 = * p_i0 +  s      * rem_nrows / nhuge_slices;
        uint32_t i1 = * p_i0 + (s + 1) * rem_nrows / nhuge_slices;

        printf("Hsl%u %ss %u+%u", s, mb->rowname, i0, i1-i0);
        fflush(stdout);
        builder_do_huge_slice(mb, H, i0, i1, scratch2size);
        builder_push_huge_slice(mb->mm, H);
        // transfer(Hq, H);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

#define xxxCOMPRESS_COMBINERS_1
#define xxxCOMPRESS_COMBINERS_2
#define xxxCOMPRESS_COMBINERS_4

#if defined(COMPRESS_COMBINERS_1) || \
    defined(COMPRESS_COMBINERS_2) || \
    defined(COMPRESS_COMBINERS_4)
#ifdef MATMUL_SUB_VSC_COMBINE_H_
#error "Please either fix or disable the assembly code for COMPRESS_COMBINERS"
/* Anyway this idea doesn't work so well and changes little to the
 * running time. */
#endif
#endif


/* {{{ staircase */
unsigned long compressed_size(unsigned long s, unsigned int defer MAYBE_UNUSED)
{
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        return iceildiv(s, 8);
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        return iceildiv(s, 4);
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        return iceildiv(s, 2);
#endif
    } else {
        return s;
    }
}

void append_compressed(vector<uint8_t>& t8, vector<uint8_t> const& S, unsigned int defer MAYBE_UNUSED)
{
    if (0) {
#ifdef  COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
    } else {
        t8.insert(t8.end(), S.begin(), S.end());
    }
}

void builder_push_vsc_slices(struct matmul_bucket_data_s * mm, vsc_slice_t * V)
{
    printf("Flushing staircase slices\n");

    mm->headers.push_back(*V->hdr);

    mm->aux.push_back(V->dispatch.size());
    mm->aux.push_back(V->steps.size());

    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        mm->aux.push_back(V->steps[l].defer);
        mm->aux.push_back(V->steps[l].tbuf_space);
    }

    /* all headers */
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        mm->headers.push_back(*V->dispatch[k].hdr);
    }
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            mm->headers.push_back(*V->dispatch[k].sub[l].hdr);
        }
    }

    /* We also add header information for the combination operations,
     * even though the payload is stored alongside with the dispatching
     * stuff for convenience. We do need the combination headers for
     * properly keeping track of where the time is spent eventually.
     *
     * The combining headers are stored en route while flushing
     * combination points.
     */

    unsigned int cidx = mm->headers.size();

    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        slice_header_t chdr[1];
        memset(chdr, 0, sizeof(slice_header_t));
        chdr->t = SLICE_TYPE_DEFER_ROW;
        chdr->i0 = V->dispatch[0].sub[l].hdr->i0;
        chdr->i1 = V->dispatch[0].sub[l].hdr->i1;
        chdr->j0 = V->hdr->j0;
        chdr->j1 = V->hdr->j1;
        chdr->nchildren = 0;
        chdr->ncoeffs = 0;
        mm->headers.push_back(*chdr);
    }

    vector<pair<pair<unsigned int, unsigned int>, pair<uint64_t, uint64_t> > > csizes;

    // unsigned int o8 = mm->t8.size();

    /* dispatch data goes in dispatching order */
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        vsc_middle_slice_t& D(V->dispatch[k]);

        /* The stream of u16 values corresponding to this strip */
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            vsc_sub_slice_t & S(D.sub[l]);
            // printf("save dispatch(%zu), sum=%x\n", S.x.size(), idiotic_sum((void*) ptrbegin(S.x), S.x.size() * sizeof(uint16_t)));
            mm->t16.insert(mm->t16.end(), S.x.begin(), S.x.end());
            S.x.clear();
        }
        /* Now the combining data */
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            if (!flush_here(k, V->dispatch.size(), V->steps[l].defer))
                continue;
            unsigned int k0 = k - (k % V->steps[l].defer);
            /* Think about the combining header ! */
            slice_header_t chdr[1];
            memset(chdr, 0, sizeof(slice_header_t));
            chdr->t = SLICE_TYPE_DEFER_CMB;
            chdr->i0 = V->dispatch[k0].sub[l].hdr->i0;
            chdr->i1 = V->dispatch[k0].sub[l].hdr->i1;
            chdr->j0 = V->dispatch[k0].sub[l].hdr->j0;
            chdr->j1 = V->dispatch[k].sub[l].hdr->j1;
            chdr->nchildren = 0;
            chdr->ncoeffs = 0;
            /* flush this one */
            vsc_sub_slice_t & S(V->dispatch[k].sub[l]);
            for( ; k0 <= k ; k0++) {
                chdr->ncoeffs += V->dispatch[k0].sub[l].hdr->ncoeffs;
            }
            unsigned int defer = V->steps[l].defer;

            ASSERT_ALWAYS(S.c.size() == chdr->ncoeffs + V->steps[l].nrows);

            csizes.push_back(make_pair(
                        make_pair(k - (k % defer), l),
                        make_pair(0, S.c.size())));
            mm->headers[cidx + l].ncoeffs += chdr->ncoeffs;
            // printf("save combine(%zu), sum=%x\n", S.c.size(), idiotic_sum((void*) ptrbegin(S.c), S.c.size()));
            // printf("m=%u\n", *max_element(S.c.begin(), S.c.end()));
            append_compressed(mm->t8, S.c, defer);
            S.c.clear();
            mm->headers.push_back(*chdr);
        }
    }

    /* This fixes the reading order w.r.t. the positioning and size of
     * the combining slices */
    uint64_t p8 = 0;
    for(unsigned int i = 0 ; i < csizes.size() ; i++) {
        // unsigned int k = csizes[i].first.first;
        unsigned int l = csizes[i].first.second;
        unsigned long s = csizes[i].second.second;
        unsigned int defer = V->steps[l].defer;
        unsigned long count = compressed_size(s, defer);
        /*
        printf("straight-order (defer %u, vstrip #%u)"
                ": @%"PRIu64" combine(%"PRIu64"), sum=%x\n",
                V->steps[csizes[i].first.second].defer, csizes[i].first.first,
                csizes[i].second.first, csizes[i].second.second,
                idiotic_sum((void*) &*(mm->t8.begin() + o8 + csizes[i].second.first), csizes[i].second.second));
                */
        csizes[i].second.first = p8;
        p8 += count;
    }
    sort(csizes.begin(), csizes.end());
    mm->aux.push_back(2 * csizes.size() + 1);
    for(unsigned int i = 0 ; i < csizes.size() ; i++) {
        /*
        printf("straight-order (defer %u, vstrip #%u)"
                ": @%"PRIu64" combine(%"PRIu64"), sum=%x\n",
                V->steps[csizes[i].first.second].defer, csizes[i].first.first,
                csizes[i].second.first, csizes[i].second.second,
                idiotic_sum((void*) &*(mm->t8.begin() + o8 + csizes[i].second.first), csizes[i].second.second));
                */
        mm->aux.push_back(csizes[i].second.first);
        mm->aux.push_back(csizes[i].second.second);
    }
    mm->aux.push_back(p8);

    mm->public_->ncoeffs += V->hdr->ncoeffs;
}
/* }}} */

void matmul_bucket_build_cache(struct matmul_bucket_data_s * mm, uint32_t * data)
{
    builder mb[1];
    builder_init(mb, mm, data);

    printf("%u rows %u cols\n", mm->public_->dim[0], mm->public_->dim[1]);

    uint32_t main_i0 = 0;

    if (mm->methods.small1 || mm->methods.small2)
        builder_do_all_small_slices(mb, &main_i0, mm->npack);

    if (mm->methods.large)
        builder_do_all_large_slices(mb, &main_i0, mm->scratch1size);

    /* Note that vsc and huge are exclusive ! */
    if (mm->methods.vsc && main_i0 < mb->nrows_t) {
        vsc_slice_t V[1];
        builder_prepare_vsc_slices(mb, V, main_i0);
        mm->scratch3size = V->tbuf_space;
        builder_push_vsc_slices(mm, V);
    }
    if (mm->methods.huge && main_i0 < mb->nrows_t) {
        builder_do_all_huge_slices(mb, &main_i0, mm->scratch2size);
    }


    /* done, at last ! */

    mm_finish_init(mm);
}

int matmul_bucket_reload_cache(struct matmul_bucket_data_s * mm)/* {{{ */
{
    FILE * f;

    f = matmul_common_reload_cache_fopen(abbytes(mm->xab,1), mm->public_, MM_MAGIC);
    if (f == NULL) { return 0; }

    for( ;; ) {
        slice_header_t hdr[1];
        MATMUL_COMMON_READ_ONE16(hdr->t, f);
        MATMUL_COMMON_READ_ONE8(hdr->nchildren, f);
        MATMUL_COMMON_READ_MANY8(hdr->pad, 5, f);
        MATMUL_COMMON_READ_ONE32(hdr->i0, f);
        MATMUL_COMMON_READ_ONE32(hdr->i1, f);
        MATMUL_COMMON_READ_ONE32(hdr->j0, f);
        MATMUL_COMMON_READ_ONE32(hdr->j1, f);
        MATMUL_COMMON_READ_ONE64(hdr->ncoeffs, f);
        if (hdr->t == SLICE_TYPE_NONE)
            break;
        mm->headers.push_back(*hdr);
    }

    MATMUL_COMMON_READ_ONE32(mm->scratch1size, f);
    MATMUL_COMMON_READ_ONE32(mm->scratch2size, f);
    MATMUL_COMMON_READ_ONE32(mm->scratch3size, f);

    size_t n16, n8, naux;
    MATMUL_COMMON_READ_ONE32(n16, f);
    MATMUL_COMMON_READ_ONE32(n8, f);
    MATMUL_COMMON_READ_ONE32(naux, f);
    mm->t16.resize(n16);
    mm->t8.resize(n8);
    mm->aux.resize(naux);
    MATMUL_COMMON_READ_MANY16(ptrbegin(mm->t16), n16, f);
    MATMUL_COMMON_READ_MANY8(ptrbegin(mm->t8), n8, f);
    MATMUL_COMMON_READ_MANY32(ptrbegin(mm->aux), naux, f);

    fclose(f);

    mm_finish_init(mm);

    return 1;
}/*}}}*/

void matmul_bucket_save_cache(struct matmul_bucket_data_s * mm)/*{{{*/
{
    FILE * f;

    f = matmul_common_save_cache_fopen(abbytes(mm->xab,1), mm->public_, MM_MAGIC);
    for(unsigned int h = 0 ; h < mm->headers.size() ; h++) {
        slice_header_t * hdr = & (mm->headers[h]);
        MATMUL_COMMON_WRITE_ONE16(hdr->t, f);
        MATMUL_COMMON_WRITE_ONE8(hdr->nchildren, f);
        MATMUL_COMMON_WRITE_MANY8(hdr->pad, 5, f);
        MATMUL_COMMON_WRITE_ONE32(hdr->i0, f);
        MATMUL_COMMON_WRITE_ONE32(hdr->i1, f);
        MATMUL_COMMON_WRITE_ONE32(hdr->j0, f);
        MATMUL_COMMON_WRITE_ONE32(hdr->j1, f);
        MATMUL_COMMON_WRITE_ONE64(hdr->ncoeffs, f);
    }
    /* useful padding */
    MATMUL_COMMON_WRITE_ONE16(0, f);
    MATMUL_COMMON_WRITE_ONE16(0, f);
    MATMUL_COMMON_WRITE_ONE32(0, f);    // nchildren+padding
    MATMUL_COMMON_WRITE_ONE32(0, f);
    MATMUL_COMMON_WRITE_ONE32(0, f);
    MATMUL_COMMON_WRITE_ONE32(0, f);
    MATMUL_COMMON_WRITE_ONE32(0, f);
    MATMUL_COMMON_WRITE_ONE64(0, f);

    MATMUL_COMMON_WRITE_ONE32(mm->scratch1size, f);
    MATMUL_COMMON_WRITE_ONE32(mm->scratch2size, f);
    MATMUL_COMMON_WRITE_ONE32(mm->scratch3size, f);
    size_t n16 = mm->t16.size();
    size_t n8 = mm->t8.size();
    size_t naux = mm->aux.size();
    MATMUL_COMMON_WRITE_ONE32(n16, f);
    MATMUL_COMMON_WRITE_ONE32(n8, f);
    MATMUL_COMMON_WRITE_ONE32(naux, f);
    MATMUL_COMMON_WRITE_MANY16(ptrbegin(mm->t16), n16, f);
    MATMUL_COMMON_WRITE_MANY8(ptrbegin(mm->t8), n8, f);
    MATMUL_COMMON_WRITE_MANY32(ptrbegin(mm->aux), naux, f);

    fclose(f);
}/*}}}*/

static inline uint32_t read32(uint16_t const * & q) /* {{{ */
{
    uint32_t res;
    res = *q++;
    res |= ((uint32_t) *q++) << 16;
    return res;
} /* }}} */

#ifndef MATMUL_SUB_SMALL1_H_
static const uint16_t * matmul_sub_small1(abobj_t x, abt * where, abt const * from, const uint16_t * q, unsigned int count)
{
        for(uint32_t c = 0 ; c < count ; c++) {
            uint32_t j = *q++;
            uint32_t di = *q++;
            // ASSERT(j < (1UL << 16));
            // ASSERT(hdr->j0 + j < hdr->j1);
            // ASSERT(hdr->i0 + di < hdr->i1);
            abadd(x, where + aboffset(x, di), from + aboffset(x, j));
        }
        return q;
}
#endif

#ifndef MATMUL_SUB_SMALL1_TR_H_
static const uint16_t * matmul_sub_small1_tr(abobj_t x, abt * where, abt const * from, const uint16_t * q, unsigned int count)
{
    for(uint32_t c = 0 ; c < count ; c++) {
        uint32_t j = *q++;
        uint32_t di = *q++;
        // ASSERT(j < (1UL << 16));
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        abadd(x, where + aboffset(x, j), from + aboffset(x, di));
    }
    return q;
}
#endif

#ifndef MATMUL_SUB_SMALL2_H_
static const uint16_t * matmul_sub_small2(abobj_t x, abt * where, abt const * from, const uint16_t * q, unsigned int count)
{
    uint32_t j = 0;
    for(uint32_t c = 0 ; c < count ; c++) {
        uint16_t h = *q++;
        j += h >> SMALL_SLICES_I_BITS;
        // ASSERT(j < pos->ncols_t);
        uint32_t di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        abadd(x, where + aboffset(x, di), from + aboffset(x, j));
    }
    return q;
}
#endif

#ifndef MATMUL_SUB_SMALL2_TR_H_
static const uint16_t * matmul_sub_small2_tr(abobj_t x, abt * where, abt const * from, const uint16_t * q, unsigned int count)
{
    uint32_t j = 0;
    for(uint32_t c = 0 ; c < count ; c++) {
        uint16_t h = *q++;
        j += h >> SMALL_SLICES_I_BITS;
        // ASSERT(j < pos->ncols_t);
        uint32_t di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        abadd(x, where + aboffset(x, j), from + aboffset(x, di));
    }
    return q;
}
#endif

struct pos_desc {
    const uint16_t * q16;
    const uint8_t * q8;
    const unsigned int * ql;
    uint32_t i;
    uint32_t nrows_t;
    uint32_t ncols_t;
};

static inline void matmul_bucket_mul_small1_vblock(struct matmul_bucket_data_s * mm, slice_header_t * hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- small1 (dense) slices"); /* {{{ */
    int usual = d == ! mm->public_->store_transposed;
    abobj_ptr x = mm->xab;
    abt * where      = dst + aboffset(x, usual ? hdr->i0 : hdr->j0);
    abt const * from = src + aboffset(x, usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        abzero(x, where, usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0));
    }

    uint32_t ncoeffs_slice = hdr->ncoeffs;
    ASSERT_ALWAYS(pos->i == hdr->i0);   // FIXME -- should disappear.

    if (usual) {
        pos->q16 = matmul_sub_small1(x, where, from, pos->q16, ncoeffs_slice);
    } else {
        pos->q16 = matmul_sub_small1_tr(x, where, from, pos->q16, ncoeffs_slice);
    }
    ASM_COMMENT("end of small1 (dense) slices"); /* }}} */
}

static inline void matmul_bucket_mul_small2(struct matmul_bucket_data_s * mm, slice_header_t * hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- small2 (dense) slices"); /* {{{ */
    int usual = d == ! mm->public_->store_transposed;
    abobj_ptr x = mm->xab;
    abt * where      = dst + aboffset(x, usual ? hdr->i0 : hdr->j0);
    abt const * from = src + aboffset(x, usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        abzero(x, where, usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0));
    }

    uint32_t ncoeffs_slice = hdr->ncoeffs;
    ASSERT_ALWAYS(pos->i == hdr->i0);   // FIXME -- should disappear.

    if (usual) {
        pos->q16 = matmul_sub_small2(x, where, from, pos->q16, ncoeffs_slice);
    } else {
        pos->q16 = matmul_sub_small2_tr(x, where, from, pos->q16, ncoeffs_slice);
    }
    /* fix alignment in any case */
    pos->q16 += (2 - 1) & - ncoeffs_slice;
    ASM_COMMENT("end of small2 (dense) slices"); /* }}} */
}

static inline void prepare_buckets_and_fences(abobj_ptr x, abt ** b, abt ** f, abt * z, const unsigned int * ql, unsigned int n)
{
    for(unsigned int k = 0 ; k < n ; k++) {
        b[k] = z;
        f[k] = b[k] + ql[k];
        z += aboffset(x, ql[k]);
    }
}

static inline void prepare_buckets(abobj_ptr x, abt ** b, abt * z, const unsigned int * ql, unsigned int n)
{
    for(unsigned int k = 0 ; k < n ; k++) {
        b[k] = z;
        z += aboffset(x, ql[k]);
    }
}

#ifndef MATMUL_SUB_LARGE_FBI_H_
static inline void matmul_sub_large_fbi(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abcopy(x, sb[*q], z, 1);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}
#endif

#ifndef MATMUL_SUB_LARGE_FBI_TR_H_
static inline void
matmul_sub_large_fbi_tr(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abadd(x, z, sb[*q]);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}
#endif

#ifndef MATMUL_SUB_LARGE_ASB_H_
// static void matmul_sub_large_asb(abobj_ptr x, abt * dst, const abt * z, const uint8_t * q, const unsigned int * ql) __attribute__((__noinline__));
static void matmul_sub_large_asb(abobj_ptr x, abt * dst, const abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* This ``applies'' the LSL_NBUCKETS_MAX small buckets whose
     * respective lengths are given by ql[0] to ql[LSL_NBUCKETS_MAX-1].
     *
     * For ql[k] <= i < ql[k+1], z[i] is added to dst[k*256+qk[i]], with
     * qk = q + ql[0] + ... + ql[k-1].
     */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abadd(x, dst + aboffset(x, *q), z);
            z+= aboffset(x, 1);
            q++;
        }
        dst += aboffset(x, 256);
    }
}
#endif

#ifndef MATMUL_SUB_LARGE_ASB_TR_H_
static inline void matmul_sub_large_asb_tr(abobj_ptr x, const abt * src, abt * z, const uint8_t * q, const unsigned int * ql)
{
    /* converse of the above */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            abcopy(x, z, src + aboffset(x, *q), 1);
            z += aboffset(x, 1);
            q++;
        }
        src += aboffset(x, 256);
    }
}
#endif

#ifndef MATMUL_SUB_LARGE_FBD_H_
static void
matmul_sub_large_fbd(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[n] such that z[i] is in array
     * pointed to by sb[q[i]]. Exactly n coefficients are expected. All
     * the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, sb[q[c]], z, 1);
        sb[q[c]]+= aboffset(x, 1);
        z += aboffset(x, 1);
    }
}
#endif

#ifndef MATMUL_SUB_LARGE_FBD_TR_H_
static inline void
matmul_sub_large_fbd_tr(abobj_ptr x, abt ** sb, abt * z, const uint8_t * q, unsigned int n)
{
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        abcopy(x, z, sb[q[c]], 1);
        sb[q[c]]+= aboffset(x, 1);
        z += aboffset(x, 1);
    }
}
#endif

static inline void matmul_bucket_mul_large(struct matmul_bucket_data_s * mm, slice_header_t * hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- large (sparse) slices"); /* {{{ */

    abobj_ptr x = mm->xab;

    int usual = d == ! mm->public_->store_transposed;

    abt * where      = dst + aboffset(x, usual ? hdr->i0 : hdr->j0);
    abt const * from = src + aboffset(x, usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        abzero(x, where, usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0));
    }

    abt * scratch = mm->scratch1;

    uint32_t j = 0;

    if (usual) {
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            uint32_t n = *pos->ql++;
            abt * bucket[LSL_NBUCKETS_MAX];
            abt const * inp = from + aboffset(x, j);
            prepare_buckets(x,bucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi(x, bucket, inp, pos->q8, n);
            matmul_sub_large_asb(x, where, scratch, pos->q8+2*n, pos->ql);
            pos->q8 += 3*n;
            // fix alignment !
            pos->q8 += n & 1;
            pos->ql += LSL_NBUCKETS_MAX;
            j = j1;
        }
    } else {
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            uint32_t n = *pos->ql++;
            abt * bucket[LSL_NBUCKETS_MAX];
            abt * outp = where + aboffset(x, j);
            prepare_buckets(x,bucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
            matmul_sub_large_asb_tr(x, from, scratch, pos->q8+2*n, pos->ql);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi_tr(x, bucket, outp, pos->q8, n);
            pos->q8 += 3 * n;
            // fix alignment !
            pos->q8 += n & 1;
            pos->ql += LSL_NBUCKETS_MAX;
            j = j1;
        }
    }
    ASM_COMMENT("end of large (sparse) slices"); /* }}} */
}

static inline void matmul_bucket_mul_huge(struct matmul_bucket_data_s * mm, slice_header_t * hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- huge (very sparse) slices"); /* {{{ */

    abobj_ptr x = mm->xab;

    int usual = d == ! mm->public_->store_transposed;

    abt * where      = dst + aboffset(x, usual ? hdr->i0 : hdr->j0);
    abt const * from = src + aboffset(x, usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        abzero(x, where, usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0));
    }

    abt * scratch = mm->scratch1;

    uint32_t j = 0;

    uint32_t di = hdr->i1 - hdr->i0;
    unsigned int nlarge = *pos->ql++;
    uint32_t di_sub = iceildiv(di, nlarge);

    if (usual) {
        /* ok -- hold your breath a second. */
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            unsigned int n = *pos->ql++;
            ASSERT_ALWAYS(n <= mm->scratch2size);
            abt * scratch2 = mm->scratch2;
            abt const * inp = src + aboffset(x, j);
            abt * bucket[HUGE_MPLEX_MAX];
            const unsigned int * Lsizes = pos->ql;
            prepare_buckets(x,bucket,scratch2,pos->ql,nlarge);
            pos->ql += nlarge;
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi(x, bucket, inp, pos->q8, n);
            pos->q8 += 2 * n;
            for(unsigned int k = 0 ; k < nlarge ; k++) {
                abt * sbucket[LSL_NBUCKETS_MAX];
                prepare_buckets(x,sbucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
                bucket[k] -= aboffset(x, Lsizes[k]);
                matmul_sub_large_fbd(x, sbucket, bucket[k], pos->q8, Lsizes[k]);
                pos->q8 += Lsizes[k];
                abt * outp = where + aboffset(x, k * di_sub);
                matmul_sub_large_asb(x, outp, scratch, pos->q8, pos->ql);
                pos->q8 += Lsizes[k];
                pos->ql += LSL_NBUCKETS_MAX;
            }
            j = j1;
        }
    } else {
        for( ; j < pos->ncols_t ; ) {
            uint32_t j1 = j + *pos->ql++;
            unsigned int n = *pos->ql++;
            abt * scratch2 = mm->scratch2;
            abt * outp = dst + aboffset(x, j);
            abt * bucket[HUGE_MPLEX_MAX];
            const unsigned int * Lsizes = pos->ql;
            prepare_buckets(x,bucket,scratch2,pos->ql,nlarge);
            pos->ql += nlarge;
            const uint8_t * q8_saved = pos->q8;
            pos->q8 += 2 * n;
            for(unsigned int k = 0 ; k < nlarge ; k++) {
                abt * sbucket[LSL_NBUCKETS_MAX];
                prepare_buckets(x,sbucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
                const uint8_t * fill = pos->q8;
                const uint8_t * apply = pos->q8 + Lsizes[k];
                const abt * inp = from + aboffset(x, k * di_sub);
                matmul_sub_large_asb_tr(x, inp, scratch, apply, pos->ql);

                matmul_sub_large_fbd_tr(x, sbucket, bucket[k], fill, Lsizes[k]);

                pos->q8 += 2 * Lsizes[k];
                pos->ql += LSL_NBUCKETS_MAX;
            }
            swap(pos->q8, q8_saved);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi_tr(x, bucket, outp, pos->q8, n);
            pos->q8 += 2 * n;
            swap(pos->q8, q8_saved);
            j = j1;
        }
    }
    ASM_COMMENT("end of huge (very sparse) slices"); /* }}} */
}

#ifndef MATMUL_SUB_VSC_DISPATCH_H
static inline void matmul_sub_vsc_dispatch(abobj_ptr x, abt * dst, abt const * src, const uint16_t * q, unsigned int count)
{
    // printf("dispatch(%u), sum=%x\n", count, idiotic_sum((void*)q, count * sizeof(uint16_t)));
    for( ; count-- ; ) {
        abcopy(x, dst, src + aboffset(x, *q++), 1);
        dst += aboffset(x, 1);
    }
}
#endif

#ifndef MATMUL_SUB_VSC_COMBINE_H_
static inline void matmul_sub_vsc_combine(abobj_ptr x, abt * dst, const abt * * mptrs, const uint8_t * q, unsigned int count, unsigned int defer MAYBE_UNUSED)
{
    // printf("combine(%u), defer %u\n", count, defer);
    // printf("combine(%u), sum=%x\n", count, idiotic_sum((void*)q, compressed_size(count, defer)));
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abadd(x, dst, mptrs[c]);
                mptrs[c] += aboffset(x, c != 0);
                dst += aboffset(x, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abadd(x, dst, mptrs[c]);
                mptrs[c] += aboffset(x, c != 0);
                dst += aboffset(x, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abadd(x, dst, mptrs[c]);
                mptrs[c] += aboffset(x, c != 0);
                dst += aboffset(x, c == 0);
            }
        }
#endif
    } else {
        for( ; count-- ; ) {
            uint8_t c = *q++;
                ASSERT(c <= defer);
            abadd(x, dst, mptrs[c]);
            mptrs[c] += aboffset(x, c != 0);
            dst += aboffset(x, c == 0);
        }
    }
}
#endif

#ifndef MATMUL_SUB_VSC_COMBINE_TR_H_
static inline void matmul_sub_vsc_combine_tr(abobj_ptr x, abt ** mptrs, const abt * qw, const uint8_t * z, unsigned int count, unsigned int defer MAYBE_UNUSED)
{
    // printf("uncombine(%u), defer %u\n", count, defer);
    // printf("uncombine(%u), sum=%x\n", count, idiotic_sum((void*)z, compressed_size(count, defer)));
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abcopy(x, mptrs[c], qw, 1);
                mptrs[c] += aboffset(x, c != 0);
                qw += aboffset(x, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abcopy(x, mptrs[c], qw, 1);
                mptrs[c] += aboffset(x, c != 0);
                qw += aboffset(x, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                abcopy(x, mptrs[c], qw, 1);
                mptrs[c] += aboffset(x, c != 0);
                qw += aboffset(x, c == 0);
            }
        }
#endif
    } else {
        for( ; count-- ; ) {
            uint8_t c = *z++;
            ASSERT(c <= defer);
            abcopy(x, mptrs[c], qw, 1);
            mptrs[c] += aboffset(x, c != 0);
            qw += aboffset(x, c == 0);
        }
    }
}
#endif

#ifndef MATMUL_SUB_VSC_DISPATCH_TR_H_
static inline void matmul_sub_vsc_dispatch_tr(abobj_ptr x, abt * qr, const abt * q, const uint16_t * z, unsigned int count)
{
    // printf("undispatch(%u), sum=%x\n", count, idiotic_sum((void*)z, count * sizeof(uint16_t)));
    for( ; count-- ; ) {
        abadd(x, qr + aboffset(x, *z++), q);
        q += aboffset(x, 1);
    }
}
#endif

/* in preparation for different steps, including the matrix
 * multiplication itself, we are rebuilding the vsc_slice_t object, or at
 * least its skeleton ; it makes it easy to run the control loops.  The
 * real data payload is of course not copied again ! */

static inline void rebuild_vsc_slice_skeleton(
        struct matmul_bucket_data_s * mm,
        vector<slice_header_t>::iterator & hdr, 
        vsc_slice_t * V, struct pos_desc * pos,
        unsigned int & nvstrips,
        unsigned int & Midx,
        unsigned int & Didx,
        unsigned int & Ridx,
        unsigned int & Cidx)
{
    *V->hdr = *hdr++;
    nvstrips = *pos->ql++;
    V->dispatch.resize(nvstrips);
    V->steps.resize(*pos->ql++);
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        V->steps[l].defer = *pos->ql++;
        V->steps[l].tbuf_space = *pos->ql++;
    }

    Midx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        vsc_middle_slice_t &D (V->dispatch[k]);
        *D.hdr = *hdr++;
    }
    Didx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        vsc_middle_slice_t &D (V->dispatch[k]);
        D.sub.resize(V->steps.size());
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            vsc_sub_slice_t & S(D.sub[l]);
            *S.hdr = *hdr++;
            V->steps[l].nrows = S.hdr->i1 - S.hdr->i0;
        }
    }
    /* Skip over the combining headers as well */
    Ridx = hdr - mm->headers.begin();
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        hdr++;
    }
    Cidx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            if (!flush_here(k, V->dispatch.size(), V->steps[l].defer))
                continue;
            hdr++;
        }
    }
}
static inline void matmul_bucket_mul_vsc(struct matmul_bucket_data_s * mm, vector<slice_header_t>::iterator & hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    abobj_ptr x = mm->xab;

    int usual = d == ! mm->public_->store_transposed;

    // abt const * from = src + aboffset(x, usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        abt * where      = dst + aboffset(x, usual ? hdr->i0 : hdr->j0);
        abzero(x, where, usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0));
    }

    abt * scratch = mm->scratch3;

    /* {{{ */

    vsc_slice_t V[1];
    unsigned int nvstrips;
    unsigned int Midx;
    unsigned int Didx;
    unsigned int Ridx;
    unsigned int Cidx;
    rebuild_vsc_slice_skeleton(mm, hdr, V, pos, nvstrips, Midx, Didx, Ridx, Cidx);

    /*}}}*/

    ASM_COMMENT("multiplication code -- vertical staircase");/*{{{*/

    /* now prepare pointers */
    vector<abt *> base_ptrs;
    vector<abt *> cptrs;
    abt * q0 = scratch;
    abt * dummy = q0;
    abzero(x, dummy, 1);
    q0++;
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        base_ptrs.push_back(q0);
        cptrs.push_back(q0);
        q0 += V->steps[l].tbuf_space;
    }
    base_ptrs.push_back(q0);

    if (usual) {
        uint32_t skipover = *pos->ql++;
        pos->ql += skipover;
        for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
            vsc_middle_slice_t const & D(V->dispatch[k]);
            const abt * qr = src + aboffset(x, D.hdr->j0);
            mm->slice_timings[Midx].t -= wct_seconds();
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_sub_slice_t const & S(D.sub[l]);
                abt * q = cptrs[l];
                unsigned int count = S.hdr->ncoeffs;
                ASSERT(base_ptrs[l] <= q);
                ASSERT(q <= base_ptrs[l+1]);
                mm->slice_timings[Didx].t -= wct_seconds();
                matmul_sub_vsc_dispatch(x, q, qr, pos->q16, count);
                q += aboffset(x, count);
                pos->q16 += count;
                mm->slice_timings[Didx].t += wct_seconds();
                Didx++;
                ASSERT(q <= base_ptrs[l+1]);
                cptrs[l] = q;
            }
            mm->slice_timings[Midx].t += wct_seconds();
            Midx++;

            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_sub_slice_t const & S(D.sub[l]);
                unsigned int defer = V->steps[l].defer;
                if (!flush_here(k,nvstrips,defer))
                    continue;

                /* our different read pointers */
                vector<abt const *> mptrs;
                mptrs.reserve(defer + 1);
                abt const * q = base_ptrs[l];
                mptrs.push_back(dummy);
                for(unsigned int k0 = k - k % defer ; k0 <= k ; k0++) {
                    mptrs.push_back(q);
                    q += V->dispatch[k0].sub[l].hdr->ncoeffs;
                }

                abt * qw = dst + aboffset(x, S.hdr->i0);
                unsigned int count = mm->headers[Cidx].ncoeffs + V->steps[l].nrows;
                ASSERT(V->steps[l].nrows == S.hdr->i1 - S.hdr->i0);
                ASSERT(q - base_ptrs[l] == (ptrdiff_t) mm->headers[Cidx].ncoeffs);

                double t = wct_seconds();
                mm->slice_timings[Cidx].t -= t;
                mm->slice_timings[Ridx+l].t -= t;
                matmul_sub_vsc_combine(x, qw, ptrbegin(mptrs), pos->q8, count, defer);
                pos->q8 += compressed_size(count, defer);
                t = wct_seconds();
                mm->slice_timings[Cidx].t += t;
                mm->slice_timings[Ridx+l].t += t;
                Cidx++;
                cptrs[l]=base_ptrs[l];
            }
        }
        ASSERT((ptrdiff_t) Cidx == hdr - mm->headers.begin());
    } else {
        /* There's quite an annoying difficulty here. Combining buffers
         * are not read in the same order here, because ``combining'' (whose
         * meaning turns into actually dispatching) occurs of course
         * _earlier_ in the operation.
         */
        pos->ql++;      // only the count, for fast skipover.
        for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
            vsc_middle_slice_t const & D(V->dispatch[k]);
            abt * qr = dst + aboffset(x, D.hdr->j0);
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_sub_slice_t const & S(D.sub[l]);
                unsigned int defer = V->steps[l].defer;

                /* Here the test is not the same. We have to fill the
                 * scratch buffers ahead of time. */
                if (k % defer)
                    continue;

                /* our different _write_ pointers */
                vector<abt *> mptrs;
                mptrs.reserve(defer + 1);
                abt * q = base_ptrs[l];
                cptrs[l] = q;
                mptrs.push_back(dummy);
                for(unsigned int k0 = k ; k0 <= when_flush(k,nvstrips,defer) ; k0++) {
                    mptrs.push_back(q);
                    q += V->dispatch[k0].sub[l].hdr->ncoeffs;
                }

                const abt * qw = src + aboffset(x, S.hdr->i0);
                // unsigned int count = mm->headers[Cidx].ncoeffs + V->steps[l].nrows;
                const uint8_t * z = pos->q8 + *pos->ql++;
                unsigned int count = *pos->ql++;

                double t = wct_seconds();
                // mm->slice_timings[Cidx].t -= t;
                mm->slice_timings[Ridx+l].t -= t;
                matmul_sub_vsc_combine_tr(x, ptrbegin(mptrs), qw, z, count, defer);
                z += compressed_size(count, defer);
                t = wct_seconds();
                // mm->slice_timings[Cidx].t += t;
                mm->slice_timings[Ridx+l].t += t;
                // Cidx++;
            }

            mm->slice_timings[Midx].t -= wct_seconds();
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_sub_slice_t const & S(D.sub[l]);
                abt * q = cptrs[l];
                unsigned int count = S.hdr->ncoeffs;
                ASSERT(base_ptrs[l] <= q);
                ASSERT(q <= base_ptrs[l+1]);
                mm->slice_timings[Didx].t -= wct_seconds();
                matmul_sub_vsc_dispatch_tr(x, qr, q, pos->q16, count);
                cptrs[l] += count;
                pos->q16 += count;
                mm->slice_timings[Didx].t += wct_seconds();
                Didx++;
                ASSERT(q <= base_ptrs[l+1]);
            }
            mm->slice_timings[Midx].t += wct_seconds();
            Midx++;
        }
        // ASSERT(Cidx == hdr - mm->headers.begin());
        pos->q8 += *pos->ql++;
    }

    ASM_COMMENT("end of vertical staircase");/*}}}*/

    /* It's an envelope type, so make sure the iterator is properly
     * placed eventually */
    hdr--;
}

///////////////////////////////////////////////////////////////////////
// just count how many times an iteration schedules a coefficient in
// the fbi/fbd/asb routines.

static inline void mm_finish_init(struct matmul_bucket_data_s * mm)
{
    struct pos_desc pos[1];

    pos->q16 = ptrbegin(mm->t16);
    pos->q8 = ptrbegin(mm->t8);
    pos->ql = ptrbegin(mm->aux);
    pos->i = 0;
    pos->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    pos->ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    /*
    for(uint16_t h = 0 ; h < mm->headers.size() ; h++) {
        slice_header_t * hdr = & (mm->headers[h]);
        printf("block %d %s [%d..%d[ x [%d..%d[, %d cld, %"PRIu64" coeffs\n",
                h, slice_name(hdr->t),
                hdr->i0, hdr->i1,
                hdr->j0, hdr->j1,
                hdr->nchildren, hdr->ncoeffs);
    }
    */
    for(vector<slice_header_t>::iterator hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        switch(hdr->t) {
            case SLICE_TYPE_SMALL1:
                break;
            case SLICE_TYPE_SMALL1_VBLOCK:
                // mm->mms1_ncoeffs += hdr->ncoeffs;
                pos->q16 += 2 * hdr->ncoeffs;
                break;
            case SLICE_TYPE_SMALL2:
                pos->q16 += hdr->ncoeffs;
                pos->q16 += (2 - 1) & - hdr->ncoeffs;
                // mm->mms2_ncoeffs += hdr->ncoeffs;
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        uint32_t n = *pos->ql++;
                        // mm->fbi_ncoeffs += n;
                        // mm->asb_ncoeffs += n;
                        pos->q8 += 3*n;
                        pos->q8 += n & 1;
                        pos->ql += LSL_NBUCKETS_MAX;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    unsigned int nlarge = *pos->ql++;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        unsigned int n = *pos->ql++;
                        // mm->fbi_ncoeffs += n;
                        // mm->fbd_ncoeffs += n;
                        // mm->asb_ncoeffs += n;
                        pos->ql += nlarge;
                        pos->ql += nlarge * LSL_NBUCKETS_MAX;
                        pos->q8 += 4*n;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_DEFER_ENVELOPE:
                {
                    vsc_slice_t V[1];
                    unsigned int nvstrips;
                    unsigned int Midx;
                    unsigned int Didx;
                    unsigned int Ridx;
                    unsigned int Cidx;
                    rebuild_vsc_slice_skeleton(mm, hdr, V, pos, nvstrips, Midx, Didx, Ridx, Cidx);
                    hdr--;
                    for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
                        unsigned int defer = V->steps[s].defer;
                        printf("Rows %"PRIu32"+%"PRIu32, mm->headers[Ridx++].i0, V->steps[s].nrows);
                        if (V->steps[s].density_upper_bound != UINT_MAX) {
                            printf(": d < %u", V->steps[s].density_upper_bound);
                        }
                        printf("; %u flushes (every %u), tbuf space %lu\n", iceildiv(nvstrips, defer), defer, V->steps[s].tbuf_space);
                    }
                }
                break;
            default:
                fprintf(stderr, "Unexpected block %s encountered\n",
                        slice_name(hdr->t));
                break;
        }
        if (hdr->j1 == pos->ncols_t) {
            pos->i = hdr->i1;
        }
    }

#if 0
    for(uint16_t h = 0 ; h < mm->headers.size() ; h++) {
        slice_header_t * hdr = & (mm->headers[h]);
        ASSERT_ALWAYS(pos->i == hdr->i0);
        switch(hdr->t) {
            case SLICE_TYPE_SMALL1_VBLOCK:
                mm->mms1_ncoeffs += hdr->ncoeffs;
                pos->q16 += 2 * hdr->ncoeffs;
                break;
            case SLICE_TYPE_SMALL2:
                pos->q16 += hdr->ncoeffs;
                pos->q16 += (2 - 1) & - hdr->ncoeffs;
                mm->mms2_ncoeffs += hdr->ncoeffs;
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        uint32_t n = *pos->ql++;
                        mm->fbi_ncoeffs += n;
                        mm->asb_ncoeffs += n;
                        // pos->q8 += 3*n;
                        pos->ql += LSL_NBUCKETS_MAX;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    unsigned int nlarge = *pos->ql++;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        unsigned int n = *pos->ql++;
                        mm->fbi_ncoeffs += n;
                        mm->fbd_ncoeffs += n;
                        mm->asb_ncoeffs += n;
                        pos->ql += nlarge;
                        pos->ql += nlarge * LSL_NBUCKETS_MAX;
                        // pos->q8 += 4*n;
                        j = j1;
                    }
                }
            case SLICE_TYPE_DEFER_ENVELOPE:
            default:
                break;
        }
        if (hdr->j1 == pos->ncols_t) {
            pos->i = hdr->i1;
        }
    }
#endif

    mm->scratch1 = (abt *) abinit(mm->xab, mm->scratch1size);
    mm->scratch2 = (abt *) abinit(mm->xab, mm->scratch2size);
    mm->scratch3 = (abt *) abinit(mm->xab, mm->scratch3size);

    mm->slice_timings.resize(mm->headers.size());
}

void matmul_bucket_zero_stats(struct matmul_bucket_data_s * mm)
{
    vector<slice_header_t>::iterator hdr;
    for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        unsigned int hidx = hdr - mm->headers.begin();
        mm->slice_timings[hidx].t = 0;
    }
}

static inline void matmul_bucket_mul_small1(struct matmul_bucket_data_s * mm, vector<slice_header_t>::iterator & hdr, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    slice_header_t & h(*hdr++);
    for(unsigned int i = 0 ; i < h.nchildren ; i++, hdr++) {
        ASSERT(hdr != mm->headers.end());
        ASSERT(hdr->t == SLICE_TYPE_SMALL1_VBLOCK);
        mm->slice_timings[hdr - mm->headers.begin()].t -= wct_seconds();
        matmul_bucket_mul_small1_vblock(mm, &*hdr, dst, src, d, pos);
        mm->slice_timings[hdr - mm->headers.begin()].t += wct_seconds();
    }
    hdr--;
}

static inline void matmul_bucket_mul_loop(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d, struct pos_desc * pos)
{
    vector<slice_header_t>::iterator hdr;

    for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        unsigned int hidx = hdr - mm->headers.begin();
        mm->slice_timings[hidx].t -= wct_seconds();
        switch(hdr->t) {
            case SLICE_TYPE_SMALL2:
                matmul_bucket_mul_small2(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_SMALL1:
                matmul_bucket_mul_small1(mm, hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                matmul_bucket_mul_large(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                matmul_bucket_mul_huge(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_DEFER_ENVELOPE:
                matmul_bucket_mul_vsc(mm, hdr, dst, src, d, pos);
                break;
            /* some slice types are ``contained'', and we should never
             * see them. NOTE that this implies in particular that the
             * timings for type "defer" (DEFER_ENVELOPE) is counted here,
             * thus including all overhead, while the timing for the
             * inner objects is counted from within the subroutines. */
            case SLICE_TYPE_SMALL1_VBLOCK:
            case SLICE_TYPE_DEFER_COLUMN:
            case SLICE_TYPE_DEFER_ROW:
            case SLICE_TYPE_DEFER_CMB:
            case SLICE_TYPE_DEFER_DIS:
            /* The current implementation no longer accepts data not
             * obeying the header structure, so the default branch also
             * aborts. */
            default:
                ASSERT_ALWAYS(0);
                break;
                break;
        }
        if (hdr->j1 == pos->ncols_t) { pos->i = hdr->i1; }
        mm->slice_timings[hidx].t += wct_seconds();
    }
}

void matmul_bucket_mul(struct matmul_bucket_data_s * mm, abt * dst, abt const * src, int d)
{
    struct pos_desc pos[1];

    pos->q16 = ptrbegin(mm->t16);
    pos->q8 = ptrbegin(mm->t8);
    pos->ql = ptrbegin(mm->aux);
    pos->i = 0;
    pos->nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    pos->ncols_t = mm->public_->dim[!mm->public_->store_transposed];

    abobj_ptr x = mm->xab;

    if (d == !mm->public_->store_transposed) {
        /* This is the ``normal'' case (matrix times vector). */
    } else {
        /* d == mm->public_->store_transposed */
        /* BEWARE, it's a priori sub-optimal ! In practice, the
         * difference isn't so striking though. */
        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }
        /* We zero out the dst area beforehand */
        abzero(x, dst, pos->ncols_t);
    }
    matmul_bucket_mul_loop(mm, dst, src, d, pos);

    mm->public_->iteration[d]++;
}

void matmul_bucket_report_vsc(struct matmul_bucket_data_s * mm, double scale, vector<slice_header_t>::iterator & hdr)
{
    uint64_t scale0;
    scale0 = (mm->public_->iteration[0] + mm->public_->iteration[1]);
    unsigned int nvstrips = hdr->nchildren;
    hdr++;
    unsigned int nsteps = hdr->nchildren;
    vector<pair<uint64_t, double> > dtime(nsteps);
    vector<pair<uint64_t, double> > ctime(nsteps);
    for(unsigned int k = 0 ; k < nvstrips ; k++) {
        ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_COLUMN);
        hdr++;
    }
    // hdr+=nvstrips;
    for(unsigned int k = 0 ; k < nvstrips ; k++) {
        for(unsigned int l = 0 ; l < nsteps ; l++) {
            ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_DIS);
            double t = mm->slice_timings[hdr - mm->headers.begin()].t;
            uint64_t nc = hdr->ncoeffs;
            dtime[l].first += nc;
            dtime[l].second += t;
            hdr++;
        }
    }
    double total_from_defer_rows = 0;
    double total_from_defer_cmbs = 0;
    for(unsigned int l = 0 ; l < nsteps ; l++) {
        ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_ROW);
        double t = mm->slice_timings[hdr - mm->headers.begin()].t;
        uint64_t nc = hdr->ncoeffs;
        ctime[l].first += nc;
        ctime[l].second += t;
        total_from_defer_rows+=t;
        hdr++;
    }
    /* Skip the combining blocks, because they're accounted for already
     * by the row blocks */
    for( ; hdr != mm->headers.end() && hdr->t == SLICE_TYPE_DEFER_CMB ; hdr++) {
        double t = mm->slice_timings[hdr - mm->headers.begin()].t;
        total_from_defer_cmbs+=t;
    }
    /* Some jitter will appear if transposed mults are performed, because
     * for the moment transposed mults don't properly store timing info
     */
    /*
    printf("jitter %.2f - %.2f = %.2f\n",
            total_from_defer_rows / scale0, total_from_defer_cmbs / scale0,
            (total_from_defer_rows - total_from_defer_cmbs) / scale0);
            */
    for(unsigned int l = 0 ; l < nsteps ; l++) {
        double t;
        uint64_t nc;
        double a;
        nc = dtime[l].first;
        t = dtime[l].second / scale0;
        a = 1.0e9 * t / nc;
        printf("defer\t%.2fs ; n=%-9"PRIu64" ; %5.2f ns/c ;"
            " scaled*%.2f : %5.2f/c\n",
            t, nc, a, scale, a * scale);
        nc = ctime[l].first;
        t = ctime[l].second / scale0;
        a = 1.0e9 * t / nc;
        printf("      + %.2fs ; n=%-9"PRIu64" ; %5.2f ns/c ;"
            " scaled*%.2f : %5.2f/c\n",
            t, nc, a, scale, a * scale);
    }
    hdr--;
}


void matmul_bucket_report(struct matmul_bucket_data_s * mm, double scale)
{
    uint64_t scale0;
    scale0 = (mm->public_->iteration[0] + mm->public_->iteration[1]);

    printf("n %"PRIu64"\n", mm->public_->ncoeffs);

    vector<slice_header_t>::iterator hdr;

    for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        if (hdr->t == SLICE_TYPE_SMALL1_VBLOCK) continue;
        if (hdr->t == SLICE_TYPE_DEFER_ENVELOPE) {
            matmul_bucket_report_vsc(mm, scale, hdr);
            continue;
        }
        double t = mm->slice_timings[hdr - mm->headers.begin()].t;
        uint64_t nc = hdr->ncoeffs;
        t /= scale0;
        double a = 1.0e9 * t / nc;
        printf("%s\t%.2fs ; n=%-9"PRIu64" ; %5.2f ns/c ;"
            " scaled*%.2f : %5.2f/c\n",
            slice_name(hdr->t), t,
            nc, a, scale, a * scale);
    }
    for(int i = 0 ; i < 40 ; i++) putchar('-');
    putchar('\n');
    for(int i = 0 ; i < SLICE_TYPE_MAX ; i++) {
        if (i == SLICE_TYPE_SMALL1_VBLOCK) continue;
        double t = 0;
        uint64_t nc = 0;
        for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
            if (hdr->t != i) continue;
            nc += hdr->ncoeffs;
            t += mm->slice_timings[hdr - mm->headers.begin()].t;
        }
        if (nc == 0) continue;
        t /= scale0;
        double a = 1.0e9 * t / nc;
        printf("%s\t%.2fs ; n=%-9"PRIu64" ; %5.2f ns/c ;"
            " scaled*%.2f : %5.2f/c\n",
            slice_name(i), t,
            nc, a, scale, a * scale);
    }
}

void matmul_bucket_auxv(struct matmul_bucket_data_s * mm MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
    if (op == MATMUL_AUX_ZERO_STATS) {
        matmul_bucket_zero_stats(mm);
        mm->public_->iteration[0] = 0;
        mm->public_->iteration[1] = 0;
    }
#if 0
    if (op == MATMUL_AUX_GET_READAHEAD) {
        unsigned int * res = va_arg(ap, unsigned int *);
        ++*res;
    }
#endif
}


void matmul_bucket_aux(struct matmul_bucket_data_s * mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_bucket_auxv (mm, op, ap);
    va_end(ap);
}
/* vim: set sw=4: */
