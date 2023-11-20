#include "cado.h"

#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"
#include "matmul_top_comm.hpp"
#include "timing.h"
#include "portability.h"

/* At some point we had this. Not sure it's still useful. */
#define ABASE_UNIVERSAL_READAHEAD_ITEMS 8

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

mmt_vec::mmt_vec(matmul_top_data & mmt, arith_generic * abase, pi_datatype_ptr pitype, int d, int flags, unsigned int n)
{
    mmt_vec_setup(*this, mmt, abase, pitype, d, flags, n);
}

/* this is for a vector which will be of interest to a group of threads
 * and jobs in direction d */
void mmt_vec_setup(mmt_vec & v, matmul_top_data & mmt, arith_generic * abase, pi_datatype_ptr pitype, int d, int flags, unsigned int n)
{
    if (abase == NULL) abase = mmt.abase;
    if (pitype == NULL) pitype = mmt.pitype;
    v.pi = mmt.pi;
    v.d = d;
    v.abase = abase;
    v.pitype = pitype;
    v.n = n;

    ASSERT_ALWAYS(n % mmt.pi->m->totalsize == 0);

    pi_comm_ptr wr = mmt.pi->wr[d];
    pi_comm_ptr xwr = mmt.pi->wr[!d];

    /* now what is the size which we are going to allocate locally */
    n /= xwr->totalsize;
    v.i0 = n * (xwr->jrank * xwr->ncores + xwr->trank);
    v.i1 = v.i0 + n;

    /* Look for readahead settings for all submatrices */
    n += ABASE_UNIVERSAL_READAHEAD_ITEMS;
    for(auto const & Mloc : mmt.matrices) {
        matmul_aux(Mloc.mm, MATMUL_AUX_GET_READAHEAD, &n);
    }

    if (flags & THREAD_SHARED_VECTOR) {
        if (wr->trank == 0) {
            v.v = abase->alloc(n, ALIGNMENT_ON_ALL_BWC_VECTORS);
            abase->vec_set_zero(v.v, n);
        }
        pi_thread_bcast(&v.v, sizeof(void*), BWC_PI_BYTE, 0, wr);
        v.siblings = NULL;
    } else {
        v.v = abase->alloc(n, ALIGNMENT_ON_ALL_BWC_VECTORS);
        abase->vec_set_zero(v.v, n);
        v.siblings = (mmt_vec **) shared_malloc(wr, wr->ncores * sizeof(mmt_vec *));
        v.siblings[wr->trank] = &v;
    }
    /* Vectors begin initialized to zero, so we have full consistency */
    v.consistency = 2;
    serialize_threads(v.pi->m);

    // pi_log_op(v.pi->m, "Hello, world");
    /* fill wrpals and mpals */
    v.wrpals[0] = (mmt_vec **) shared_malloc(v.pi->wr[0], v.pi->wr[0]->ncores * sizeof(mmt_vec *));
    v.wrpals[0][v.pi->wr[0]->trank] = &v;
    serialize_threads(v.pi->m);
    v.wrpals[1] = (mmt_vec **) shared_malloc(v.pi->wr[1], v.pi->wr[1]->ncores * sizeof(mmt_vec *));
    v.wrpals[1][v.pi->wr[1]->trank] = &v;
    serialize_threads(v.pi->m);
    v.mpals = (mmt_vec **) shared_malloc(v.pi->m, v.pi->m->ncores * sizeof(mmt_vec *));
    v.mpals[v.pi->m->trank] = &v;
    serialize_threads(v.pi->m);

}

mmt_vec::~mmt_vec()
{
    if (d == -1) return;
    pi_comm_ptr wr = pi->wr[d];
    serialize_threads(wr);
    if (rsbuf[0]) abase->free(rsbuf[0]);
    if (rsbuf[1]) abase->free(rsbuf[1]);
    if (siblings) {
        abase->free(v);
        shared_free(wr, siblings);
    } else {
        if (wr->trank == 0)
            abase->free(v);
    }
    shared_free(pi->wr[0], wrpals[0]);
    shared_free(pi->wr[1], wrpals[1]);
    shared_free(pi->m, mpals);
}
/* }}} */

/* This is a companion call to mmt_my_own_offset_in_items. Here, we do
 * not give the offset of data we're particularly the "owner" of, but
 * data that is owned by job jrank _and the same trank as ours_.
 */
size_t mmt_my_own_offset_in_items(mmt_vec const & v, unsigned int jrank)
{
    pi_comm_ptr wr = v.pi->wr[v.d];
    size_t eblock = (v.i1 - v.i0) /  wr->totalsize;
    int pos = jrank * wr->ncores + wr->trank;
    return pos * eblock;
}

/* my "own" offset is the added offset within my locally stored data area
 * which represents the data range I am the owner of. This data range
 * correspond to the index range v.i0 + offset to v.i0 + offset + size
 */
size_t mmt_my_own_offset_in_items(mmt_vec const & v)
{
    pi_comm_ptr wr = v.pi->wr[v.d];
    return mmt_my_own_offset_in_items(v, wr->jrank);
}

size_t mmt_my_own_offset_in_bytes(mmt_vec const & v)
{
    return v.abase->vec_elt_stride(mmt_my_own_offset_in_items(v));
}

arith_generic::elt * mmt_my_own_subvec(mmt_vec & v, unsigned int jrank)
{
    return v.abase->vec_subvec(v.v, mmt_my_own_offset_in_items(v, jrank));
}

arith_generic::elt * mmt_my_own_subvec(mmt_vec & v)
{
    return v.abase->vec_subvec(v.v, mmt_my_own_offset_in_items(v));
}

arith_generic::elt const * mmt_my_own_subvec(mmt_vec const & v)
{
    return v.abase->vec_subvec(v.v, mmt_my_own_offset_in_items(v));
}

size_t mmt_my_own_size_in_items(mmt_vec const & v)
{
    pi_comm_ptr wr = v.pi->wr[v.d];
    size_t eblock = (v.i1 - v.i0) /  wr->totalsize;
    return eblock;
}

size_t mmt_my_own_size_in_bytes(mmt_vec const & v)
{
    return v.abase->vec_elt_stride(mmt_my_own_size_in_items(v));
}

/* This copies **ONLY** the data we are supposed to own from v to w.
 * mmt_own_vec_set2 is slightly special, since in this case we might in
 * this way be picking data from vector areas owned by other blocks (or
 * writing there). Therefore we convey the info about which vector piece
 * we care about with the first argument z.
 *
 */
void mmt_own_vec_set2(mmt_vec const & z, mmt_vec & w, mmt_vec const & v)
{
    if (&v == &w) return;
    ASSERT_ALWAYS(z.d == v.d);
    ASSERT_ALWAYS(z.d == w.d);
    size_t off = mmt_my_own_offset_in_items(z);
    size_t sz = mmt_my_own_size_in_items(z);
    ASSERT_ALWAYS(sz == mmt_my_own_size_in_items(v));
    ASSERT_ALWAYS(sz == mmt_my_own_size_in_items(w));
    v.abase->vec_set(
            w.abase->vec_subvec(w.v, off),
            v.abase->vec_subvec(v.v, off),
            sz);
}
void mmt_own_vec_set(mmt_vec & w, mmt_vec const & v)
{
    ASSERT_ALWAYS(v.abase == w.abase);
    mmt_own_vec_set2(v, w, v);
    w.consistency = 1;
}
/*
void mmt_vec_swap(mmt_vec & w, mmt_vec & v)
{
    mmt_vec foo;
    memcpy(foo,v,sizeof(mmt_vec));
    memcpy(v,w,sizeof(mmt_vec));
    memcpy(w,foo,sizeof(mmt_vec));
}
*/


/* This is a no-op if storage is shared across threads in direction v.d.
 * If it's not, then each thread has its own copy. We take each
 * thread as authoritative with respect to one specific data range, and
 * all threads replicate this authoritative source to their own data
 *
 * The post-condition is that the data range is the same across all
 * threads in direction v.d. We don't care whether it's consistent across
 * nodes.  (In particular, this call is not an mmt_vec_broadcast, since
 * we only communicate between threads)
 */
void mmt_vec_share_across_threads(mmt_vec & v)
{
    if (mmt_vec_is_shared(v)) return;
    pi_comm_ptr wr = v.pi->wr[v.d];
    for(unsigned int t = 0 ; t < wr->ncores ; t++) {
        if (t == wr->trank)
            continue;
        mmt_vec const & w = mmt_vec_sibling(v, t);

        for(unsigned int j = 0 ; j < wr->njobs ; j++) {
            size_t off = mmt_my_own_offset_in_items(w, j);
            size_t sz  = mmt_my_own_size_in_items(w);

            v.abase->vec_set(
                    v.abase->vec_subvec(v.v, off),
                    v.abase->vec_subvec(w.v, off),
                    sz);
        }
    }
}

void mmt_full_vec_set(mmt_vec & w, mmt_vec const & v)
{
    /* DO **NOT** early-quit when v==w, because we might be calling this
     * with v and w being siblings, maybe equal for one of the threads.
     */
    // same remark as above
    // ASSERT_ALWAYS(v.abase == w.abase);
    ASSERT_ALWAYS(v.d == w.d);
    ASSERT_ALWAYS(mmt_my_own_size_in_items(v) == mmt_my_own_size_in_items(w));
    pi_comm_ptr wr = w.pi->wr[v.d];
    if (w.siblings) {
        /* w is not shared. Maybe v is, we don't really care */
        if (w.v != v.v) {
            v.abase->vec_set(w.v, v.v, v.i1 - v.i0);
        }
    } else if (mmt_vec_is_shared(v) || v.consistency == 2) {
        /* We make do with a single data source in v (either because
         * there _is_ only one, or because all are assumed consistent
         * anyway
         */
        if (w.v != v.v) {
            if (wr->trank == 0)
                w.abase->vec_set(w.v, v.v, v.i1 - v.i0);
            serialize_threads(wr);
        }
    } else {
        /* w is shared, but we need to read from multiple data sources.
         *
         * It is not entirely clear if we need to read only the data that
         * is relevant to the current node, or to all nodes. Given that
         * in the other cases, we copy the full range between i0 and i1,
         * it seems reasonable to attempt to do the same.
         */
        mmt_vec const & vt = mmt_vec_sibling(v, wr->trank);

        for(unsigned int j = 0 ; j < wr->njobs ; j++) {
            size_t off = mmt_my_own_offset_in_items(vt, j);
            size_t sz  = mmt_my_own_size_in_items(vt);

            w.abase->vec_set(
                    v.abase->vec_subvec(w.v, off),
                    v.abase->vec_subvec(vt.v, off),
                    sz);
        }
        serialize_threads(wr);
    }
    w.consistency = v.consistency;
}

void mmt_full_vec_set_zero(mmt_vec & v)
{
    if (v.siblings) {
        v.abase->vec_set_zero(v.v, v.i1 - v.i0);
    } else {
        serialize_threads(v.pi->wr[v.d]);
        if (v.pi->wr[v.d]->trank == 0)
            v.abase->vec_set_zero(v.v, v.i1 - v.i0);
    }
    v.consistency = 2;
    serialize_threads(v.pi->wr[v.d]);
}

void mmt_vec_set_basis_vector_at(mmt_vec & v, int k, unsigned int j)
{
    mmt_full_vec_set_zero(v);
    mmt_vec_add_basis_vector_at(v,k,j);
}

void mmt_vec_add_basis_vector_at(mmt_vec & v, int k, unsigned int j)
{
    if (v.i0 <= j && j < v.i1) {
        if (v.siblings) {
            v.abase->simd_set_ui_at(v.abase->vec_item(v.v, j - v.i0), k, 1);
        } else {
            serialize_threads(v.pi->wr[v.d]);
            if (v.pi->wr[v.d]->trank == 0)
                v.abase->simd_set_ui_at(v.abase->vec_item(v.v, j - v.i0), k, 1);
            serialize_threads(v.pi->wr[v.d]);
        }
    }
}

void mmt_vec_add_basis_vector(mmt_vec & v, unsigned int j)
{
    mmt_vec_add_basis_vector_at(v, 0, j);
}

void mmt_vec_set_basis_vector(mmt_vec & v, unsigned int j)
{
    mmt_vec_set_basis_vector_at(v, 0, j);
}

void mmt_vec_downgrade_consistency(mmt_vec & v)
{
    ASSERT_ALWAYS(v.consistency == 2);
    size_t erase[2][2];
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
        serialize_threads(v.pi->wr[v.d]);
    if (v.siblings) {
        erase[0][0] = 0;
        erase[0][1] = off;
        erase[1][0] = off + sz;
        erase[1][1] = v.i1 - v.i0;
    } else {
        /* There are no siblings, which means that this vector is shared
         * across all threads in this direction. Let only one thread do
         * the job.
         */
        if (v.pi->wr[v.d]->trank == 0) {
            erase[0][0] = 0;
            /* because we are rank 0, this is the minimal offset for this set
             * of threads */
            erase[0][1] = off;
            erase[1][0] = off + sz * v.pi->wr[v.d]->ncores;
            erase[1][1] = v.i1 - v.i0;
        } else {
            erase[0][0] = 0;
            erase[0][1] = 0;
            erase[1][0] = 0;
            erase[1][1] = 0;
        }
    }
    for(int i = 0 ; i < 2 ; i++) {
        if (erase[i][1] != erase[i][0]) {
            v.abase->vec_set_zero(
                    v.abase->vec_subvec(v.v, erase[i][0]),
                    erase[i][1] - erase[i][0]);
        }
    }
    v.consistency = 1;
    serialize_threads(v.pi->wr[v.d]);
}


#if 0
// looks utterly bogus, so...
/* On a shared vector which is assumed to be commonly known by all
 * nodes/threads, select only one portion to remain, and zero out the
 * rest (for later use by e.g.  matmul_top_mul_comm or other integrated
 * function).
 */
static void mmt_own_vec_clear_complement(matmul_top_data & mmt, int d)
{
    mmt_comm_ptr mdst = mmt.wr[d];
    pi_comm_ptr pidst = mmt.pi->wr[d];
    if (mdst->v.flags & THREAD_SHARED_VECTOR)
        serialize_threads(pidst);
    if (pidst->trank == 0 || !(mdst->v.flags & THREAD_SHARED_VECTOR)) {
        if (pidst->jrank == 0 && pidst->trank == 0) {
            /* ok, we keep the data */
        } else {
            mdst->v.abase->vec_set_zero(mdst->v.abase, mdst->v.v, mdst->i1 - mdst->i0);
        }
    }
}
#endif
void mmt_vec_clear_padding(mmt_vec & v, size_t unpadded, size_t padded)
{
    /* This can be applied no matter what the consistency argument says
     * */
    serialize(v.pi->m);
    if (unpadded >= padded) return;

    size_t s0 = unpadded >= v.i0 ? (unpadded - v.i0) : 0;
    size_t s1 = padded >= v.i0 ? (padded - v.i0) : 0;
    s0 = MIN(s0, v.i1 - v.i0);
    s1 = MIN(s1, v.i1 - v.i0);

    if (s1 - s0)
        v.abase->vec_set_zero(
                v.abase->vec_subvec(v.v, s0), s1-s0);

    serialize(v.pi->m);
}

mmt_vec & mmt_vec_sibling(mmt_vec & v, unsigned int i)
{
    if (v.siblings) {
        return *v.siblings[i];
    } else {
        return v;
    }
}

mmt_vec const & mmt_vec_sibling(mmt_vec const & v, unsigned int i)
{
    if (v.siblings) {
        return *v.siblings[i];
    } else {
        return v;
    }
}

/* {{{ generic interfaces for load/save */
/* {{{ load */
int mmt_vec_load(mmt_vec & v, const char * filename_pattern, unsigned int itemsondisk, unsigned int block_position)
{
    serialize(v.pi->m);
    int tcan_print = v.pi->m->trank == 0 && v.pi->m->jrank == 0;

    ASSERT_ALWAYS(strstr(filename_pattern, "%u-%u") != NULL);

    int char2 = v.abase->is_characteristic_two();
    int splitwidth = char2 ? 64 : 1;
    unsigned int Adisk_width = splitwidth;
    unsigned int Adisk_multiplex = v.abase->simd_groupsize() / Adisk_width;

    size_t sizeondisk = v.abase->vec_elt_stride(itemsondisk);
    arith_generic::elt * mychunk = mmt_my_own_subvec(v);
    size_t mysize = mmt_my_own_size_in_bytes(v);
    size_t bigstride = v.abase->vec_elt_stride(1);
    size_t smallstride = bigstride / Adisk_multiplex;

    int global_ok = 1;

    for(unsigned int b = 0 ; global_ok && b < Adisk_multiplex ; b++) {
        unsigned int b0 = block_position + b * Adisk_width;
        char * filename;
        int rc;
        rc = asprintf(&filename, filename_pattern, b0, b0 + splitwidth);
        ASSERT_ALWAYS(rc >= 0);
        double tt = -wct_seconds();
        if (tcan_print) {
            printf("Loading %s ...", filename);
            fflush(stdout);
        }
        pi_file_handle f;
        int ok = pi_file_open(f, v.pi, v.d, filename, "rb");
        /* "ok" is globally consistent after pi_file_open */
        if (!ok) {
            if (v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                fprintf(stderr, "ERROR: failed to load %s: %s\n", filename, strerror(errno));
            }
        } else {
            serialize(v.pi->m);
            ssize_t s = pi_file_read_chunk(f, mychunk, mysize, sizeondisk,
                    bigstride, b * smallstride, (b+1) * smallstride);
            int ok = s >= 0 && (size_t) s == sizeondisk / Adisk_multiplex;
            /* "ok" is globally consistent after pi_file_read_chunk */
            if (!ok) {
                if (v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                    fprintf(stderr, "ERROR: failed to load %s: short read, %s\n", filename, errno ? strerror(errno) : "no error reported via errno");
                }
            }
            /* Always reduce mod p after load */
            for(size_t i = 0 ; i < mmt_my_own_size_in_items(v) ; i++) {
                v.abase->reduce(v.abase->vec_item(mychunk, i));
            }
            v.consistency = ok;
            /* not clear it's useful, but well. */
            if (ok) mmt_vec_broadcast(v);
            serialize_threads(v.pi->m);
            pi_file_close(f);
        }
        free(filename);
        tt += wct_seconds();
        if (ok && tcan_print) {
            char buf[20], buf2[20];
            printf(" done [%s in %.2fs, %s/s]\n",
                    size_disp(sizeondisk / Adisk_multiplex, buf),
                    tt,
                    size_disp(sizeondisk / Adisk_multiplex / tt, buf2));
        }
        global_ok = global_ok && ok;
    }

    serialize_threads(v.pi->m);
    return global_ok;
}
/* }}} */
/* {{{ save */
int mmt_vec_save(mmt_vec & v, const char * filename_pattern, unsigned int itemsondisk, unsigned int block_position)
{
    serialize_threads(v.pi->m);
    int tcan_print = v.pi->m->trank == 0 && v.pi->m->jrank == 0;

    ASSERT_ALWAYS(strstr(filename_pattern, "%u-%u") != NULL);

    int char2 = v.abase->is_characteristic_two();
    int splitwidth = char2 ? 64 : 1;
    unsigned int Adisk_width = splitwidth;
    unsigned int Adisk_multiplex = v.abase->simd_groupsize() / Adisk_width;

    size_t sizeondisk = v.abase->vec_elt_stride(itemsondisk);
    arith_generic::elt * mychunk = mmt_my_own_subvec(v);
    size_t mysize = mmt_my_own_size_in_bytes(v);
    size_t bigstride = v.abase->vec_elt_stride(1);
    size_t smallstride = bigstride / Adisk_multiplex;

    int global_ok = 1;

    for(unsigned int b = 0 ; b < Adisk_multiplex ; b++) {
        unsigned int b0 = block_position + b * Adisk_width;
        char * filename;
        int rc = asprintf(&filename, filename_pattern, b0, b0 + splitwidth);
        ASSERT_ALWAYS(rc >= 0);
        char * tmpfilename;
        rc = asprintf(&tmpfilename, "%s.tmp", filename);
        ASSERT_ALWAYS(rc >= 0);
        double tt = -wct_seconds();
        if (tcan_print) {
            printf("Saving %s ...", filename);
            fflush(stdout);
        }
        pi_file_handle f;
        int ok = pi_file_open(f, v.pi, v.d, tmpfilename, "wb");
        /* "ok" is globally consistent after pi_file_open */
        if (!ok) {
            if (v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                fprintf(stderr, "WARNING: failed to save %s: %s\n", filename, strerror(errno));
                unlink(tmpfilename);    // just in case
            }
        } else {
            ASSERT_ALWAYS(v.consistency == 2);
            serialize_threads(v.pi->m);
            ssize_t s = pi_file_write_chunk(f, mychunk, mysize, sizeondisk,
                    bigstride, b * smallstride, (b+1) * smallstride);
            serialize_threads(v.pi->m);
            ok = s >= 0 && (size_t) s == sizeondisk / Adisk_multiplex;
            /* "ok" is globally consistent after pi_file_write_chunk */
            if (!ok && v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                fprintf(stderr, "ERROR: failed to save %s: short write, %s\n", filename, errno ? strerror(errno) : "no error reported via errno");
            }
            ok = pi_file_close(f);
            if (!ok && v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                fprintf(stderr, "ERROR: failed to save %s: failed fclose, %s\n", filename, errno ? strerror(errno) : "no error reported via errno");
            }
            if (v.pi->m->trank == 0 && v.pi->m->jrank == 0) {
                ok = rename(tmpfilename, filename) == 0;
                if (!ok) {
                    fprintf(stderr, "ERROR: failed to save %s: failed rename, %s\n", filename, errno ? strerror(errno) : "no error reported via errno");
                }
            }
            pi_bcast(&ok, 1, BWC_PI_INT, 0, 0, v.pi->m);
        }
        free(filename);
        free(tmpfilename);
        tt += wct_seconds();
        if (tcan_print) {
            char buf[20], buf2[20];
            printf(" done [%s in %.2fs, %s/s]\n",
                    size_disp(sizeondisk / Adisk_multiplex, buf),
                    tt,
                    size_disp(sizeondisk / Adisk_multiplex / tt, buf2));
        }
        global_ok = global_ok && ok;
    }

    serialize_threads(v.pi->m);
    return global_ok;
}
/* }}} */
/* }}} */

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

void mmt_vec_set_random_through_file(mmt_vec & v, const char * filename_pattern, unsigned int itemsondisk, gmp_randstate_t rstate, unsigned int block_position)
{
    /* FIXME: this generates the complete vector on rank 0, saves it, and
     * loads it again. But I'm a bit puzzled by the choice of saving a
     * number of items which is n0[d]. Seems to me that this is in fact
     * incorrect, we want n0[!d] here.
     */
    
    arith_generic * A = v.abase;
    parallelizing_info_ptr pi = v.pi;
    int tcan_print = v.pi->m->trank == 0 && v.pi->m->jrank == 0;

    int char2 = v.abase->is_characteristic_two();
    int splitwidth = char2 ? 64 : 1;
    unsigned int Adisk_width = splitwidth;
    unsigned int Adisk_multiplex = v.abase->simd_groupsize() / Adisk_width;

    ASSERT_ALWAYS(itemsondisk % Adisk_multiplex == 0);
    unsigned int loc_itemsondisk = itemsondisk / Adisk_multiplex;

    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for(unsigned int b = 0 ; b < Adisk_multiplex ; b++) {
            unsigned int b0 = block_position + b * Adisk_width;
            char * filename;
            int rc;
            rc = asprintf(&filename, filename_pattern, b0, b0 + splitwidth);
            ASSERT_ALWAYS(rc >= 0);

            /* we want to create v.n / Adisk_multiplex entries --
             * but we can't do that with access to just A. So we
             * generate slightly more, and rely on itemsondisk to do
             * the job of properly cutting the overflowing data.
             */

            size_t nitems = iceildiv(v.n, Adisk_multiplex);
            arith_generic::elt * y;
            y = A->alloc(nitems);
            A->vec_set_zero(y, nitems);
            A->vec_set_random(y, nitems, rstate);
            double tt = -wct_seconds();
            if (tcan_print) {
                printf("Creating random vector %s...", filename);
                fflush(stdout);
            }
            FILE * f = fopen(filename, "wb");
            ASSERT_ALWAYS(f);
            rc = fwrite(y, A->vec_elt_stride(1), loc_itemsondisk, f);
            ASSERT_ALWAYS(rc == (int) loc_itemsondisk);
            fclose(f);
            tt += wct_seconds();
            if (tcan_print) {
                size_t sizeondisk = A->vec_elt_stride(loc_itemsondisk);
                size_t fraction = sizeondisk / Adisk_multiplex;
                char buf[20], buf2[20];
                printf(" done [%s in %.2fs, %s/s]\n",
                        size_disp(fraction, buf),
                        tt,
                        size_disp(fraction / tt, buf2));
            }
            A->free(y);

            free(filename);
        }
    }
    int ok = mmt_vec_load(v, filename_pattern, itemsondisk, block_position);
    ASSERT_ALWAYS(ok);
}

unsigned long mmt_vec_hamming_weight(mmt_vec const & y) {
    ASSERT_ALWAYS(y.consistency == 2);
    unsigned long w = y.abase->vec_simd_hamming_weight(y.v, y.i1 - y.i0);
    /* all threads / cores in wiring wr[y.d] share the same data and
     * thus deduce the same count */
    pi_allreduce(NULL, &w, 1, BWC_PI_UNSIGNED_LONG, BWC_PI_SUM, y.pi->wr[!y.d]);
    return w;
}

/* this is inconsistent in the sense that it's balancing-dependent */
void mmt_vec_set_random_inconsistent(mmt_vec & v, gmp_randstate_t rstate)
{
    mmt_full_vec_set_zero(v);
    v.abase->vec_set_random(mmt_my_own_subvec(v), mmt_my_own_size_in_items(v), rstate);
    v.consistency=1;
    mmt_vec_allreduce(v);
}

void mmt_vec_set_x_indices(mmt_vec & y, uint32_t * gxvecs, int m, unsigned int nx)
{
    int shared = mmt_vec_is_shared(y);
    arith_generic * A = y.abase;
    mmt_full_vec_set_zero(y);
    if (!shared || y.pi->wr[y.d]->trank == 0) {
        for(int j = 0 ; j < m ; j++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                uint32_t i = gxvecs[j*nx+k];
                // set bit j of entry i to 1.
                if (i < y.i0 || i >= y.i1)
                    continue;
                {
                    arith_generic::elt & x = A->vec_item(y.v, i - y.i0);
                    A->simd_add_ui_at(x, j, 1);
                }
            }
        }
    }
    y.consistency=2;
    if (shared)
        serialize_threads(y.pi->wr[y.d]);
}

/* Set to the zero vector, except for the first n entries that are taken
 * from the vector v
 */
void mmt_vec_set_expanded_copy_of_local_data(mmt_vec & y, const arith_generic::elt * v, unsigned int n)
{
    int shared = !y.siblings;
    arith_generic * A = y.abase;
    mmt_full_vec_set_zero(y);
    if (!shared || y.pi->wr[y.d]->trank == 0) {
        for(unsigned int i = y.i0 ; i < y.i1 && i < n ; i++) {
            A->set( A->vec_item(y.v, i - y.i0), A->vec_item(v, i));
        }
    }
    y.consistency=2;
    if (shared)
        serialize_threads(y.pi->wr[y.d]);
}



