#include "cado.h"

#include <stdio.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "async.h"
#include "xdotprod.h"
#include "rolling.h"
#include "matops.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"

int exit_code = 0;

struct blstate {
    matmul_top_data mmt;
    mmt_vec y, my;
    mpfq_vbase A;
    gmp_randstate_t rstate;

    /* We'll need several intermediary n*n matrices. These will be
     * allocated everywhere (we set flags to 0, so as to avoid having
     * shared vectors) */
    mmt_vec V[3];
    mat64 L[3];
    bit_vector D[3];

    /* Here are the semantics of the data fields above.
     *
     * For iteration n, we let n0 = n%3, n1 = (n-1)%3, n2 = (n-2)%3.
     *
     * V[n0] is the V vector which is the input to iteration n. It does
     *       not consist of independent  vectors.
     * D[n0] is the extracted set of columns, computed from V[n0] (and
     *       also from D[n1]. Together, V[n0] and D[n0] allow to compute
     *       the basis of the n-th sub vector space W, although no
     *       explicit data is reserved to its storage.
     * L[n0] is computed from V[n0] and D[n0], and is some sort of local
     *       inverse (iirc, it's W_i^{inv} in most accounts).
     *
     * of course *[n1] and *[n2] are the same for the previous steps.
     *
     * In order to compute the vector for step n + 1, the data to be used
     * is V[*], L[*], and D[n0, n1]
     */
};

/* given a 64 by 64 symmetric matrix A (given as 64 uint64_t's), and an input
 * bitmap S given in the form of one uint64_t, compute B and T so that
 * the following hold:
 *
 *    (i)    B = T * B = B * T
 *    (ii)   B * A * T = T
 *    (iii)  rank(B) = rank(A)
 *    (iv)   (1-S)*T = 1-S
 *
 * where S is identified with the diagonal matrix whose entries are given by S.
 *
 * In other words, this inverts (ii) a maximal minor (iii) in the matrix A,
 * selecting in priority (iv) for defining this minor the row indices
 * which are zero in S.
 *
 * The matrix A is clobbered.
 *
 * This routine does some allocation on the stack. Speed is not critical.
 */
uint64_t extraction_step(uint64_t * B, uint64_t * A, uint64_t S)
{
    int order[64], reorder[64];
    uint64_t B0[64];
    uint64_t T = 0;
    /* convert to a priority list, in a "save trees" style.  */
    for(int o=64,z=0,i=64;i-->0;) order[S&(UINT64_C(1)<<i)?--o:z++]=i;
    for(int i = 0 ; i < 64 ; i++) B0[i] = UINT64_C(1)<<i;
    for(int i = 0 ; i < 64 ; i++) reorder[i]=-1;
    for(int j = 0 ; j < 64 ; j++) {
        int oj = order[j];
        uint64_t mj = UINT64_C(1)<<oj;
        int p = -1;
        for(int i = 0 ; i < 64 ; i++) {
            int oi = order[i];
            uint64_t mi = UINT64_C(1)<<oi;
            if (T & mi) continue;
            if (A[oi] & mj) {
                p = i;
                break;
            }
        }
        if (p < 0) continue;

        int op = order[p];
        /* Of course it's important to use indices op and oj here ! */
        reorder[op] = oj;
        uint64_t mp = UINT64_C(1) << op;
        /* We have a pivot, great. */
        ASSERT_ALWAYS(!(T & mp));
        T |= mp;
        /* add row op to all rows except itself */
        for(int i = 0 ; i < 64 ; i++) {
            if (i == p) continue;
            int oi = order[i];
            uint64_t x = ~-!(A[oi] & mj);
            B0[oi] ^= B0[op] & x;
            A[oi] ^= A[op] & x;
        }
    }
    /* Now at this point, we have some more work to do.
     *
     * A*T is now almost the identity matrix -- its square is diagonal
     * with zeros and ones, so A*T is just an involution. The reorder[]
     * array just has this involution.
     *
     * B is such that B*original_A = new_A.
     *
     * The matrix we want to return is new_A*T*b*T. We use reorder[] to
     * actually copy this into A.
     */
    memset(B, 0, sizeof(B0));
    for(int i = 0 ; i < 64 ; i++) {
        if (reorder[i] >= 0)
            B[reorder[i]]=B0[i]&T;
    }
    return T;
}


void blstate_init(struct blstate * bl, parallelizing_info_ptr pi, param_list_ptr pl)
{
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;
    /* Note that THREAD_SHARED_VECTOR can't work in a block Lanczos
     * context, since both ways are used for input and output.
     *
     * Well, at least it does not seem to be as easy as it is in the
     * block Wiedemann context. Would be neat to find a way, though,
     * since otherwise this represents quite a significant memory
     * footprint in the end.
     */

    gmp_randinit_default(bl->rstate);

    mpfq_vbase_oo_field_init_byfeatures(A,
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, bw->ys[1]-bw->ys[0],
            MPFQ_DONE);

    matmul_top_init(mmt, A, pi, pl, bw->dir);

    /* it's not really in the plans yet */
    ASSERT_ALWAYS(mmt->nmatrices == 1);

    mmt_vec_init(mmt,0,0, bl->y,   bw->dir, 0, mmt->n[bw->dir]);
    mmt_vec_init(mmt,0,0, bl->my, !bw->dir, 0, mmt->n[!bw->dir]);

    for(int i = 0 ; i < 3 ; i++) {
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_init(bl->D[i], bw->n);
        /* We need as well the two previous vectors. For these, distributed
         * storage will be ok. */
        mmt_vec_init(mmt,0,0, bl->V[i], bw->dir, 0, mmt->n[bw->dir]);
    }

}

void blstate_clear(struct blstate * bl)
{
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;

    serialize(mmt->pi->m);
    for(int i = 0 ; i < 3 ; i++) {
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_clear(bl->D[i]);
        /* We need as well the two previous vectors. For these, distributed
         * storage will be ok. */
        mmt_vec_clear(mmt, bl->V[i]);
    }
    mmt_vec_clear(mmt, bl->y);
    mmt_vec_clear(mmt, bl->my);
    matmul_top_clear(bl->mmt);
    A->oo_field_clear(bl->A);
    gmp_randclear(bl->rstate);
}

void blstate_set_start(struct blstate * bl)
{
    matmul_top_data_ptr mmt = bl->mmt;

    /* D = identity, L too, and V = 0 */
    for(int i = 0 ; i < 3 ; i++) {
        bit_vector_set(bl->D[i], 1);
        mat64_set_identity(bl->L[i]);
    }
    /* matmul_top_vec_init has already set V to zero */
    /* for bw->dir=0, mmt->n0[0] is the number of rows. */
    mmt_vec_set_random_through_file(bl->V[0], "blstart.0", mmt->n0[bw->dir], bl->rstate);
    mmt_own_vec_set(bl->y, bl->V[0]);
}

void blstate_load(struct blstate * bl, unsigned int iter)
{
    unsigned int i0 = iter % 3;
    unsigned int i1 = (iter+3-1) % 3;
    unsigned int i2 = (iter+3-2) % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    parallelizing_info_ptr pi = mmt->pi;

    char * filename;
    int rc = asprintf(&filename, "blstate.%u", iter);
    ASSERT_ALWAYS(rc >= 0);
    pi_file_handle f;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Loading %s...", filename); fflush(stdout); }

    pi_file_open(f, mmt->pi, bw->dir, filename, "rb");
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        bit_vector_read_from_stream(bl->D[i1], f->f);
        size_t rc;
        rc = fread(bl->L[i1], sizeof(mat64), 1, f->f);
        ASSERT_ALWAYS(rc == (size_t) 1);
        rc = fread(bl->L[i2], sizeof(mat64), 1, f->f);
        ASSERT_ALWAYS(rc == (size_t) 1);
    }
    mmt_vec_load_stream(f, bl->V[i0], mmt->n0[bw->dir]);
    mmt_vec_load_stream(f, bl->V[i1], mmt->n0[bw->dir]);
    mmt_vec_load_stream(f, bl->V[i2], mmt->n0[bw->dir]);

    pi_file_close(f);
    if (tcan_print) { printf("done\n"); fflush(stdout); }
    free(filename);
}

void blstate_save(struct blstate * bl, unsigned int iter)
{
    unsigned int i0 = iter % 3;
    unsigned int i1 = (iter+3-1) % 3;
    unsigned int i2 = (iter+3-2) % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    parallelizing_info_ptr pi = mmt->pi;

    char * filename;
    int rc = asprintf(&filename, "blstate.%u", iter);
    ASSERT_ALWAYS(rc >= 0);
    pi_file_handle f;
    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Saving %s...", filename); fflush(stdout); }

    pi_file_open(f, mmt->pi, bw->dir, filename, "wb");
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        bit_vector_write_to_stream(bl->D[i1], f->f);
        size_t rc;
        rc = fwrite(bl->L[i1], sizeof(mat64), 1, f->f);
        ASSERT_ALWAYS(rc == (size_t) 1);
    }
    mmt_vec_save_stream(f, bl->V[i0], mmt->n0[bw->dir]);
    mmt_vec_save_stream(f, bl->V[i1], mmt->n0[bw->dir]);
    mmt_vec_save_stream(f, bl->V[i2], mmt->n0[bw->dir]);

    pi_file_close(f);
    if (tcan_print) { printf("done\n"); fflush(stdout); }
    free(filename);
}

/* Given a *BINARY* vector block of *EXACTLY* 64 vectors (that is, we
 * force bw->n == 64), compute a 64*64 matrix such that m * v is in row
 * reduced echelon form. Also return the rank.
 *
 * Vector v is modified by this process (and put in RREF).
 *
 * This is a collective operation.
 */
int mmt_vec_echelon(mat64_ptr m, mmt_vec_ptr v0)
{
    mat64_set_identity(m);
    uint64_t * v = mmt_my_own_subvec(v0);
    size_t eblock = mmt_my_own_size_in_items(v0);
    /* This is the total number of non-zero coordinates of the vector v */
    size_t n = v0->n;
    /* In all what follows, we'll talk about v being a 64*n matrix, with
     * [v[i]&1] being "the first row", and so on.  */
    uint64_t usedrows = 0;
    int rank = 0;
    for(int i = 0 ; i < 64 ; i++) {
        uint64_t mi = UINT64_C(1) << i;
        /* Find the earliest column which has non-zero in the i-th row */
        unsigned int j;
        for(j = 0 ; j < eblock ; j++) {
            if (v[j] & mi) break;
        }
        if (j == eblock) j = n;
        else j += v0->i0 + mmt_my_own_offset_in_items(v0);
        unsigned int jmin;
        pi_allreduce(&j, &jmin, 1, BWC_PI_UNSIGNED, BWC_PI_MIN, v0->pi->m);
        if (jmin == n) {
            /* zero row */
            continue;
        }
        usedrows |= mi;
        rank++;
        /* zero out the other coefficients of this column. This means
         * that we need to collect the control data from the thread which
         * owns this column, and then act accordingly */
        uint64_t control = 0;
        if (jmin == j) {
            control = v[j - (v0->i0 + mmt_my_own_offset_in_items(v0))];
            ASSERT_ALWAYS(control & mi);
            control ^= mi;
        }
        /* TODO: once we require mpi-3.0, use MPI_UINT64_T instead */
        ASSERT_ALWAYS(sizeof(unsigned long long) == sizeof(uint64_t));
        pi_allreduce(NULL, &control, 1, BWC_PI_UNSIGNED_LONG_LONG, BWC_PI_MAX, v0->pi->m);
        /* add row i to all rows where we had a coeff in column j */
        /* we'll do that for all coefficients in the block, but on m this
         * is just one single operation */
        /* in this notation (as elsewhere with mat64 data), m[i] is
         * understood as row i of m (pay attention to the fact that this
         * differs a bit from what happens with y, where y[0] is in fact
         * the first coordinate of the row vector block, which can be
         * understood as the first column.
         */
        addmul_To64_o64(m, control, m[i]);
        for(unsigned int k = 0 ; k < eblock ; k++) {
            v[k] ^= control & -!!(v[k] & mi);
        }
    }
    /* put in RREF -- well, almost, since the only thing we do here is
     * that non-zero rows are before zero rows. */
    mat64 Z, N;
    int nZ = 0, nN = 0;
    memset(Z, 0, sizeof(mat64));
    memset(N, 0, sizeof(mat64));
    for(int i = 0 ; i < 64 ; i++) {
        uint64_t mi = UINT64_C(1) << i;
        if (usedrows & mi) {
            N[nN++] = m[i];
        } else {
            Z[nZ++] = m[i];
        }
    }
    ASSERT_ALWAYS(nN == rank);
    ASSERT_ALWAYS(nZ == 64 - rank);
    memcpy(m, N, rank * sizeof(uint64_t));
    memcpy(m + rank, Z, (64 - rank) * sizeof(uint64_t));

    return rank;
}



/* This saves first the vector V, and then the vector V*A. By
 * construction, the two are orthogonal. It is the responsibility of the
 * external program to check exactly where we have kernel vectors in V,
 * using the comparison with V*A.
 */
void blstate_save_result(struct blstate * bl, unsigned int iter)
{
    mat64 m0, m1, m2;
    int r;
    unsigned int i0 = iter % 3;
    matmul_top_data_ptr mmt = bl->mmt;
    parallelizing_info_ptr pi = mmt->pi;

    /* bw->dir=0: mmt->n0[bw->dir] = number of rows */
    /* bw->dir=1: mmt->n0[bw->dir] = number of columns */

    pi_file_handle f;
    const char * filename = "bl-auxiliary.bin";
    pi_file_open(f, mmt->pi, bw->dir, filename, "wb");

    int tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { printf("Saving %s...\n", filename); fflush(stdout); }

    mmt_full_vec_set(bl->y, bl->V[i0]);

    /* save raw V because it's conceivably useful */
    ASSERT_ALWAYS(bl->y->n == mmt->n[bw->dir]);
    mmt_vec_save_stream(f, bl->y, mmt->n0[bw->dir]);

    mmt_vec_twist(mmt, bl->y);
    matmul_top_mul_cpu(mmt, 0, bl->y->d, bl->my, bl->y);
    mmt_vec_allreduce(bl->my);
    mmt_vec_untwist(mmt, bl->my);

    /* save V*M as well because it's conceivably useful
     * XXX TODO: save in the other direction !
     *
     * pi_file_open should not have the "inner" flag. This should be a
     * property of the _write and _read calls. But before we do that, we
     * should clean up the mess done in mksol & gather.
     */
    ASSERT_ALWAYS(bl->my->n == mmt->n[!bw->dir]);
    mmt_vec_save_stream(f, bl->my, mmt->n0[!bw->dir]);

    /* Do some rank calculations */
    /* m0 is just temp stuff. */
    mmt_full_vec_set(bl->y, bl->V[i0]);
    r = mmt_vec_echelon(m0, bl->y);
    if (tcan_print) printf("\trank(V) == %d\n", r);
    /* m1 will be largest rank matrix such that m1*V*M == 0 */
    r = mmt_vec_echelon(m1,bl-> my);
    if (tcan_print) printf("\trank(V*M) == %d\n", r);

    /* good, so now let's look for real nullspace elements. Since
     * we've put V*M in RREF, we know what transformation of V
     * creates zeros in V*M. First, take out the combinations which
     * yield independent rows in V*M -- these are uninteresting */
    for(int i = 0 ; i < r ; i++) {
        m1[i] = 0;
    }
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        size_t rc = fwrite(m1, sizeof(mat64), 1, f->f);
        ASSERT_ALWAYS(rc == (size_t) 1);
    }
    /* Now set y = m1 * V */
    mmt_full_vec_set(bl->y, bl->V[i0]);
    mul_N64_T6464(mmt_my_own_subvec(bl->y),
            mmt_my_own_subvec(bl->y),
            m1,
            mmt_my_own_size_in_items(bl->y));
    bl->y->consistency = 1;
    /* Compute m2 such that m2 * m1 * V is in RREF. We discard the
     * combinations which lead to zero, so that m2*m1 should have rank
     * precisely the rank of the nullspace */
    r = mmt_vec_echelon(m2, bl->y);
    if (tcan_print) printf("\trank(V cap nullspace(M)) == %d\n", r);
    /* Now we have the really interesting basis of the nullspace */
    for(int i = r ; i < 64 ; i++) {
        m2[i] = 0;
    }
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        size_t rc = fwrite(m2, sizeof(mat64), 1, f->f);
        ASSERT_ALWAYS(rc == (size_t) 1);
    }
    serialize(mmt->pi->m);
    /* Now apply m2*m1 to v, for real */
    mmt_full_vec_set(bl->y, bl->V[i0]);
    mul_N64_T6464(mmt_my_own_subvec(bl->y),
            mmt_my_own_subvec(bl->y),
            m1,
            mmt_my_own_size_in_items(bl->y));
    mul_N64_T6464(mmt_my_own_subvec(bl->y),
            mmt_my_own_subvec(bl->y),
            m2,
            mmt_my_own_size_in_items(bl->y));


    /* Now save the reduced kernel basis */
    ASSERT_ALWAYS(bl->y->n == mmt->n[bw->dir]);
    mmt_vec_save_stream(f, bl->y, mmt->n0[bw->dir]);

    pi_file_close(f);
    if (tcan_print) { printf("Saving %s...done\n", filename); fflush(stdout); }

    mmt_vec_save(bl->y, "blsolution.0", mmt->n0[bw->dir]);
}


void * bl_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int tcan_print = bw->can_print && pi->m->trank == 0;
    struct timing_data timing[1];

    struct blstate bl[1];
    /* shortcuts */
    matmul_top_data_ptr mmt = bl->mmt;
    mpfq_vbase_ptr A = bl->A;

    /* so many features we do not support ! */
    ASSERT_ALWAYS(bw->m == bw->n);
    ASSERT_ALWAYS(bw->ys[0] == 0);
    ASSERT_ALWAYS(bw->ys[1] == bw->n);

    /* Don't think we stand any chance with interleaving with block
     * Lanczos... */
    ASSERT_ALWAYS(!pi->interleaved);

    // int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    block_control_signals();

    blstate_init(bl, pi, pl);

    mpfq_vbase_tmpl AxA;
    mpfq_vbase_oo_init_templates(AxA, A, A);

    size_t nelts_for_nnmat = bw->n * (bw->n / bl->A->groupsize(bl->A));

    serialize(pi->m);

    /* note that we don't know how to do checking for BL. Sure, we can
     * check for mutual orthogonality of vector blocks at the
     * checkpoints, but that is rather lame.
     */

    matmul_top_comm_bench(mmt, bw->dir);

    if (bw->start == 0) {
        blstate_set_start(bl);
        blstate_save(bl, 0);
    } else {
        blstate_load(bl, bw->start);
    }

    int auto_end = 0;

    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        length = mmt->n[bw->dir] / (bw->n - 0.76) + 10;
        /* allow some deviation */
        length += 2*integer_sqrt(length);
        /* Because bw is a global variable, we protect its use */
        if (serialize_threads(pi->m)) {
            bw->end = length;
        }
        auto_end = length;
        serialize(pi->m);
    }

    if (tcan_print) {
        printf ("Target iteration is %u ; going to %u\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }


    timing_init(timing, 4 * mmt->nmatrices, bw->start, bw->interval * iceildiv(bw->end, bw->interval));
    for(int i = 0 ; i < mmt->nmatrices ; i++) {
        timing_set_timer_name(timing, 4*i, "CPU%d", i);
        timing_set_timer_items(timing, 4*i, mmt->matrices[i]->mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%d", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%d", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%d", i);
    }

    void * vav;
    void * vaav;

    cheating_vec_init(A, &vav, nelts_for_nnmat);
    cheating_vec_init(A, &vaav, nelts_for_nnmat);

    /* TODO: Put that in the state file. */
    int sum_Ni = 0;

    for(int s = bw->start ; s < bw->end ; s += bw->interval ) {
        /* XXX For BL, how much do we care about twisting ? */
        /* The BW iteration is such that for a stored matrix M', we do
         * the iteration with P*M', for some matrix P. More precisely, we
         * do the following, where w is a row vector.
         *      Transpose(w) <- M' * v
         *      v <- P * Transpose(w) = P * M' * v
         *
         * (in a nutshell, P transforms row index N*x+y to N*y+x -- but
         * it's slightly more complicated than this. This is because we
         * do communications with all-gather and reduce-scatter).
         *
         * So in order to have an equivalent of multiplying by M, we want
         * P*M' = S^-1 * M * S, i.e. M' = (SP)^-1 M S.
         *
         * Now if in fact, we are doing the iteration followed by its
         * transpose, this may end up somewhat different.
         *
         * Note though that it's a bit incorrect to think like this. In
         * the BL iteration, this matrix P has little bearing anyway.
         * What we really do is
         *      Transpose(w) <- M' * v, followed by
         *      Transpose(v) <- w * M',
         * so that v <- Transpose(M')*M' * v, no matter what P is:
         * All communications are done with allreduce, so there is no
         * inherent shuffling in the parallelized matrix multiplication.
         *
         * When we have a stored matrix M', we do multiply by M'. Now let
         * us say that we have run the balancing process just like it
         * would have been for BW.  So we have M' = (SP)^-1 M S. If we
         * multiply by M'^T afterwards, we end up multiplying by S^-1 *
         * M^T * M * S. So twisting a vector means computing S^-1*V --
         * that should be essentially the same as for bwc.
         *
         * [description above is for solving M*w=0, while for BL we are
         * really in the factoring case where we care about w*M=0. So
         * then, twisting would mean v --> v*S*P, which is what
         * matmul_top_twist_vector does for bw->dir==0 and nullspace=left]
         *
         * Bottom line: it's harmless to do the balancing the very same
         * way we do with BW.
         */

        for(int k = 0 ; k < 3 ; k++) {
            mmt_vec_twist(mmt, bl->V[k]);
        }
        mmt_own_vec_set(bl->y, bl->V[s % 3]);

        /* I *think* that this is necessary */
        mmt_vec_broadcast(bl->y);

        /* FIXME: for BL, we use 8 timers while only 4 are defined in the
         * structure.
         */
        serialize(pi->m);
        int i, i0, i1, i2;

        for(i = 0 ; i < bw->interval ; i++) {
            i0 = (s + i) % 3;
            i1 = (s + i + 3 - 1) % 3;
            i2 = (s + i + 3 - 2) % 3;

            mmt_full_vec_set_zero(bl->my);

            mmt_vec_ptr yy[2] = {bl->y, bl->my};

            for(int d = 0 ; d < 2 ; d++) {
                /* The first part of this loop must be guaranteed to be free
                 * of any mpi calls */
                /* at this point, timer is [0] (CPU) */

                {
                    matmul_top_mul_cpu(mmt, 0, yy[d]->d, yy[!d], yy[d]);
                    timing_next_timer(timing);  /* now timer is [1] (cpu-wait) */
                }
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [2] (COMM) */
                mmt_vec_allreduce(yy[!d]);

                timing_next_timer(timing);  /* now timer is [3] (comm-wait) */
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [0] (CPU) */
            }

            /* We have now:
             *  V_{n}           bl->V[0]
             *  V_{n-1}         bl->V[1]
             *  V_{n-2}         bl->V[2]
             *  A * V_n         mcol->v
             */

            AxA->dotprod(A, A, vav,
                    mmt_my_own_subvec(bl->V[i0]),
                    mmt_my_own_subvec(bl->y),
                    mmt_my_own_size_in_items(bl->y));

            pi_allreduce(NULL, vav,
                    nelts_for_nnmat, bl->mmt->pitype, BWC_PI_SUM, pi->m);

            AxA->dotprod(A, A, vaav,
                    mmt_my_own_subvec(bl->y),
                    mmt_my_own_subvec(bl->y),
                    mmt_my_own_size_in_items(bl->y));

            pi_allreduce(NULL, vaav, nelts_for_nnmat,
                    bl->mmt->pitype, BWC_PI_SUM, pi->m);


            ASSERT_ALWAYS(bl->D[i0]->n == 64);

            {
                size_t eblock = mmt_my_own_size_in_items(bl->y);
                ASSERT_ALWAYS(bl->y->abase->vec_elt_stride(bl->y->abase, 1) == sizeof(uint64_t));

                // Here are the operations we will now perform
                //
                // X := (IN/*-Dn0*/)*Vn + Dn0*VA;
                // C0 := ((IN/*-Dn0*/)*vav + Dn0*vaav) * Ln0 * Vn;
                // C1 := Dn0 * vav * Ln1 * Vn1;
                // C2 := Dn0 * vav * (IN-Dn1) * Ln2 * Vn2;
                //
                // all are performed locally.

                // we set up X in mcol->v

                uint64_t * V0 = mmt_my_own_subvec(bl->V[i0]);
                uint64_t * V1 = mmt_my_own_subvec(bl->V[i1]);
                uint64_t * V2 = mmt_my_own_subvec(bl->V[i2]);
                uint64_t * VA  = mmt_my_own_subvec(bl->y);
                uint64_t * X   = mmt_my_own_subvec(bl->y);
                uint64_t D0;
                uint64_t D1 = bl->D[i1]->p[0];
                mat64_ptr mvav = (mat64_ptr) vav;
                mat64_ptr mvaav = (mat64_ptr) vaav;
                mat64_ptr mL0 = bl->L[i0];
                mat64_ptr mL1 = bl->L[i1];
                mat64_ptr mL2 = bl->L[i2];
                mat64 m0, m1, m2, t;

                /* We need to save vav for use a wee bit later in this loop. */
                memcpy(t, mvav, sizeof(mat64));

                D0 = bl->D[i0]->p[0] = extraction_step(mL0, t, D1);

                int Ni = bit_vector_popcount(bl->D[i0]);
                sum_Ni += Ni;
                // int defect = bw->n - Ni;
                // printf("step %d, dim=%d\n", s+i, Ni);
                if (Ni == 0) {
                    break;
                }

            
                for(size_t i = 0 ; i < eblock ; i++) {
                    X[i] = V0[i] ^ (VA[i] & D0);
                }

                for(int i = 0 ; i < 64 ; i++) t[i] = mvav[i] ^ (mvaav[i] & D0);
                mul_6464_6464(m0, mL0, t);
                for(int i = 0 ; i < 64 ; i++) mvav[i] &= D0;
                mul_6464_6464(m1, mL1, mvav);
                for(int i = 0 ; i < 64 ; i++) mL2[i] &= ~D1;
                mul_6464_6464(m2, mL2, mvav);

                addmul_N64_6464(X, V0, m0, eblock);
                addmul_N64_6464(X, V1, m1, eblock);
                addmul_N64_6464(X, V2, m2, eblock);
                bl->y->consistency = 1;

                /* So. We need to get ready for the next step, don't we ? */
                mmt_vec_broadcast(bl->y);
                mmt_own_vec_set(bl->V[i2], bl->y);
                bl->V[i2]->consistency = 1;
            }
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        for(int k = 0 ; k < 3 ; k++) {
            if (bl->V[k]->consistency < 2)
                mmt_vec_broadcast(bl->V[k]);
            mmt_vec_untwist(mmt, bl->V[k]);
        }

        serialize(pi->m);
        if (i < bw->interval) {
            if (tcan_print) printf("Finished at iteration %d, sum_dim=%d=N*(%u-%.3f)\n", s+i, sum_Ni, bw->n, bw->n-(double)sum_Ni/(s+i));

            blstate_save_result(bl, s+i);
            /* We need to cheat somewhat */
            if (serialize_threads(pi->m)) {
                bw->end = s + i;
            }
            serialize_threads(pi->m);
            timing->end_mark = s + i;
            break;
        }


        blstate_save(bl, s+bw->interval);
        serialize(pi->m);

        if (tcan_print) printf("N=%d ; sum_dim=%d=N*(%u-%.3f)\n", s+i, sum_Ni, bw->n, bw->n-(double)sum_Ni/(s+i));
        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "blocklanczos");
    }

    // we can't do as we do with BW, since the process really stops in
    // the end, continuing the last block to the checkpoint mark does not
    // make sense.
    timing_final_tally(pi, timing, tcan_print, "blocklanczos");

    cheating_vec_clear(A, &vav, nelts_for_nnmat);
    cheating_vec_clear(A, &vaav, nelts_for_nnmat);

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (auto_end && bw->end == auto_end) {
        if (tcan_print) {
            printf("FAILED blocklanczos (no collapse found for inner product).\n");
        }
        if (serialize_threads(pi->m)) {
            exit_code = 1;
        }
        serialize_threads(pi->m);
    } else {
        if (tcan_print) {
            printf("Done blocklanczos.\n");
        }
    }
    serialize(pi->m);

    blstate_clear(bl);

    timing_clear(timing);

    return NULL;
}

/* The unit test for extraction_step, which is a static function here,
 * actually includes the source file. We just want to make sure we don't
 * expose main()
 */
#ifndef BL_TESTING
int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }
    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(bl_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return exit_code;
}
#endif
