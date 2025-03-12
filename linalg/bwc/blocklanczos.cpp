#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstdint>              // for uint64_t, UINT64_C
#include <cstring>              // for memcpy, memset

#include <memory>

#include <gmp.h>                 // for gmp_randclear, gmp_randinit_default
#include "fmt/format.h"
#include "fmt/base.h"

#include "gmp_aux.h"
#include "matmul.hpp"              // for matmul_public_s
#include "parallelizing_info.hpp"
#include "matmul_top.hpp"
#include "select_mpi.h"
#include "params.h"
#include "misc.h"
#include "bw-common.h"
#include "async.hpp"
#include "bblas.hpp"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "bit_vector.h"
#include "macros.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "matmul_top_comm.hpp"
#include "blocklanczos.hpp"
#include "blocklanczos_extraction.hpp"
#include "bwc_filenames.hpp"
#include "utils_cxx.hpp"


blstate::blstate(parallelizing_info_ptr pi, cxx_param_list & pl)
    : A(arith_generic::instance(bw->p, bw->ys[1]-bw->ys[0]))
    , mmt(A.get(), pi, pl, bw->dir)
    , AxA(arith_cross_generic::instance(A.get(), A.get()))
    , y(mmt, nullptr, nullptr,  bw->dir, 0, mmt.n[bw->dir])
    , my(mmt, nullptr, nullptr, !bw->dir, 0, mmt.n[!bw->dir])
    , V {
        mmt_vec(mmt, nullptr, nullptr, bw->dir, 0, mmt.n[bw->dir]),
        mmt_vec(mmt, nullptr, nullptr, bw->dir, 0, mmt.n[bw->dir]),
        mmt_vec(mmt, nullptr, nullptr, bw->dir, 0, mmt.n[bw->dir]),
    }
{
    /* Note that THREAD_SHARED_VECTOR can't work in a block Lanczos
     * context, since both ways are used for input and output.
     *
     * Well, at least it does not seem to be as easy as it is in the
     * block Wiedemann context. Would be neat to find a way, though,
     * since otherwise this represents quite a significant memory
     * footprint in the end.
     */

    /* it's not really in the plans yet */
    ASSERT_ALWAYS(mmt.matrices.size() == 1);

    for(int i = 0 ; i < 3 ; i++) {
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_init(D[i], bw->n);
        /* We need as well the two previous vectors. For these, distributed
         * storage will be ok. */
    }
}

blstate::~blstate()
{
    serialize(mmt.pi->m);
    for(int i = 0 ; i < 3 ; i++) {
        /* We also need D_n, D_{n-1}, D_{n-2}. Those are in fact bitmaps.
         * Not clear that the bitmap type is really the one we want, though. */
        bit_vector_clear(D[i]);
    }
}

void blstate::set_start()
{
    /* D = identity, L too, and V = 0 */
    for(int i = 0 ; i < 3 ; i++) {
        bit_vector_set(D[i], 1);
        L[i] = 1;
    }
    /* matmul_top_vec_init has already set V to zero */
    /* for bw->dir=0, mmt.n0[0] is the number of rows. */
    mmt_vec_set_random_through_file(V[0], "blstart.0.V{}-{}.i0", mmt.n0[bw->dir], rstate, 0);
    mmt_own_vec_set(y, V[0]);
}

void blstate::load( unsigned int iter)
{
    unsigned int const i0 = iter % 3;
    unsigned int const i1 = (iter+3-1) % 3;
    unsigned int const i2 = (iter+3-2) % 3;
    parallelizing_info_ptr pi = mmt.pi;

    auto filename_base = fmt::format("blstate.{}.", iter);
    int const tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { fmt::print("Loading {}* ...", filename_base); fflush(stdout); }

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        auto tmp = filename_base + "control";
        auto f = fopen_helper(tmp, "rb");
        bit_vector_read_from_stream(D[i1], f.get());
        size_t rc;
        rc = fread(L[i1].data(), sizeof(mat64), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
        rc = fread(L[i2].data(), sizeof(mat64), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
    }

    std::string tmp;
    int rc;

    tmp = filename_base + bwc_V_file::pattern(iter);
    rc = mmt_vec_load(V[i0], tmp, mmt.n0[bw->dir], 0);
    ASSERT_ALWAYS(rc >= 0);

    tmp = filename_base + bwc_V_file::pattern(iter + 1);
    rc = mmt_vec_load(V[i1], tmp, mmt.n0[bw->dir], 0);
    ASSERT_ALWAYS(rc >= 0);

    tmp = filename_base + bwc_V_file::pattern(iter + 2);
    ASSERT_ALWAYS(rc >= 0);
    rc = mmt_vec_load(V[i2], tmp, mmt.n0[bw->dir], 0);
    ASSERT_ALWAYS(rc >= 0);

    if (tcan_print) { fmt::print("done\n"); fflush(stdout); }
}

void blstate::save(unsigned int iter)
{
    unsigned int const i0 = iter % 3;
    unsigned int const i1 = (iter+3-1) % 3;
    unsigned int const i2 = (iter+3-2) % 3;
    parallelizing_info_ptr pi = mmt.pi;

    auto filename_base = fmt::format("blstate.{}", iter);
    int const tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { fmt::print("Saving {}.* ...", filename_base); fflush(stdout); }

    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        auto tmp = filename_base + ".control";
        auto f = fopen_helper(tmp, "wb");
        bit_vector_write_to_stream(D[i1], f.get());
        const size_t rc = fwrite(L[i1].data(), sizeof(mat64), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
    }

    std::string tmp;

    tmp = filename_base + bwc_V_file::pattern(iter);
    mmt_vec_save(V[i0], tmp, mmt.n0[bw->dir], 0);

    tmp = filename_base + bwc_V_file::pattern(iter + 1);
    mmt_vec_save(V[i1], tmp, mmt.n0[bw->dir], 0);

    tmp = filename_base + bwc_V_file::pattern(iter + 2);
    mmt_vec_save(V[i2], tmp, mmt.n0[bw->dir], 0);

    if (tcan_print) { fmt::print("done\n"); fflush(stdout); }
}

/* Given a *BINARY* vector block of *EXACTLY* 64 vectors (that is, we
 * force bw->n == 64), compute a 64*64 matrix such that m * v is in row
 * reduced echelon form. Also return the rank.
 *
 * Vector v is modified by this process (and put in RREF).
 *
 * This is a collective operation.
 */
static int mmt_vec_echelon(mat64 & m, mmt_vec const & v0)
{
    m = 1;
    auto * v = (uint64_t *) mmt_my_own_subvec(v0);
    size_t const eblock = mmt_my_own_size_in_items(v0);
    /* This is the total number of non-zero coordinates of the vector v */
    size_t const n = v0.n;
    /* In all what follows, we'll talk about v being a 64*n matrix, with
     * [v[i]&1] being "the first row", and so on.  */
    uint64_t usedrows = 0;
    int rank = 0;
    for(int i = 0 ; i < 64 ; i++) {
        uint64_t const mi = UINT64_C(1) << i;
        /* Find the earliest column which has non-zero in the i-th row */
        unsigned int j;
        for(j = 0 ; j < eblock ; j++) {
            if (v[j] & mi) break;
        }
        if (j == eblock) j = n;
        else j += v0.i0 + mmt_my_own_offset_in_items(v0);
        unsigned int jmin;
        pi_allreduce(&j, &jmin, 1, BWC_PI_UNSIGNED, BWC_PI_MIN, v0.pi->m);
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
            control = v[j - (v0.i0 + mmt_my_own_offset_in_items(v0))];
            ASSERT_ALWAYS(control & mi);
            control ^= mi;
        }
        /* TODO: once we require mpi-3.0, use MPI_UINT64_T instead */
        ASSERT_ALWAYS(sizeof(unsigned long long) == sizeof(uint64_t));
        pi_allreduce(nullptr, &control, 1, BWC_PI_UNSIGNED_LONG_LONG, BWC_PI_MAX, v0.pi->m);
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
    Z = 0;
    N = 0;
    for(int i = 0 ; i < 64 ; i++) {
        uint64_t const mi = UINT64_C(1) << i;
        if (usedrows & mi) {
            N[nN++] = m[i];
        } else {
            Z[nZ++] = m[i];
        }
    }
    ASSERT_ALWAYS(nN == rank);
    ASSERT_ALWAYS(nZ == 64 - rank);
    memcpy(m.data(), N.data(), rank * sizeof(uint64_t));
    memcpy(m.data() + rank, Z.data(), (64 - rank) * sizeof(uint64_t));

    return rank;
}



/* This saves first the vector V, and then the vector V*A. By
 * construction, the two are orthogonal. It is the responsibility of the
 * external program to check exactly where we have kernel vectors in V,
 * using the comparison with V*A.
 */
void blstate::save_result( unsigned int iter)
{
    mat64 m0, m1, m2;
    int r;
    unsigned int const i0 = iter % 3;
    parallelizing_info_ptr pi = mmt.pi;

    /* bw->dir=0: mmt.n0[bw->dir] = number of rows */
    /* bw->dir=1: mmt.n0[bw->dir] = number of columns */

    std::string filename_base = "blaux.";


    int const tcan_print = bw->can_print && pi->m->trank == 0;
    if (tcan_print) { fmt::print("Saving {}* ...\n", filename_base); fflush(stdout); }

    mmt_full_vec_set(y, V[i0]);

    /* save raw V because it's conceivably useful */
    ASSERT_ALWAYS(y.n == mmt.n[bw->dir]);

    mmt_vec_save(y, filename_base + "Y{}-{}", mmt.n0[bw->dir], 0);

    mmt_vec_twist(mmt, y);
    matmul_top_mul_cpu(mmt, 0, y.d, my, y);
    mmt_vec_allreduce(my);
    mmt_vec_untwist(mmt, my);

    /* save V*M as well because it's conceivably useful
     * XXX TODO: save in the other direction !
     *
     * pi_file_open should not have the "inner" flag. This should be a
     * property of the _write and _read calls. But before we do that, we
     * should clean up the mess done in mksol & gather.
     */
    ASSERT_ALWAYS(my.n == mmt.n[!bw->dir]);
    ASSERT_ALWAYS(y.n == mmt.n[bw->dir]);
    mmt_vec_save(my, filename_base + "MY{}-{}", mmt.n0[!bw->dir], 0);

    /* Do some rank calculations */
    /* m0 is just temp stuff. */
    mmt_full_vec_set(y, V[i0]);
    r = mmt_vec_echelon(m0, y);
    if (tcan_print) fmt::print("\trank(V) == {}\n", r);
    /* m1 will be largest rank matrix such that m1*V*M == 0 */
    r = mmt_vec_echelon(m1,  my);
    if (tcan_print) fmt::print("\trank(V*M) == {}\n", r);

    /* good, so now let's look for real nullspace elements. Since
     * we've put V*M in RREF, we know what transformation of V
     * creates zeros in V*M. First, take out the combinations which
     * yield independent rows in V*M -- these are uninteresting */
    for(int i = 0 ; i < r ; i++) {
        m1[i] = 0;
    }
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        auto f = fopen_helper(filename_base + "M1", "wb");
        size_t const rc = fwrite(m1.data(), sizeof(mat64), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
    }
    /* Now set y = m1 * V */
    mmt_full_vec_set(y, V[i0]);
    mul_N64_T6464((uint64_t *) mmt_my_own_subvec(y),
            (uint64_t *) mmt_my_own_subvec(y),
            m1,
            mmt_my_own_size_in_items(y));
    y.consistency = 1;
    /* Compute m2 such that m2 * m1 * V is in RREF. We discard the
     * combinations which lead to zero, so that m2*m1 should have rank
     * precisely the rank of the nullspace */
    r = mmt_vec_echelon(m2, y);
    if (tcan_print) fmt::print("\trank(V cap nullspace(M)) == {}\n", r);
    /* Now we have the really interesting basis of the nullspace */
    for(int i = r ; i < 64 ; i++) {
        m2[i] = 0;
    }
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        auto f = fopen_helper(filename_base + "M2", "wb");
        size_t const rc = fwrite(m2.data(), sizeof(mat64), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
    }
    serialize(mmt.pi->m);
    /* Now apply m2*m1 to v, for real */
    mmt_full_vec_set(y, V[i0]);
    mul_N64_T6464((uint64_t *) mmt_my_own_subvec(y),
            (uint64_t *) mmt_my_own_subvec(y),
            m1,
            mmt_my_own_size_in_items(y));
    mul_N64_T6464((uint64_t *) mmt_my_own_subvec(y),
            (uint64_t *) mmt_my_own_subvec(y),
            m2,
            mmt_my_own_size_in_items(y));


    /* Now save the reduced kernel basis */
    ASSERT_ALWAYS(y.n == mmt.n[bw->dir]);
    mmt_vec_save(y, filename_base + "YR{}-{}", mmt.n0[bw->dir], 0);

    if (tcan_print) { fmt::print("Saving {}* ...done\n", filename_base); fflush(stdout); }

    mmt_vec_save(y, "blsolution{}-{}.0", mmt.n0[bw->dir], 0);
}


int blstate::operator()(parallelizing_info_ptr pi)
{
    int exit_code = 0;

    int const tcan_print = bw->can_print && pi->m->trank == 0;
    struct timing_data timing[1];

    size_t const nelts_for_nnmat = bw->n * (bw->n / A->simd_groupsize());

    serialize(pi->m);

    /* note that we don't know how to do checking for BL. Sure, we can
     * check for mutual orthogonality of vector blocks at the
     * checkpoints, but that is rather lame.
     */

    matmul_top_comm_bench(mmt, bw->dir);

    if (bw->start == 0) {
        set_start();
        save(0);
    } else {
        load(bw->start);
    }

    int auto_end = 0;

    if (bw->end == 0) {
        /* Decide on an automatic ending value */
        unsigned int length;
        length = mmt.n[bw->dir] / (bw->n - 0.76) + 10;
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
        fmt::print ("Target iteration is {} ; going to {}\n", bw->end,
                bw->interval * iceildiv(bw->end, bw->interval));
    }


    timing_init(timing, 4 * mmt.matrices.size(), bw->start, bw->interval * iceildiv(bw->end, bw->interval));
    for(size_t i = 0 ; i < mmt.matrices.size() ; i++) {
        timing_set_timer_name(timing, 4*i, "CPU%zu", i);
        timing_set_timer_items(timing, 4*i, mmt.matrices[i].mm->ncoeffs);
        timing_set_timer_name(timing, 4*i+1, "cpu-wait%zu", i);
        timing_set_timer_name(timing, 4*i+2, "COMM%zu", i);
        timing_set_timer_name(timing, 4*i+3, "comm-wait%zu", i);
    }

    arith_generic::elt * vav;
    arith_generic::elt * vaav;

    vav = A->alloc(nelts_for_nnmat, mat64::alignment);
    vaav = A->alloc(nelts_for_nnmat, mat64::alignment);

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
            mmt_vec_twist(mmt, V[k]);
        }
        mmt_own_vec_set(y, V[s % 3]);

        /* I *think* that this is necessary */
        mmt_vec_broadcast(y);

        /* FIXME: for BL, we use 8 timers while only 4 are defined in the
         * structure.
         */
        serialize(pi->m);
        int i, i0, i1, i2;

        for(i = 0 ; i < bw->interval ; i++) {
            i0 = (s + i) % 3;
            i1 = (s + i + 3 - 1) % 3;
            i2 = (s + i + 3 - 2) % 3;

            mmt_full_vec_set_zero(my);

            mmt_vec * yy[2] = {&y, &my};

            for(int d = 0 ; d < 2 ; d++) {
                /* The first part of this loop must be guaranteed to be free
                 * of any mpi calls */
                /* at this point, timer is [0] (CPU) */

                {
                    matmul_top_mul_cpu(mmt, 0, yy[d]->d, *yy[!d], *yy[d]);
                    timing_next_timer(timing);  /* now timer is [1] (cpu-wait) */
                }
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [2] (COMM) */
                mmt_vec_allreduce(*yy[!d]);

                timing_next_timer(timing);  /* now timer is [3] (comm-wait) */
                serialize(pi->m);           /* for measuring waits only */

                timing_next_timer(timing);  /* now timer is [0] (CPU) */
            }

            /* We have now:
             *  V_{n}           V[0]
             *  V_{n-1}         V[1]
             *  V_{n-2}         V[2]
             *  A * V_n         mcol->v
             */

            A->vec_set_zero(vav, nelts_for_nnmat);

            AxA->add_dotprod(vav,
                    mmt_my_own_subvec(V[i0]),
                    mmt_my_own_subvec(y),
                    mmt_my_own_size_in_items(y));

            pi_allreduce(nullptr, vav,
                    nelts_for_nnmat, mmt.pitype, BWC_PI_SUM, pi->m);

            A->vec_set_zero(vaav, nelts_for_nnmat);

            AxA->add_dotprod(vaav,
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(y),
                    mmt_my_own_size_in_items(y));

            pi_allreduce(nullptr, vaav, nelts_for_nnmat,
                    mmt.pitype, BWC_PI_SUM, pi->m);


            ASSERT_ALWAYS(D[i0]->n == 64);

            {
                size_t const eblock = mmt_my_own_size_in_items(y);
                ASSERT_ALWAYS(y.abase->elt_stride() == sizeof(uint64_t));

                // Here are the operations we will now perform
                //
                // X := (IN/*-Dn0*/)*Vn + Dn0*VA;
                // C0 := ((IN/*-Dn0*/)*vav + Dn0*vaav) * Ln0 * Vn;
                // C1 := Dn0 * vav * Ln1 * Vn1;
                // C2 := Dn0 * vav * (IN-Dn1) * Ln2 * Vn2;
                //
                // all are performed locally.

                // we set up X in mcol->v

                uint64_t * V0 = (uint64_t *) mmt_my_own_subvec(V[i0]);
                uint64_t * V1 = (uint64_t *) mmt_my_own_subvec(V[i1]);
                uint64_t * V2 = (uint64_t *) mmt_my_own_subvec(V[i2]);
                uint64_t * VA  = (uint64_t *) mmt_my_own_subvec(y);
                uint64_t * X   = (uint64_t *) mmt_my_own_subvec(y);
                uint64_t D0;
                uint64_t const D1 = D[i1]->p[0];
                mat64 & mvav = *(mat64*) vav;
                mat64 & mvaav = *(mat64*) vaav;
                mat64 & mL0 = L[i0];
                mat64  const& mL1 = L[i1];
                mat64 & mL2 = L[i2];
                mat64 m0, m1, m2, t;

                /* We need to save vav for use a wee bit later in this loop. */
                t = mvav;

                D0 = D[i0]->p[0] = extraction_step(mL0.data(), t.data(), D1);

                int const Ni = bit_vector_popcount(D[i0]);
                sum_Ni += Ni;
                // int defect = bw->n - Ni;
                // fmt::print("step {}, dim={}\n", s+i, Ni);
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
                y.consistency = 1;

                /* So. We need to get ready for the next step, don't we ? */
                mmt_vec_broadcast(y);
                mmt_own_vec_set(V[i2], y);
                V[i2].consistency = 1;
            }
            timing_check(pi, timing, s+i+1, tcan_print);
        }
        serialize(pi->m);

        for(int k = 0 ; k < 3 ; k++) {
            if (V[k].consistency < 2)
                mmt_vec_broadcast(V[k]);
            mmt_vec_untwist(mmt, V[k]);
        }

        serialize(pi->m);
        if (i < bw->interval) {
            if (tcan_print) fmt::print("Finished at iteration {}, sum_dim={}=N*({}-{:.3f})\n", s+i, sum_Ni, bw->n, bw->n-(double)sum_Ni/(s+i));

            save_result(s+i);
            /* We need to cheat somewhat */
            if (serialize_threads(pi->m)) {
                bw->end = s + i;
            }
            serialize_threads(pi->m);
            timing->end_mark = s + i;
            break;
        }


        save(s+bw->interval);
        serialize(pi->m);

        if (tcan_print) fmt::print("N={} ; sum_dim={}=N*({}-{:.3f})\n", s+i, sum_Ni, bw->n, bw->n-(double)sum_Ni/(s+i));
        // reached s + bw->interval. Count our time on cpu, and compute the sum.
        timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "blocklanczos");
    }

    // we can't do as we do with BW, since the process really stops in
    // the end, continuing the last block to the checkpoint mark does not
    // make sense.
    timing_final_tally(pi, timing, tcan_print, "blocklanczos");

    A->free(vav);
    A->free(vaav);

#if 0
    pi_log_clear(pi->m);
    pi_log_clear(pi->wr[0]);
    pi_log_clear(pi->wr[1]);
#endif

    if (auto_end && bw->end == auto_end) {
        if (tcan_print) {
            fmt::print("FAILED blocklanczos (no collapse found for inner product).\n");
        }
        if (serialize_threads(pi->m)) {
            exit_code = 1;
        }
        serialize_threads(pi->m);
    } else {
        if (tcan_print) {
            fmt::print("Done blocklanczos.\n");
        }
    }
    serialize(pi->m);

    timing_clear(timing);

    return exit_code;
}

static int exit_code = 0;

static void * bl_prog(parallelizing_info_ptr pi, cxx_param_list & pl, void * arg MAYBE_UNUSED)
{
    /* so many features we do not support ! */
    ASSERT_ALWAYS(bw->m == bw->n);
    ASSERT_ALWAYS(bw->ys[0] == 0);
    ASSERT_ALWAYS(bw->ys[1] == bw->n);

    /* Don't think we stand any chance with interleaving with block
     * Lanczos... */
    ASSERT_ALWAYS(!pi->interleaved);

    // int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;

    blstate bl(pi, pl);

    bl(pi);

    return nullptr;
}


// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

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
    if (bw->ys[0] < 0) { fmt::print(stderr, "no ys value set\n"); exit(1); }

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(bl_prog, pl, 0);

    parallelizing_info_finish();

    bw_common_clear(bw);

    return exit_code;
}
