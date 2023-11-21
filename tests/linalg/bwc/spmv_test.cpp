#include "cado.h" // IWYU pragma: keep
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>                 // for gmp_randclear, gmp_randinit_default
#include <stdint.h>              // for uint32_t, uint64_t
#include <string.h>              // for memset
#include "balancing.hpp"           // for balancing_clear, balancing_init, bal...
#include "parallelizing_info.hpp"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "params.h"
#include "bw-common.h"
#include "arith-generic.hpp"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "macros.h"

int verbose = 0;

void mmt_vec_set_0n(mmt_vec & v, size_t items)
{
    serialize(v.pi->m);
    /* For debug: set to vector [0, 1, ..., n[
     *
     * Here we do something a bit fishy. We don't really have a way
     * to set to the integer i, when in fact we're talking a simd
     * thing: IOW, we have set_ui_at, but no set_ui. So let's do a
     * dirty cast */
    // ASSERT_ALWAYS((size_t) v.abase->elt_stride() <= sizeof(uint64_t));
    ASSERT_ALWAYS(v.abase->elt_stride() >= sizeof(unsigned int));
    // size_t nwords = (size_t) v.abase->elt_stride(v.abase) / sizeof(uint64_t);
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
    arith_generic::elt * data = mmt_my_own_subvec(v);
    /* Put 0's everywhere, and put i at other places (just with a dirty
     * cast) */
    v.abase->vec_set_zero(mmt_my_own_subvec(v), mmt_my_own_size_in_items(v));
    for(size_t s = 0 ; s < sz ; s++) {
        arith_generic::elt & u = v.abase->vec_item(data, s);
        /* yes, this is awful */
        auto ptr = reinterpret_cast<unsigned int *>(&u);
        *ptr = (v.i0 + off + s < items) ? v.i0 + off + s : 0;
    }
    v.consistency = 1;
    serialize(v.pi->m);
    mmt_vec_broadcast(v);
    serialize(v.pi->m);
}

/* check that v[i] == p[i] */
void mmt_vec_check_equal_0n(mmt_vec const & v, size_t items)
{
    serialize(v.pi->m);
    ASSERT_ALWAYS(v.consistency == 2);
    // ASSERT_ALWAYS((size_t) v.abase->elt_stride() <= sizeof(uint64_t));
    ASSERT_ALWAYS(v.abase->elt_stride() >= sizeof(unsigned int));
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
    arith_generic::elt const * data = mmt_my_own_subvec(v);
    arith_generic::elt * temp_alloc = v.abase->alloc(1);
    arith_generic::elt & temp = *temp_alloc;
    auto ptr = reinterpret_cast<unsigned int *>(temp_alloc);
    for(size_t s = 0 ; s < sz ; s++) {
        arith_generic::elt const & u = v.abase->vec_item(data, s);
        v.abase->set(temp, u);
        /* yes, this is awful */
        *ptr ^= (v.i0 + off + s < items ? v.i0 + off + s : 0);
        ASSERT_ALWAYS(v.abase->is_zero(temp));
    }
    v.abase->free(temp_alloc);
}

/* check that v[i] == p[i] */
void mmt_vec_check_equal_0n_permuted(mmt_vec const & v, size_t items, uint32_t * p)
{
    serialize(v.pi->m);
    ASSERT_ALWAYS(v.consistency == 2);
    // ASSERT_ALWAYS((size_t) v.abase->elt_stride() <= sizeof(uint64_t));
    ASSERT_ALWAYS(v.abase->elt_stride() >= sizeof(unsigned int));
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
    arith_generic::elt const * data = mmt_my_own_subvec(v);
    arith_generic::elt * temp_alloc = v.abase->alloc(1);
    arith_generic::elt & temp = *temp_alloc;
    auto ptr = reinterpret_cast<unsigned int *>(temp_alloc);
    for(size_t s = 0 ; s < sz ; s++) {
        arith_generic::elt const & u = v.abase->vec_item(data, s);
        v.abase->set(temp, u);
        /* yes, this is awful */
        *ptr ^= (v.i0 + off + s < items ? p[v.i0 + off + s] : 0);
        ASSERT_ALWAYS(v.abase->is_zero(temp));
    }
    v.abase->free(temp_alloc);
}

/* check that v[i] == p^-1[i] */
void mmt_vec_check_equal_0n_inv_permuted(mmt_vec const & v, size_t items, uint32_t * p)
{
    serialize(v.pi->m);
    // ASSERT_ALWAYS((size_t) v.abase->elt_stride() <= sizeof(uint64_t));
    ASSERT_ALWAYS(v.abase->elt_stride() >= sizeof(unsigned int));
    size_t off = mmt_my_own_offset_in_items(v);
    size_t sz = mmt_my_own_size_in_items(v);
    arith_generic::elt const * data = mmt_my_own_subvec(v);
    arith_generic::elt * temp_alloc = v.abase->alloc(1);
    arith_generic::elt & temp = *temp_alloc;
    auto ptr = reinterpret_cast<unsigned int *>(temp_alloc);
    for(size_t s = 0 ; s < sz ; s++) {
        arith_generic::elt const & u = v.abase->vec_item(data, s);
        v.abase->set(temp, u);
        /* yes, this is awful */
        if (v.i0 + off + s < items) {
            ASSERT_ALWAYS(*ptr < items);
            ASSERT_ALWAYS(p[*ptr] == v.i0 + off + s);
            *ptr = 0;
        }
        ASSERT_ALWAYS(v.abase->is_zero(temp));
    }
    serialize_threads(v.pi->wr[v.d]);
    v.abase->free(temp_alloc);
}

/* This only does a multiplication */

void * tst_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    ASSERT_ALWAYS(!pi->interleaved);

    if (verbose) {
        pi_log_init(pi->m);
        pi_log_init(pi->wr[0]);
        pi_log_init(pi->wr[1]);
    }

    gmp_randstate_t rstate;

    // int tcan_print = bw->can_print && pi->m->trank == 0;

    // int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    // int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, bw->n));

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, bw->seed);

    /* recall that d here is only for one optimized direction. The other
     * direction should still work. Note that changing the optimized
     * direction is going to trigger a different choice of cache files. */
    /* TODO check this with the various combinations of d */
    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    /* check that mmt_apply_identity does what we expect. This function
     * goes from a vector in one direction to a vector in another
     * direction (but leaves to us the task of doing the reduction,
     * however).
     */
    {
        mmt_vec v  (mmt,0,0, 0, /* shared ! */ 1, mmt.n[0]);
        mmt_vec vi (mmt,0,0, 1,                0, mmt.n[1]);
        mmt_vec vii(mmt,0,0, 0,                0, mmt.n[0]);
        serialize(pi->m);
        mmt_vec_set_0n(v, mmt.n0[0]);
        /* We save the Z files, although it's useless, the checking done
         * here is allright */
        mmt_vec_save(v, "Z%u-%u.0", mmt.n0[0], 0);
        mmt_apply_identity(vi, v);
        mmt_vec_allreduce(vi);
        mmt_vec_clear_padding(vi, mmt.n0[1], mmt.n0[0]);
        mmt_vec_check_equal_0n(vi, mmt.n0[1]);
        mmt_vec_save(vi, "ZI%u-%u.0", mmt.n0[1], 0);
        mmt_apply_identity(vii, vi);
        mmt_vec_allreduce(vii);
        mmt_vec_check_equal_0n(vii, MIN(mmt.n0[0], mmt.n0[1]));
        mmt_vec_save(vi, "ZII%u-%u.0", mmt.n0[0], 0);
    }


    /* Now check that mmt_apply_S and mmt_unapply_S do what they are
     * supposed to.  */
    {
        balancing bb;
        const char * bname = param_list_lookup_string(pl, "balancing");
        balancing_init(bb);

        if (mmt.pi->m->jrank == 0 && mmt.pi->m->trank == 0) {
            balancing_read(bb, bname);
        }
        pi_bcast(&bb, sizeof(balancing), BWC_PI_BYTE, 0, 0, mmt.pi->m);
        /* fix the mmt.rowperm and mmt.colperm */
        if (bb.rowperm) {
            if (mmt.pi->m->jrank || mmt.pi->m->trank) {
                bb.rowperm = (uint32_t *) malloc(bb.trows * sizeof(uint32_t));
            }
            pi_bcast(bb.rowperm, bb.trows * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, mmt.pi->m);
        }
        if (bb.colperm) {
            if (mmt.pi->m->jrank || mmt.pi->m->trank) {
                bb.colperm = (uint32_t *) malloc(bb.tcols * sizeof(uint32_t));
            }
            pi_bcast(bb.colperm, bb.tcols * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, mmt.pi->m);
        }



        for(int test_shared = 0 ; test_shared < 2 ; test_shared++) {
            serialize(pi->m);
            uint32_t * xr = bb.rowperm;
            uint32_t * xc = bb.colperm;
            uint32_t *freeme[2] = {NULL,NULL};
            if ((mmt.matrices[0].bal.flags & FLAG_REPLICATE) && (xr || xc)) {
                ASSERT_ALWAYS(xc || xr);
                ASSERT_ALWAYS(mmt.matrices[0].bal.trows == mmt.matrices[0].bal.tcols);
                /* P is the permutation which sends
                 * sub-block nv*i+j to sub-block nh*j+i
                 */
                unsigned int nh = mmt.matrices[0].bal.nh;
                unsigned int nv = mmt.matrices[0].bal.nv;
                size_t z = mmt.matrices[0].bal.trows / (nh * nv);
                if (!xr) {
                    /* implicit Sr is P * Sc */
                    xr = (uint32_t *) malloc(mmt.matrices[0].bal.trows * sizeof(uint32_t));
                    freeme[0] = xr;
                    /* The image of i is Sc(P(i)) */
                    for(size_t i = 0 ; i < nh ; i++) {
                        for(size_t j = 0 ; j < nv ; j++) {
                            for(size_t k = 0 ; k < z ; k++) {
                                /* The image of (nv*i+j)*z+k by P is
                                 * (nh*j+i)*z+k */
                                /* And its image by P*Sc is therefore
                                 * xc[(nh*j+i)*z+k];
                                 */
                                xr[(nv*i+j)*z+k] = xc[(nh*j+i)*z+k];
                            }
                        }
                    }
                }
                if (!xc) {
                    /* implicit Sc is P^-1 * Sr */
                    xc = (uint32_t *) malloc(mmt.matrices[0].bal.tcols * sizeof(uint32_t));
                    freeme[1] = xc;
                    /* The image of i is Sr(P^-1(i)) */
                    for(size_t i = 0 ; i < nh ; i++) {
                        for(size_t j = 0 ; j < nv ; j++) {
                            for(size_t k = 0 ; k < z ; k++) {
                                /* The image of (nh*j+i)*z+k by P^-1 is
                                 * (nv*i+j)*z+k */
                                /* And its image by P^-1*Sc is therefore
                                 * xr[(nv*i+j)*z+k];
                                 */
                                xc[(nh*j+i)*z+k] = xr[(nv*i+j)*z+k];
                            }
                        }
                    }
                }
            }
            if (xc) {
                mmt_vec v(mmt,0,0, 1, test_shared, mmt.n[1]);

                /* Because we want to test the permutation, we need to
                 * fill our dummy vector with the full range [0..n[, not
                 * just [0.n0[ padded with zeroes. Two reasons for this:
                 *  - the permutation wors on the index set [0..n[ anyway
                 *    (not just [0..n0[),
                 *  - in order to tell whether we're doing the right
                 *    thing, we value uniqueness of data in the vector.
                 */
                mmt_vec_set_0n(v, mmt.n[1]);
                /* trsp(v) <- trsp(v) * S: coefficient i goes to position S(i).
                 * Hence we must check that in position j, we have what
                 * was beforehand in position S^{-1}(i), right. Which is
                 * precisely S^{-1}(j).
                 */
                /* apply == 1 d == 1 Sc defined:   v <- v * Sc
                 * apply == 1 d == 1 Sc implicit:  v <- v * Sr
                 */
                mmt_vec_apply_S(mmt, 0, v);
                mmt_vec_check_equal_0n_inv_permuted(v, mmt.n[1], xc);
                /* apply == 0 d == 1 Sc defined:   v <- v * Sc^-1
                 * apply == 0 d == 1 Sc implicit:  v <- v * Sr^-1
                 */
                mmt_vec_unapply_S(mmt, 0, v);
                mmt_vec_check_equal_0n(v, mmt.n[1]);

                /* Do the same check in the other direction as well. Most
                 * probably redundant, though. */
                mmt_vec_set_0n(v, mmt.n[1]);
                mmt_vec_unapply_S(mmt, 0, v);
                mmt_vec_check_equal_0n_permuted(v, mmt.n[1], xc);
                mmt_vec_apply_S(mmt, 0, v);
                mmt_vec_check_equal_0n(v, mmt.n[1]);
            }
            if (xr) {
                /* same for row vectors */
                mmt_vec v(mmt,0,0, 0, test_shared, mmt.n[0]);

                mmt_vec_set_0n(v, mmt.n[0]);
                /* v <- v * S:  coefficient i goes to position S(i).  */
                /* apply == 1 d == 0 Sr defined:   v <- v * Sr
                 * apply == 1 d == 0 Sr implicit:  v <- v * Sc
                 */
                mmt_vec_apply_S(mmt, 0, v);
                mmt_vec_check_equal_0n_inv_permuted(v, mmt.n[0], xr);
                /*
                 * apply == 0 d == 0 Sr defined:   v <- v * Sr^-1
                 * apply == 0 d == 0 Sr implicit:  v <- v * Sc^-1
                 */
                mmt_vec_unapply_S(mmt, 0, v);
                mmt_vec_check_equal_0n(v, mmt.n[0]);

                mmt_vec_set_0n(v, mmt.n[0]);
                mmt_vec_unapply_S(mmt, 0, v);
                mmt_vec_check_equal_0n_permuted(v, mmt.n[0], xr);
                mmt_vec_apply_S(mmt, 0, v);
                mmt_vec_check_equal_0n(v, mmt.n[1]);
            }
            if (freeme[0]) free(freeme[0]);
            if (freeme[1]) free(freeme[1]);
        }
        balancing_clear(bb);
    }


    /* matrix times vector product */
    /* The file pair (Y.0, MY.0) will be checked for correctness against the
     * result obtained by short_matmul */
    {
        mmt_vec y( mmt,0,0, 1, /* shared ! */ 1, mmt.n[1]);
        mmt_vec my(mmt,0,0, 0,                0, mmt.n[0]);
        serialize(pi->m);
        mmt_vec_set_random_through_file(y, "Y%u-%u.0", mmt.n0[1], rstate, 0);
        /* recall that for all purposes, bwc operates with M*T^-1 and not M
         */
        mmt_vec_apply_T(mmt, y);
        mmt_vec_twist(mmt, y);
        matmul_top_mul_cpu(mmt, 0, y.d, my, y);
        /* watch out -- at this point we are *NOT* doing
         * matmul_top_mul_comm, so there's no multiplication by P^-1 !
         */
        mmt_vec_allreduce(my);
        mmt_vec_untwist(mmt, my);
        mmt_vec_save(my, "MY%u-%u.0", mmt.n0[0], 0);
    }

    /* vector times matrix product */
    /* The file pair (W, WM) will be checked for correctness against the
     * result obtained by short_matmul */
    {
        mmt_vec w(  mmt,0,0, 0, /* shared ! */ 1, mmt.n[0]);
        mmt_vec wm( mmt,0,0, 1,                0, mmt.n[1]);
        serialize(pi->m);
        mmt_vec_set_random_through_file(w, "W%u-%u.0", mmt.n0[0], rstate, 0);
        mmt_vec_twist(mmt, w);
        matmul_top_mul_cpu(mmt, 0, w.d, wm, w);
        /* It's not the same if we do allreduce+untwist, thus stay on the
         * image side, or reduce (changing side) + untwist, which brings
         * us back to the original side.
         *
         * So this test is probably most relevant for rectangular
         * matrices, but we may do it as well for square matrices. The
         * only catch is that in the latter case, we should pay attention
         * to what it does exactly.
         */
        mmt_vec_allreduce(wm);
        mmt_vec_untwist(mmt, wm);
        mmt_vec_unapply_T(mmt, wm);
        mmt_vec_save(wm, "WM%u-%u.0", mmt.n0[1], 0);
    }

    /* indices_twist */
    {
        uint32_t * xx;
        unsigned int m = 2;
        unsigned int nx = 2;
        xx = (uint32_t *) malloc(m * nx * sizeof(uint32_t));
        for(unsigned int k = 0 ; k < m * nx ; k++) {
            xx[k] = gmp_urandomm_ui(rstate, mmt.n0[0]);
        }
        pi_bcast(xx, m * nx * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);
        for(int d = 0 ; d < 2 ; d++)  {
            char * tmp;
            int rc;
            mmt_vec x(mmt,0,0, d, /* shared ! */ 1, mmt.n[d]);

            /* prepare a first vector */
            mmt_vec_set_x_indices(x, xx, m, nx);
            rc = asprintf(&tmp, "Xa%s.%d", "%u-%u", d);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(x, tmp, mmt.n[d], 0);
            free(tmp);

            /* and then a second vector, which should be equal */
            indices_twist(mmt, xx, nx * m, d);
            mmt_vec_set_x_indices(x, xx, m, nx);
            mmt_vec_untwist(mmt, x);
            rc = asprintf(&tmp, "Xb%s.%d", "%u-%u", d);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(x, tmp, mmt.n[d], 0);
            free(tmp);

            /* now for the twisted versions, too */
            mmt_vec_set_x_indices(x, xx, m, nx);
            mmt_vec_twist(mmt, x);
            rc = asprintf(&tmp, "XTa%s.%d", "%u-%u", d);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(x, tmp, mmt.n[d], 0);
            free(tmp);

            indices_twist(mmt, xx, nx * m, d);
            mmt_vec_set_x_indices(x, xx, m, nx);
            rc = asprintf(&tmp, "XTb%s.%d", "%u-%u", d);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(x, tmp, mmt.n[d], 0);
            free(tmp);
        }
        free(xx);
    }

    /* could also be tested:
     * - simply the permutation P
     */

#if 0
    matmul_top_apply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xA", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xB", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_apply_P_apply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xC", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_S_unapply_P(mmt, d);
    mmt_vec_save(mmt, NULL, "xD", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_apply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xE", d, 0);

    mmt_vec_load(mmt, NULL, "Y", d, 0);
    matmul_top_unapply_S(mmt, d);
    mmt_vec_save(mmt, NULL, "xF", d, 0);
#endif

#if 0
    // load to mcol->v, apply the permutation from the balancing.
    // Communicate, save mcol->v
    // rr, r => zv eq Pr*Sr^-1*Pr^-1*zu
    // rr, m->fw_colperm[rr], => zv eq Pr*Sr*zu
    //
    // other direction: zv - Pr^-1*Sr^-1*zu
    // twist
    // bw->dir==1 d==1 Pr^-1*Sr^-1*zu
    // bw->dir==1 d==0 Pr*Sr*zu
    mmt_vec_load(mmt, NULL, "U", d, 0);
    matmul_top_twist_vector(mmt, d);
    mmt_vec_save(mmt, NULL, "V", d, 0);
#endif

#if 0
    // untwist
    // rr, m->fw_colperm[rr], => zv - Sr^-1*Pr^-1*zu;
    // other direction: Sr*Pr*zu
    // bw->dir==1 d==0 Sr^-1*Pr^-1*zu;
    // bw->dir==1 d==1 Sr*Pr*zu;
    mmt_vec_load(mmt, NULL, "U", d, 0);
    matmul_top_untwist_vector(mmt, d);
    mmt_vec_save(mmt, NULL, "V", d, 0);
#endif
    // apply_identity(mmt, d);
    // mmt_vec_save(mmt, NULL, "V", !d, 0);

    serialize(pi->m);

#if 0
    mmt_vec oldv[2];
    memcpy(oldv[0], mmt.wr[0]->v, sizeof(mmt_vec));
    memcpy(oldv[1], mmt.wr[1]->v, sizeof(mmt_vec));
    memset(mmt.wr[0]->v, 0, sizeof(mmt_vec));
    memset(mmt.wr[1]->v, 0, sizeof(mmt_vec));
    mmt_vec_init(mmt, NULL, NULL, NULL, 0, flags[1 ^ 0]);
    mmt_vec_init(mmt, NULL, NULL, NULL, 1, flags[1 ^ 1]);

    serialize_threads(mmt.pi->m);

    size_t stride = abbytes(mmt.abase, 1);
    ASSERT_ALWAYS((mrow->v.flags & THREAD_SHARED_VECTOR) == 0);
    if (pirow->trank == 0 && pirow->jrank == 0) {
        mpfq_generic_copy(stride, mmt.wr[!bw->dir]->v.v, oldv[!bw->dir]->v, mmt.wr[!bw->dir]->i1 - mmt.wr[!bw->dir]->i0);
    }
    serialize_threads(mmt.pi->m);

    matmul_top_mul_comm(mmt, bw->dir);
    mmt_vec_save(mmt, NULL, CHECK_FILE_BASE, bw->dir, bw->interval);

    mmt_vec_clear(mmt, NULL, 0);
    mmt_vec_clear(mmt, NULL, 1);
    memcpy(mmt.wr[0]->v, oldv[0], sizeof(mmt_vec));
    memcpy(mmt.wr[1]->v, oldv[1], sizeof(mmt_vec));
    memset(oldv[0], 0, sizeof(mmt_vec));
    memset(oldv[1], 0, sizeof(mmt_vec));
#endif

    // A->oo_field_clear(A);
    gmp_randclear(rstate);

    if (verbose) {
        pi_log_print_all(pi);

        pi_log_clear(pi->m);
        pi_log_clear(pi->wr[0]);
        pi_log_clear(pi->wr[1]);
    }

    return NULL;
}


void usage()
{
    exit(1);
}

int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    param_list_decl_usage(pl, "v", "turn on some demo logging");
    matmul_top_decl_usage(pl);

    /* declare local parameters and switches. */
    param_list_configure_switch(pl, "v", &verbose);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    pi_go(tst_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}

