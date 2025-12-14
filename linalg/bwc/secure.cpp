#include "cado.h" // IWYU pragma: keep
#include <cerrno>               // for ENOENT, errno
#include <climits>              // for INT_MIN
#include <cstdint>              // for uint32_t
#include <cstdlib>
#include <cstdio>

#include <string>                // for string
#include <memory>
#include <algorithm>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "gmp_aux.h"
#include "bw-common.h"
#include "macros.h"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "matmul_top_vec.hpp"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "select_mpi.h"
#include "xvectors.hpp"
#include "mmt_vector_pair.hpp"
#include "utils_cxx.hpp"
#include "bwc_filenames.hpp"
#include "misc.h"

/* We create the check data based on:
 *
 *  - the random seed
 *  - the X vector that was created by the prep program.
 *
 * The check data is made of four distinct sets of files.
 *
 * The first two are singletons, and depend on the random seed only.
 *
 * Ct0-$nchecks.0-$m (also referred to as T): a random matrix of size bw->m*nchecks (yes it's named like this because the data is written row-major, and the first interval in the name customarily denotes the number of items in a major division).
 *
 * Cr0-$nchecks.0-$nchecks (also referred to as R): a sequence of random matrices of size nchecks * nchecks
 *
 * Multiple instances of the other two can exist, depending on various
 * check distances.
 *
 * Cv0-$nchecks.$s (also referred to as C) : check vector for distance $s. Depends on X.
 * Cd0-$nchecks.$s (also referred to as D) : check vector for distance $s. Depends on X, T, and R.
 *
 * Cv0-<splitwidth>.<j> == trsp(M)^j * X * Ct (See note (T))
 * Cd0-<splitwidth>.<j> == \sum_{0<=i<j} trsp(M)^i * X * Ct * Cr[i] (See note (T))
 * 
 * (T): This assumes that nullspace=right. If nullspace=left, replace M by trsp(M).
 */

static void * sec_prog(parallelizing_info_ptr pi, cxx_param_list & pl, void * arg MAYBE_UNUSED)
{

    int const fake = param_list_lookup_string(pl, "random_matrix") != nullptr;

    ASSERT_ALWAYS(!pi->interleaved);

    int const tcan_print = bw->can_print && pi->m->trank == 0;

    int const withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, nchecks));

    /* We need that in order to do matrix products */
    std::unique_ptr<arith_cross_generic> AxA(arith_cross_generic::instance(A.get(), A.get()));

    matmul_top_data mmt(A.get(), pi, pl, bw->dir);
    pi_datatype_ptr A_pi = mmt.pitype;

    /* we work in the opposite direction compared to other programs */
    mmt_vector_pair myy(mmt, !bw->dir);
    mmt_vec & my = myy[0];

    mmt_vec dvec(mmt, nullptr, nullptr, !bw->dir, /* shared ! */ 1, mmt.n[!bw->dir]);

    unsigned int const unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    /* Because we're a special case, we _expect_ to work opposite to
     * optimized direction. So we pass bw->dir even though _we_ are going
     * to call mmt_mul with !bw->dir.
     */

    serialize(pi->m);

    /* To fill Cr and Ct, we need a random state */
    cxx_gmp_randstate rstate;

#if 0
    /* After all, a zero seed is fine, too */
    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = time(NULL);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }
#endif
    gmp_randseed_ui(rstate, bw->seed);
    if (tcan_print) {
        fmt::print("// Random generator seeded with {}\n", bw->seed);
    }

    /* Ct is a constant projection matrix of size bw->m * nchecks */
    /* It depends only on the random seed. We create it if start==0, or
     * reload it otherwise. */
    std::string const Tfilename = fmt::format("Ct0-{}.0-{}", nchecks, bw->m);
    size_t const T_coeff_size = A->vec_elt_stride(bw->m);
    arith_generic::elt * Tdata;
    Tdata = A->alloc(bw->m);

    /* Cr is a list of matrices of size nchecks * nchecks */
    /* It depends only on the random seed */
    std::string const Rfilename = fmt::format("Cr0-{}.0-{}", nchecks, nchecks);

    std::unique_ptr<FILE, delete_FILE> Rfile;
    size_t const R_coeff_size = A->vec_elt_stride(nchecks);


    /* {{{ First check consistency of existing files with the bw->start
     * value. We wish to abort early (and not touch any existing file!)
     * if an inconsistency is detected.
     */
    int consistency = 1;
    if (pi->m->jrank == 0 && pi->m->trank == 0) {
        Rfile = fopen_helper(Rfilename, "ab");
        size_t Rsz = file_bytes(Rfile.get());

        if (bw->start == 0) {
            if (Rsz) {
                fmt::print(stderr, "Refusing to overwrite {} with new random data\n", Rfilename);
                consistency = 0;
            }
        } else {
            if (!Rfile) {
                fmt::print(stderr,
                        "Cannot expand non-existing {}"
                        " with new random data\n", Rfilename);
                consistency = 0;
            } else if (Rsz != (size_t) bw->start * R_coeff_size) {
                fmt::print(stderr,
                        "Cannot expand {} ({} entries)"
                        " starting at position {}\n",
                        Rfilename,
                        (unsigned int) (Rsz / R_coeff_size), bw->start);
                consistency = 0;
            }
        }

        /* Non-destructively open for writing */
        auto T = fopen_helper(Tfilename, bw->start == 0 ? "wb" : "rb");
        size_t Tsz = file_bytes(T.get());

        if (!T) {
            fmt::print(stderr, "File {} cannot be opened\n",
                    Tfilename);
            consistency = 0;
        } else if (bw->start == 0 && Tsz) {
            fmt::print(stderr, "Refusing to overwrite {}"
                    " with new random data\n", Tfilename);
            consistency = 0;
        } else if (bw->start && Tsz != T_coeff_size) {
            fmt::print(stderr, "File {} has wrong size ({} != {}),"
                    " cannot expand check data\n", Tfilename,
                    Tsz, T_coeff_size);
            consistency = 0;
        }

        /* The non-master branch does exactly the same! */
        pi_bcast(&consistency, 1, BWC_PI_INT, 0, 0, pi->m);
        if (!consistency) pi_abort(EXIT_FAILURE, pi->m);

        /* {{{ create or load T, based on the random seed. */

        /* When start > 0, the call below does not care about the data it
         * generates, it only cares about the side effect to the random
         * state. We do it just in order to keep the random state
         * synchronized compared to what would have happened if we
         * started with start=0. It's cheap enough anyway.
         *
         * (also, random generation matters only at the leader node)
         */
        A->vec_set_random(Tdata, bw->m, rstate);
        if (bw->start == 0) {
            size_t rc = fwrite(Tdata, A->vec_elt_stride(bw->m), 1, T.get());
            ASSERT_ALWAYS(rc == 1);
            if (tcan_print) fmt::print("Saved {}\n", Tfilename);
        } else {
            /* We should be reading the same data, unless we changed
             * the seed.
             */
            size_t rc = fread(Tdata, A->vec_elt_stride(bw->m), 1, T.get());
            ASSERT_ALWAYS(rc == 1);
            if (tcan_print) fmt::print("Loaded {}\n", Tfilename);
        }
        /* }}} */

        ASSERT_ALWAYS(Rfile);

        pi_bcast(Tdata, bw->m, A_pi, 0, 0, pi->m);
    } else {
        pi_bcast(&consistency, 1, BWC_PI_INT, 0, 0, pi->m);
        if (!consistency) pi_abort(EXIT_FAILURE, pi->m);
        pi_bcast(Tdata, bw->m, A_pi, 0, 0, pi->m);
    }
    /* }}} */

    /* same remark as above. We want the random state to be in sync even
     * if bw->start>0 and we don't _really_ have stuff to generate for
     * this data that's already there.
     */

    if (bw->start) {
        arith_generic::elt * Rdata = A->alloc(nchecks, ALIGNMENT_ON_ALL_BWC_VECTORS);
        for(int k = 0 ; k < bw->start ; k++) {
            /* same remark as above */
            A->vec_set_random(Rdata, nchecks, rstate);
        }

        A->free(Rdata);
    }

    /* {{{ create initial Cv and Cd, or load them if start>0 */
    if (bw->start == 0) {
        if (tcan_print)
            fmt::print("We have start=0: creating Cv0-{}.0 as an expanded copy of X*T\n", nchecks);

        std::vector<uint32_t> gxvecs;
        unsigned int nx = 0;

        if (!fake) {
            gxvecs = load_x(bw->m, nx, pi);
        } else {
            gxvecs = set_x_fake(bw->m, nx, pi);
        }

        mmt_full_vec_set_zero(my);
        ASSERT_ALWAYS(bw->m % nchecks == 0);

        for(int c = 0 ; c < bw->m ; c += nchecks) {
            mmt_vec_set_x_indices(dvec, gxvecs, c, c + nchecks, nx);
            AxA->addmul_tiny(
                    mmt_my_own_subvec(my),
                    mmt_my_own_subvec(dvec),
                    A->vec_subvec(Tdata, c),
                    mmt_my_own_size_in_items(my));
        }
        /* addmul_tiny degrades consistency ! */
        my.consistency = 1;
        mmt_vec_broadcast(my);
        mmt_vec_save(my, bwc_Cv_file::pattern(0), unpadded, 0);

        mmt_full_vec_set_zero(dvec);
    } else {
        int ok;

        ok = mmt_vec_load(my,   bwc_Cv_file::pattern(bw->start), unpadded, 0);
        ASSERT_ALWAYS(ok);
        ok = mmt_vec_load(dvec, bwc_Cd_file::pattern(bw->start), unpadded, 0);
        ASSERT_ALWAYS(ok);
    }
    /* }}} */

    /* {{{ adjust the list of check stops according to the check_stops
     * and interval parameters
     */
    serialize_threads(pi->m);
    if (pi->m->trank == 0) {
        /* the bw object is global ! */
        bw_set_length_and_interval_krylov(bw, mmt.n0);
    }
    serialize_threads(pi->m);
    ASSERT_ALWAYS(bw->end % bw->interval == 0);

    /* easier to deal with a copy on the stack */
    std::vector<int> check_stops(bw->check_stops, bw->check_stops + bw->number_of_check_stops);

    /* see #30025 -- we want a sensible default */
    if (check_stops.empty()) {
        check_stops.push_back(0);
        check_stops.push_back(bw->interval);
        int const a = std::min(16, bw->interval / 2);
        if (a) {
            /* if interval == 1, don't bother */
            check_stops.push_back(a);
            check_stops.push_back(bw->interval + a);
        }
    }

    /* if something was provided, make really sure that at least the
     * interval value is within the list of check stops */
    if (std::find(check_stops.begin(), check_stops.end(), bw->interval) == check_stops.end()) {
        check_stops.push_back(bw->interval);
    }
    std::ranges::sort(check_stops);
    serialize_threads(pi->m);

    if (tcan_print) {
        fmt::print("Computing trsp(x)*M^k for check stops k={}\n",
                join(check_stops, ","));
    }
    serialize(pi->m);
    /* }}} */

    // {{{ kill the warning about wrong spmv direction
    for(auto const & Mloc : mmt.matrices) {
        Mloc.mm->iteration[!bw->dir] = INT_MIN;
    }
    // }}}

    int k = bw->start;
    for(int const next : check_stops) {
        serialize(pi->m);
        if (next == 0) {
            /* if 0 is in check_stops, we don't want to create files such
             * as Cd0-64.0 */
            continue;
        }
        if (next < k) {
            /* This may happen when start is passed and is beyond the
             * first check stop.
             */
            continue;
        }
        mmt_vec_twist(mmt, my);
        mmt_vec_twist(mmt, dvec);

        /* Allocate an area in memory where we store the random
         * coefficient. We do so in order to avoid the situation where an
         * interrupted execution of the binary has already committed to
         * disk the random coefficients (possibly cut halfway because of
         * stdio buffering), while no Cd file has been written yet.
         * Therefore this is mostly a matter of consistency.
         */
        arith_generic::elt * Rdata_stream = nullptr;
        int const k0 = k;
        if (next > k0) {
            Rdata_stream = A->alloc(nchecks * (next - k0), ALIGNMENT_ON_ALL_BWC_VECTORS);
            A->vec_set_zero(Rdata_stream, nchecks * (next - k0));
            A->vec_set_random(Rdata_stream, nchecks * (next - k0), rstate);
            pi_bcast(Rdata_stream, nchecks * (next - k0), A_pi, 0, 0, pi->m);
        }
    
        for( ; k < next ; k++) {
            /* new random coefficient in R */
            arith_generic::elt * Rdata = A->vec_subvec(Rdata_stream, (k-k0) * nchecks);
            /* At this point Rdata should be consistent across all
             * threads */
            AxA->addmul_tiny(
                    mmt_my_own_subvec(dvec),
                    mmt_my_own_subvec(my),
                    Rdata,
                    mmt_my_own_size_in_items(my));
            /* addmul_tiny degrades consistency ! */
            dvec.consistency = 1;
            mmt_vec_broadcast(dvec);
            pi_log_op(mmt.pi->m, "iteration %d", k);
            matmul_top_mul(mmt, myy.vectors(), nullptr);

            if (tcan_print) {
                putchar('.');
                fflush(stdout);
            }
        }
        serialize(pi->m);
        serialize_threads(mmt.pi->m);
        mmt_vec_untwist(mmt, my);
        mmt_vec_untwist(mmt, dvec);

        mmt_vec_save(my,   bwc_Cv_file::pattern(k), unpadded, 0);
        mmt_vec_save(dvec, bwc_Cd_file::pattern(k), unpadded, 0);
        if (pi->m->trank == 0 && pi->m->jrank == 0 && next > k0) {
            const size_t rc = fwrite(Rdata_stream, A->vec_elt_stride(nchecks), next - k0, Rfile.get());
            ASSERT_ALWAYS(rc == (size_t) (next - k0));
            fflush(Rfile.get());
        }
        A->free(Rdata_stream);
    }

    A->free(Tdata);

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

    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(sec_prog, pl, nullptr);

    parallelizing_info_finish();

    bw_common_clear(bw);
    return 0;
}
