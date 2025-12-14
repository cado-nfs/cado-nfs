#include "cado.h" // IWYU pragma: keep

#include <cerrno>              // for ENOENT, errno
#include <climits>             // for UINT_MAX
#include <cstdlib>             // for exit, EXIT_FAILURE
#include <cstring>             // for strncmp, strlen, strrchr
#include <cstdio>               // for size_t, sscanf, fclose, fopen, fread

#include <algorithm>            // for sort, min
#include <map>                  // for map<>::mapped_type, _Rb_tree_iterator
#include <memory>
#include <stdexcept>            // for runtime_error
#include <string>               // for string, basic_string, operator<<, cha...
#include <utility>              // for pair, make_pair
#include <vector>               // for vector

#include <sys/stat.h>           // for stat

#include <gmp.h>                // for mpz_cmp_ui
#include "fmt/base.h"
#include "fmt/format.h"

#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "bw-common.h"          // for bw, bw_common_clear, bw_common_decl_u...
#include "macros.h"             // for ASSERT_ALWAYS, MAYBE_UNUSED
#include "misc.h"               // ok_NOKNOK
#include "params.h"             // for param_list_clear, param_list_init
#include "portability.h" // asprintf // IWYU pragma: keep
#include "select_mpi.h"         // for MPI_Abort, MPI_Comm_rank, MPI_COMM_WORLD
#include "utils_cxx.hpp"
#include "bwc_filenames.hpp"

// arguments we need
// m
// n
// prime
//
// wdir only as a way to shorten file names.
// nullspace ??

static void vec_alloc(arith_generic * A, arith_generic::elt *& z, size_t vsize)
{
    z = A->alloc(vsize, ALIGNMENT_ON_ALL_BWC_VECTORS);
    A->vec_set_zero(z, vsize);
}

static void vec_free(arith_generic * A, arith_generic::elt *& z, size_t vsize MAYBE_UNUSED)
{
    A->free(z);
}

static int vec_read(arith_generic * A, void * z, bwc_file_base const & v, size_t vsize, const char * prefix = nullptr)
{
    fmt::print("{} {} ...", prefix, v);
    auto f = fopen_helper(std::string(v), "rb", true);
    if (f) {
        size_t const rc = fread(z, A->elt_stride(), vsize, f.get());
        if (rc == vsize) {
            fmt::print("{}", " done\n");
            return 1;
        }
    }
    fmt::print("{}", " failed\n");
    return 0;
}

static size_t vec_items(arith_generic * A, bwc_file_base const & v)
{
    struct stat sbuf[1];
    const std::string s(v);
    int const rc = stat(s.c_str(), sbuf);
    if (rc < 0 && errno == ENOENT)
        return 0;
    ASSERT_ALWAYS(rc == 0);
    return sbuf->st_size / A->elt_stride();
}

template<typename T>
static size_t common_size(arith_generic * Ac, std::vector<T> const & Cfiles, const char * name)
{
    size_t vsize = 0;
    T vsize_first;

    for(auto & C : Cfiles) {
        size_t items = vec_items(Ac, C);
        if (items == 0) {
            fmt::print("{} has disappeared\n", C);
            continue;
        }
        if (vsize == 0) {
            vsize = items;
            vsize_first = C;
        } else if (vsize != items) {
            fmt::print(stderr,
                    "File sizes disagree for {} ({} items) and {} ({} items)\n",
                    vsize_first, vsize, C, items);
            exit(EXIT_FAILURE);
        }
    }
    if (vsize) fmt::print("{} files have {} coordinates\n", name, vsize);
    return vsize;
}

using vseq_t = std::map<bwc_column_range, std::vector<bwc_V_file> >;

static void check_V_files(arith_generic * Ac, vseq_t & Vsequences, std::vector<bwc_Cv_file> & Cfiles, int & nfailed)/*{{{*/
{
    if (Cfiles.empty()) return;

    size_t const nchecks = Ac->simd_groupsize();
    size_t vsize = common_size(Ac, Cfiles, "Cv");

    for(unsigned int i0 = 0 ; i0 + 1 < Cfiles.size() ; i0++) {
        bwc_Cv_file const & C_i0(Cfiles[i0]);
        arith_generic::elt * Cv_i0;
        vec_alloc(Ac, Cv_i0, vsize);

        int has_read_Cv_i0 = 0;

        for(unsigned int i1 = i0 + 1 ; i1 < Cfiles.size() ; i1++) {
            bwc_Cv_file const & C_i1(Cfiles[i1]);

            fmt::print("Doing checks for distance {} using {} and {}\n",
                    C_i1.stretch - C_i0.stretch,
                    C_i0, C_i1);

            arith_generic::elt * Cv_i1;
            vec_alloc(Ac, Cv_i1, vsize);

            int has_read_Cv_i1 = 0;

            /* {{{ check all V files together */
            for(auto const & it : Vsequences) {
                unsigned int j0 = it.first[0];
                unsigned int j1 = it.first[1];
                auto const & Vs(it.second);

                fmt::print(" checks on V files for sequence {}-{}\n", j0, j1);

                std::unique_ptr<arith_generic> Av(arith_generic::instance(bw->p,j1-j0));

                std::unique_ptr<arith_cross_generic> AcxAv(arith_cross_generic::instance(Ac, Av.get()));


                /* {{{ Check that all V files here have the proper size */
                for(auto const & V : Vs) {
                    size_t items = vec_items(Av.get(), V);
                    if (items != vsize)
                        throw std::invalid_argument(fmt::format(
                                    "{} has {} coordinates,"
                                    " different from expected {}\n",
                                    V, items, vsize));
                }
                /* }}} */

                arith_generic::elt * Vv;
                vec_alloc(Av.get(), Vv, vsize);
                unsigned int Vv_iter = UINT_MAX;

                arith_generic::elt * dotprod_scratch[2];
                vec_alloc(Av.get(), dotprod_scratch[0], nchecks);
                vec_alloc(Av.get(), dotprod_scratch[1], nchecks);

                unsigned int j = 0;
                for(unsigned int i = 0 ; i < Vs.size() ; i++) {
                    for(j = i + 1 ; j < Vs.size() ; j++) {
                        if (Vs[j].n + C_i0.stretch == Vs[i].n + C_i1.stretch)
                            break;
                    }
                    if (j == Vs.size()) continue;

                    if (!has_read_Cv_i0++)
                        if (!vec_read(Ac, Cv_i0, C_i0, vsize, " "))
                            continue;

                    if (!has_read_Cv_i1++)
                        if (!vec_read(Ac, Cv_i1, C_i1, vsize, " "))
                            continue;

                    // We might want to store how many checks we had for
                    // each vector. This can be done in a global map.
                    // Vs[j].checks++;

                    std::string vi { Vs[i] };
                    if (Vs[i].dirname() == C_i1.dirname())
                        vi = Vs[i].basename();

                    std::string vj { Vs[j] };
                    if (Vs[j].dirname() == C_i1.dirname())
                        vj = Vs[j].basename();

                    fmt::print("  check {} against {}\n", vi, vj);
                    if (Vs[i].n != Vv_iter) {
                        if (!vec_read(Ac, Vv, Vs[i], vsize, "   "))
                            continue;
                        Vv_iter = Vs[i].n;
                    }

                    Av->vec_set_zero(dotprod_scratch[0], nchecks);

                    /* compute the dot product */
                    AcxAv->add_dotprod(
                            dotprod_scratch[0],
                            Cv_i1,
                            Vv,
                            vsize);

                    if (!vec_read(Ac, Vv, Vs[j], vsize, "   "))
                        continue;

                    Vv_iter = Vs[j].n;

                    Av->vec_set_zero(dotprod_scratch[1], nchecks);
                    AcxAv->add_dotprod(
                            dotprod_scratch[1],
                            Cv_i0,
                            Vv,
                            vsize);

                    int const cmp = Av->vec_cmp(dotprod_scratch[0], dotprod_scratch[1], nchecks);

                    std::string diag = fmt::format(
                            "  check {} against {} -> {}\n",
                            vi, vj, ok_NOKNOK(cmp == 0));

                    fmt::print("{}", diag);

                    if (cmp != 0) {
                        nfailed++;
                        fmt::print(stderr, "{}", diag);
                    }
                }
                vec_free(Av.get(), dotprod_scratch[0], nchecks);
                vec_free(Av.get(), dotprod_scratch[1], nchecks);
                vec_free(Ac, Vv, vsize);

            }
            /* }}} */

            vec_free(Ac, Cv_i1, vsize);
        }
        vec_free(Ac, Cv_i0, vsize);
    }
}/*}}}*/

static void check_A_files(arith_generic * Ac, std::vector<bwc_V_file> const & Vfiles, std::vector<bwc_A_file> const & Afiles, std::vector<bwc_Cd_file> const & Dfiles, bwc_Cr_file & R, bwc_Ct_file & T, int & nfailed)
{
    if (Dfiles.empty())
        return;
    arith_generic::elt * Dv = nullptr;
    size_t const vsize = common_size(Ac, Dfiles, "Cd");
    size_t const nchecks = Ac->simd_groupsize();

    size_t rsize = vec_items(Ac, R) / nchecks;
    fmt::print("Cr file has {} coordinates\n", rsize);

    ASSERT_ALWAYS(vec_items(Ac, T) ==  (size_t) bw->m);

    vec_alloc(Ac, Dv, vsize);

    arith_generic::elt * Tdata = nullptr;
    Tdata = Ac->alloc(bw->m, ALIGNMENT_ON_ALL_BWC_VECTORS);
    {
        auto f = fopen_helper(T, "rb");
        const size_t rc = fread(Tdata, Ac->vec_elt_stride(bw->m), 1, f.get());
        ASSERT_ALWAYS(rc == 1);
    }

    arith_generic::elt * Rdata = nullptr;
    Rdata = Ac->alloc(Ac->vec_elt_stride(nchecks) * rsize, ALIGNMENT_ON_ALL_BWC_VECTORS);
    {
        auto f = fopen_helper(R, "rb");
        const size_t rc = fread(Rdata, Ac->vec_elt_stride(nchecks), rsize, f.get());
        ASSERT_ALWAYS(rc == rsize);
    }

    for(auto const & D : Dfiles) {
        if (D.stretch == 0)
            continue;
        if (D.stretch > rsize) {
            fmt::print(stderr, "Cannot do checks using {}, too few items in R file\n", R);
            continue;
        }
        fmt::print("Doing A file checks for distance {} using {}\n"
                "  (as well as {} and {})\n",
                D.stretch, D, T, R);
        int has_read_D = 0;
        /* first scan potential base files V, and the restrict to cases
         * where we have all the required A files to do a check...
         */
        for(auto const & V0 : Vfiles) {
            /* Look for "A%u-%u.%u-%u" % (V0.j0, V0.j1, V0.n,
             * V0.n+D.stretch), or any collection of files that would
             * concatenate to exactly this file */
            std::vector<bwc_A_file> a_list;
            size_t n_reach = V0.n;
            /* A_files are sorted */
            for(auto const & A : Afiles) {
                if (A.j1 <= V0.j0) continue;
                if (A.j0 > V0.j0) break;
                if (A.n0 > n_reach) break;
                if (A.n1 <= n_reach) continue;
                if (A.n0 <= n_reach) {
                    ASSERT_ALWAYS(n_reach == V0.n || n_reach == A.n0);
                    /* Since we (still) collate A files, it's important
                     * to be able to do this check even if V starts in
                     * the middle of the A file -- this accounts for the
                     * n_reach == V0.n possibility. Otherwise, we must
                     * go from one A file to the next, hence the
                     * requirement that n_reach == A.n0 .*/
                    n_reach = A.n1;
                    a_list.push_back(A);
                }
                if (n_reach >= V0.n + D.stretch)
                    break;
            }

            if (n_reach < V0.n + D.stretch)
                continue;

            fmt::print("  check {} against {} entries of{}\n",
                    V0, D.stretch, join(a_list, " "));

            std::unique_ptr<arith_generic> Av(arith_generic::instance(bw->p, V0.j1 - V0.j0));
            std::unique_ptr<arith_cross_generic> AcxAv(arith_cross_generic::instance(Ac, Av.get()));

            arith_generic::elt * dotprod_scratch[3];
            vec_alloc(Av.get(), dotprod_scratch[0], nchecks);
            vec_alloc(Av.get(), dotprod_scratch[1], nchecks);
            vec_alloc(Av.get(), dotprod_scratch[2], nchecks);

            /* read data from the A files. We redo the detection loop
             * that we had above */
            n_reach = V0.n;
            for(auto const & A : Afiles) {
                if (A.j1 <= V0.j0) continue;
                if (A.j0 > V0.j0) break;
                if (A.n0 > n_reach) break;
                if (A.n1 <= n_reach) continue;
                if (A.n0 <= n_reach) {
                    ASSERT_ALWAYS(n_reach == V0.n || n_reach == A.n0);
                    fmt::print("   read {} small {}*{} matrices from {}\n",
                            std::min(A.n1, V0.n + D.stretch) - n_reach,
                            bw->m, bw->n, A);
                    auto a = fopen_helper(A, "rb");
                    for(unsigned int p = n_reach ; p < A.n1 && p < V0.n + D.stretch ; p++) {
                        Av->vec_set_zero(dotprod_scratch[1], nchecks);
                        for(int c = 0 ; c < bw->m ; c += (int) nchecks) {
                            size_t rc;
                            for(int r = 0 ; r < (int) nchecks ; r++) {
                                size_t const simd = Av->simd_groupsize();
                                size_t const rowsize = (A.j1 - A.j0) / simd;
                                size_t const matsize = bw->m * rowsize;
                                size_t const nmats = p - A.n0;

                                rc = (size_t) fseek(a.get(),
                                            nmats * Av->vec_elt_stride(matsize) +
                                            (c + r) * Av->vec_elt_stride(rowsize) +
                                            Av->vec_elt_stride((V0.j0 - A.j0) / simd),
                                        SEEK_SET);
                                ASSERT_ALWAYS(rc == 0);
                                rc = fread(Av->vec_subvec(dotprod_scratch[2], r), Av->elt_stride(), 1, a.get());
                                ASSERT_ALWAYS(rc == 1);
                            }
                            AcxAv->add_dotprod(
                                    dotprod_scratch[1],
                                    Ac->vec_subvec(Tdata, c),
                                    dotprod_scratch[2],
                                    nchecks);
                        }
                        AcxAv->add_dotprod(
                                dotprod_scratch[0],
                                Ac->vec_subvec(Rdata, (p - V0.n) * nchecks),
                                dotprod_scratch[1],
                                nchecks);
                    }
                    n_reach = A.n1;
                }
                if (n_reach >= V0.n + D.stretch)
                    break;
            }

            int can_check = 1;

            if (!has_read_D++)
                if (!vec_read(Ac, Dv, D, vsize, "   "))
                    can_check = 0;

            if (can_check) {
                arith_generic::elt * Vv;

                vec_alloc(Av.get(), Vv, vsize);
                if (!vec_read(Av.get(), Vv, V0, vsize, "   "))
                    can_check = 0;

                if (can_check) {
                    Av->vec_set_zero(dotprod_scratch[1], nchecks);
                    AcxAv->add_dotprod(
                            dotprod_scratch[1],
                            Dv,
                            Vv,
                            vsize);
                }

                vec_free(Av.get(), Vv, vsize);
            }

            if (can_check) {
                int const cmp = Av->vec_cmp(dotprod_scratch[0], dotprod_scratch[1], nchecks);

                std::string diag = fmt::format(
                        "  check {} against {} entries of{} -> {}\n",
                        V0, D.stretch, join(a_list, " "),
                        ok_NOKNOK(cmp == 0));
                fmt::print("{}", diag);

                if (cmp != 0) {
                    nfailed++;
                    fmt::print(stderr, "{}", diag);
                }
            }

            if (!can_check)
                fmt::print("{}", "  (check aborted because of missing files)\n");

            vec_free(Av.get(), dotprod_scratch[2], nchecks);
            vec_free(Av.get(), dotprod_scratch[1], nchecks);
            vec_free(Av.get(), dotprod_scratch[0], nchecks);
        }
    }
    Ac->free(Rdata);
    Ac->free(Tdata);
    vec_free(Ac, Dv, vsize);
}

/* This is *not* a parallel program, so we depart significantly from the
 * way programs such as krylov or mksol are written.
 *
 */
static void * check_prog(cxx_param_list & pl MAYBE_UNUSED, int argc, char const * argv[])
{
    int const withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int const nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    std::unique_ptr<arith_generic> const Ac(arith_generic::instance(bw->p, nchecks));

    std::vector<bwc_Cv_file> Cfiles;
    std::vector<bwc_Cd_file> Dfiles;
    std::vector<bwc_Cr_file> Rfiles;
    std::vector<bwc_Ct_file> Tfiles;
    std::vector<bwc_V_file> Vfiles;
    std::vector<bwc_A_file> Afiles;
    /* we classify the files below, but we don't use them (yet) */
    std::vector<bwc_F_file> Ffiles;
    std::vector<bwc_S_file> Sfiles;

    for(int i = 0 ; i < argc ; i++) {
        std::string s(argv[i]);
        if (s.back() == '~')
            continue;
        if (decomposed_path(s).extension() == ".rhs")
            continue;
        if (bwc_file_base::match(Cfiles, s)
                || bwc_file_base::match(Dfiles, s)
                || bwc_file_base::match(Rfiles, s)
                || bwc_file_base::match(Tfiles, s)
                || bwc_file_base::match(Afiles, s)
                || bwc_file_base::match(Ffiles, s)
                || bwc_file_base::match(Sfiles, s)
                || bwc_file_base::match(Vfiles, s))
            continue;
        else
            throw std::invalid_argument(fmt::format(
                        "Invalid file name: {}\n", s));
    }

    /* {{{ some consistency checks */
    for(auto & C : Cfiles) {
        ASSERT_ALWAYS(C.j0 == 0);
        ASSERT_ALWAYS(C.j1 == (unsigned int) nchecks);
    }
    for(auto & D : Dfiles) {
        ASSERT_ALWAYS(D.j0 == 0);
        ASSERT_ALWAYS(D.j1 == (unsigned int) nchecks);
    }
    ASSERT_ALWAYS(Rfiles.size() <= 1);
    ASSERT_ALWAYS(Tfiles.size() <= 1);
    for(auto & R : Rfiles) {
        ASSERT_ALWAYS(R.nchecks == (unsigned int) nchecks);
    }
    for(auto & T : Tfiles) {
        ASSERT_ALWAYS(T.nchecks == (unsigned int) nchecks);
        /* T files for different m's could maybe coexist, at least in
         * theory. However since both C and D depend on T, this does not
         * seem very viable */
        ASSERT_ALWAYS(T.m == (unsigned int) bw->m);
    }
    /* }}} */

    /* {{{ sort */
    std::ranges::sort(Cfiles);
    std::ranges::sort(Dfiles);
    std::ranges::sort(Rfiles);
    std::ranges::sort(Afiles);
    std::ranges::sort(Vfiles);
    // std::ranges::sort(Ffiles);
    // std::ranges::sort(Sfiles);
    /* }}} */

    /* {{{ split V files in sequences -- what for ? */
    /* How many (j0-j1) ranges do we have for V files */
    vseq_t Vsequences;
    for(auto const & V : Vfiles)
        Vsequences[V].push_back(V);
    for(auto & seq: Vsequences)
        std::ranges::sort(seq.second);
    /* }}} */

    int nfailed = 0;

    check_V_files(Ac.get(), Vsequences, Cfiles, nfailed);

    /* Check A files using V, D, T, and R */
    if ((Tfiles.empty() || Rfiles.empty()) && !Dfiles.empty()) {
        fmt::print(stderr, "{}", "It makes no sense to provide Cd files and no Cr and Ct file\n");
        exit(EXIT_FAILURE);
    } else if (!Tfiles.empty() && !Rfiles.empty() && !Dfiles.empty()) {
        check_A_files(Ac.get(), Vfiles, Afiles, Dfiles, Rfiles.front(), Tfiles.front(), nfailed);
    }

    if (nfailed) {
        std::string diag = fmt::format("{} checks FAILED !!!!!!!!!!!!!!!!!\n", nfailed);
        fmt::print("{}", diag);
        fmt::print(stderr, "{}", diag);
        exit(EXIT_FAILURE);
    }


#if 0
    uint32_t * gxvecs = nullptr;
    unsigned int nx = 0;
    if (!fake) {
        load_x(&gxvecs, bw->m, &nx, pi);
    } else {
        set_x_fake(&gxvecs, bw->m, &nx, pi);
    }

    if (!bw->skip_online_checks) {
        /* We do the dot product by working on the local vector chunks.
         * Therefore, we must really understand the check vector as
         * playing a role in the very same direction of the y vector!
         */
        mmt_vec_init(mmt, Ac, Ac_pi,
                check_vector, bw->dir, THREAD_SHARED_VECTOR, mmt->n[bw->dir]);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        mmt_vec_load(check_vector, CHECK_FILE_BASE, bw->interval,  mmt->n0[bw->dir]);
        if (tcan_print) { printf("done\n"); }
    }

    if (!bw->skip_online_checks) {
        ahead = Ac->alloc(nchecks, ALIGNMENT_ON_ALL_BWC_VECTORS);
    }

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    void * xymats;

    if (tcan_print) {
        printf("Each thread allocates %zd kb for the Ac matrices\n",
                Ac->vec_elt_stride(Ac, bw->m*bw->interval) >> 10);
    }
    if (!bw->skip_online_checks) {
        Ac->vec_set_zero(Ac, ahead, nchecks);
        AxAc->dotprod(Ac, Ac, ahead,
                mmt_my_own_subvec(check_vector),
                mmt_my_own_subvec(ymy[0]),
                mmt_my_own_size_in_items(ymy[0]));
    }

    for(int i = 0 ; i < bw->interval ; i++) {
        /* Compute the product by x */
        x_dotprod(Ac->vec_subvec(Ac, xymats, i * bw->m),
                gxvecs, bw->m, nx, ymy[0], 1);

        matmul_top_mul(mmt, ymy, timing);

        timing_check(pi, timing, s+i+1, tcan_print);
    }
    serialize(pi->m);

    /* See remark above. */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    if (!bw->skip_online_checks) {
        /* Last dot product. This must cancel ! */
        x_dotprod(ahead, gxvecs, nchecks, nx, ymy[0], -1);

        pi_allreduce(nullptr, ahead, nchecks, mmt->pitype, BWC_PI_SUM, pi->m);
        if (!Ac->vec_is_zero(Ac, ahead, nchecks)) {
            printf("Failed check at iteration %d\n", s + bw->interval);
            exit(1);
        }
    }

    mmt_vec_untwist(mmt, ymy[0]);

    /* Now (and only now) collect the xy matrices */
    pi_allreduce(nullptr, xymats,
            bw->m * bw->interval,
            mmt->pitype, BWC_PI_SUM, pi->m);

    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        char * tmp;
        int rc;
        rc = asprintf(&tmp, A_FILE_PATTERN, ys[0], ys[1], s, s+bw->interval);
        FILE * f = fopen(tmp, "wb");
        rc = fwrite(xymats, Ac->vec_elt_stride(Ac, 1), bw->m*bw->interval, f);
        if (rc != bw->m*bw->interval) {
            fprintf(stderr, "Ayee -- short write\n");
            // make sure our input data won't be deleted -- this
            // chunk will have to be redone later, maybe the disk
            // failure is temporary (?)
        }
        fclose(f);
        free(tmp);
    }

    mmt_vec_save(ymy[0], v_name, s + bw->interval, unpadded);

    if (pi->m->trank == 0 && pi->m->jrank == 0)
        keep_rolling_checkpoints(v_name, s + bw->interval);

    serialize(pi->m);

    // reached s + bw->interval. Count our time on cpu, and compute the sum.
    timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "check");

    free(gxvecs);
    free(v_name);

    for(int i = 0 ; i < mmt->nmatrices + nmats_odd ; i++) {
        mmt_vec_clear(mmt, ymy[i]);
    }
    free(ymy);
#endif

    return nullptr;
}


int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

    param_list_usage_header(pl,
            "File names are checked together (all relevant checks tying files to eachother within the provided argument list are done). The program tells which checks have been done\n"
            "Options are as follows. Note that not all are relevant to this program specifically:\n");

    bw_common_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    check_prog(pl, argc, argv);

    bw_common_clear(bw);

    return 0;
}

