/* Schirokauer maps

Input:

* A list of the (npurged) a,b pairs. This is obtained from the
  purgedfile.
* A matrix of (small_nrows) rows and (npurged) cols, which indicates
  the contents of each relation-set. This is obtained from the
  indexfile.
* The sub-group order (ell) such that ell | p-1
  Note: All computations are done mod ell^2.
* (eps): the exponent used in the computation of the Schirokauer maps.
  Note: eps = ppcm(eps_i), where eps_i = ell^(deg(f_i)) - 1 and f = f_1 ... f_k
mod ell

Output

* A matrix of (small_nrows) rows and (nmaps)=deg(f) cols (mpz_t).  For each
  relation (rel) the (nmaps) Schirokauer maps are computed as the second
  least-significant digit of the ell-adic representation of the polynomial
  equal to (rel^eps - 1) / ell.

  In case of two algebraic sides, SM's are computed for sides 0..1 in that
  order.

*/

#include "cado.h" // IWYU pragma: keep

#include "cado_poly.h"
#include "filter_io.h"
#include "gmp_aux.h"
#include "gzip.h"
#include "macros.h"
#include "mpz_poly.h"
#include "omp_proxy.h"
#include "params.h"

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <gmp.h>

#include "purgedfile.h"
#include "select_mpi.h"
#include "sm_utils.hpp"
#include "stats.h"
#include "timing.h"
#include "verbose.h"

static stats_data_t stats; /* struct for printing progress */

static void * thread_sm(void * context_data, earlyparsed_relation_ptr rel)
{
    auto & ps = *static_cast<std::vector<pair_and_sides>*>(context_data);
    ps[rel->num] = pair_and_sides(rel->a, rel->b, rel->active_sides[0],
                                  rel->active_sides[1]);

    return nullptr;
}

static std::vector<sm_relset>
build_rel_sets(char const * purgedname, char const * indexname,
               uint64_t * small_nrows, std::vector<mpz_poly_srcptr> const & F,
               cxx_mpz const & ell2)
{
    uint64_t nrows, ncols, len_relset;
    uint64_t r[MAX_LEN_RELSET];
    int64_t e[MAX_LEN_RELSET];
    int ret;

    /* array of (a,b) pairs from (purgedname) file */

    purgedfile_read_firstline(purgedname, &nrows, &ncols);

    std::vector<pair_and_sides> pairs(nrows);

    /* For each rel, read the a,b-pair and init the corresponding poly pairs[]
     */
    fprintf(stdout, "\n# Reading %" PRIu64 " (a,b) pairs\n", nrows);
    fflush(stdout);
    char const * fic[2] = {purgedname, nullptr};
    filter_rels(fic, (filter_rels_callback_t)thread_sm, &pairs,
                EARLYPARSE_NEED_AB_HEXA, nullptr, nullptr);

    /* Array of (small_nrows) relation-sets built from array (pairs) and
       (indexname) file  */
    FILE * ix = fopen_maybe_compressed(indexname, "r");
    ASSERT_ALWAYS(ix != nullptr);

    ret = fscanf(ix, "%" SCNu64 "\n", small_nrows);
    ASSERT_ALWAYS(ret == 1);

    std::vector<sm_relset> rels;
    rels.reserve(*small_nrows);

    fprintf(stdout, "\n# Building %" PRIu64 " relation-sets\n", *small_nrows);
    fflush(stdout);
    uint64_t i;
    stats_init(stats, stdout, &i, nbits(*small_nrows) - 5, "Computed",
               "relation-sets", "", "relsets");
    for (i = 0; i < *small_nrows; i++) {
        ret = fscanf(ix, "%" SCNu64 "", &len_relset);
        ASSERT_ALWAYS(ret == 1 && len_relset < MAX_LEN_RELSET);

        for (uint64_t k = 0; k < len_relset; k++) {
            ret = fscanf(ix, " %" SCNx64 ":%" SCNd64 "", &r[k], &e[k]);
            ASSERT_ALWAYS(ret == 2);
        }

        rels.emplace_back(
            sm_build_one_relset(r, e, len_relset, pairs, F, ell2));

        if (stats_test_progress(stats))
            stats_print_progress(stats, i, 0, 0, 0);
    }
    stats_print_progress(stats, *small_nrows, 0, 0, 1);
    fclose_maybe_compressed(ix, indexname);

    return rels;
}

static void print_all_sm(FILE * out, std::vector<sm_side_info> const & sm_info,
                         int nb_polys, std::vector<cxx_mpz_poly> const & sm)
{
    ASSERT_ALWAYS(sm_info.size() == (size_t)nb_polys);
    ASSERT_ALWAYS(sm.size() == (size_t)nb_polys);
    for (int side = 0, c = 0; side < nb_polys; side++) {
        if (sm_info[side].nsm == 0)
            continue;
        if (c++)
            fprintf(out, " ");
        print_sm(out, sm_info[side], sm[side]);
    }
    fprintf(out, "\n");
}

// Basic MPI communications

static void MPI_Send_mpz(mpz_srcptr z, int dst)
{
    mp_size_t nlimbs = mpz_size(z);
    MPI_Send(&nlimbs, 1, CADO_MPI_MP_SIZE_T, dst, 0, MPI_COMM_WORLD);
    MPI_Send(&z->_mp_size, 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, dst, 0,
             MPI_COMM_WORLD);
    MPI_Send(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, dst, 0, MPI_COMM_WORLD);
}

static void MPI_Recv_mpz(mpz_ptr z, int src)
{
    mp_size_t nlimbs;
    MPI_Recv(&nlimbs, 1, CADO_MPI_MP_SIZE_T, src, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    _mpz_realloc(z, nlimbs);
    MPI_Recv(&z->_mp_size, 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, src, 0,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, src, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
}

static void MPI_Send_mpz_poly(mpz_poly_srcptr poly, int dst)
{
    MPI_Send(&poly->deg, 1, MPI_INT, dst, 0, MPI_COMM_WORLD);
    for (int i = 0; i <= poly->deg; ++i)
        MPI_Send_mpz(mpz_poly_coeff_const(poly, i), dst);
}

static void MPI_Recv_mpz_poly(mpz_poly_ptr poly, int src)
{
    MPI_Recv(&poly->deg, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    ASSERT_ALWAYS(poly->deg + 1 >= 0);
    /* FIXME -- what about already allocated coefficients ? */
    ASSERT_ALWAYS(poly->alloc == 0);
    mpz_poly_realloc(poly, poly->deg + 1);
    for (int i = 0; i <= poly->deg; ++i)
        MPI_Recv_mpz(mpz_poly_coeff(poly, i), src);
}

static void MPI_Send_relset(sm_relset const & relset, int dst, int nb_polys)
{
    ASSERT_ALWAYS(relset.nb_polys() == nb_polys);
    for (int i = 0; i < nb_polys; ++i) {
        MPI_Send_mpz_poly(relset.num[i], dst);
        MPI_Send_mpz_poly(relset.denom[i], dst);
    }
}

static void MPI_Recv_relset(sm_relset & relset, int src, int nb_polys)
{
    relset = sm_relset(nb_polys);
    for (int i = 0; i < nb_polys; ++i) {
        MPI_Recv_mpz_poly(relset.num[i], src);
        MPI_Recv_mpz_poly(relset.denom[i], src);
    }
}

static void MPI_Send_res(std::vector<cxx_mpz_poly> const & res, int dst,
                         std::vector<sm_side_info> const & sm_info,
                         int nb_polys)
{
    ASSERT_ALWAYS(sm_info.size() == (size_t)nb_polys);
    ASSERT_ALWAYS(res.size() == (size_t)nb_polys);
    for (int side = 0; side < nb_polys; side++) {
        if (sm_info[side].nsm == 0)
            continue;
        MPI_Send_mpz_poly(res[side], dst);
    }
}

static void MPI_Recv_res(std::vector<cxx_mpz_poly> & res, int src,
                         std::vector<sm_side_info> const & sm_info,
                         int nb_polys)
{
    res.assign(nb_polys, cxx_mpz_poly());
    for (int side = 0; side < nb_polys; side++) {
        if (sm_info[side].nsm == 0)
            continue;
        MPI_Recv_mpz_poly(res[side], src);
    }
}

// Pthread part: on each node, we use shared memory instead of mpi

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "purged", "(required) purged file");
    param_list_decl_usage(pl, "index", "(required) index file");
    param_list_decl_usage(pl, "out", "output file (stdout if not given)");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
    param_list_decl_usage(pl, "nsms",
                          "number of SM on side 0,1,... (default is "
                          "computed by the program)");
    param_list_decl_usage(pl, "t",
                          "number of threads on each mpi job (default 1)");
    verbose_decl_usage(pl);
}

static void usage(char const * argv, char const * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit(EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

int main(int argc, char const ** argv)
{
    MPI_Init(&argc, (char ***)&argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int const idoio = (rank == 0); // Am I the job allowed to do I/O ?
    double t0;

    char const * argv0 = argv[0];

    char const * polyfile = nullptr;
    char const * purgedfile = nullptr;
    char const * indexfile = nullptr;
    char const * outfile = nullptr;

    cxx_param_list pl;
    cxx_cado_poly cpoly;

    uint64_t nb_relsets;
    cxx_mpz ell, ell2;

    /* read params */
    declare_usage(pl);

    if (argc == 1)
        usage(argv[0], nullptr, pl);

    argc--, argv++;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        if (idoio) {
            fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
            usage(argv0, nullptr, pl);
        }
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == nullptr) {
        if (idoio) {
            fprintf(stderr, "Error: parameter -poly is mandatory\n");
            param_list_print_usage(pl, argv0, stderr);
        }
        exit(EXIT_FAILURE);
    }

    /* Read purged filename from command line */
    if ((purgedfile = param_list_lookup_string(pl, "purged")) == nullptr) {
        if (idoio) {
            fprintf(stderr, "Error: parameter -purged is mandatory\n");
            param_list_print_usage(pl, argv0, stderr);
        }
        exit(EXIT_FAILURE);
    }

    /* Read index filename from command line */
    if ((indexfile = param_list_lookup_string(pl, "index")) == nullptr) {
        if (idoio) {
            fprintf(stderr, "Error: parameter -index is mandatory\n");
            param_list_print_usage(pl, argv0, stderr);
        }
        exit(EXIT_FAILURE);
    }

    /* Read outfile filename from command line ; defaults to stdout. */
    outfile = param_list_lookup_string(pl, "out");

    /* Read ell from command line (assuming radix 10) */
    if (!param_list_parse(pl, "ell", ell)) {
        if (idoio) {
            fprintf(stderr, "Error: parameter -ell is mandatory\n");
            param_list_print_usage(pl, argv0, stderr);
        }
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    if (!cado_poly_read(cpoly, polyfile)) {
        if (idoio) {
            fprintf(stderr, "Error reading polynomial file\n");
        }
        exit(EXIT_FAILURE);
    }

    /* negative value means that the value that will be used is the value
     * computed later by sm_side_info_init */
    std::vector<int> nsm_arg(cpoly->nb_polys, -1);
    /* Read number of sm to be printed from command line */
    param_list_parse_int_args_per_side(pl, "nsm", nsm_arg.data(),
                                       cpoly->nb_polys,
                                       ARGS_PER_SIDE_DEFAULT_AS_IS);

    std::vector<mpz_poly_srcptr> F(cpoly->nb_polys);

    for (int side = 0; side < cpoly->nb_polys; side++) {
        F[side] = cpoly->pols[side];
        if (nsm_arg[side] > F[side]->deg) {
            if (idoio) {
                fprintf(stderr,
                        "Error: nsm%d=%d can not exceed the degree=%d\n", side,
                        nsm_arg[side], F[side]->deg);
            }
            exit(EXIT_FAILURE);
        }
    }

    char const * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

    if (param_list_warn_unused(pl)) {
        if (idoio) {
            usage(argv0, nullptr, pl);
        } else {
            exit(EXIT_FAILURE);
        }
    }

    /* Print ell and ell^2 */
    mpz_mul(ell2, ell, ell);
    if (idoio)
        fmt::print("# Sub-group order:\nell = {}\n# Computation is done "
                   "modulo ell2 = ell^2:\nell2 = {}\n",
                   ell, ell2);

    std::vector<sm_side_info> sm_info;

    for (int side = 0; side < cpoly->nb_polys; side++) {
        sm_info.emplace_back(F[side], ell, 0);
        sm_info[side].set_mode(sm_mode_string);
    }

    for (int side = 0; side < cpoly->nb_polys; side++) {
        if (idoio) {
            fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side,
                    side);
            mpz_poly_fprintf(stdout, F[side]);
            printf("# SM info on side %d:\n", side);
            sm_info[side].print(stdout);
        }
        if (nsm_arg[side] >= 0) {
            if (idoio && nsm_arg[side] < sm_info[side].nsm)
                printf("# Warning: as per the command line parameters, we will "
                       "only compute %d SMs instead of %d on side %d. Be sure "
                       "to know what you are doing!\n",
                       nsm_arg[side], sm_info[side].nsm, side);
            sm_info[side].nsm = nsm_arg[side]; /* command line wins */
        }
        if (idoio)
            printf("# Will compute %d SMs on side %d\n", sm_info[side].nsm,
                   side);

        /* do some consistency checks */
        if (sm_info[side].unit_rank != sm_info[side].nsm) {
            if (idoio)
                fprintf(stderr,
                        "# On side %d, unit rank is %d, computing %d SMs ; "
                        "weird.\n",
                        side, sm_info[side].unit_rank, sm_info[side].nsm);
            /* for the 0 case, we haven't computed anything: prevent the
             * user from asking SM data anyway */
            ASSERT_ALWAYS(sm_info[side].unit_rank != 0);
        }
    }
    fflush(stdout);

    // If nsm (or nsm_arg) is 0 on one side, then set F[side] to nullptr to
    // desactivate the corresponding computations.
    // (note that it does not seem to be a good idea to forcibly deactivate
    // the SMs on a given side, I think).
    for (int side = 0; side < cpoly->nb_polys; ++side) {
        if (sm_info[side].nsm == 0) {
            F[side] = nullptr;
        }
    }

    ///////////////////////
    // Only process 0 constructs the relation sets.
#ifdef HAVE_OPENMP
    unsigned int thmax = omp_get_max_threads();
#else
    unsigned int const thmax = 1;
#endif
    std::vector<sm_relset> rels;
    if (rank == 0) {
        rels = build_rel_sets(purgedfile, indexfile, &nb_relsets, F, ell2);
        fprintf(stdout,
                "\n# Computing Schirokauer maps for %" PRIu64
                " relation-sets, using %d threads and %d jobs.\n",
                nb_relsets, thmax, size);
        fflush(stdout);
    }
    MPI_Bcast(&nb_relsets, 1, CADO_MPI_UINT64_T, 0, MPI_COMM_WORLD);

    ///////////////////////
    // Send a share of the rel sets to each process (round Robin)
    uint64_t const nb_parts = (nb_relsets - 1) / size + 1; // ceiling
    std::vector<sm_relset> part_rels(nb_parts);
    if (rank == 0) {
        for (uint64_t i = 0; i < nb_parts; ++i) {
            part_rels[i] = rels[i * size];
            for (int j = 1; j < size; ++j) {
                if (i * size + j < nb_relsets)
                    MPI_Send_relset(rels[i * size + j], j, cpoly->nb_polys);
            }
        }
    } else {
        for (uint64_t i = 0; i < nb_parts; ++i) {
            if (i * size + rank < nb_relsets)
                MPI_Recv_relset(part_rels[i], 0, cpoly->nb_polys);
        }
    }

    ///////////////////////
    // Process the relsets.
    t0 = seconds();

    // nb_parts vectors of nb_polys polynomials
    std::vector<std::vector<cxx_mpz_poly>> dst;
    dst.assign(nb_parts,
               std::vector<cxx_mpz_poly>(cpoly->nb_polys, cxx_mpz_poly()));

    /* updated only by thread 0 */
    uint64_t count_processed_sm = 0;
    stats_init(stats, stdout, &count_processed_sm, nbits(nb_relsets) - 5,
               "Computed", "SMs", "", "SMs");

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
#ifdef HAVE_OPENMP
        unsigned int thid = omp_get_thread_num();
        unsigned int thnb = omp_get_num_threads();
#else
        unsigned int const thid = 0;
        unsigned int const thnb = 1;
#endif
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
        for (uint64_t i = 0; i < nb_relsets; i++) {
            for (int side = 0; side < cpoly->nb_polys; side++) {
                if (sm_info[side].nsm == 0)
                    continue;
                mpz_poly_reduce_frac_mod_f_mod_mpz(
                    part_rels[i].num[side], part_rels[i].denom[side],
                    sm_info[side].f0, sm_info[side].ell2);
                cxx_mpz_poly SM;
                cxx_mpz_poly u;
                mpz_poly_set(u, part_rels[i].num[side]);
                sm_info[side].compute_piecewise(SM, u);
                mpz_poly_swap(dst[i][side], SM);
            }
            if (thid == 0) {
                /* Static schedule, all threads can reasonably be expected
                 * to progress at the same speed */
                count_processed_sm += thnb;
                if (stats_test_progress(stats))
                    stats_print_progress(stats, count_processed_sm, 0, 0, 0);
            }
        }
    }
    stats_print_progress(stats, nb_parts, 0, 0, 1);

    fprintf(stderr, "Job %d: processed all relsets in %f s\n", rank,
            seconds() - t0);

    // Send back results and print
    if (rank != 0) { // sender
        for (uint64_t i = 0; i < nb_parts; ++i) {
            if (i * size + rank < nb_relsets)
                MPI_Send_res(dst[i], 0, sm_info, cpoly->nb_polys);
        }
    } else { // rank 0 receives and prints. (round Robin again)
        FILE * out = outfile ? fopen(outfile, "w") : stdout;
        ASSERT_ALWAYS(out != nullptr);
        int nsm_total = 0;
        for (int side = 0; side < cpoly->nb_polys; side++) {
            nsm_total += sm_info[side].nsm;
        }
        fmt::print(out, "{} {} {}\n", nb_relsets, nsm_total, ell);
        std::vector<cxx_mpz_poly> res(cpoly->nb_polys);
        for (uint64_t i = 0; i < nb_parts; ++i) {
            print_all_sm(out, sm_info, cpoly->nb_polys, dst[i]);
            for (int j = 1; j < size; ++j) {
                if (i * size + j < nb_relsets) {
                    MPI_Recv_res(res, j, sm_info, cpoly->nb_polys);
                    print_all_sm(out, sm_info, cpoly->nb_polys, res);
                }
            }
        }
        if (outfile)
            fclose(out);
    }

    MPI_Finalize();

    return 0;
}
