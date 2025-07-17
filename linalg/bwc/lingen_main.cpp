#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <sys/utsname.h>

#include "fmt/base.h"

#include "gmp_aux.h"
#include "bw-common.h"
#include "arith-hard.hpp"
#include "lingen.hpp"
#include "lingen_bmstatus.hpp"
#include "lingen_bw_dimensions.hpp"
#include "lingen_checkpoints.hpp"
#include "lingen_expected_pi_length.hpp"
#include "lingen_io_matpoly.hpp"
#include "lingen_io_wrappers.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_tuning.hpp"
#include "logline.hpp"
#include "memusage.h"
#include "misc.h"
#include "params.h"
#include "portability.h"

/* If non-zero, then reading from A is actually replaced by reading from
 * a random generator */
static unsigned int input_length = 0;
static unsigned int random_input_length = 0;

static int allow_zero_on_rhs = 0;

static int global_flag_ascii = 0;
static int global_flag_tune = 0;
static int split_input_file = 0;  /* unsupported ; do acollect by ourselves */
static int split_output_file = 0; /* do split by ourselves */
static int rank0_exit_code = EXIT_SUCCESS;

static cxx_gmp_randstate rstate;


static void lingen_decl_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "ascii",
            "read and write data in ascii");
    param_list_decl_usage(pl, "timings",
            "provide timings on all output lines");
    param_list_decl_usage(pl, "tune",
            "activate tuning mode");
    param_list_decl_usage(pl, "allow_zero_on_rhs",
            "do not cry if the generator corresponds to a zero contribution on the RHS vectors");

    /* we must be square ! */
    param_list_decl_usage(pl, "mpi", "number of MPI nodes across which the execution will span, with mesh dimensions");
    param_list_decl_usage(pl, "thr", "number of threads (on each node) for the program, with mesh dimensions");

    param_list_decl_usage(pl, "nrhs",
            "number of columns that correspond to rhs vectors");
    param_list_decl_usage(pl, "rhs",
            "file with rhs vectors (only the header is read)");

    param_list_decl_usage(pl, "afile",
            "input sequence file");
    param_list_decl_usage(pl, "input_length",
            "input sequence length (defaults to auto-detect)");
    param_list_decl_usage(pl, "random-input-with-length",
            "use surrogate for input");
    param_list_decl_usage(pl, "split-input-file",
            "work with split files on input");
    param_list_decl_usage(pl, "split-output-file",
            "work with split files on output");
    param_list_decl_usage(pl, "random_seed",
            "seed the random generator");
    param_list_decl_usage(pl, "ffile",
            "output generator file");

#if 0
    param_list_decl_usage(pl, "lingen_mpi_threshold",
            "use MPI matrix operations above this size");
    param_list_decl_usage(pl, "lingen_threshold",
            "use recursive algorithm above this size");
#endif

    param_list_configure_switch(pl, "--tune", &global_flag_tune);
    param_list_configure_switch(pl, "--ascii", &global_flag_ascii);
    param_list_configure_alias(pl, "seed", "random_seed");

}/*}}}*/

template<bool is_binary>
static unsigned int count_lucky_columns(bmstatus<is_binary> & bm)/*{{{*/
{
    auto & d = bm.d;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const b = m + n;
    int const luck_mini = expected_pi_length(d);
    MPI_Bcast(bm.lucky.data(), b, MPI_UNSIGNED, 0, bm.com[0]);
    unsigned int nlucky = 0;
    for(unsigned int j = 0 ; j < b ; nlucky += bm.lucky[j++] >= luck_mini) ;
    return nlucky;
}/*}}}*/

    template<bool is_binary>
static int check_luck_condition(bmstatus<is_binary> & bm)/*{{{*/
{
    auto & d = bm.d;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const nlucky = count_lucky_columns(bm);

    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    if (!rank) {
        fmt::print("Number of lucky columns: {} ({} wanted)\n", nlucky, n);
    }

    if (nlucky == n)
        return 1;

    if (!rank) {
        fmt::print(stderr, "Could not find the required set of solutions (nlucky={})\n", nlucky);
    }
    if (random_input_length) {
        static int once=0;
        if (once++) {
            if (!rank) {
                fmt::print(stderr, "Solution-faking loop crashed\n");
            }
            MPI_Abort(bm.com[0], EXIT_FAILURE);
        }
        if (!rank) {
            fmt::print("Random input: faking successful computation\n");
        }
        for(unsigned int j = 0 ; j < n ; j++) {
            unsigned int const s = (j * 1009) % (m+n);
            bm.lucky[s]  = expected_pi_length(d);
            bm.delta[s] -= expected_pi_length(d);
        }
        return check_luck_condition(bm);
    }

    return 0;
}/*}}}*/

static void print_node_assignment(MPI_Comm comm)/*{{{*/
{
    int rank;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    struct utsname me[1];
    int const rc = uname(me);
    if (rc < 0) { perror("uname"); MPI_Abort(comm, 1); }
    size_t const sz = 1 + sizeof(me->nodename);
    std::unique_ptr<char[]> global(new char[size * sz]);
    std::fill_n(global.get(), size * sz, 0);
    memcpy(global.get() + rank * sz, me->nodename, sizeof(me->nodename));

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            global.get(), sz, MPI_BYTE, comm);
    if (rank == 0) {
        char name[80];
        int len=80;
        MPI_Comm_get_name(comm, name, &len);
        name[79]=0;
        for(int i = 0 ; i < size ; i++) {
            fmt::print("# {} rank {}: {}\n", name, i, global.get() + i * sz);
        }
    }
}/*}}}*/

/* Counting memory usage in the recursive algorithm.
 *
 * The recursive algorithm is designed to allow the allocated memory for
 * the input to be reused for placing the output. Some memory might have
 * been saved by upper layers. We also have some local allocation.
 *
 * Notations: The algorithm starts at depth 0 with an
 * input length L, and the notation \ell_i denotes L/2^(i+1). We have
 * \ell_i=2\ell_{i+1}. The notation \alpha denotes m/(m+n). Note that the
 * input has size \alpha*(1-\alpha)*L times (m+n)^2*\log_2(p) (divided by
 * r^2 if relevant).
 *
 * We define five quantities. All are understood as multiples of
 * (m+n)^2*\log_2(p).
 *
 * MP(i) is the extra storage needed for the MP operation at depth i.
 *
 * MUL(i) is the extra storage needed for the MUL operation at depth i.
 *
 * IO(i) is the common size of the input and output data of the call at
 *       depth i. We have
 *              IO(i) = 2\alpha\ell_i
 *
 * ST(i) is the storage *at all levels above the current one* (i.e. with
 *    depth strictly less than i) for the data that is still live and
 *    need to exist until after we return. This count is maximized in the
 *    leftmost branch, where chopped E at all levels must be kept.
 *    chopped E at depth i (not counted in ST(i) !) is:
 *          \alpha(1+\alpha) \ell_i
 *    (counted as the degree it takes to make the necessary data that
 *    we want to use to compute E_right),
 *    so the cumulated cost above is twice the depth 0 value, minus the
 *    depth i value, i.e.
 *              ST(i) = \alpha(1+\alpha)(L-2\ell_i).
 * SP(i) is the "spike" at depth i: not counting allocation that is
 *    already reserved for IO or ST, this is the amount of extra memory
 *    that is required by the call at depth i. We have:
 *      SP(i) = max {
 *              \alpha\ell_i,
 *              \alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+SP(i+1),
 *             2\alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *             2\alpha\ell_i-2\alpha\ell_i+SP(i+1)
 *             4\alpha\ell_i-2\alpha\ell_i+MUL(i)
 *             }
 *            = max {
 *              \alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+SP(i+1),
 *             2\alpha\ell_i-2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *             4\alpha\ell_i-2\alpha\ell_i+MUL(i)
 *                           }
 * 
 * Combining this together, and using
 * ST(i)+\alpha(1+\alpha)\ell_i=ST(i+1), we have:
 *
 * IO(i)+ST(i)+SP(i) = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              ST(i) + 2\alpha\ell_i+\alpha(1+\alpha)\ell_i+MP(i),
 *              ST(i) + 4\alpha\ell_i+MUL(i)
 *                      }
 *                   = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              \alpha(1+\alpha)(L-\ell_i)  + 2\alpha\ell_i + MP(i),
 *              \alpha(1+\alpha)(L-2\ell_i) + 4\alpha\ell_i + MUL(i)
 *                      }
 *                   = max {
 *              IO(i+1)+ST(i+1)+SP(i+1),
 *              \alpha((1+\alpha)L+(1-\alpha)\ell_i) + MP(i),
 *              \alpha((1+\alpha)L+2(1-\alpha)\ell_i) + MUL(i),
 *                      }
 *
 * Let RMP(i) be the amount of memory that is reserved while we are doing
 * the MP operation, and define RMUL similarly. We have:
 *      RMP(i)  = \alpha(1+\alpha)(L-\ell_i)  + 2\alpha\ell_i
 *      RMUL(i) = \alpha(1+\alpha)(L-2\ell_i) + 4\alpha\ell_i
 * whence:
 *      RMP(i) = \alpha((1+\alpha)L+(1-\alpha)\ell_i)
 *      RMUL(i) = \alpha((1+\alpha)L+2(1-\alpha)\ell_i)
 *
 * We have RMP(i) <= RMUL(i) <= RMP(0) <= RMUL(0) = 2\alpha*L. We'll use
 * the un-simplified expression later.
 *
 * Furthermore IO(infinity)=SP(infinity)=0, and ST(infinity)=\alpha(1+\alpha)L
 *
 * So that eventually, the amount of reserved memory for the whole
 * algorithm is RMUL(0)=2\alpha*L (which is 2/(1-\alpha)=2*(1+m/n) times
 * the input size). On top of that we have the memory required
 * for the transforms.
 *
 *
 * When going MPI, matrices may be rounded with some inaccuracy.
 * Splitting in two a 3x3 matrix leads to a 2x2 chunk, which is 1.77
 * times more than the simplistic proportionality rule.
 *
 * Therefore it makes sense to distinguish between matrices of size
 * m*(m+n) and (m+n)*(m+n). If we recompute RMUL(i) by taking this into
 * account, we obtain:
 *      [m/r][(m+n)/r][(1+\alpha)(L-2\ell_i)] + [(m+n)/r]^2*[4\alpha\ell_i]
 * where we only paid attention to the rounding issues with dimensions,
 * as those are more important than for degrees. Bottom line, the max is
 * expected to be for i=0, and that will be made only of pi matrices.
 */

/* Some of the early reading must be done before we even start, since
 * the code that we run depends on the input size.
 */

/* We don't have a header file for this one */
extern "C" void check_for_mpi_problems();

template<typename arith_hard_type, bool is_binary = arith_hard_type::is_binary>
static size_t K_elts_to_bytes(arith_hard_type const & ab, size_t x)
{
    if constexpr (is_binary)
        return iceildiv((x),ULONG_BITS) * sizeof(unsigned long);
    else
        return ab.vec_elt_stride(x);
}



template<bool is_binary>
struct blanket_allowances;
#ifndef LINGEN_BINARY
template<>
struct blanket_allowances<false> {
    matpoly_ft<fft_transform_info>::memory_guard blanket_ft { SIZE_MAX };
};
#else
template<>
struct blanket_allowances<true> {
    matpoly_ft<gf2x_cantor_fft_info>::memory_guard blanket_ft { SIZE_MAX };
    matpoly_ft<gf2x_ternary_fft_info>::memory_guard blanket_ft2 { SIZE_MAX };
};
#endif

template<bool is_binary>
static int wrapped_main(int argc, char const *argv[])
{
    cxx_param_list pl;

    bw_common_decl_usage(pl);
    lingen_decl_usage(pl);
    logline_decl_usage(pl);
    lingen_tuning_decl_usage(pl);
    lingen_checkpoint<is_binary>::decl_usage(pl);
    lingen_io_matpoly<is_binary>::decl_usage(pl);
    tree_stats::declare_usage(pl);

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);
    /* {{{ interpret our parameters */

    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    param_list_parse_int(pl, "allow_zero_on_rhs", &allow_zero_on_rhs);
    param_list_parse_uint(pl, "random-input-with-length", &random_input_length);
    param_list_parse_uint(pl, "input-length", &input_length);
    param_list_parse_int(pl, "split-output-file", &split_output_file);
    param_list_parse_int(pl, "split-input-file", &split_input_file);

    auto const afile = param_list_parse<std::string>(pl, "afile");

    if (bw->m == -1) {
	fmt::print(stderr, "no m value set\n");
	exit(EXIT_FAILURE);
    }
    if (bw->n == -1) {
	fmt::print(stderr, "no n value set\n");
	exit(EXIT_FAILURE);
    }
    if (!global_flag_tune && afile.empty() && !random_input_length) {
        fmt::print(stderr, "No afile provided\n");
        exit(EXIT_FAILURE);
    }

    /* we allow ffile and ffile to be both NULL */
    auto ffile = param_list_parse<std::string>(pl, "ffile");
    if (ffile.empty() && !afile.empty())
        ffile = afile + ".gen";
    ASSERT_ALWAYS(afile.empty() == ffile.empty());

    bmstatus<is_binary> bm(bw->m, bw->n, bw->p);

    const char * rhs_name = param_list_lookup_string(pl, "rhs");
    if (!global_flag_tune && !random_input_length) {
        if (!rhs_name) {
            fmt::print(stderr, "# When using lingen, you must either supply --random-input-with-length, or provide a rhs, or possibly provide rhs=none\n");
        } else if (strcmp(rhs_name, "none") == 0) {
            rhs_name = nullptr;
        }
    }
    if (param_list_parse_uint(pl, "nrhs", &(bm.d.nrhs)) && rhs_name) {
        fmt::print(stderr, "# the command line arguments rhs= and nrhs= are incompatible\n");
        exit(EXIT_FAILURE);
    }
    if (rhs_name && strcmp(rhs_name, "none") != 0) {
        if (!rank)
            get_rhs_file_header(rhs_name, nullptr, &(bm.d.nrhs), nullptr);
        MPI_Bcast(&bm.d.nrhs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    gmp_randseed_ui(rstate, bw->seed);

#if 0
    bm.lingen_threshold = 10;
    bm.lingen_mpi_threshold = 1000;
    param_list_parse_uint(pl, "lingen_threshold", &(bm.lingen_threshold));
    param_list_parse_uint(pl, "lingen_mpi_threshold", &(bm.lingen_mpi_threshold));
    if (bm.lingen_mpi_threshold < bm.lingen_threshold) {
        bm.lingen_mpi_threshold = bm.lingen_threshold;
        fmt::print(stderr, "Argument fixing: setting lingen_mpi_threshold={} (because lingen_threshold={})\n",
                bm.lingen_mpi_threshold, bm.lingen_threshold);
    }


#if defined(CADO_FAKEMPI_H)
    bm.lingen_mpi_threshold = UINT_MAX;
#endif
#endif

    /* }}} */

    /* TODO: we should rather use lingen_platform.
     */
    /* {{{ Parse MPI args. Make bm.com[0] a better mpi communicator */
    bm.mpi_dims[0] = 1;
    bm.mpi_dims[1] = 1;
    param_list_parse_intxint(pl, "mpi", bm.mpi_dims);
    {
        /* Display node index wrt MPI_COMM_WORLD */
        print_node_assignment(MPI_COMM_WORLD);

        /* Reorder all mpi nodes so that each node gets the given number
         * of jobs, but close together.
         */
        int const mpi[2] = { bm.mpi_dims[0], bm.mpi_dims[1], };
        int thr[2] = {1,1};
#ifdef  HAVE_OPENMP
        if (param_list_parse_intxint(pl, "thr", thr)) {
            if (omp_get_max_threads() >= thr[0] * thr[1]) {
                if (!rank)
                    fmt::print("# Limiting number of openmp threads to {}\n",
                            thr[0] * thr[1]);
                omp_set_num_threads(thr[0] * thr[1]);
            } else {
                if (!rank)
                    fmt::print("# Number of openmp threads is capped at {}"
                            ", which is below thr={}{}. Keeping as it is\n",
                            omp_get_max_threads(),
                            thr[0], thr[1]);
            }
        }
#else
        if (param_list_parse_intxint(pl, "thr", thr)) {
            if (thr[0]*thr[1] != 1) {
                if (!rank) {
                    fmt::print(stderr, "This program only wants openmp for multithreading. Ignoring thr argument.\n");
                }
                param_list_add_key(pl, "thr", "1x1", PARAMETER_FROM_CMDLINE);
            }
        }
#endif

#ifdef  CADO_FAKEMPI_H
        if (mpi[0]*mpi[1] > 1) {
            fmt::print(stderr, "non-trivial option mpi= can't be used with fakempi. Please do an MPI-enabled build (MPI=1)\n");
            exit(EXIT_FAILURE);
        }
#endif
        if (!rank)
            fmt::print("# size={} mpi={}{} thr={}{}\n", size, mpi[0], mpi[1], thr[0], thr[1]);
        ASSERT_ALWAYS(size == mpi[0] * mpi[1]);
        if (bm.mpi_dims[0] != bm.mpi_dims[1]) {
            if (!rank)
                fmt::print(stderr, "The current lingen code is limited to square splits ; here, we received a {} x {} split, which will not work\n",
                    bm.mpi_dims[0], bm.mpi_dims[1]);
            abort();
        }
        int const irank = rank / mpi[1];
        int const jrank = rank % mpi[1];
        bm.com[0] = MPI_COMM_WORLD;
        /* MPI Api has some very deprecated prototypes */
        MPI_Comm_set_name(bm.com[0], (char*) "world");

        char commname[32];
        snprintf(commname, sizeof(commname), "row%d\n", irank);
        MPI_Comm_split(MPI_COMM_WORLD, irank, jrank, &(bm.com[1]));
        MPI_Comm_set_name(bm.com[1], commname);

        snprintf(commname, sizeof(commname), "col%d\n", jrank);
        MPI_Comm_split(MPI_COMM_WORLD, jrank, irank, &(bm.com[2]));
        MPI_Comm_set_name(bm.com[2], commname);

        print_node_assignment(bm.com[0]);

        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        if ((bm.d.m + bm.d.n) / simd < (unsigned int) mpi[0]) {
            fmt::print("########################################################\n");
            fmt::print("# Warning: this run will leave some resources idle:\n"
                   "# the matrices of size {}*{} and {}*{} can be split into\n"
                   "# chunks of minimal size {}, whence an mpi split over {}*{} is useless\n",
                   bm.d.m,
                   bm.d.m + bm.d.n,
                   bm.d.m + bm.d.n,
                   bm.d.m + bm.d.n,
                   simd,
                   mpi[0],
                   mpi[0]);
            fmt::print("########################################################\n");
        }
    }
    /* }}} */

    /* lingen tuning accepts some arguments. We look them up so as to
     * avoid failures down the line */
    lingen_tuning_lookup_parameters(pl);
    
    tree_stats::interpret_parameters(pl);
    logline_interpret_parameters(pl);
    lingen_checkpoint<is_binary>::interpret_parameters(pl);
    lingen_io_matpoly<is_binary>::interpret_parameters(pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* TODO: read the a files in scattered mode */

    /* Don't police memory right now, we don't care */
    typename matpoly<is_binary>::memory_guard const main_memory(SIZE_MAX);

    std::unique_ptr<lingen_input_wrapper_base<is_binary>> A_series;

    if (random_input_length) {
        A_series.reset(new lingen_random_input<is_binary>(&bm.d.ab, bm.d.m, bm.d.n, rstate, random_input_length));
    } else {
        A_series.reset(new lingen_file_input<is_binary>(&bm.d.ab, bm.d.m, bm.d.n, afile, global_flag_ascii, input_length));
    }

    /* run the mpi problem detection only if we're certain that we're at
     * least close to the ballpark where this sort of checks make sense.
     */
    if (K_elts_to_bytes(bm.d.ab, (size_t) A_series->guessed_length() * (size_t) (bm.d.m + bm.d.n)) >= (1 << 28)) {
        check_for_mpi_problems();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* This will cause the initial read */
    std::unique_ptr<lingen_E_from_A<is_binary>> E_series(new lingen_E_from_A<is_binary>(bm.d, *A_series));

    bm.t = E_series->t0;

    size_t const L = E_series->guessed_length();

    {
        typename matpoly<is_binary>::memory_guard const blanket(SIZE_MAX);
        blanket_allowances<is_binary> const dummy;
        try {
            bm.hints = lingen_tuning(bm.d, L - bm.t, bm.com[0], pl);
        } catch (std::overflow_error const & e) {
            fputs(e.what(), stderr);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    if (global_flag_tune) {
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }

    size_t const safe_guess = global_flag_ascii ? ceil(1.05 * L) : L;

    /* c0 is (1+m/n) times the input size */
    size_t const c0 = K_elts_to_bytes(bm.d.ab,
                iceildiv(bm.d.m + bm.d.n, bm.mpi_dims[0]) *
                iceildiv(bm.d.m + bm.d.n, bm.mpi_dims[1]) *
                iceildiv(bm.d.m*safe_guess, bm.d.m+bm.d.n));
    typename matpoly<is_binary>::memory_guard const main(2*c0);

    if (!rank) {
        char buf[20];
        fmt::print("# Estimated memory for JUST transforms (per node): {}\n",
                size_disp(2*c0, buf));
        fmt::print("# Estimated peak total memory (per node): max at depth {}: {}\n",
                bm.hints.ipeak,
                size_disp(bm.hints.peak, buf));
    }

    int const go_mpi = bm.companion(0, L).go_mpi();

    if (go_mpi) {
        if (!rank) {
            if (size > 1) {
                fmt::print("Expected length {} exceeds MPI threshold,"
                       " going MPI now.\n",
                       L);
            } else {
                fmt::print("Expected length {} exceeds MPI threshold, "
                       "but the process is not running in an MPI context.\n",
                       L);
            }
        }
        MPI_Barrier(bm.com[0]);
    }

    std::unique_ptr<lingen_output_wrapper_base<is_binary>> Fdst;
    std::unique_ptr<lingen_output_wrapper_base<is_binary>> Fdst_rhs;
    
    if (random_input_length) {
        Fdst.reset(new lingen_output_to_sha1sum<is_binary>(& bm.d.ab, bm.d.n, bm.d.n, "F"));
        Fdst_rhs.reset(new lingen_output_to_sha1sum<is_binary>(& bm.d.ab, bm.d.nrhs, bm.d.n, "Frhs"));
    } else if (split_output_file) {
        std::string const & pattern = ffile;
        Fdst.reset(new lingen_output_to_splitfile<is_binary>(& bm.d.ab, bm.d.n, bm.d.n, pattern + ".sols{2}-{3}.{0}-{1}", global_flag_ascii));
        Fdst_rhs.reset(new lingen_output_to_splitfile<is_binary>(& bm.d.ab, bm.d.nrhs, bm.d.n, pattern + ".sols{2}-{3}.{0}-{1}.rhs", global_flag_ascii));
    } else {
        Fdst.reset(new lingen_output_to_singlefile<is_binary>(& bm.d.ab, bm.d.n, bm.d.n, ffile, global_flag_ascii));
        Fdst_rhs.reset(new lingen_output_to_singlefile<is_binary>(& bm.d.ab, bm.d.nrhs, bm.d.n, ffile + ".rhs", global_flag_ascii));
    }

    if (go_mpi && size > 1) {
        bigmatpoly_model const model(bm.com, bm.mpi_dims[0], bm.mpi_dims[1]);
        bigmatpoly<is_binary> E(& bm.d.ab, model, bm.d.m, bm.d.m + bm.d.n, safe_guess);
        lingen_scatter<bigmatpoly<is_binary>> fill_E(E);
        lingen_F0<is_binary> const & F0 = *E_series;
        pipe(*E_series, fill_E, "Read");
        bm.delta.assign(bm.d.m + bm.d.n, F0.t0);
        logline_init_timer();
        auto pi = bw_biglingen_collective(bm, E);
        bm.stats.final_print();
        bm.display_deltas();
        if (!rank) fmt::print("(pi.alloc = {})\n", pi.my_cell().capacity());
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        pi.zero_pad(simd * iceildiv(pi.get_size(), simd));
        if (check_luck_condition(bm)) {
            lingen_gather_reverse<bigmatpoly<is_binary>> read_PI(pi);
            lingen_F_from_PI Fsrc(bm, read_PI, F0);
            pipe(Fsrc, *Fdst, "Written", true);
            Fsrc.write_rhs(*Fdst_rhs);
        }
    } else if (!rank) {
        /* We do this only in the rank==0 case, since we have really
         * nothing to do at the other ranks.
         */

        /* We don't want to bother with memory problems in the non-mpi
         * case when the tuning was done for MPI: this is because the
         * per-transform ram was computed in the perspective of an MPI
         * run, and not for a plain run.
         */
        blanket_allowances<is_binary> const dummy;
        matpoly<is_binary> E(& bm.d.ab, bm.d.m, bm.d.m + bm.d.n, safe_guess);
        lingen_scatter<matpoly<is_binary>> fill_E(E);
        lingen_F0<is_binary> const & F0 = *E_series;
        pipe(*E_series, fill_E, "Read");
        bm.delta.assign(bm.d.m + bm.d.n, F0.t0);
        logline_init_timer();
        auto pi = bw_lingen_single(bm, E);
        bm.stats.final_print();
        bm.display_deltas();
        if (!rank) fmt::print("(pi.alloc = {})\n", pi.capacity());
        constexpr const unsigned int simd = is_binary ? ULONG_BITS : 1;
        pi.zero_pad(simd * iceildiv(pi.get_size(), simd));
        if (check_luck_condition(bm)) {
            lingen_gather_reverse<matpoly<is_binary>> read_PI(pi);
            lingen_F_from_PI<is_binary> Fsrc(bm, read_PI, F0);
            pipe(Fsrc, *Fdst, "Written", true);
            Fsrc.write_rhs(*Fdst_rhs);
        }
    }

    if (!rank && random_input_length) {
        fmt::print("t_basecase = {:.2f}\n", bm.t_basecase);
        fmt::print("t_mp = {:.2f}\n", bm.t_mp);
        fmt::print("t_mul = {:.2f}\n", bm.t_mul);
        fmt::print("t_cp_io = {:.2f}\n", bm.t_cp_io);
        auto const peakmem = PeakMemusage();
        if (peakmem > 0)
            fmt::print("# PeakMemusage (MB) = {} (VmPeak: can be misleading)\n", peakmem >> 10);
    }

    return 0;   // ignored.
}

/* We do this so that the dtors of the data that gets allocated within
 * main are allowed to use MPI_Comm_rank.
 */
// coverity[root_function]
int main(int argc, char const * argv[])
{
#ifdef  HAVE_OPENMP
    if (getenv("OMP_DYNAMIC") == NULL) {
        /* Change the default behavior with respect to dynamic thread
         * allocation, but do it with *lower* priority than the
         * environment variable (so that the possibility of changing the
         * behavior at runtime is retained).
         */
        omp_set_dynamic(true);
    }
#endif

    bw_common_init(bw, &argc, &argv);

#ifdef LINGEN_BINARY
    constexpr bool is_binary = true;
#else
    constexpr bool is_binary = false;
#endif
    wrapped_main<is_binary>(argc, argv);

    bw_common_clear(bw);
    return rank0_exit_code;
}

/* vim:set sw=4 sta et: */
