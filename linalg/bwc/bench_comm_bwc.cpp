#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstdint>              // for uint32_t
#include <memory>
#include <string>                // for string, operator+
#include <gmp.h>                 // for gmp_randclear, gmp_randinit_default
#ifdef HAVE_RESOURCE_H
#include <sys/resource.h>	/* for cputime */
#endif
#include <sys/time.h>	/* for gettimeofday */
#include "matmul.hpp"              // for matmul_public_s
#include "parallelizing_info.hpp"
#include "matmul_top.hpp"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.hpp"
#include "bw-common.h"
#include "async.hpp"
#include "xdotprod.hpp"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "fmt/core.h"            // for check_format_string
#include "fmt/printf.h" // fmt::fprintf // IWYU pragma: keep
#include "fmt/format.h"
#include "macros.h"
#include "timing.h"

using namespace fmt::literals;

void * bench_comm_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != nullptr;
    fake = fake || param_list_lookup_string(pl, "static_random_matrix") != nullptr;
    if (fake) bw->skip_online_checks = 1;
    int const tcan_print = bw->can_print && pi->m->trank == 0;

    int const ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        fprintf(stderr, "bench_bwc does not work in the interleaved setting\n");
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<arith_generic> const A(arith_generic::instance(bw->p, ys[1]-ys[0]));
    block_control_signals();

    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    serialize(pi->m);
    matmul_top_comm_bench(mmt, bw->dir);

    if (tcan_print) {
        printf("Done bench.\n");
    }
    serialize(pi->m);

    return nullptr;
}

// coverity[root_function]
int main(int argc, char * argv[])
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
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    ASSERT_ALWAYS(param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(!param_list_lookup_string(pl, "solutions"));

    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(bench_comm_prog, pl, nullptr);

    parallelizing_info_finish();
    bw_common_clear(bw);

    return 0;
}

