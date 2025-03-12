#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>

#include <gmp.h>

#include "arith-generic.hpp"
#include "bw-common.h"
#include "gmp_aux.h"
#include "macros.h"
#include "matmul_top.hpp"
#include "matmul_top_vec.hpp"
#include "mmt_vector_pair.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "runtime_numeric_cast.hpp"
#include "select_mpi.h"
#include "timing.h"

static void * bench_cpu_prog(parallelizing_info_ptr pi, cxx_param_list & pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != nullptr;
    fake = fake || param_list_lookup_string(pl, "static_random_matrix") != nullptr;
    if (fake) bw->skip_online_checks = 1;
    int const tcan_print = bw->can_print && pi->m->trank == 0;

    int const ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        fprintf(stderr,
                "bench_cpu_bwc does not work in the interleaved setting\n");
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<arith_generic> const A(arith_generic::instance(bw->p, ys[1]-ys[0]));

    matmul_top_data mmt(A.get(), pi, pl, bw->dir);

    /* we allocate as many vectors as we have matrices, plus one if the
     * number of matrices is odd (so we always have an even number of
     * vectors). If the number of matrices is odd, then
     * the first vector may be shared.  Otherwise, I believe it cannot
     * (but I'm not really sure)
     *
     * Storage for vectors need actually not be present at all times.
     * This could be improved.
     */

    mmt_vector_pair ymy(mmt, bw->dir);

    /* I have absolutely no idea why, but the two --apparently useless--
     * serializing calls around the next block seem to have a beneficial
     * impact on the SEGv's we see every now and then with --mca
     * mpi_leave_pinned 1
     */
    serialize(pi->m);
    {
        cxx_gmp_randstate rstate;

        unsigned long const g = pi->m->jrank * pi->m->ncores + pi->m->trank;
        gmp_randseed_ui(rstate, bw->seed + g);
        mmt_vec_set_random_inconsistent(ymy[0], rstate);
        mmt_vec_truncate(mmt, ymy[0]);
    }
    serialize(pi->m);

    serialize_threads(pi->m);


    for(int streak = 1 ; streak < bw->interval ; streak <<= 1) {
        if (tcan_print)
            printf("Measuring timings for %d multiplications (cpu-bound)\n", streak);
        mmt_vec_twist(mmt, ymy[0]);

        serialize(pi->m);
        double timers[3] = {0,};
        timers[2] = -wct_seconds();
        thread_seconds_user_sys(timers);
        timers[0] *= -1;
        timers[1] *= -1;

        for(int i = 0 ; i < streak ; i++) {
            // matmul_top_mul(mmt, ymy, timing);
            {
                int const d = ymy[0].d;
                int const NM = runtime_numeric_cast<int>(mmt.matrices.size());
                int const nmats_odd = NM & 1;
                int midx = (d ? (NM - 1) : 0);
                for(int l = 0 ; l < NM ; l++) {
                    mmt_vec & src = ymy[l];
                    int const last = l == (NM - 1);
                    int const lnext = last && !nmats_odd ? 0 : (l+1);
                    mmt_vec & dst = ymy[lnext];

                    src.consistency = 2;
                    matmul_top_mul_cpu(mmt, midx, d, dst, src);


                    midx += d ? -1 : 1;
                }
            }
        }
        thread_seconds_user_sys(timers);
        timers[2] += wct_seconds();
        serialize(pi->m);

        char buf[40];
        snprintf(buf, 40, "%.2f@%.1f%% ", timers[2]/streak, 100.0*(timers[0]+timers[1])/timers[2]);
        grid_print(pi, buf, strlen(buf) + 1, tcan_print);

        double slowest = timers[2]/streak;
        double fastest = timers[2]/streak;
        pi_allreduce(nullptr, &slowest, 1, BWC_PI_DOUBLE, BWC_PI_MAX, pi->m);
        pi_allreduce(nullptr, &fastest, 1, BWC_PI_DOUBLE, BWC_PI_MIN, pi->m);
        if (tcan_print)
            printf("Cpu time spread: %.2f ... %.2f (delta = %.2f)\n",
                    fastest, slowest, slowest - fastest);

        const double max_measurement_time = 120;
        if (fastest * streak >= max_measurement_time) {
            if (tcan_print) {
                printf("Stopping measurement now, as we've spent more than %.0f s on timing %d rounds\n", max_measurement_time, streak);
            }
            break;
        }

        mmt_vec_untwist(mmt, ymy[0]);

        serialize(pi->m);
    }

    if (tcan_print) {
        printf("Done bench.\n");
    }
    serialize(pi->m);

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
    if (bw->ys[0] < 0) { fprintf(stderr, "no ys value set\n"); exit(1); }

    ASSERT_ALWAYS(param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(!param_list_lookup_string(pl, "solutions"));

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(bench_cpu_prog, pl, nullptr);

    parallelizing_info_finish();
    
    bw_common_clear(bw);

    return 0;
}

