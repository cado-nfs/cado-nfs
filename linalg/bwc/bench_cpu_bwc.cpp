#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include "rolling.h"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "fmt/core.h"            // for check_format_string
#include "fmt/printf.h" // fmt::fprintf // IWYU pragma: keep
#include "fmt/format.h"
#include "macros.h"
#include "matmul_top_vec.hpp"
using namespace fmt::literals;

double
wct_seconds (void)
{
    struct timeval tv[1];
    gettimeofday (tv, NULL);
    return (double)tv->tv_sec + (double)tv->tv_usec*1.0e-6;
}

void thread_seconds_user_sys(double * res, double m)
{
    struct rusage ru[1];
#ifdef HAVE_RUSAGE_THREAD
    getrusage(RUSAGE_THREAD, ru);
#else
#error "implement me"
#endif
    res[0] += m * ((double)ru->ru_utime.tv_sec + (double)ru->ru_utime.tv_usec/1.0e6);
    res[1] += m * ((double)ru->ru_stime.tv_sec + (double)ru->ru_stime.tv_usec/1.0e6);
}

void * bench_cpu_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    int fake = param_list_lookup_string(pl, "random_matrix") != NULL;
    fake = fake || param_list_lookup_string(pl, "static_random_matrix") != NULL;
    if (fake) bw->skip_online_checks = 1;
    int tcan_print = bw->can_print && pi->m->trank == 0;

    int ys[2] = { bw->ys[0], bw->ys[1], };
    if (pi->interleaved) {
        fprintf(stderr,
                "bench_cpu_bwc does not work in the interleaved setting\n");
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<arith_generic> A(arith_generic::instance(bw->p, ys[1]-ys[0]));
    block_control_signals();

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

    int nmats_odd = mmt.nmatrices & 1;

    mmt_vec * ymy = new mmt_vec[mmt.nmatrices + nmats_odd];
    matmul_top_matrix_ptr mptr;
    mptr = (matmul_top_matrix_ptr) mmt.matrices + (bw->dir ? (mmt.nmatrices - 1) : 0);
    for(int i = 0 ; i < mmt.nmatrices ; i++) {
        int shared = (i == 0) & nmats_odd;
        mmt_vec_setup(ymy[i], mmt,0,0, bw->dir ^ (i&1), shared, mptr->n[bw->dir]);
        mmt_full_vec_set_zero(ymy[i]);

        mptr += bw->dir ? -1 : 1;
    }
    if (nmats_odd) {
        mmt_vec_setup(ymy[mmt.nmatrices], mmt,0,0, !bw->dir, 0, mmt.matrices[0]->n[bw->dir]);
        mmt_full_vec_set_zero(ymy[mmt.nmatrices]);
    }

    /* I have absolutely no idea why, but the two --apparently useless--
     * serializing calls around the next block seem to have a beneficial
     * impact on the SEGv's we see every now and then with --mca
     * mpi_leave_pinned 1
     */
    serialize(pi->m);
    {
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
        unsigned long g = pi->m->jrank * pi->m->ncores + pi->m->trank;
        gmp_randseed_ui(rstate, bw->seed + g);
        mmt_vec_set_random_inconsistent(ymy[0], rstate);
        mmt_vec_truncate(mmt, ymy[0]);
        gmp_randclear(rstate);
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
        thread_seconds_user_sys(timers, -1);

        for(int i = 0 ; i < streak ; i++) {
            // matmul_top_mul(mmt, ymy, timing);
            {
                int d = ymy[0].d;
                int nmats_odd = mmt.nmatrices & 1;
                int midx = (d ? (mmt.nmatrices - 1) : 0);
                for(int l = 0 ; l < mmt.nmatrices ; l++) {
                    mmt_vec & src = ymy[l];
                    int last = l == (mmt.nmatrices - 1);
                    int lnext = last && !nmats_odd ? 0 : (l+1);
                    mmt_vec & dst = ymy[lnext];

                    src.consistency = 2;
                    matmul_top_mul_cpu(mmt, midx, d, dst, src);


                    midx += d ? -1 : 1;
                }
            }
        }
        thread_seconds_user_sys(timers, 1);
        timers[2] += wct_seconds();
        serialize(pi->m);

        char buf[40];
        snprintf(buf, 40, "%.2f@%.1f%% ", timers[2]/streak, 100.0*(timers[0]+timers[1])/timers[2]);
        grid_print(pi, buf, strlen(buf) + 1, tcan_print);

        double slowest = timers[2]/streak;
        double fastest = timers[2]/streak;
        pi_allreduce(NULL, &slowest, 1, BWC_PI_DOUBLE, BWC_PI_MAX, pi->m);
        pi_allreduce(NULL, &fastest, 1, BWC_PI_DOUBLE, BWC_PI_MIN, pi->m);
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

    delete[] ymy;

    return NULL;
}

// coverity[root_function]
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

    ASSERT_ALWAYS(param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(!param_list_lookup_string(pl, "solutions"));

    catch_control_signals();

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(bench_cpu_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}

