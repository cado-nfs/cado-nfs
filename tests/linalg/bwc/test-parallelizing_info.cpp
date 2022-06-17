#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "misc.h"
#include "select_mpi.h"
#include "parallelizing_info.hpp"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"

int verbose = 0;

void * test_code(parallelizing_info_ptr pi, param_list_ptr pl MAYBE_UNUSED, void * dummy MAYBE_UNUSED)
{
    serialize(pi->m);
    char * report_string;
    int rc = asprintf(&report_string, "J%uT%u\n", pi->m->jrank, pi->m->trank);
    ASSERT_ALWAYS(rc >= 0);
    size_t report_string_size = strlen(report_string) + 1;

        size_t max_report_size = 0;
        pi_allreduce(&report_string_size, &max_report_size, 1, BWC_PI_SIZE_T, BWC_PI_MAX, pi->m);
        void * all_reports = malloc(pi->m->totalsize * max_report_size);
        memset(all_reports, 0, pi->m->totalsize * max_report_size);
        memcpy(pointer_arith(all_reports, max_report_size * (pi->m->jrank * pi->m->ncores + pi->m->trank)), report_string, report_string_size);
        pi_allgather(NULL, 0, 0,
                all_reports, max_report_size, BWC_PI_BYTE, pi->m);

        if (max_report_size > 1 && pi->m->jrank == 0 && pi->m->trank == 0) {
            for(unsigned int j = 0 ; j < pi->m->njobs ; j++) {
                for(unsigned int t = 0 ; t < pi->m->ncores ; t++) {
                    char * locreport = (char *) pointer_arith(all_reports, max_report_size * (j * pi->m->ncores + t));
                    if (verbose) printf("##### J%uT%u timing report:\n%s",
                            j, t, locreport);
                    char * their;
                    rc = asprintf(&their, "J%uT%u\n", j, t);
                    ASSERT_ALWAYS(rc >= 0);
                    ASSERT_ALWAYS(strcmp(locreport, their) == 0);
                    free(their);
                }
            }
        }
        serialize(pi->m);
        free(all_reports);

    free(report_string);
    return NULL;
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    char ** argv0 = argv;

    MPI_Init(&argc, &argv);

    cxx_param_list pl;

    parallelizing_info_init();
    parallelizing_info_decl_usage(pl);
    param_list_decl_usage(pl, "v", "turn on some more logging");

    param_list_configure_switch(pl, "v", &verbose);

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    parallelizing_info_lookup_parameters(pl);

    param_list_warn_unused(pl);

    pi_go(test_code, pl, NULL);

    parallelizing_info_finish();
    MPI_Finalize();
}

