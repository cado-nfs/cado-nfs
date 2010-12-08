#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "parallelizing_info.h"
#include "select_mpi.h"

#include "utils.h"
#include "mf.h"
#include "balancing.h"
#include "balancing_workhorse.h"

// FIXME -- the command line currently does not allow changing this.
int transposing = 1;

void usage()
{
    fprintf(stderr, "Usage: ./mf-dobal [options] --matrix <mfile> <bfile>\n"
    "This program is an MPI program. Collectively, nodes build a split\n"
    "version of the matrix found in file mfile, the balancing being\n"
    "computed according to the balancing file bfile.\n"
    "Options recognized:\n"
    "\t(none)\n");
    exit(1);
}

void * all(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    matrix_u32 mat;
    mat->transpose=1;
    mat->bfile = param_list_lookup_string(pl, "balancing");
    mat->mfile = param_list_lookup_string(pl, "matrix");
    balancing_get_matrix_u32(pi, pl, mat);
    return NULL;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    param_list pl;
    param_list_init(pl);
    int wild = 0;
    argv++, argc--;
    param_list_configure_knob(pl, "--transpose", &transposing);
    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv))
	    continue;
	if (argv[0][0] != '-' && wild == 0) {
            param_list_add_key(pl, "balancing", argv[0], PARAMETER_FROM_CMDLINE);
	    wild++, argv++, argc--;
	    continue;
	}
	fprintf(stderr, "Unknown option %s\n", argv[0]);
	exit(1);
    }

    // just as a reminder -- this is looked up from balancing.c
    param_list_lookup_string(pl, "balancing_use_auxfile");

    if (!param_list_lookup_string(pl, "balancing")) usage();
    pi_go(all, pl, 0);

    param_list_clear(pl);
    MPI_Finalize();
}
