#include "cado.h"

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "bwc_config.h"
#include "select_mpi.h"
#include "bw-common-mpi.h"

int bw_common_init_mpi(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv)
{
    bw_common_init_defaults(bw);

#ifdef  MPI_LIBRARY_MT_CAPABLE
    int req = MPI_THREAD_MULTIPLE;
    int prov;
    MPI_Init_thread(p_argc, p_argv, req, &prov);
    if (req != prov) {
        fprintf(stderr, "Cannot init mpi with MPI_THREAD_MULTIPLE ;"
                " got %d != req %d\n",
                prov, req);
        exit(1);
    }
#else
    MPI_Init(p_argc, p_argv);
#endif
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bw->can_print = rank == 0 || getenv("CAN_PRINT");

    return bw_common_init_shared(bw, pl, p_argc, p_argv);
}
int bw_common_clear_mpi(struct bw_params * bw)
{
    bw_common_clear(bw);
    MPI_Finalize();
    return 0;
}
