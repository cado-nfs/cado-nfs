#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "select_mpi.h"
#include "misc.h"

long failed_offset = 0;

int check_if_large_allgather_is_ok(MPI_Comm comm0)
{
    /* This check is here to provoke early failure under openmpi, with
     * some collection of parameters.
     *
     * When the psm2 mtl is used, we have a limit on the maximum size of
     * messages, it seems. It seems totally backwards to me. Short of a
     * way to check ahead of time what the limit is, we'll simply try and
     * see if this is indeed the case.
     * 
     * Here is an excerpt of the error message
     Message size 5326413824 bigger than supported by PSM2 API. Max = 4294967296
     [node-1:04815] *** An error occurred in MPI_Allgather
     [node-1:04815] *** reported by process [3170238465,0]
     [node-1:04815] *** on communicator MPI_COMM_WORLD
     [node-1:04815] *** MPI_ERR_OTHER: known error not in list
     [node-1:04815] *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
     [node-1:04815] ***    and potentially your MPI job)
     [node-1.nancy.grid5000.fr:04793] 1 more process has sent help message help-mtl-psm2.txt / message too big

     * Are we safe with the OFI mtl ? not really. As it happens, this
     * same check fails as well, with the rather unpleasant specifics
     * that the MPI_Allgather call returns fine -- only the copied data
     * doesn't quite meet the expectations... :-((((
     *
     * It seems that to avoir this error, one has to use --mca pml ob1 .
     * At least on omnipath. I haven't checked on infiniband.
     */
    /*
     * This function returns:
     *  0 on success.
     *  a non-zero MPI Error code if MPI_Allgather returned one.
     *  -1 if no MPI Error code was returned, but the result of Allgather
     *  was wrong.
     *  -2 if memory allocation failed.
     *
     * (note that the MPI document guarantees that MPI error codes are
     * positive integers)
     */

    size_t s = 1 << 16;
    unsigned int items_at_a_time = (1 << 16) + 1;
    MPI_Datatype mpi_ft;
    MPI_Type_contiguous(s, MPI_BYTE, &mpi_ft);
    MPI_Type_commit(&mpi_ft);
    int comm_size0, comm_rank0;
    MPI_Comm_size(comm0, &comm_size0);
    MPI_Comm_rank(comm0, &comm_rank0);
    MPI_Comm comm;

    /* The test makes more sense for larger communicators, but as a
     * matter of fact it works already for only 2 jobs, and that gives us
     * more control on the allocated size.
     */
    const int number_of_nodes_for_test = 2;

    MPI_Comm_split(comm0, comm_rank0 < number_of_nodes_for_test, comm_rank0, &comm);
    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    int rc = 0;

    if (comm_rank0 < number_of_nodes_for_test) {
        MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
        void * data = malloc(items_at_a_time * comm_size * s);
        memset(data, 0, items_at_a_time * comm_size * s);
        int alloc_ok = data != NULL;
        MPI_Allreduce(MPI_IN_PLACE, &alloc_ok, 1, MPI_INT, MPI_MIN, comm);
        if (alloc_ok) {
            memset(((char*)data) + items_at_a_time * s * comm_rank, 0x42, items_at_a_time * s);
            rc = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                    data, items_at_a_time,
                    mpi_ft, comm);
            if (rc == 0) {
                void * p = memchr(data, 0, items_at_a_time * comm_size * s);
                if (p != NULL) {
                    /* We found a zero, we shouldn't ! */
                    rc = -1;
                    failed_offset = ((char*)p)-(char*)data;
                }
            }
            MPI_Barrier(comm);
        } else {
            rc = -2;
        }
        if (data) free(data);
        MPI_Type_free(&mpi_ft);
        MPI_Comm_free(&comm);
    }

    MPI_Barrier(comm0);
    return rc;
}

/* Theoretically ucx on omnipath should work, but it doesn't seem
 * obvious...
 *
 * At least this issue was purportedly fixed:
    https://github.com/openucx/ucx/issues/750
*/
const char * error_explanation_mpi[] = {
    "As a workaround, you might want to use another (mtl-type) layer.",
    "With OpenMPI, here is what can be done (it does not seem useful to compile an",
    "ucx-enabled version):",
    /*
       "Compile ucx-1.5.1 (I guess v1.5 can do):",
       "  sudo apt install libfabric-dev pkg-config librdmacm-dev",
       "  git clone https://github.com/openucx/ucx.git",
       "  cd ucx ; git co v1.5.1",
       "  ./autogen.sh",
       "  ./configure --prefix=/opt/ucx-1.5.1 && make -j64 && make install",
       "  cd ../",
       */
    "Compile openmpi-4.0.1 (I only tested this version):",
    "  tar xf openmpi-4.0.1.tar.bz2",
    "  cd openmpi-4.0.1",
    "  ./configure --prefix=/opt",
        // " --with-ucx=/opt/ucx-1.5.1",
    " \\",
    "    --disable-mpi-fortran --without-cuda --disable-opencl \\",
    "    && make -j64 && make install",
    "Compile using the software above (set all PATHs correctly first).",
    "Run as follows, depending on the hardware:",
    "  On Mellanox infiniband ConnectX-3:",
    "    mpiexec -n 4 --prefix /opt/openmpi-4.0.1/ \\",
    "      --mca btl_openib_allow_ib true \\",
    "      --mca btl openib \\",
    "      --mca mtl ^ofi \\",
    "      --mca pml ob1 \\",
    "      ./a.out",
    "  On Intel OmniPath:",
    "    mpiexec -n 4 --prefix /opt/openmpi-4.0.1/ \\",
    "      --mca btl_openib_allow_ib true \\",
    "      --mca btl openib \\",
    "      --mca mtl ^psm2 \\",
    "      --mca pml ob1 \\",
    "      ./a.out",
};


void check_large_allgather_or_abort(const char * prefix)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("%schecking if large MPI_Allgather (>4GB) works: ...\n", prefix);
    int err = check_if_large_allgather_is_ok(MPI_COMM_WORLD);
    if (!rank) {
        printf("%schecking if large MPI_Allgather (>4GB) works: %s\n", prefix,
                ok_NOK(err == 0));
        if (err != 0)
            fprintf(stderr, "%schecking if large MPI_Allgather (>4GB) works: %s\n", prefix,
                    ok_NOK(err == 0));
    }
    if (err == 0) return;

    if (err == -2) {
        fprintf(stderr, "%sCould not allocate memory buffer."
                " Proceeding anyway, since this is just a test,"
                " but the program will probably abort sooner"
                " or later anyway.\n", prefix);
        return;
    }

    int someone_has_minusone = (err == -1);
    MPI_Allreduce(MPI_IN_PLACE, &someone_has_minusone, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (someone_has_minusone) {
        long * offsets = malloc(size * sizeof(long));
        offsets[rank] = failed_offset;
        MPI_Gather(&failed_offset, 1, MPI_LONG,
                offsets, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        if (!rank) {
            for(int i = 0 ; i < size ; i++) {
                fprintf(stderr, "node %d failed_offset = 0x%lx\n", i, offsets[i]);
        }
        }
        free(offsets);
    }

    if (!rank) {
        if (err > 0) { /* return an MPI Error if we've got one. */
            /* MPI_ERR_OTHER... mostly useless */
            char error[1024];
            int errorlen = sizeof(error);
            MPI_Error_string(err, error, &errorlen);
            fprintf(stderr, "%sMPI error returned:\n%s%s\n",
                    prefix, prefix, error);
        }
        size_t s = sizeof(error_explanation_mpi)/sizeof(error_explanation_mpi[0]);
        fprintf(stderr, "%s%s\n%s%s\n",
                prefix,
                "A basic test of an MPI_Allgather with large messages (>4GB) failed.",
                prefix,
                "This could be due to the PSM2 layer having an API limit at 4GB per message."
               );
        for(size_t i = 0 ; i < s ; ++i) {
            fprintf(stderr, "%s%s\n", prefix, error_explanation_mpi[i]);
        }
    }
    abort();
}

int check_if_large_mpi_send_is_ok(MPI_Comm comm0)
{
    /* The test makes more sense for larger communicators, but as a
     * matter of fact it works already for only 2 jobs, and that gives us
     * more control on the allocated size.
     */
    const int number_of_nodes_for_test = 2;

    int comm_rank0;
    MPI_Comm_rank(comm0, &comm_rank0);

    MPI_Comm comm;
    MPI_Comm_split(comm0, comm_rank0 < number_of_nodes_for_test, comm_rank0, &comm);

    int err = 0;

    if (comm_rank0 < number_of_nodes_for_test) {
        int comm_rank;
        MPI_Comm_rank(comm, &comm_rank);

        size_t chunk = 3<<29;

        void * data = malloc(chunk);

        int alloc_ok = data != NULL;
        MPI_Allreduce(MPI_IN_PLACE, &alloc_ok, 1, MPI_INT, MPI_MIN, comm);
        if (alloc_ok) {
            memset(data, 0x42, chunk);
            if (comm_rank == 1) {
                err = MPI_Send(data, chunk, MPI_BYTE, 0, 0xbeef, comm);
            } else if (comm_rank == 0) {
                err = MPI_Recv(data, chunk, MPI_BYTE, 1, 0xbeef, comm, MPI_STATUS_IGNORE);
            }
            free(data);
        } else {
            err = -2;
        }
        MPI_Barrier(comm);
    }
    MPI_Barrier(comm0);

    return err;
}

void check_large_mpi_send_or_abort(const char * prefix)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("%schecking if large MPI_Send (>1GB) works: ...\n", prefix);
    int err = check_if_large_mpi_send_is_ok(MPI_COMM_WORLD);
    if (!rank) {
        printf("%schecking if large MPI_Send (>1GB) works: %s\n", prefix,
                ok_NOK(err == 0));
        if (err != 0)
            fprintf(stderr, "%schecking if large MPI_Send (>1GB) works: %s\n", prefix,
                    ok_NOK(err == 0));
    }
    if (err == 0) return;

    if (err == -2) {
        fprintf(stderr, "%sCould not allocate memory buffer."
                " Proceeding anyway, since this is just a test,"
                " but the program will probably abort sooner"
                " or later anyway.\n", prefix);
        return;
    }

    if (!rank) {
        /* MPI_ERR_OTHER... mostly useless */
        char error[1024];
        int errorlen = sizeof(error);
        MPI_Error_string(err, error, &errorlen);
        fprintf(stderr, "%sMPI error returned:\n%s%s\n",
                prefix, prefix, error);
        fprintf(stderr, "%s%s\n", prefix,
                "A basic test of an MPI_Send with large messages (>1GB) failed."
               );
        size_t s = sizeof(error_explanation_mpi)/sizeof(error_explanation_mpi[0]);
        for(size_t i = 0 ; i < s ; ++i) {
            fprintf(stderr, "%s%s\n", prefix, error_explanation_mpi[i]);
        }
    }
    abort();
}

void check_for_mpi_problems()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        printf("# Skipping MPI checks (only %d node)\n", size);
    } else {
        check_large_mpi_send_or_abort("# ");
        // check_large_allgather_or_abort("# ");
        check_large_allgather_or_abort("# ");
    }
}

#ifdef WANT_MAIN
int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    check_for_mpi_problems();
    MPI_Finalize();
}
#endif

