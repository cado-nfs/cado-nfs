/* This is (almost) standalone */
#define _POSIX_C_SOURCE 200112L

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_UTSNAME_H
#include <sys/utsname.h>
#endif
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#ifndef STANDALONE
#include "select_mpi.h"
#else
#include <mpi.h>
#define ASSERT_ALWAYS(x) do { if (!(x)) abort(); } while (0)
#endif

int rank;
int size;

#define BUFFERSIZE      (1 << 28)

void share_file(const char * fname, int root, size_t total, MPI_Comm comm)
{
    FILE * f;
    f = fopen(fname, rank == root ? "r" : "w");
    if (f == NULL)
        MPI_Abort(comm,1);

    char * buf = malloc(BUFFERSIZE);
    size_t fsz = 0;
    time_t t0 = time(NULL);
    time_t t1 = t0 + 1;
    for( ; ; ) {
        int n = 0;
        if (rank == root)
            n = fread(buf, 1, BUFFERSIZE, f);

        MPI_Bcast(&n, 1, MPI_INT, root, comm);
        fsz += n;
        if (n == 0)
            break;
        MPI_Bcast(buf, n, MPI_BYTE, root, comm);
        if (rank != root) {
            int m=fwrite(buf, 1, n, f);
            if (m != n)
                abort();
        } else {
            time_t t = time(NULL);
            if (t >= t1) {
                printf("%.1f MB in %d s (%.2f MB/s) [%.1f%%]\n",
                        fsz * 1.0e-6,
                        (int) (t-t0), (double) fsz * 1.0e-6 / (t-t0),
                        100.0 * (double) fsz/ total
                        );
                t1 = t + 1;
            }
        }
    }
    fclose(f);
    int t = time(NULL) - t0;
    if (rank == 0)
        printf(" broadcasted in %d s (%.2f MB/s)\n",
                t, (double) fsz * 1.0e-6 / t);
    free(buf);
}

int main(int argc, char * argv[])
{
    struct stat sbuf[1];
#ifdef HAVE_UTSNAME_H
    struct utsname u[1];
    char *nodename = u[0].nodename;
#else
    char *nodename = "Unknown node name";
#endif

#ifdef HAVE_UTSNAME_H
    uname(u);
#endif
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("node %d/%d on %s\n", rank,size,nodename);

    int duplicate=0;
    size_t minname=strlen(nodename) + 1;
    size_t maxname=strlen(nodename) + 1;
    MPI_Allreduce(MPI_IN_PLACE, &maxname, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &minname, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    ASSERT_ALWAYS(minname == maxname);
    char * allnames = malloc(size * maxname);
    MPI_Allgather(nodename, maxname, MPI_BYTE, allnames, maxname, MPI_BYTE, MPI_COMM_WORLD);
    for(int i = 0 ; i < rank ; i++) {
        if (memcmp(allnames + i * maxname, nodename, maxname) == 0) {
            fprintf(stderr, "%s on node %d/%d duplicates node %d/%d\n",
                    nodename, rank, size, i, size);
            duplicate=1;
        }
    }

    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, duplicate, rank, &comm);

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if (!duplicate)
    printf("reduced node %d/%d on %s\n", rank,size,nodename);

    for(int i = 1 ; i < argc && !duplicate; i++) {
        int rc;
        memset(sbuf, 0, sizeof(struct stat));
        rc = stat(argv[i], sbuf);
        if (rc < 0 && errno != ENOENT)
            abort();
        int ok = rc == 0;
        size_t * allsizes = malloc(size * sizeof(size_t));
        size_t me = sbuf->st_size;
        MPI_Allgather(&me, sizeof(size_t), MPI_BYTE,
                allsizes, sizeof(size_t), MPI_BYTE, comm);
        size_t szmax = 0;
        for(int k = 0 ; k < size ; k++) {
            ok = ok && allsizes[k] == me;
            if (allsizes[k] > szmax)
                szmax = allsizes[k];
        }
        MPI_Allreduce(MPI_IN_PLACE, &ok,    1, MPI_INT, MPI_SUM, comm);
        if (ok) {
            if (rank == 0)
                printf("%s ok everywhere (%zu MB)\n", argv[i], szmax >> 20);
            continue;
        } else {
            if (rc == 0 && me < szmax) {
                printf("node %d/%d, %s is only %zu < %zu. Removed\n",
                        rank,size,argv[i],me,szmax);
                unlink(argv[i]);
                rc = -1;
            }
            ok = rc == 0;
            MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_SUM, comm);
            if (ok != 1) {
                if (rank == 0)
                    fprintf(stderr, "warning: %d<%d full files found for %s\n", ok, size, argv[i]);
                // MPI_Abort(comm,1);
            }
        }
        int node=0;
        if (rc == 0) node = rank;
        MPI_Allreduce(MPI_IN_PLACE, &node, 1, MPI_INT, MPI_MAX, comm);
        if (rank == 0) {
            printf("%s (%zu MB, node %d)\n", argv[i], szmax >> 20, node);
            fflush(stdout);
        }
        share_file(argv[i], node, szmax, comm);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_free(&comm);

    MPI_Finalize();
}

