#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include "bwc_config.h"
#include "select_mpi.h"
#include "parallelizing_info.h"
#include "portability.h"
#include "macros.h"
#include "misc.h"
#include "verbose.h"


#include <sys/time.h>   // gettimeofday
#ifdef  HAVE_UTSNAME_H
#include <sys/utsname.h>
#endif

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
#include "cpubinding.h"
#endif  /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */

static inline void pi_comm_init_pthread_things(pi_comm_ptr w, const char * desc)
{
    struct pthread_things * res;

    res = malloc(sizeof(struct pthread_things));

    barrier_init(res->bh, w->ncores);
    my_pthread_barrier_init(res->b, NULL, w->ncores);
    my_pthread_mutex_init(res->m, NULL);
    res->desc = strdup(desc);

    w->th = res;
}

static inline void pi_comm_destroy_pthread_things(pi_comm_ptr w)
{
    barrier_destroy(w->th->bh);
    my_pthread_barrier_destroy(w->th->b);

    my_pthread_mutex_destroy(w->th->m);
    /* Beware ! Freeing mustn't happen more than once ! */
    free(w->th->desc);
    free(w->th);
    w->th = NULL;
}

static void print_several(unsigned int n1, unsigned int n2, char a, char b, unsigned int w)
{
    char * toto;
    toto = malloc(w+2);
    printf("%c",b);
    memset(toto, a, w+1);
    toto[w+1]='\0';
    for(unsigned int i = 0 ; i < n1 ; i++) {
        toto[w]=a;
        for(unsigned int j = 0 ; j < n2-1 ; j++) {
            printf("%s", toto);
        }
        toto[w]=b;
        printf("%s", toto);
    }
    printf("\n");
    free(toto);
}

struct pi_go_helper_s {
    void *(*fcn)(parallelizing_info_ptr, param_list pl, void*);
    parallelizing_info_ptr p;
    param_list_ptr pl;
    void * arg;
};

void * pi_go_helper_func(struct pi_go_helper_s * s)
{
    pi_interleaving_enter(s->p);
#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    cpubinding_do_pinning(s->p->cpubinding_info, s->p->wr[1]->trank, s->p->wr[0]->trank);
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */
    int debug = 0;
    param_list_parse_int(s->pl, "debug-parallel-bwc", &debug);

    if (debug) {
        pi_log_init(s->p->m);
        pi_log_init(s->p->wr[0]);
        pi_log_init(s->p->wr[1]);
    }

    void * ret = (s->fcn)(s->p, s->pl, s->arg);

    if (debug) {
        pi_log_clear(s->p->m);
        pi_log_clear(s->p->wr[0]);
        pi_log_clear(s->p->wr[1]);
    }

    pi_interleaving_flip(s->p);
    pi_interleaving_leave(s->p);
    return ret;
}

typedef void * (*pthread_callee_t)(void*);

/*
   void pi_errhandler(MPI_Comm * comm, int * err, ...)
   {
   int size;
   int rank;
   MPI_Comm_size(*comm, &size);
   MPI_Comm_rank(*comm, &rank);
   fprintf(stderr, "Fatal MPI error ;\n");
   fprintf(stderr, " job %d in a comm. of size %d ;\n", rank, size);
   if (err == NULL) {
   fprintf(stderr, " no error code.\n");
   } else {
   char buf[MPI_MAX_ERROR_STRING];
   int len = sizeof(buf);
   MPI_Error_string(*err, buf, &len);
   fprintf(stderr, " %s\n", buf);
   }
   abort();
   }
   */

void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print)
{
    pi_comm_ptr wr = pi->m;
    unsigned long maxsize = siz;
    pi_allreduce(NULL, &maxsize, 1, BWC_PI_UNSIGNED_LONG, BWC_PI_MAX, wr);

    char * strings = malloc(wr->totalsize * maxsize);
    memset(strings, 0, wr->totalsize * maxsize);

    int me = wr->jrank * wr->ncores + wr->trank;

    /* instead of doing memcpy, we align the stuff. */
    char * fmt;
    int rc = asprintf(&fmt, "%%-%lus", maxsize-1);
    ASSERT_ALWAYS(rc >= 0);
    snprintf(strings + me * maxsize, maxsize, fmt, buf);
    free(fmt);
    pi_allreduce(NULL, strings, wr->totalsize * maxsize, BWC_PI_BYTE, BWC_PI_BXOR, wr);

    char * ptr = strings;

    serialize(pi->m);

    if (print) {
        /* There's also the null byte counted in siz. So that makes one
         * less.
         */
        unsigned int nj0 = pi->wr[0]->njobs;
        unsigned int nj1 = pi->wr[1]->njobs;
        unsigned int nt0 = pi->wr[0]->ncores;
        unsigned int nt1 = pi->wr[1]->ncores;
        print_several(nj0, nt0, '-', '+', maxsize-1);
        for(unsigned int j1 = 0 ; j1 < nj1 ; j1++) {
            for(unsigned int t1 = 0 ; t1 < nt1 ; t1++) {
                printf("|");
                for(unsigned int j0 = 0 ; j0 < nj0 ; j0++) {
                    for(unsigned int t0 = 0 ; t0 < nt0 ; t0++) {
                        int fence = t0 == nt0 - 1;
                        printf("%s%c", ptr, fence ? '|' : ' ');
                        ptr += maxsize;
                    }
                }
                printf("\n");
            }
            print_several(nj0, nt0, '-', '+', maxsize-1);
        }
    }
    serialize(wr);
    free(strings);
}

static void get_node_number_and_prefix(parallelizing_info_ptr pi)
{
    int len = strlen(pi->nodename);
    int minlen = len;
    int maxlen = len;
    MPI_Allreduce(MPI_IN_PLACE, &minlen, 1, MPI_INT, MPI_MIN, pi->m->pals);
    MPI_Allreduce(MPI_IN_PLACE, &maxlen, 1, MPI_INT, MPI_MAX, pi->m->pals);
    int j;
    for(j = 0 ; j < minlen ; j++) {
        int c = pi->nodename[j];
        if (isdigit(c)) {
            c = -1;
        }
        int cmin = c;
        int cmax = c;
        MPI_Allreduce(MPI_IN_PLACE, &cmin, 1, MPI_INT, MPI_MIN, pi->m->pals);
        MPI_Allreduce(MPI_IN_PLACE, &cmax, 1, MPI_INT, MPI_MAX, pi->m->pals);
        if (cmin != cmax || cmin < 0)
            break;
    }
    memcpy(pi->nodeprefix, pi->nodename, j);
    pi->nodeprefix[j]='\0';
    if (j && pi->nodeprefix[j-1] == '-')
        pi->nodeprefix[j-1] = '\0';
    memset(pi->nodenumber_s, ' ', maxlen - j);
    pi->nodenumber_s[maxlen-j]='\0';
    memcpy(pi->nodenumber_s + maxlen-len, pi->nodename + j, len-j);
    if (j == maxlen) {
        // all nodes have the same name. Most probably because we're only
        // on one node. Arrange for the display to look vaguely right.
        memcpy(pi->nodeprefix, pi->nodename, PI_NAMELEN);
        memcpy(pi->nodenumber_s, pi->nodename, PI_NAMELEN);
    }
}

static void display_process_grid(parallelizing_info_ptr pi)
{
    char * all_node_ids = malloc(PI_NAMELEN * pi->m->njobs);
    char * all_node_ids2 = malloc(PI_NAMELEN * pi->m->njobs);
    typedef int cpair[2];
    cpair my_coords = { pi->wr[0]->jrank, pi->wr[1]->jrank, };
    cpair * all_cpairs = malloc(sizeof(cpair) * pi->m->njobs);

    MPI_Allgather(pi->nodenumber_s, PI_NAMELEN, MPI_BYTE,
            all_node_ids, PI_NAMELEN, MPI_BYTE, pi->m->pals);
    MPI_Allgather(my_coords, 2, MPI_INT, all_cpairs, 2, MPI_INT, pi->m->pals);

    for(unsigned int i = 0 ; i < pi->m->njobs ; i++) {
        int x = all_cpairs[i][0];
        int y = all_cpairs[i][1];
        unsigned int j = y * pi->wr[0]->njobs + x;
        if (j >= pi->m->njobs) {
            fprintf(stderr, "uh ?\n");
            MPI_Abort(pi->m->pals, 42);
        }
        memcpy(all_node_ids2 + j * PI_NAMELEN,
                all_node_ids + i * PI_NAMELEN,
                PI_NAMELEN);
    }

    if (pi->m->jrank == 0) {
        if (strlen(pi->nodeprefix)) {
            printf("%d nodes on %s\n", pi->m->njobs, pi->nodeprefix);
        } else {
            printf("%d nodes, no common name prefix\n", pi->m->njobs);
        }
        char * pad = malloc(PI_NAMELEN);
        memset(pad, ' ', PI_NAMELEN);
        int node_id_len = strlen(pi->nodenumber_s);
        pad[node_id_len]='\0';
        if (node_id_len) {
            if (pi->thr_orig[0]) {
                pad[node_id_len/2]='\'';
            } else {
                pad[node_id_len/2]='.';
            }
        }
        int mapping_error = 0;
        for(unsigned int i = 0 ; i < pi->wr[1]->njobs ; i++) {       // i == y
        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            for(unsigned int j = 0 ; j < pi->wr[0]->njobs ; j++) {       // j == x
            for(unsigned int jt = 0 ; jt < pi->wr[0]->ncores ; jt++) {       // jt == x
                char * what = pad;
                if (!it && !jt) {
                    what = all_node_ids2 + PI_NAMELEN * (i * pi->wr[0]->njobs + j); 
                }
                if (pi->thr_orig[0]) {
                    /* then we have it==jt==0, but the real leader is not
                     * this one. Check at least that we have the proper
                     * identity (same leader)
                     */
                    ASSERT_ALWAYS(!it && !jt);
                    unsigned int i0 = i - (i % pi->thr_orig[0]);
                    unsigned int j0 = j - (j % pi->thr_orig[1]);
                    const char * ref = all_node_ids2 + PI_NAMELEN * (i0 * pi->wr[0]->njobs + j0); 
                    if (strcmp(ref, what) != 0) {
                        mapping_error=1;
                    }
                    if (i != i0 || j != j0)
                        what = pad;
                }
                printf(" %s", what);
            }
            }
            printf("\n");
        }
        }
        if (mapping_error) {
            fprintf(stderr, "Error: we have not obtained the expected node mapping; the only_mpi=1 setting is VERY sensible to the process mapping chosen by the MPI implementation. Chances are you got it wrong. A working setup with openmpi-1.8.2 is a uniqueified host file, with --mca rmaps_base_mapping_policy slot (which is the default setting).\n");
            exit(1);
        }
        free(pad);
    }
    free(all_node_ids2);
    free(all_node_ids);
    free(all_cpairs);
    MPI_Barrier(pi->m->pals);
}

static void pi_init_mpilevel(parallelizing_info_ptr pi, param_list pl)
{
    int err;

    memset(pi, 0, sizeof(parallelizing_info));

#ifdef HAVE_UTSNAME_H
    {
        struct utsname u[1];
        uname(u);
        memcpy(pi->nodename, u->nodename, MIN(PI_NAMELEN, sizeof(u->nodename)));
        pi->nodename[PI_NAMELEN-1]='\0';
        char * p = strchr(pi->nodename, '.');
        if (p) *p='\0';
    }
#else
    strncpy(pi->nodename, "unknown", sizeof(pi->nodename));
#endif

#ifndef NDEBUG
    /* Must make sure that we have proper ASSERTS after MPI calls then. */
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif

    int mpi[2] = {1,1};
    int thr[2] = {1,1};
    int only_mpi = 0;

    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);
    param_list_parse_int(pl, "only_mpi", &only_mpi);

    unsigned int nhj = mpi[0];
    unsigned int nvj = mpi[1];
    unsigned int nhc = thr[0];
    unsigned int nvc = thr[1];

    pi->m->njobs = nhj * nvj;
    pi->m->ncores = nhc * nvc;
    pi->m->totalsize = pi->m->njobs * pi->m->ncores;
    // MPI_Errhandler my_eh;
    // MPI_Comm_create_errhandler(&pi_errhandler, &my_eh);
    // MPI_Comm_set_errhandler(MPI_COMM_WORLD, my_eh);

    if (only_mpi) {
        int grank;
        MPI_Comm_rank(MPI_COMM_WORLD, & grank);
        if (grank == 0)
            printf("Making %d independent single-threaded MPI jobs, instead of %d %d-thread jobs\n",
                pi->m->totalsize, pi->m->njobs, pi->m->ncores);
        int tt = pi->m->ncores;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int nnum = rank / tt;
        int nidx = rank % tt;
        int ni = nnum / nvj;
        int nj = nnum % nvj;
        int ti = nidx / nvc;
        int tj = nidx % nvc;
        int i = ni * nhc + ti;
        int j = nj * nvc + tj;
        int nrank = i * (nvj * nvc) + j;

        MPI_Comm_split(MPI_COMM_WORLD, 0, nrank, & pi->m->pals);
        pi->m->ncores = 1;
        pi->m->njobs = pi->m->totalsize;
        memcpy(pi->thr_orig, thr, 2*sizeof(int));
        mpi[0] *= thr[0]; thr[0]=1;
        mpi[1] *= thr[1]; thr[1]=1;
        nhj = mpi[0];
        nvj = mpi[1];
        nhc = thr[0];
        nvc = thr[1];
    } else {
        pi->m->pals = MPI_COMM_WORLD;
    }

    MPI_Comm_rank(pi->m->pals, (int*) & pi->m->jrank);

    int size;
    MPI_Comm_size(pi->m->pals, & size);
    if ((unsigned int) size != pi->m->njobs) {
        if (pi->m->jrank == 0) {
            fprintf(stderr, "Inconsistency -- exactly %u == %u * %u"
                    " MPI jobs are needed -- got %u\n",
                    pi->m->njobs, nhj, nvj, size);
        }
        exit(1);
    }

    pi->m->trank = 0;

    /* Init to the simple case */
    pi->wr[0]->njobs = nvj;
    pi->wr[0]->ncores = nvc;
    pi->wr[1]->njobs = nhj;
    pi->wr[1]->ncores = nhc;
    pi->wr[0]->totalsize = pi->wr[0]->njobs * pi->wr[0]->ncores;
    pi->wr[1]->totalsize = pi->wr[1]->njobs * pi->wr[1]->ncores;

    // Here we want the _common row number_ ; not to be confused with the
    // rank !
    unsigned int jcommon[2];
    jcommon[0] = pi->m->jrank / nvj;
    jcommon[1] = pi->m->jrank % nvj;

    pi->wr[0]->trank = 0;
    pi->wr[1]->trank = 0;

    for(int d = 0 ; d < 2 ; d++) {
        // A subgroup contains all process having the same jcommon value.
        // Therefore, we have as many horizontal MPI barriers set up as
        // one finds MPI job numbers across a column.
        err = MPI_Comm_split(pi->m->pals, jcommon[d],
                pi->m->jrank, &pi->wr[d]->pals);
        ASSERT_ALWAYS(!err);
        MPI_Comm_rank(pi->wr[d]->pals, (int*) &pi->wr[d]->jrank);
    }

    get_node_number_and_prefix(pi);
    display_process_grid(pi);

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    /* prepare the cpu binding messages, and print the unique messages we
     * receive */
    if (!verbose_enabled(CADO_VERBOSE_PRINT_BWC_CPUBINDING)) {
        pi->cpubinding_info = cpubinding_get_info(NULL, pl, thr);
    } else {
        char * cpubinding_messages;
        pi->cpubinding_info = cpubinding_get_info(&cpubinding_messages, pl, thr);
        int msgsize = 0;
        if (cpubinding_messages)
            msgsize = strlen(cpubinding_messages);
        MPI_Allreduce(MPI_IN_PLACE, &msgsize, 1, MPI_INT, MPI_MAX, pi->m->pals);
        if (msgsize == 0) {
            if (cpubinding_messages)
                free(cpubinding_messages);
        }
        msgsize++;
        int chunksize = PI_NAMELEN + msgsize;
        char * big_pool = malloc(pi->m->njobs * chunksize);
        memset(big_pool, 0, pi->m->njobs * chunksize);
        if (cpubinding_messages) {
            int rc;
            rc = strlcpy(big_pool + pi->m->jrank * chunksize, pi->nodename, PI_NAMELEN);
            ASSERT_ALWAYS(rc == (int) strlen(pi->nodename));
            rc = strlcpy(big_pool + pi->m->jrank * chunksize + PI_NAMELEN, cpubinding_messages, msgsize);
            ASSERT_ALWAYS(rc == (int) strlen(cpubinding_messages));
            free(cpubinding_messages);
        }

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                big_pool, chunksize, MPI_BYTE, pi->m->pals);

        if (pi->m->jrank == 0) {
            // const char * refnode = NULL;
            const char * ref = NULL;
            for(unsigned int i = 0 ; i < pi->m->njobs ; i++) {
                const char * node = big_pool + i * chunksize;
                const char * msg = node + PI_NAMELEN;
                if (!ref || strcmp(msg, ref) != 0) {
                    printf("cpubinding messages on %s:\n%s", node, msg);
                    ref = msg;
                    /*
                       refnode = node;
                       } else {
                       printf("cpubinding messages on %s: same as %s\n", node, refnode);
                       */
                }
            }
        }
        free(big_pool);
    }
#else
    if (param_list_lookup_string(pl, "cpubinding")) {
        if (pi->m->jrank == 0) {
#ifdef HAVE_HWLOC
                    printf("cpubinding: parameter ignored (no C++11)\n");
#elif defined(HAVE_CXX11)
                    printf("cpubinding: parameter ignored (no hwloc)\n");
#else
                    printf("cpubinding: parameter ignored (no hwloc, no C++11)\n");
#endif
        }
    }
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */
}

/* How do we build a rank in pi->m from two rank in pi->wr[inner] and
 * pi->wr->outer ?
 */
static int mrank_from_tworanks(parallelizing_info_ptr pi, int d, int ji, int jo)
{
    int jj[2];
    jj[d] = ji;
    jj[!d] = jo;
    return jj[0] + jj[1] * pi->wr[0]->njobs;
    /* d == 0:
     * ji == rank when reading horizontally
     * jo == rank when reading vertically
     */
}


static parallelizing_info *
pi_grid_init(parallelizing_info_ptr pi)
{
    // unsigned int nvj = pi->wr[0]->njobs;
    // unsigned int nhj = pi->wr[1]->njobs;
    unsigned int nvc = pi->wr[0]->ncores;
    unsigned int nhc = pi->wr[1]->ncores;

    /* used in several places for doing snprintf */
    char buf[20];

    // OK. So far we have something reasonable except for per-thread
    // setup. We have to first duplicate our structure so as to
    // specialize the pi things, and then set up agreed barriers.

    // the global barrier is not too much work.
    pi_comm_init_pthread_things(pi->m, "main");

    // column and row barriers are more tricky.
    parallelizing_info * grid;
    grid = (parallelizing_info *) malloc(nhc * nvc * sizeof(parallelizing_info));
    // setting to zero is important, since several fields are ensured as
    // null.
    memset(grid, 0, nhc * nvc * sizeof(parallelizing_info));

    char commname[64];
#ifndef MPI_LIBRARY_MT_CAPABLE
    MPI_Comm_set_name(pi->m->pals, "main");
    snprintf(commname, sizeof(commname), "rows%u-%u", 
            pi->wr[1]->jrank * pi->wr[1]->ncores,
            (pi->wr[1]->jrank + 1) * pi->wr[1]->ncores - 1);
    MPI_Comm_set_name(pi->wr[0]->pals, commname);
    snprintf(commname, sizeof(commname), "cols%u-%u", 
            pi->wr[0]->jrank * pi->wr[0]->ncores,
            (pi->wr[0]->jrank + 1) * pi->wr[0]->ncores - 1);
    MPI_Comm_set_name(pi->wr[1]->pals, commname);
#endif

    for(unsigned int k = 0 ; k < nhc * nvc ; k++) {
        parallelizing_info_ptr e = grid[k];
        memcpy(e, pi, sizeof(parallelizing_info));
        unsigned int k0 = k % nvc;
        unsigned int k1 = k / nvc;
        e->m->trank = k;
        e->wr[0]->trank = k0;
        e->wr[1]->trank = k1;
#ifdef  MPI_LIBRARY_MT_CAPABLE
        MPI_Comm_dup(pi->m->pals, &(e->m->pals));
        MPI_Comm_dup(pi->wr[0]->pals, &(e->wr[0]->pals));
        MPI_Comm_dup(pi->wr[1]->pals, &(e->wr[1]->pals));

        /* give a name */
        snprintf(commname, sizeof(commname), "main.%u", k);
        MPI_Comm_set_name(e->m->pals, commname);
        snprintf(commname, sizeof(commname), "row%u.%u", 
                e->wr[1]->jrank * e->wr[1]->ncores + e->wr[1]->trank,
                e->wr[0]->trank);
        MPI_Comm_set_name(e->wr[0]->pals, commname);
        snprintf(commname, sizeof(commname), "col%u.%u", 
                e->wr[0]->jrank * e->wr[0]->ncores + e->wr[0]->trank,
                e->wr[1]->trank);
        MPI_Comm_set_name(e->wr[1]->pals, commname);
#endif  /* MPI_LIBRARY_MT_CAPABLE */
    }

    // row barriers.
    // we've got pi->wr[1]->ncores cores working on separate rows (that's
    // the number of cores we encounter when walking down one column). So
    // each row leader is at index c * pi->wr[0]->ncores, because
    // pi->wr[0]->ncores is the number of cores working on separate
    // columns.
    //
    for(unsigned int c = 0 ; c < pi->wr[1]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "r%u", 
                pi->wr[1]->jrank * pi->wr[1]->ncores + c);
        unsigned int Nc = c * pi->wr[0]->ncores;
        pi_comm_ptr leader = grid[Nc]->wr[0];
        pi_comm_init_pthread_things(leader, buf);
        // replicate.
        for(unsigned int k = 0 ; k < pi->wr[0]->ncores ; k++) {
            grid[k + Nc]->wr[0]->th = leader->th;
        }
    }

    // column barriers.
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "c%u", 
                pi->wr[0]->jrank * pi->wr[0]->ncores + c);
        pi_comm_ptr leader = grid[c]->wr[1];
        pi_comm_init_pthread_things(leader, buf);
        // replicate.
        for(unsigned int k = 0 ; k < pi->wr[1]->ncores ; k++) {
            unsigned int Nk = k * pi->wr[0]->ncores;
            grid[c + Nk]->wr[1]->th = leader->th;
        }
    }

    return grid;
}

#if 0
/* TODO: rewrite this using grid_print instead */
static void
pi_grid_print_sketch(parallelizing_info_ptr pi, parallelizing_info * grid)
{
    char buf[20];
    for(unsigned int ij = 0 ; ij < pi->wr[1]->njobs ; ij++) {
        // note that at this point, we're not yet MT.
        int err = MPI_Barrier(pi->m->pals);
        ASSERT_ALWAYS(!err);

        if (pi->wr[1]->jrank != ij) continue;

        if (pi->m->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);

        for(unsigned int it = 0 ; it < pi->wr[1]->ncores ; it++) {
            if (pi->wr[0]->jrank == 0) printf("|");
            for(unsigned int jj = 0 ; jj < pi->wr[0]->njobs ; jj++) {
                err = MPI_Barrier(pi->wr[0]->pals);
                ASSERT_ALWAYS(!err);
                for(unsigned int jt = 0 ; jt < pi->wr[0]->ncores ; jt++) {
                    parallelizing_info_srcptr tpi;
                    tpi = grid[it*pi->wr[0]->ncores+jt];
                    snprintf(buf, sizeof(buf), "J%uT%u%s%s",
                            tpi->m->jrank,
                            tpi->m->trank,
                            tpi->wr[0]->th->desc,
                            tpi->wr[1]->th->desc);
                    err = MPI_Bcast(buf, sizeof(buf),
                            MPI_BYTE, jj, pi->wr[0]->pals);
                    ASSERT_ALWAYS(!err);
                    int fence = jt == pi->wr[0]->ncores - 1;
                    if (pi->wr[0]->jrank == 0)
                        printf("%-12s%c", buf, fence ? '|' : ' ');
                }
            }
            if (pi->wr[0]->jrank == 0) printf("\n");
        }
        if (pi->wr[1]->jrank == ij && pi->wr[0]->jrank == 0)
            print_several(pi->wr[0]->njobs,pi->wr[0]->ncores, '-', '+', 12);
    }
}
#endif

/* given a grid of processes which has beed set up with pi_grid_init,
 * start the processes.
 *
 * The ngrids argument is here to accomodate the interleaved case.
 */
static void pi_go_mt_now(
        parallelizing_info ** grids,
        unsigned int ngrids,
        unsigned int n,
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    my_pthread_t * tgrid;
    struct pi_go_helper_s * helper_grid;

    tgrid = (my_pthread_t *) malloc(n * ngrids * sizeof(my_pthread_t));
    helper_grid = malloc(n * ngrids * sizeof(struct pi_go_helper_s));

    for(unsigned int k = 0 ; k < n * ngrids ; k++) {
        helper_grid[k].fcn = fcn;
        helper_grid[k].pl = pl;
        helper_grid[k].arg = arg;
        helper_grid[k].p = (parallelizing_info_ptr) grids[k/n] + (k%n);
        my_pthread_create(tgrid+k, NULL,
                (pthread_callee_t) pi_go_helper_func, helper_grid+k);
    }

    // nothing done here. We're the main job.
    for(unsigned int k = 0 ; k < n * ngrids ; k++) {
        my_pthread_join(tgrid[k], NULL);
    }
    free(helper_grid);
    free(tgrid);
}

static void pi_grid_clear(parallelizing_info_ptr pi, parallelizing_info * grid)
{
    // destroy row barriers.
    for(unsigned int c = 0 ; c < pi->wr[1]->ncores ; c++) {
        unsigned int Nc = c * pi->wr[0]->ncores;
        pi_comm_destroy_pthread_things(grid[Nc]->wr[0]);
    }
    // destroy column barriers
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        pi_comm_destroy_pthread_things(grid[c]->wr[1]);
    }
#ifdef  MPI_LIBRARY_MT_CAPABLE
    for(unsigned int k = 0 ; k < pi->m->ncores ; k++) {
        MPI_Comm_free(&grid[k]->m->pals);
        MPI_Comm_free(&grid[k]->wr[0]->pals);
        MPI_Comm_free(&grid[k]->wr[1]->pals);
    }
#endif  /* MPI_LIBRARY_MT_CAPABLE */

    /* Don't do pi_log_clear, just like we haven't done pi_log_init. The
     * sub-threads may define it, in which case it's their responsibility
     * to clear the thing.
     */

    free(grid);
}

static void pi_clear_mpilevel(parallelizing_info_ptr pi)
{
#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    int thr[2] = { pi->wr[0]->ncores, pi->wr[1]->ncores };
    cpubinding_free_info(pi->cpubinding_info, thr);
#endif /* defined(HAVE_HWLOC) && defined(HAVE_CXX11) */

    pi_comm_destroy_pthread_things(pi->m);

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_free(&pi->wr[d]->pals);
    }
    if (pi->m->pals != MPI_COMM_WORLD)
        MPI_Comm_free(&pi->m->pals);

    // MPI_Errhandler_free(&my_eh);
}

#if 0
static void shout_going_mt(pi_comm_ptr m)
{
    int err = MPI_Barrier(m->pals);
    ASSERT_ALWAYS(!err);
    if (m->jrank == m->njobs - 1)
        printf("going multithread now\n");
    err = MPI_Barrier(m->pals);
    ASSERT_ALWAYS(!err);
}
#endif

static void pi_go_inner_not_interleaved(
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    parallelizing_info pi;
    parallelizing_info * grid;
    pi_init_mpilevel(pi, pl);
    grid = pi_grid_init(pi);
    // shout_going_mt(pi->m);
    // pi_grid_print_sketch(pi, grids[0]);
    pi_go_mt_now(&grid, 1, pi->m->ncores, fcn, pl, arg);
    pi_grid_clear(pi, grid);
    pi_clear_mpilevel(pi);
}

static void pi_go_inner_interleaved(
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    parallelizing_info pi[2];
    parallelizing_info * grids[2];

    pi_init_mpilevel(pi[0], pl);
    memcpy(pi[1], pi[0], sizeof(parallelizing_info));

    pi[0]->interleaved = malloc(sizeof(pi_interleaving));
    pi[0]->interleaved->idx = 0;

    pi[1]->interleaved = malloc(sizeof(pi_interleaving));
    pi[1]->interleaved->idx = 1;

    /* Now the whole point is that it's the _same_ barrier ! */
    my_pthread_barrier_t * b = malloc(sizeof(my_pthread_barrier_t));
    my_pthread_barrier_init(b, NULL, 2 * pi[0]->m->ncores);
    pi[0]->interleaved->b = b;
    pi[1]->interleaved->b = b;

    pi_dictionary d;
    pi_dictionary_init(d);
    pi[0]->dict = d;
    pi[1]->dict = d;

    grids[0] = pi_grid_init(pi[0]);
    grids[1] = pi_grid_init(pi[1]);

    // shout_going_mt(pi[0]->m);
    // pi_grid_print_sketch(pi, grids[0]);

    pi_go_mt_now(grids, 2, pi[0]->m->ncores, fcn, pl, arg);

    pi_grid_clear(pi[0], grids[0]);
    pi_grid_clear(pi[1], grids[1]);

    my_pthread_barrier_destroy(b);
    pi_dictionary_clear(d);

    free(pi[0]->interleaved);
    free(pi[1]->interleaved);

    pi_clear_mpilevel(pi[0]);
    /* pi[1] is a simple copy, it doesn't have to be cleared. */
}

void pi_dictionary_init(pi_dictionary_ptr d)
{
    my_pthread_rwlock_init(d->m, NULL);
    d->e = NULL;
}

void pi_dictionary_clear(pi_dictionary_ptr d)
{
    /* last grab of the lock -- altough it would be preferred if we
     * weren't mt of course ! */
    my_pthread_rwlock_wrlock(d->m);
    pi_dictionary_entry_ptr ne = d->e;
    for( ; ne ; ) {
        pi_dictionary_entry_ptr nne = ne->next;
        free(ne);
        ne = nne;
    }
    d->e = NULL;
    my_pthread_rwlock_unlock(d->m);
    my_pthread_rwlock_destroy(d->m);
}

void pi_store_generic(parallelizing_info_ptr pi, unsigned long key, unsigned long who, void * value)
{
    pi_dictionary_ptr d = pi->dict;
    ASSERT_ALWAYS(d != NULL);
    pi_dictionary_entry_ptr n = malloc(sizeof(pi_dictionary_entry));
    n->key = key;
    n->who = who;
    n->value = value;
    my_pthread_rwlock_wrlock(d->m);
    n->next = d->e;
    d->e = n;
    my_pthread_rwlock_unlock(d->m);
}

void * pi_load_generic(parallelizing_info_ptr pi, unsigned long key, unsigned long who)
{
    pi_dictionary_ptr d = pi->dict;
    ASSERT_ALWAYS(d != NULL);
    void * r = NULL;
    my_pthread_rwlock_rdlock(d->m);
    pi_dictionary_entry_ptr e = d->e;
    for( ; e ; e = e->next) {
        if (e->key == key && e->who == who) {
            r = e->value;
            break;
        }
    }
    my_pthread_rwlock_unlock(d->m);
    return r;
}

void parallelizing_info_decl_usage(param_list pl)/*{{{*/
{
    param_list_decl_usage(pl, "mpi", "number of MPI nodes across which the execution will span, with mesh dimensions");
    param_list_decl_usage(pl, "thr", "number of threads (on each node) for the program, with mesh dimensions");
    param_list_decl_usage(pl, "interleaving", "whether we should start two interleaved sets of threads. Only supported by some programs, has no effect on others.");
    param_list_decl_usage(pl, "only_mpi", "replace threads by distinct MPI jobs");
    param_list_decl_usage(pl, "debug-parallel-bwc", "enable heavy debugging messages for desperate cases");

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    cpubinding_decl_usage(pl);
#else
    param_list_decl_usage(pl, "cpubinding", "(ignored, no sufficient sotftware support)");
#endif
}
/*}}}*/

void parallelizing_info_lookup_parameters(param_list_ptr pl)/*{{{*/
{
    /* These will all be looked later (within this file, though).
     */
    param_list_lookup_string(pl, "mpi");
    param_list_lookup_string(pl, "thr");
    param_list_lookup_string(pl, "interleaving");
    param_list_lookup_string(pl, "only_mpi");
    param_list_lookup_string(pl, "debug-parallel-bwc");

#if defined(HAVE_HWLOC) && defined(HAVE_CXX11)
    cpubinding_lookup_parameters(pl);
#else
    param_list_lookup_string(pl, "cpubinding");
#endif
}
/*}}}*/


void pi_go(
        void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
        param_list pl,
        void * arg)
{
    int interleaving = 0;
    param_list_parse_int(pl, "interleaving", &interleaving);
    if (interleaving) {
        pi_go_inner_interleaved(fcn, pl, arg);
    } else {
        pi_go_inner_not_interleaved(fcn, pl, arg);
    }
}


void pi_log_init(pi_comm_ptr wr)
{
    struct pi_log_book * lb = malloc(sizeof(struct pi_log_book));
    memset(lb, 0, sizeof(struct pi_log_book));
    wr->log_book = lb;

    char commname[MPI_MAX_OBJECT_NAME];
    int namelen = sizeof(commname);
    MPI_Comm_get_name(wr->pals, commname, &namelen);

    if (wr->jrank == 0 && wr->trank == 0) {
        printf("Enabled logging for %s:%s\n", commname, wr->th->desc);
    }
}

void pi_log_clear(pi_comm_ptr wr)
{
    if (wr->log_book)
        free(wr->log_book);
    wr->log_book = (void*) 0xdeadbeef;
}

void pi_log_op(pi_comm_ptr wr, const char * fmt, ...)
{
    va_list ap;
    struct pi_log_book * lb = wr->log_book;
    if (!lb)
        return;

    va_start(ap, fmt);

    ASSERT_ALWAYS(lb->next < PI_LOG_BOOK_ENTRIES);

    struct pi_log_entry * e = &(lb->t[lb->next]);

    gettimeofday(e->tv, NULL);

    vsnprintf(e->what, sizeof(e->what), fmt, ap);
    if (wr->ncores > 1)
        fprintf(stderr, "%s:%d:%d %s\n", wr->th->desc, wr->jrank, wr->trank, e->what);
    va_end(ap);

    lb->next++;
    lb->next %= PI_LOG_BOOK_ENTRIES;
    lb->hsize++;
}

static void pi_log_print_backend(pi_comm_ptr wr, const char * myname, char ** strings, int * n, int alloc)
{
    struct pi_log_book * lb = wr->log_book;
    if (!lb) return;

    ASSERT_ALWAYS(lb->next < PI_LOG_BOOK_ENTRIES);

    int h = lb->hsize;
    if (h >= PI_LOG_BOOK_ENTRIES) {
        h = PI_LOG_BOOK_ENTRIES;
    }
    int i0 = lb->next - h;
    if (i0 < 0)
        i0 += PI_LOG_BOOK_ENTRIES;
    char commname[MPI_MAX_OBJECT_NAME];
    int namelen = sizeof(commname);
    MPI_Comm_get_name(wr->pals, commname, &namelen);

    int rc;

    for(int i = i0 ; h-- ; ) {
        struct pi_log_entry * e = &(lb->t[i]);

        ASSERT_ALWAYS(*n < alloc);
        rc = asprintf(&(strings[*n]), "%" PRIu64 ".%06u (%s) %s %s", 
                (uint64_t) e->tv->tv_sec,
                (unsigned int) e->tv->tv_usec,
                myname,
                e->what,
                commname);
        ASSERT_ALWAYS(rc != -1);
        ++*n;

        i++;
        if (i == PI_LOG_BOOK_ENTRIES)
            i = 0;
    }
}

void pi_log_print(pi_comm_ptr wr)
{
    /* FIXME stdout or stderr ? */
    int alloc = PI_LOG_BOOK_ENTRIES;
    char ** strings = malloc(alloc * sizeof(char*));
    int n = 0;
    pi_log_print_backend(wr, "", strings, &n, alloc);
    for(int i = 0 ; i < n ; i++) {
        puts(strings[i]);
        free(strings[i]);
    }
    fflush(stdout);
    free(strings);
}

typedef int (*sortfunc_t)(const void *, const void *);

int p_strcmp(char const * const * a, char const * const * b)
{
    /* Input is formatted so that sorting makes sense. */
    return strcmp(*a, *b);
}

void pi_log_print_all(parallelizing_info_ptr pi)
{
    int alloc = 3 * PI_LOG_BOOK_ENTRIES;
    char ** strings = malloc(alloc * sizeof(char*));
    int n = 0;
    char * myname;
    int rc;
    rc = asprintf(&myname, "%s%s", pi->wr[0]->th->desc, pi->wr[1]->th->desc);
    ASSERT_ALWAYS(rc != -1);

    pi_log_print_backend(pi->m, myname, strings, &n, alloc);
    pi_log_print_backend(pi->wr[0], myname, strings, &n, alloc);
    pi_log_print_backend(pi->wr[1], myname, strings, &n, alloc);
    qsort(strings, n, sizeof(char*), (sortfunc_t) &p_strcmp);

    for(int i = 0 ; i < n ; i++) {
        puts(strings[i]);
        free(strings[i]);
    }
    free(myname);
    free(strings);
}

/* {{{ predefined types and operations */
static struct pi_datatype_s pi_predefined_types[] = {
    { MPI_INT,	                NULL, sizeof(int) },
    { MPI_DOUBLE,	        NULL, sizeof(double) },
    { MPI_BYTE,	                NULL, 1 },
    { MPI_UNSIGNED,	        NULL, sizeof(unsigned) },
    { MPI_UNSIGNED_LONG,	NULL, sizeof(unsigned long) },
    { MPI_UNSIGNED_LONG_LONG,	NULL, sizeof(unsigned long long) },
    { MPI_LONG,	                NULL, sizeof(long) },
};

/* TODO: mpi-3.0 introduced (at last) MPI_UINT64_T and friends. There are
 * several places where this can come in handy. So once we set our mind
 * on requiring mpi-3.0, we can do some simplifications here and there.
 */
pi_datatype_ptr BWC_PI_INT              = pi_predefined_types + 0;
pi_datatype_ptr BWC_PI_DOUBLE           = pi_predefined_types + 1;
pi_datatype_ptr BWC_PI_BYTE             = pi_predefined_types + 2;
pi_datatype_ptr BWC_PI_UNSIGNED         = pi_predefined_types + 3;
pi_datatype_ptr BWC_PI_UNSIGNED_LONG    = pi_predefined_types + 4;
pi_datatype_ptr BWC_PI_UNSIGNED_LONG_LONG    = pi_predefined_types + 5;
pi_datatype_ptr BWC_PI_LONG             = pi_predefined_types + 6;


struct pi_op_s BWC_PI_MIN[1]	= { { .stock = MPI_MIN, } };
struct pi_op_s BWC_PI_MAX[1]	= { { .stock = MPI_MAX, } };

static void pi_dispatch_op_add_stock(void *, void *, int *, MPI_Datatype *);
static void pi_dispatch_op_add_custom(void *, void *, size_t, pi_datatype_ptr);

struct pi_op_s BWC_PI_SUM[1]	= { {
        .stock = MPI_SUM,
        .f_stock = pi_dispatch_op_add_stock,
        .f_custom = pi_dispatch_op_add_custom,
    } };
struct pi_op_s BWC_PI_BXOR[1]	= { { .stock = MPI_BXOR, } };
struct pi_op_s BWC_PI_BAND[1]	= { { .stock = MPI_BAND, } };
struct pi_op_s BWC_PI_BOR[1]	= { { .stock = MPI_BOR, } };

struct reduction_function {
    MPI_Datatype datatype;
    MPI_Op op;
    thread_reducer_t f;
};

/* {{{ predefined reduction functions */

/* always: a == my value, b == the other value */
static void reducer_byte_min(const unsigned char * b, unsigned char * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_byte_max(const unsigned char * b, unsigned char * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_byte_sum(const unsigned char * b, unsigned char * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_byte_bxor(const unsigned char * b, unsigned char * a, int s) { for( ; s-- ; a++, b++) *a ^= *b; }
static void reducer_int_min(const int * b, int * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_int_max(const int * b, int * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_int_sum(const int * b, int * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_int_band(const int * b, int * a, int s) { for( ; s-- ; a++, b++) *a &= *b; }
static void reducer_int_bor(const int * b, int * a, int s) { for( ; s-- ; a++, b++) *a |= *b; }
static void reducer_uint_min(const unsigned int * b, unsigned int * a, unsigned int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_uint_max(const unsigned int * b, unsigned int * a, unsigned int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_uint_sum(const unsigned int * b, unsigned int * a, unsigned int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_uint_band(const unsigned int * b, unsigned int * a, unsigned int s) { for( ; s-- ; a++, b++) *a &= *b; }
static void reducer_uint_bor(const unsigned int * b, unsigned int * a, unsigned int s) { for( ; s-- ; a++, b++) *a |= *b; }
static void reducer_double_min(const double * b, double * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_double_max(const double * b, double * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_double_sum(const double * b, double * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_long_min(const long * b, long * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_long_max(const long * b, long * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_long_sum(const long * b, long * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_longlong_min(const long long * b, long long * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_longlong_max(const long long * b, long long * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_longlong_sum(const long long * b, long long * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_ulong_min(const unsigned long * b, unsigned long * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_ulong_max(const unsigned long * b, unsigned long * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_ulong_sum(const unsigned long * b, unsigned long * a, int s) { for( ; s-- ; a++, b++) *a += *b; }
static void reducer_ulonglong_min(const unsigned long long * b, unsigned long long * a, int s) { for( ; s-- ; a++, b++) if (*b < *a) *a = *b; }
static void reducer_ulonglong_max(const unsigned long long * b, unsigned long long * a, int s) { for( ; s-- ; a++, b++) if (*b > *a) *a = *b; }
static void reducer_ulonglong_sum(const unsigned long long * b, unsigned long long * a, int s) { for( ; s-- ; a++, b++) *a += *b; }

struct reduction_function predefined_functions[] = {
    { MPI_BYTE,          MPI_MIN,  (thread_reducer_t) reducer_byte_min, },
    { MPI_BYTE,          MPI_MAX,  (thread_reducer_t) reducer_byte_max, },
    { MPI_BYTE,          MPI_SUM,  (thread_reducer_t) reducer_byte_sum, },
    { MPI_BYTE,          MPI_BXOR, (thread_reducer_t) reducer_byte_bxor, },
    { MPI_INT,           MPI_MIN,  (thread_reducer_t) reducer_int_min, },
    { MPI_INT,           MPI_MAX,  (thread_reducer_t) reducer_int_max, },
    { MPI_INT,           MPI_SUM,  (thread_reducer_t) reducer_int_sum, },
    { MPI_INT,           MPI_BAND, (thread_reducer_t) reducer_int_band, },
    { MPI_INT,           MPI_BOR,  (thread_reducer_t) reducer_int_bor, },
    { MPI_UNSIGNED,      MPI_MIN,  (thread_reducer_t) reducer_uint_min, },
    { MPI_UNSIGNED,      MPI_MAX,  (thread_reducer_t) reducer_uint_max, },
    { MPI_UNSIGNED,      MPI_SUM,  (thread_reducer_t) reducer_uint_sum, },
    { MPI_UNSIGNED,      MPI_BAND, (thread_reducer_t) reducer_uint_band, },
    { MPI_UNSIGNED,      MPI_BOR,  (thread_reducer_t) reducer_uint_bor, },
    { MPI_DOUBLE,        MPI_MIN,  (thread_reducer_t) reducer_double_min, },
    { MPI_DOUBLE,        MPI_MAX,  (thread_reducer_t) reducer_double_max, },
    { MPI_DOUBLE,        MPI_SUM,  (thread_reducer_t) reducer_double_sum, },
    { MPI_LONG,          MPI_MIN,  (thread_reducer_t) reducer_long_min, },
    { MPI_LONG,          MPI_MAX,  (thread_reducer_t) reducer_long_max, },
    { MPI_LONG,          MPI_SUM,  (thread_reducer_t) reducer_long_sum, },
    { MPI_UNSIGNED_LONG, MPI_MIN,  (thread_reducer_t) reducer_ulong_min, },
    { MPI_UNSIGNED_LONG, MPI_MAX,  (thread_reducer_t) reducer_ulong_max, },
    { MPI_UNSIGNED_LONG, MPI_SUM,  (thread_reducer_t) reducer_ulong_sum, },
    { MPI_LONG_LONG,          MPI_MIN,  (thread_reducer_t) reducer_longlong_min, },
    { MPI_LONG_LONG,          MPI_MAX,  (thread_reducer_t) reducer_longlong_max, },
    { MPI_LONG_LONG,          MPI_SUM,  (thread_reducer_t) reducer_longlong_sum, },
    { MPI_UNSIGNED_LONG_LONG, MPI_MIN,  (thread_reducer_t) reducer_ulonglong_min, },
    { MPI_UNSIGNED_LONG_LONG, MPI_MAX,  (thread_reducer_t) reducer_ulonglong_max, },
    { MPI_UNSIGNED_LONG_LONG, MPI_SUM,  (thread_reducer_t) reducer_ulonglong_sum, },
    { 0, 0, NULL, },
};

static thread_reducer_t lookup_reduction_function(MPI_Datatype datatype, MPI_Op op)
{
    struct reduction_function * x = predefined_functions;
    for( ; x->f && ! (x->datatype == datatype && x->op == op) ; x++) ;
    if (!x->f) {
        fprintf(stderr, "bad reduction function\n");
        abort();
    }
    return x->f;
}

/* }}} */

/* }}} */

/* {{{ user-defined (mpfq) types and operations */

static int pi_mpi_attribute_key;

/* we need to present the mpi implementation with proper user functions
 * for operating on the mpfq types. Do so only for addition, since it's
 * the only one we use eventually.
 *
 * Of course we could conceivably do so for more operations.
 */
static void pi_dispatch_op_add_stock(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
    /* get the mpfq_vbase_ptr associated to the datatype */
    int got_it;
    mpfq_vbase_ptr abase;
    MPI_Type_get_attr(*datatype, pi_mpi_attribute_key, (void*) &abase, &got_it);
    assert(got_it);
    abase->vec_add(abase, inoutvec, inoutvec, invec, *len);
}
static void pi_dispatch_op_add_custom(void *invec, void *inoutvec, size_t len, pi_datatype_ptr datatype)
{
    /* FIXME: mpfq's vec_add should really take size_t arguments */
    ASSERT_ALWAYS(len < (size_t) UINT_MAX);
    datatype->abase->vec_add(datatype->abase, inoutvec, inoutvec, invec, len);
}

/* XXX must be called in single threaded context  */
static void pi_init_attribute_things()
{
    MPI_Type_create_keyval(MPI_TYPE_DUP_FN, MPI_TYPE_NULL_DELETE_FN, &pi_mpi_attribute_key, NULL);
    MPI_Op_create(BWC_PI_SUM->f_stock, 1, &BWC_PI_SUM->custom);
}

/* XXX must be called in single threaded context  */
static void pi_clear_attribute_things()
{
    MPI_Op_free(&BWC_PI_SUM->custom);
    MPI_Type_free_keyval(&pi_mpi_attribute_key);
}

pi_datatype_ptr pi_alloc_mpfq_datatype(parallelizing_info_ptr pi, mpfq_vbase_ptr abase)
{
    pi_datatype_ptr ptr = shared_malloc_set_zero(pi->m, sizeof(struct pi_datatype_s));
    if (pi->m->trank == 0) {
        ptr->abase = abase;
        ptr->item_size = abase->vec_elt_stride(abase, 1);
        MPI_Type_contiguous(ptr->item_size, MPI_BYTE, &ptr->datatype);
        MPI_Type_commit(&ptr->datatype);
        MPI_Type_set_attr(ptr->datatype, pi_mpi_attribute_key, abase);
    }
    serialize_threads(pi->m);
    return ptr;
}

void pi_free_mpfq_datatype(parallelizing_info_ptr pi, pi_datatype_ptr ptr)
{
    serialize_threads(pi->m);
    if (pi->m->trank == 0) {
        MPI_Type_delete_attr(ptr->datatype, pi_mpi_attribute_key);
        MPI_Type_free(&ptr->datatype);
    }
    shared_free(pi->m, ptr);
}

/* This may *not* use MPI_Reduce_local, because we might not have an MPI
 * implementation in the first place, just a set of placeholers. So it is
 * important that we find the operation to perform by ourselves
 */
void pi_reduce_local(void *inbuf, void *inoutbuf, size_t count,
                    pi_datatype_ptr datatype, pi_op_ptr op)
{
    if (datatype->abase) {
        if (op->f_custom) {
            op->f_custom(inbuf, inoutbuf, count, datatype);
        } else {
            fprintf(stderr, "undefined operation on mpfq type\n");
            abort();
        }
    } else {
        thread_reducer_t f = lookup_reduction_function(datatype->datatype, op->stock);
        (*f)(inbuf, inoutbuf, count);
    }
}

/* }}} */

/* {{{ collective operations with pi_comm_ptrs
 *
 * we provide a calling interface which is similar to mpi.
 *
 * we inherit the same crappy notion of datatypes, and use the same
 * MPI_Datatype data for it. When not building with MPI, our "fake MPI"
 * headers provides magic constants which are sufficient to our needs.
 *
 * We define the following functions:
 *
 *      {pi,pi_thread}_{bcast,allreduce,data_eq}
 *
 * where data_eq is a collective operation returning true or false
 * depending on whether the pointed data area is consistent across all
 * callers.
 *
 * Being thread-level collective functions, these routines assume that
 * the provided pointers are different on the calling threads (otherwise
 * there's no point).
 *
 * Notice also that as long as the underlying mpi implementation cannot
 * be considered thread-safe in the MPI_THREAD_MULTIPLE meaning, it is
 * not possible to call any of the pi collective functions simultaneously
 * from two distinct parallel communicators. That is, either calls are with
 * pi->m, or with pi->wr[d] provided that pi->wr[!d]->trank == 0. The
 * SEVERAL_THREADS_PLAY_MPI_BEGIN construct may be useful, as follows:
 *      SEVERAL_THREADS_PLAY_MPI_BEGIN(pi->wr[!d]) {
 *              pi-collective operation on pi->wr[d]
 *      }
 *      SEVERAL_THREADS_PLAY_MPI_END();
 */

struct pi_collective_arg {
    pi_comm_ptr wr;
    /* this pointer is written from all threads. Its contents are
     * meaningful only for the root thread */
    void * ptr;
    unsigned int root;
    pi_datatype_ptr datatype;
    size_t count;
    /* last two only for reduction */
    pi_op_ptr op;
    void * dptr;
};

/* {{{ broadcast
 *
 * thread-level strategy:
 * upon entering the barrier, the "root" threads fills in
 * a->wr->th->utility_ptr its own pointer. Upon leaving the barrier, all
 * threads copy data from there to their own buffer.
 *
 * mpi-level strategy:
 * broadcast on all other nodes from the root first, then on all threads
 * separately on each node.
 */
static void pi_thread_bcast_in(int s MAYBE_UNUSED, struct pi_collective_arg * a)
{
    /* here we have a mutex locked */
    if (a->wr->trank == a->root) {
        a->wr->th->utility_ptr = a->ptr;
    }
}

static void pi_thread_bcast_out(int s MAYBE_UNUSED, struct pi_collective_arg * a)
{
    /* here we know that exactly one thread has just filled the
     * a->wr->th->utility_ptr, so that a memcpy is possible, as long as
     * the owner of the data area has not touched it again */
    if (a->ptr != a->wr->th->utility_ptr) {       /* placate valgrind */
        memcpy(a->ptr, a->wr->th->utility_ptr, a->count * a->datatype->item_size);
    }
}

/* broadcast the data area pointed to by [ptr, ptr+size[ to all threads.
 * Pointers at calling threads must differ. */
void pi_thread_bcast(void * ptr, size_t count, pi_datatype_ptr datatype, unsigned int root, pi_comm_ptr wr)
{
    struct pi_collective_arg a[1];
    a->wr = wr;
    a->ptr = ptr;
    a->root = root;
    a->datatype = datatype;
    a->count = count;
    barrier_wait(wr->th->bh,
            (void(*)(int,void*)) &pi_thread_bcast_in,
            (void(*)(int,void*)) &pi_thread_bcast_out,
            a);
    serialize_threads(wr);
}

void pi_bcast(void * ptr, size_t count, pi_datatype_ptr datatype, unsigned int jroot, unsigned int troot, pi_comm_ptr wr)
{
    int err;

    ASSERT(jroot < wr->njobs);
    ASSERT(troot < wr->ncores);
    if (wr->trank == troot) {
        err = MPI_Bcast(ptr, count, datatype->datatype, jroot, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    pi_thread_bcast(ptr, count, datatype, troot, wr);
}

/* }}} */

/* {{{ allreduce
 *
 * thread-level strategy (dumb): use the provided sendbuf's, one after another,
 * as "last worked on" areas. The last being reduced in is then used to
 * communicate the data to all other threads.
 *
 * mpi-level strategy: thread-level allreduce, then in-place mpi-level allreduce
 * from the leader threads, the thread-level broadcast of the result.
 */
/* this gets called with s = (count-1), (count-2), ..., 0 */
static void pi_thread_allreduce_in(int s, struct pi_collective_arg * a)
{
    /* utility_ptr points to the "last worked on" area.
     * first call: set the "last worked on" area to our pointer.
     * further calls: reduce our data with the "lwo" data, and set the
     * "lwo" pointer to our modified data in our pointer.
     *
     * last call: reduction done, copy to dptr.
     */
    if (s < (int) a->wr->ncores - 1) {
        pi_reduce_local(a->wr->th->utility_ptr, a->ptr, a->count, a->datatype, a->op);
    }
    a->wr->th->utility_ptr = a->ptr;
}

static void pi_thread_allreduce_out(int s MAYBE_UNUSED, struct pi_collective_arg * a)
{
    /* almost the same as pi_thread_bcast_out, except that we work on
     * dptr...
     */
    if (a->dptr != a->wr->th->utility_ptr) {       /* placate valgrind */
        memcpy(a->dptr, a->wr->th->utility_ptr, a->count * a->datatype->item_size);
    }
}

/* ptr == dptr is supported */
void pi_thread_allreduce(void *ptr, void *dptr, size_t count,
                pi_datatype_ptr datatype, pi_op_ptr op, pi_comm_ptr wr)
{
    /* this code is buggy in the not-in-place case. let's fall back
     * to the in-place situation unconditionally.
     *
     * nature of the bug: if we don't do the fallback below, then clearly
     * a->ptr is being written to. We should clean this up, and use const
     * when possible.
     */
    if (ptr && ptr != dptr) {
        memcpy(dptr, ptr, count * datatype->item_size);
        ptr = NULL;
    }
    struct pi_collective_arg a[1];
    a->wr = wr;
    a->ptr = ptr ? ptr : dptr;  /* handle in-place */
    a->root = 0;        /* meaningless for allreduce */
    a->datatype = datatype;
    a->count = count;
    a->op = op;
    a->dptr = dptr;
    barrier_wait(wr->th->bh,
            (void(*)(int,void*)) &pi_thread_allreduce_in,
            (void(*)(int,void*)) &pi_thread_allreduce_out,
            a);
    serialize_threads(wr);
}


/* This intentionally has the same prototype as MPI_Allreduce */
/* Only a few operations and datatypes are supported.
 *
 * NULL at sendbuf means in-place.
 */
void pi_allreduce(void *sendbuf, void *recvbuf, size_t count,
        pi_datatype_ptr datatype, pi_op_ptr op, pi_comm_ptr wr)
{
    ASSERT_ALWAYS(count <= (size_t) INT_MAX);
    pi_thread_allreduce(sendbuf, recvbuf, count, datatype, op, wr);
    if (wr->trank == 0) {
        if (datatype->abase) {
            /* Then it's a type for which we supposedly have written an
             * overloaded operation. */
            MPI_Allreduce(MPI_IN_PLACE, recvbuf, count, datatype->datatype, op->custom, wr->pals);
        } else {
            MPI_Allreduce(MPI_IN_PLACE, recvbuf, count, datatype->datatype, op->stock, wr->pals);
        }
    }
    /* now it's just a matter of broadcasting to all threads */
    pi_thread_bcast(recvbuf, count, datatype, 0, wr);
}

/* }}} */

/* {{{ data_eq */

int pi_thread_data_eq(void *buffer, size_t count, pi_datatype_ptr datatype, pi_comm_ptr wr)
{
    void ** bufs = shared_malloc(wr, wr->ncores * sizeof(void*));
    int * oks = shared_malloc(wr, wr->ncores * sizeof(int));

    bufs[wr->trank] = buffer;
    oks[wr->trank] = 1;

    /* given 5 threads, here is for example what is done at each step.
     *
     * stride=1
     *  ok0 = (x0==x1)
     *  ok2 = (x2==x3)
     * stride=2
     *  ok0 = x0==x1==x2==x3
     * stride 4
     *  ok0 = x0==x1==x2==x3==x4
     */
    for(unsigned int stride = 1 ; stride < wr->ncores ; stride <<= 1) {
        unsigned int i = wr->trank;
        serialize_threads(wr);       /* because of this, we don't
                                           leave the loop with break */
        if (i + stride >= wr->ncores) continue;
        /* many threads are doing nothing here */
        if (i & stride) continue;
        oks[i] = oks[i] && oks[i + stride] && memcmp(bufs[i], bufs[i + stride], datatype->item_size * count) == 0;
    }
    serialize_threads(wr);
    int ok = oks[0];
    shared_free(wr, oks);
    shared_free(wr, bufs);
    return ok;
}

int pi_data_eq(void *buffer, size_t count, pi_datatype_ptr datatype, pi_comm_ptr wr)
{
    int ok = pi_thread_data_eq(buffer, count, datatype, wr);

    ASSERT_ALWAYS(count <= (size_t) INT_MAX);
    if (wr->trank == 0) {
        MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_BAND, wr->pals);
    }
    pi_thread_bcast(&ok, 1, BWC_PI_INT, 0, wr);
    return ok;
}

/* }}} */

/* }}} */

/*{{{ collective malloc and free (thread-level of course) */
void * shared_malloc(pi_comm_ptr wr, size_t size)
{
    void * ptr = NULL;
    if (wr->trank == 0) ptr = malloc(size);
    pi_thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0, wr);
    return ptr;
}

void * shared_malloc_set_zero(pi_comm_ptr wr, size_t size)
{
    void * ptr = NULL;
    if (wr->trank == 0) {
        ptr = malloc(size);
        memset(ptr, 0, size);
    }
    pi_thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0, wr);
    return ptr;
}

void shared_free(pi_comm_ptr wr, void * ptr)
{
    serialize_threads(wr);
    if (wr->trank == 0) free(ptr);
}
/*}}}*/

int serialize__(pi_comm_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
{
    int err;
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(%s:%u) barrier #%d, %u/%u on %s [%p]\n", s, l, w->th_count,
            w->trank, w->ncores,
            w->th->desc, w->th->b);
#endif
    my_pthread_barrier_wait(w->th->b);
    if (w->trank == 0) {
        err = MPI_Barrier(w->pals);
        ASSERT_ALWAYS(!err);
    }
    // struct timeval tv[1];
    // gettimeofday(tv, NULL);
    // printf("%.2f\n", tv->tv_sec + (double) tv->tv_usec / 1.0e6);
    // sleep(1);
    return w->jrank == 0 && w->trank == 0;
}

int serialize_threads__(pi_comm_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
{
#ifdef  CONCURRENCY_DEBUG
    // note how w->th_count is normally a thread-private thing !
    w->th_count++;
    printf("(%s:%u) tbarrier #%d, %u/%u on %s [%p]\n", s, l, w->th_count,
            w->trank, w->ncores,
            w->th->desc, w->th->b);
#endif
    my_pthread_barrier_wait(w->th->b);
    // struct timeval tv[1];
    // gettimeofday(tv, NULL);
    // printf("%.2f\n", tv->tv_sec + (double) tv->tv_usec / 1.0e6);
    // sleep(1);
    return w->trank == 0;
}

// XXX Many serializing calls in the two routines below have no
// functional use, apart from the fact that they exert some pressure on
// the mpi+pthreads way of doing things.

static void say_hello(pi_comm_ptr w, parallelizing_info_ptr pi MAYBE_UNUSED)
{
    serialize(w);
    for(unsigned int j = 0 ; j < w->njobs ; j++) {
        serialize(w);
        if (w->jrank != j)
            continue;
        my_pthread_mutex_lock(w->th->m);
#ifdef CONCURRENCY_DEBUG
        /* Make it less verbose -- if it ever hangs in there, then
         * we can re-enable it */
        printf("(%s) J%uT%u ; %s%s ; (%s:j%ut%u) (%s:j%ut%u)\n",
                w->th->desc,
                pi->m->jrank,
                pi->m->trank,
                pi->wr[0]->th->desc,
                pi->wr[1]->th->desc,
                pi->wr[0]->th->desc, pi->wr[0]->jrank, pi->wr[0]->trank,
                pi->wr[1]->th->desc, pi->wr[1]->jrank, pi->wr[1]->trank
              );
#endif
        my_pthread_mutex_unlock(w->th->m);
    }
}

void pi_hello(parallelizing_info_ptr pi)
{
    if (serialize(pi->m)) {
#ifdef  CONCURRENCY_DEBUG
        printf("Doing hello world loop\n");
#endif
    }
    say_hello(pi->m, pi);

    // before changing the grain of the barriers, we have to serialize at
    // the thread level (i.e. wait for the mpi serialization calls to
    // complete). Otherwise, inter-job communication occuring within
    // say_hello might clash with the inter-job calls done within
    // serialize().
    // Rule is: never let mpi-level synchronisations on a communicator
    // wander outside an execution period for which all threads sharing
    // this communicator are doing the same.
    for(unsigned int i = 0 ; i < pi->wr[0]->totalsize ; i++) {
        serialize(pi->m);
        serialize_threads(pi->m);
        if (i == pi->wr[0]->jrank * pi->wr[0]->ncores + pi->wr[0]->trank) {
            say_hello(pi->wr[1], pi);
        }
    }

    for(unsigned int i = 0 ; i < pi->wr[1]->totalsize ; i++) {
        serialize(pi->m);
        serialize_threads(pi->m);
        if (i == pi->wr[1]->jrank * pi->wr[1]->ncores + pi->wr[1]->trank) {
            say_hello(pi->wr[0], pi);
        }
    }

    if (serialize(pi->m)) {
#ifdef CONCURRENCY_DEBUG
        printf("OK: Finished hello world loop\n");
#endif
    }
    serialize(pi->m);
}
/*{{{ new i/o layer */
/* Open a file handle on the leader node. The two pi_comm_ptr's
 * represent the inner and outer levels; they're mostly unused for this
 * function call, but requested for consistency with the pi_file_read and
 * pi_file_write calls.
 *
 * returns true on success, false on error.
 *
 * In multithreaded context, f must represent a different data area on
 * all threads.
 */
int pi_file_open(pi_file_handle_ptr f, parallelizing_info_ptr pi, int inner, const char * name, const char * mode)
{
    memset(f, 0, sizeof(pi_file_handle));
    f->pi = pi;
    f->inner = inner;
    f->outer = !inner;
    f->f = NULL;
    f->name = strdup(name);
    f->mode = strdup(mode);
    int rc = 1;
    if (f->pi->m->jrank == 0 && f->pi->m->trank == 0) {
        f->f = fopen(name, mode);
        if (f->f == NULL)
            rc = 0;
    }
    pi_bcast(&rc, 1, BWC_PI_INT, 0, 0, f->pi->m);
    return rc;
}

/* close the file handle */
void pi_file_close(pi_file_handle_ptr f)
{
    free(f->name);
    free(f->mode);
    if (f->pi->m->jrank == 0 && f->pi->m->trank == 0)
        fclose(f->f);
    memset(f, 0, sizeof(pi_file_handle));
}

int area_is_zero(const void * src, ptrdiff_t offset0, ptrdiff_t offset1)
{
    const char * b0 = pointer_arith_const(src, offset0);
    const char * b1 = pointer_arith_const(src, offset1);
    for( ; b0 < b1 ; b0++) {
        if (*b0) return 0;
    }
    return 1;
}

ssize_t pi_file_write_leader(pi_file_handle_ptr f, void * buf, size_t size)
{
    unsigned long ret;
    if (f->pi->m->jrank == 0 && f->pi->m->trank == 0)
        ret = fwrite(buf, 1, size, f->f);
    pi_bcast(&ret, 1, BWC_PI_UNSIGNED_LONG, 0, 0, f->pi->m);
    return (ssize_t) ret;
}

ssize_t pi_file_read_leader(pi_file_handle_ptr f, void * buf, size_t size)
{
    unsigned long ret;
    if (f->pi->m->jrank == 0 && f->pi->m->trank == 0)
        ret = fread(buf, 1, size, f->f);
    pi_bcast(&ret, 1, BWC_PI_UNSIGNED_LONG, 0, 0, f->pi->m);
    return (ssize_t) ret;
}


ssize_t pi_file_write(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize)
{
    ASSERT_ALWAYS(size <= ULONG_MAX);

    /* Some sanity checking. This is not technically required: of
     * course we would be able to cope with uneven data sets.
     * However, having even sizes greatly simplifies the pictures, as
     * it alleviates the need to collect the sizes for all threads */
    unsigned long size_ul = size;
    ASSERT_ALWAYS(pi_data_eq(&size_ul, 1, BWC_PI_UNSIGNED_LONG, f->pi->m));

    unsigned int nci = f->pi->wr[f->inner]->ncores;
    unsigned int nco = f->pi->wr[f->outer]->ncores;
    unsigned int nji = f->pi->wr[f->inner]->njobs;
    unsigned int njo = f->pi->wr[f->outer]->njobs;
    unsigned int ci  = f->pi->wr[f->inner]->trank;
    unsigned int co  = f->pi->wr[f->outer]->trank;
    unsigned int ji  = f->pi->wr[f->inner]->jrank;
    unsigned int jo  = f->pi->wr[f->outer]->jrank;
    unsigned int jm  = f->pi->m->jrank;

    /* the inner threads will do writes and MPI sends in one go. For
     * outer threads, it is not possible to concatenate similarly,
     * because we have writes relative to the mpi-level interleaved.
     */
    void * sbuf = shared_malloc_set_zero(f->pi->m, size * nci);

    long res = 0;
    int failed = 0;

    for(unsigned int xjo = 0 ; xjo < njo ; xjo++) {
        for(unsigned int xco = 0 ; xco < nco ; xco++) {
            for(unsigned int xji = 0 ; xji < nji ; xji++) {
                unsigned int xj = mrank_from_tworanks(f->pi, f->inner, xji, xjo);
                if (jm == xj) {
                    ASSERT_ALWAYS (ji == xji && jo == xjo);
                    /* fill the buffer with the contents relevant to the
                     * inner communicator (nci threads at work here) */
                    serialize_threads(f->pi->m);
                    if (xco == co) memcpy(sbuf + ci * size, buf, size);
                    serialize_threads(f->pi->m);
                }
                if (f->pi->m->trank == 0) {
                    /* only one thread per node does something. */
                    if (!jm) {
                        /* leader node to receive data */
                        if (xj) {       /* except its own... */
                            MPI_Recv(sbuf, size * nci, MPI_BYTE,
                                    xj, xj, f->pi->m->pals,
                                    MPI_STATUS_IGNORE);
                        }
                        if (!failed) {
                            size_t wanna_write = MIN(size * nci, totalsize - res);
                            if (wanna_write < size * nci) {
                                ASSERT_ALWAYS(area_is_zero(sbuf, wanna_write, size * nci));
                            }
                            ssize_t x = fwrite(sbuf, 1, wanna_write, f->f);
                            res += x;
                            if (x < (ssize_t) wanna_write) {
                                failed=1;
                            }
                        }
                    } else if (jm == xj) {
                        /* my turn to send data */
                        MPI_Send(sbuf, size * nci, MPI_BYTE,
                                0, xj, f->pi->m->pals);
                    }
                }
            }
        }
    }
    shared_free(f->pi->m, sbuf);
    pi_bcast(&res, 1, BWC_PI_LONG, 0, 0, f->pi->m);
    return res;
}

ssize_t pi_file_read(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize)
{
    ASSERT_ALWAYS(size <= ULONG_MAX);

    /* Some sanity checking. This is not technically required: of
     * course we would be able to cope with uneven data sets.
     * However, having even sizes greatly simplifies the pictures, as
     * it alleviates the need to collect the sizes for all threads */
    unsigned long size_ul = size;
    ASSERT_ALWAYS(pi_data_eq(&size_ul, 1, BWC_PI_UNSIGNED_LONG, f->pi->m));

    unsigned int nci = f->pi->wr[f->inner]->ncores;
    unsigned int nco = f->pi->wr[f->outer]->ncores;
    unsigned int nji = f->pi->wr[f->inner]->njobs;
    unsigned int njo = f->pi->wr[f->outer]->njobs;
    unsigned int ci  = f->pi->wr[f->inner]->trank;
    unsigned int co  = f->pi->wr[f->outer]->trank;
    unsigned int ji  = f->pi->wr[f->inner]->jrank;
    unsigned int jo  = f->pi->wr[f->outer]->jrank;
    unsigned int jm  = f->pi->m->jrank;

    void * sbuf = shared_malloc_set_zero(f->pi->m, size * nci);

    long res = 0;
    int failed = 0;

    for(unsigned int xjo = 0 ; xjo < njo ; xjo++) {
        for(unsigned int xco = 0 ; xco < nco ; xco++) {
            for(unsigned int xji = 0 ; xji < nji ; xji++) {
                unsigned int xj = mrank_from_tworanks(f->pi, f->inner, xji, xjo);
                /* for job (ji,jo), nci*nco threads enter here. Unless
                 * (ji == xji && jo == xjo), this loop becomes trivial.
                 */
                if (f->pi->m->trank == 0) {
                    /* only one thread per node does something. */
                    if (!jm) {
                        /* leader node to read then send data */
                        memset(sbuf, 0, nci * size);
                        if (!failed) {
                            size_t wanna_read = MIN(size * nci, totalsize - res);
                            ssize_t x = fread(sbuf, 1, wanna_read, f->f);
                            res += x;
                            if (x < (ssize_t) wanna_read) {
                                failed=1;
                            }
                        }
                        if (xj) {       /* except to itself... */
                            MPI_Send(sbuf, size * nci, MPI_BYTE,
                                    xj, xj, f->pi->m->pals);
                        }
                    } else if (jm == xj) {
                        /* my turn to receive data */
                        MPI_Recv(sbuf, size * nci, MPI_BYTE,
                                0, xj, f->pi->m->pals,
                                MPI_STATUS_IGNORE);
                    }
                }
                if (jm == xj) {
                    ASSERT_ALWAYS (ji == xji && jo == xjo);
                    serialize_threads(f->pi->m);
                    /* fill the buffer with the contents relevant to the
                     * inner communicator (nci threads at work here) */
                    if (xco == co) memcpy(buf, sbuf + ci * size, size);
                    serialize_threads(f->pi->m);
                }
            }
        }
    }
    shared_free(f->pi->m, sbuf);
    pi_bcast(&res, 1, BWC_PI_LONG, 0, 0, f->pi->m);
    return res;
}
/*}}}*/

void pi_interleaving_flip(parallelizing_info_ptr pi)
{
    if (!pi->interleaved)
        return;
    my_pthread_barrier_wait(pi->interleaved->b);
}

void pi_interleaving_enter(parallelizing_info_ptr pi)
{
    if (!pi->interleaved || pi->interleaved->idx == 0)
        return;

    my_pthread_barrier_wait(pi->interleaved->b);
}

void pi_interleaving_leave(parallelizing_info_ptr pi)
{
    if (!pi->interleaved || pi->interleaved->idx == 1)
        return;

    my_pthread_barrier_wait(pi->interleaved->b);
}

#if 0
/* XXX the name is not very informative, notably about the important
 * point that pointer must be identical.
 */
/* Considering a data area with a pointer identical on all calling
 * threads, tell whether the contents agree on all MPI jobs+threads
 * (returns 1 if this is the case)
 *
 * More exactly, this function only cares about the pointer on thread 0.
 */
int mpi_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    /* TODO: think about allocating less */
    int cmp = 0;
    if (pi->m->trank == 0) {
        void * b_and = malloc(sz);
        void * b_or = malloc(sz);
        MPI_Allreduce(buffer, b_and, sz, MPI_BYTE, MPI_BAND, pi->m->pals);
        MPI_Allreduce(buffer, b_or, sz, MPI_BYTE, MPI_BOR, pi->m->pals);
        cmp = memcmp(b_and, b_or, sz);
        free(b_and);
        free(b_or);
    }
    thread_bcast(pi->m, &cmp, sizeof(int), 0);
    return cmp == 0;
}
#endif

/* XXX must be called in single threaded context  */
void parallelizing_info_init()
{
    pi_init_attribute_things();
}
/* XXX must be called in single threaded context  */
void parallelizing_info_finish()
{
    pi_clear_attribute_things();
}



