#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// we're doing open close mmap truncate...
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>
#include "bwc_config.h"
#include "select_mpi.h"
#include "parallelizing_info.h"
#include "macros.h"

#include <sys/time.h>   // gettimeofday

static inline void pi_wiring_init_pthread_things(pi_wiring_ptr w, const char * desc)
{
    struct pthread_things * res;

    res = malloc(sizeof(struct pthread_things));

    barrier_init(res->bh, w->ncores);
    my_pthread_barrier_init(res->b, NULL, w->ncores);
    my_pthread_mutex_init(res->m, NULL);
    res->desc = strdup(desc);

    w->th = res;
}

static inline void pi_wiring_destroy_pthread_things(pi_wiring_ptr w)
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
    void * ret = (s->fcn)(s->p, s->pl, s->arg);
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

/* TODO: There would be a fairly easy generalization to that, and quite
 * handy. Modify the meaning of the ``siz'' argument to mean something
 * which doesn't have to be equal on all threads. Instead, pick the
 * maximum and use it for formatting.
 */

void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print)
{
    pi_wiring_ptr wr = pi->m;

    char * strings;
    strings = malloc(wr->njobs * wr->ncores * siz);

    int me = wr->jrank * wr->ncores + wr->trank;

    /* instead of doing memcpy, we align the stuff. */
    char * fmt;
    int rc = asprintf(&fmt, "%%-%zus", siz-1);
    ASSERT_ALWAYS(rc >= 0);
    snprintf(strings + me * siz, siz, fmt, buf);
    free(fmt);

    /* ceinture et bretelles */

    serialize(wr);
    for(unsigned int j = 0 ; j < wr->njobs ; j++) {
        for(unsigned int t = 0 ; t < wr->ncores ; t++) {
            /* we serialize because of the thread_broadcast bug */
            // serialize_threads(wr);
            global_broadcast(wr, strings + (j * wr->ncores + t) * siz,
                    siz, j, t);
        }
    }
    serialize(wr);

    char * ptr = strings;

    if (print) {
        /* There's also the null byte counted in siz. So that makes one
         * less.
         */
        unsigned int nj0 = pi->wr[0]->njobs;
        unsigned int nj1 = pi->wr[1]->njobs;
        unsigned int nt0 = pi->wr[0]->ncores;
        unsigned int nt1 = pi->wr[1]->ncores;
        print_several(nj0, nt0, '-', '+', siz-1);
        for(unsigned int j1 = 0 ; j1 < nj1 ; j1++) {
            for(unsigned int t1 = 0 ; t1 < nt1 ; t1++) {
                printf("|");
                for(unsigned int j0 = 0 ; j0 < nj0 ; j0++) {
                    for(unsigned int t0 = 0 ; t0 < nt0 ; t0++) {
                        int fence = t0 == nt0 - 1;
                        printf("%s%c", ptr, fence ? '|' : ' ');
                        ptr += siz;
                    }
                }
                printf("\n");
            }
            print_several(nj0, nt0, '-', '+', siz-1);
        }
    }
    serialize(wr);
    free(strings);
}

static void pi_init_mpilevel(parallelizing_info_ptr pi, param_list pl)
{
    int err;
    int mpi[2] = { 1, 1, };
    int thr[2] = { 1, 1, };

    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);

    unsigned int nhj = mpi[0];
    unsigned int nvj = mpi[1];
    unsigned int nhc = thr[0];
    unsigned int nvc = thr[1];

    memset(pi, 0, sizeof(parallelizing_info));

    pi->m->njobs = nhj * nvj;
    pi->m->ncores = nhc * nvc;
    pi->m->totalsize = pi->m->njobs * pi->m->ncores;

#ifndef NDEBUG
    /* Must make sure that we have proper ASSERTS after MPI calls then. */
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif

    // MPI_Errhandler my_eh;
    // MPI_Comm_create_errhandler(&pi_errhandler, &my_eh);
    // MPI_Comm_set_errhandler(MPI_COMM_WORLD, my_eh);

    MPI_Comm_dup(MPI_COMM_WORLD, & pi->m->pals);
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
    pi_wiring_init_pthread_things(pi->m, "main");

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
        pi_wiring_ptr leader = grid[Nc]->wr[0];
        pi_wiring_init_pthread_things(leader, buf);
        // replicate.
        for(unsigned int k = 0 ; k < pi->wr[0]->ncores ; k++) {
            grid[k + Nc]->wr[0]->th = leader->th;
        }
    }

    // column barriers.
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        snprintf(buf, sizeof(buf), "c%u", 
                pi->wr[0]->jrank * pi->wr[0]->ncores + c);
        pi_wiring_ptr leader = grid[c]->wr[1];
        pi_wiring_init_pthread_things(leader, buf);
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
        pi_wiring_destroy_pthread_things(grid[Nc]->wr[0]);
    }
    // destroy column barriers
    for(unsigned int c = 0 ; c < pi->wr[0]->ncores ; c++) {
        pi_wiring_destroy_pthread_things(grid[c]->wr[1]);
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
    pi_wiring_destroy_pthread_things(pi->m);

    for(int d = 0 ; d < 2 ; d++) {
        MPI_Comm_free(&pi->wr[d]->pals);
    }
    MPI_Comm_free(&pi->m->pals);

    // MPI_Errhandler_free(&my_eh);
}

#if 0
static void shout_going_mt(pi_wiring_ptr m)
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

void pi_go(void *(*fcn)(parallelizing_info_ptr, param_list pl, void *),
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


void pi_log_init(pi_wiring_ptr wr)
{
    struct pi_log_book * lb = malloc(sizeof(struct pi_log_book));
    memset(lb, 0, sizeof(struct pi_log_book));
    wr->log_book = lb;

    char commname[MPI_MAX_OBJECT_NAME];
    int namelen = sizeof(commname);
    MPI_Comm_get_name(wr->pals, commname, &namelen);
    if (wr->jrank == 0 && wr->trank == 0) {
        printf("Enabled logging for %s\n", commname);
    }
}

void pi_log_clear(pi_wiring_ptr wr)
{
    if (wr->log_book)
        free(wr->log_book);
    wr->log_book = (void*) 0xdeadbeef;
}

void pi_log_op(pi_wiring_ptr wr, const char * fmt, ...)
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
    va_end(ap);

    lb->next++;
    lb->next %= PI_LOG_BOOK_ENTRIES;
    lb->hsize++;
}

static void pi_log_print_backend(pi_wiring_ptr wr, const char * myname, char ** strings, int * n, int alloc)
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
        rc = asprintf(&(strings[*n]), "%"PRIu64".%06u (%s) %s %s", 
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

void pi_log_print(pi_wiring_ptr wr)
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

struct thread_broadcast_arg {
    pi_wiring_ptr wr;
    void ** ptr;
    unsigned int i;
};

void thread_broadcast_in(int s MAYBE_UNUSED, struct thread_broadcast_arg * a)
{
    if (a->wr->trank == a->i) {
        a->wr->th->utility_ptr = *a->ptr;
    }
}

void thread_broadcast_out(int s MAYBE_UNUSED, struct thread_broadcast_arg * a)
{
    *a->ptr = a->wr->th->utility_ptr;
}

void thread_broadcast(pi_wiring_ptr wr, void ** ptr, unsigned int i)
{
    struct thread_broadcast_arg a[1];
    a->wr = wr;
    a->ptr = ptr;
    a->i = i;
    barrier_wait(wr->th->bh,
            (void(*)(int,void*)) &thread_broadcast_in,
            (void(*)(int,void*)) &thread_broadcast_out,
            a);
#if 0
    // FIXME: this design is flawed. Correct reentrancy is not achieved.
    // All the stuff around the barrier here should be done with the
    // mutex locked, otherwise disaster may occur.
    // A fix would incur either basically rewriting barrier code, so as
    // to place stuff at the right spot. Note that there _is_ some
    // barrier code in bw/barrier.[ch] which could easily be expanded for
    // this purpose.
    if (wr->trank == i) {
        wr->th->utility_ptr = *ptr;
    }
    my_pthread_barrier_wait(wr->th->b);
    *ptr = wr->th->utility_ptr;
#endif
}

void global_broadcast(pi_wiring_ptr wr, void * ptr, size_t size, unsigned int j, unsigned int t)
{
    int err;

    ASSERT(j < wr->njobs);
    ASSERT(t < wr->ncores);
    if (wr->trank == t) {
        err = MPI_Bcast(ptr, size, MPI_BYTE, j, wr->pals);
        ASSERT_ALWAYS(!err);
    }
    void * leader_ptr = ptr;
    thread_broadcast(wr, &leader_ptr, t);
    if (ptr != leader_ptr)
        memcpy(ptr, leader_ptr, size);
}


int serialize__(pi_wiring_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
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

int serialize_threads__(pi_wiring_ptr w, const char * s MAYBE_UNUSED, unsigned int l MAYBE_UNUSED)
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

void say_hello(pi_wiring_ptr w, parallelizing_info_ptr pi MAYBE_UNUSED)
{
    serialize(w);
    for(unsigned int j = 0 ; j < w->njobs ; j++) {
        serialize(w);
        if ((unsigned int) w->jrank != j)
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

void hello(parallelizing_info_ptr pi)
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

static int get_counts_and_displacements(pi_wiring_ptr w, int my_size,
        int * displs, int * counts)
{
    counts[w->jrank] = my_size;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, counts, 1, MPI_INT, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* Now all threads know the size argument of their friend threads. We
     * need to know the displacement values. These are obtained by
     * looking at the progressive sum over threads.
     *
     * We're serialization-happy here, since anyway I/O s occur rarely
     */
    /* The _void thingy makes gcc happy, and I'm not too comfortable
     * about deciding whether those pesky aliasing warnings are false
     * positives or not. Better shut them up for good. */
    void * allcounts_void;
    if (w->trank == 0) {
        allcounts_void = malloc(w->njobs*w->ncores*sizeof(int));
    }
    thread_broadcast(w, &allcounts_void, 0);
    int * allcounts = (int *) allcounts_void;
    for(unsigned int k = 0 ; k < w->njobs ; k++) {
        allcounts[k * w->ncores + w->trank] = counts[k];
    }
    serialize_threads(w);
    int dis = 0;
    for(unsigned int k = 0 ; k < w->njobs ; k++) {
        unsigned int t;
        for(t = 0 ; t < w->trank ; t++) {
            dis += allcounts[k * w->ncores + t];
        }
        displs[k] = dis;
        dis += allcounts[k * w->ncores + t];
        for(t++ ; t < w->ncores ; t++) {
            dis += allcounts[k * w->ncores + t];
        }
    }
    serialize_threads(w);
    if (w->trank == 0) {
        free(allcounts);
    }
    return dis;
}

/* d indicates whether we read the ranks row by row (d==0) or column by
 * column */
int get_counts_and_displacements_2d(parallelizing_info_ptr pi, int d,
        int my_size, int * displs, int * counts)
{
    pi_wiring_ptr w = pi->m;

    memset(displs, -1, w->njobs * sizeof(int));

    counts[w->jrank] = my_size;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, counts, 1, MPI_INT, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    /* Now all threads know the size argument of their friend threads. We
     * need to know the displacement values. These are obtained by
     * looking at the progressive sum over threads.
     *
     * We're serialization-happy here, since anyway I/O s occur rarely
     */
    void * allcounts_void;
    if (w->trank == 0) {
        allcounts_void = malloc(w->njobs*w->ncores*sizeof(int));
    }
    thread_broadcast(w, &allcounts_void, 0);
    int * allcounts = (int *) allcounts_void;

    // Given an orientation value d, the numbering for the thread (it,jt)
    // on job (ij,jj) is:
    //
    // it+pi->wr[d]->ncores*(ij+pi->wr[d]->njobs*(jt+pi->wr[!d]->ncores*jj))
    // where ij ranges over pi->wr[d]->njobs and jj over pi->wr[!d]->njobs
    // (it and jt similar of course)
    //
    // Since we're filling the array in a multithreaded way, we're only
    // interested by the relevant value of it,jt ; the contribution
    // depending on ij,jj can be added later on. So we split the sum:
    //
    // it+pi->wr[d]->ncores*pi->wr[d]->njobs*jt
    // +
    // pi->wr[d]->ncores*(ij+pi->wr[d]->njobs*pi->wr[!d]->ncores*jj)
    // == 
    // it+pi->wr[d]->ncores*pi->wr[d]->njobs*jt
    // +
    // pi->wr[d]->ncores*ij + pi->m->ncores*pi->wr[d]->njobs*jj

    unsigned int v0;

    v0 = pi->wr[d]->trank+pi->wr[!d]->trank*pi->wr[d]->ncores*pi->wr[d]->njobs;
    for(unsigned int jj = 0 ; jj < pi->wr[!d]->njobs ; jj++) {
        unsigned int v = v0;
        /* d == 1 : we'll order offsets column by column */
        /* d == 0 : we'll order row column by row */
        for(unsigned int ij = 0 ; ij < pi->wr[d]->njobs ; ij++) {
            /* relative to pi->m, which job exactly corresponds to this
             * count ?
             * it's row_index * ncols + col_index. But we have to know
             * which is the row index.
             * d == 1 : ij is the row index
             * d == 0 : jj is the row index
             */
            unsigned int which = 0;
            which += ij * (d == 1 ? pi->wr[0]->njobs : 1);
            which += jj * (d == 0 ? pi->wr[0]->njobs : 1);
            ASSERT(which < pi->m->njobs);
            allcounts[v] = counts[which];
            v += pi->wr[d]->ncores;
        }
        v0 += pi->m->ncores*pi->wr[d]->njobs;
    }

    serialize_threads(w);

    int * alldisps = (int*) malloc(w->njobs*w->ncores*sizeof(int));
    int dis = 0;
    for(unsigned int k = 0 ; k < w->njobs * w->ncores ; k++) {
        alldisps[k] = dis;
        dis += allcounts[k];
    }

    v0 = pi->wr[d]->trank+pi->wr[!d]->trank*pi->wr[d]->ncores*pi->wr[d]->njobs;
    for(unsigned int jj = 0 ; jj < pi->wr[!d]->njobs ; jj++) {
        unsigned int v = v0;
        for(unsigned int ij = 0 ; ij < pi->wr[d]->njobs ; ij++) {
            unsigned int which = 0;
            which += ij * (d == 1 ? pi->wr[0]->njobs : 1);
            which += jj * (d == 0 ? pi->wr[0]->njobs : 1);
            ASSERT(which < pi->m->njobs);
            displs[which] = alldisps[v];
            v += pi->wr[d]->ncores;
        }
        v0 += pi->m->ncores*pi->wr[d]->njobs;
    }

    free(alldisps);

#if 0
    for(unsigned int i = 0 ; i < w->njobs ; i++) {
        my_pthread_mutex_lock(w->th->m);
        printf("J%uT%u %u @ %u\n", i, w->trank, counts[i], displs[i]);
        my_pthread_mutex_unlock(w->th->m);
    }
#endif

    serialize_threads(w);

#if 0
    ASSERT_ALWAYS(displs[w->jrank] >= 0);
    if (displs[w->jrank] == 0)
        ASSERT_ALWAYS(w->trank == 0);
#endif

    if (w->trank == 0) {
        free(allcounts);
    }

    return dis;
}

/* That's a workalike to MPI_File_write_ordered, except that we work with
 * the finer-grained pi_wiring_ptr's.
 */
int pi_save_file(pi_wiring_ptr w, const char * name, void * buf, size_t mysize)
{
    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * recvcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements(w, mysize, displs, recvcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * recvbuf = NULL;
    int fd = -1;        // only used by leader
    int rc;
    int err;

    if (leader) {
        fd = open(name, O_RDWR|O_CREAT, 0666);
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n", name, strerror(errno));
            goto pi_save_file_leader_init_done;
        }

        rc = ftruncate(fd, wsiz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(%s): %s\n",
                    name, strerror(errno));
            close(fd);
            goto pi_save_file_leader_init_done;
        }
        recvbuf = mmap(NULL, wsiz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (recvbuf == MAP_FAILED) {
            fprintf(stderr, "mmap(%s): %s\n",
                    name, strerror(errno));
            recvbuf = NULL;
            close(fd);
            goto pi_save_file_leader_init_done;
        }
pi_save_file_leader_init_done:
        ;
    }

    /* Now all threads from job zero see the area mmaped by their leader
     * thread. It does not make sense on other jobs of course, since for
     * these, recvbuf is unused. */
    if (w->jrank == 0) {
        thread_broadcast(w, &recvbuf, 0);
    }

    /* Rather unfortunate, but error checking requires some checking. As
     * mentioned earlier, we don't feel concerned a lot by this, since
     * inour context I/Os are rare enough
     */
    int ok = recvbuf != NULL;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Bcast(&ok, 1, MPI_INT, 0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
    if (!ok) return 0;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Gatherv(buf, mysize, MPI_BYTE,
                recvbuf, recvcounts, displs, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    // what a pain in the Xss. Because we don't exit with the mutex
    // locked, we can't do much... The main thread might end up calling
    // munmap before we complete the call to MPI_Gatherv.

    serialize_threads(w);

    free(recvcounts);
    free(displs);

    if (leader) {
        munmap(recvbuf, wsiz);
        rc = ftruncate(fd, siz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(): %s\n", strerror(errno));
            /* If only ftruncate failed, don't return an error */
        }
        close(fd);
    }

    return 1;
}

/* d indicates whether we read the ranks row by row (d==0) or column by
 * column */
int pi_save_file_2d(parallelizing_info_ptr pi, int d, const char * name, void * buf, size_t mysize)
{
    pi_wiring_ptr w = pi->m;

    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * recvcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements_2d(pi, d, mysize, displs, recvcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * recvbuf = NULL;
    int fd = -1;        // only used by leader
    int rc;
    int err;

    if (leader) {
        fd = open(name, O_RDWR|O_CREAT, 0666);
        if (fd < 0) {
            fprintf(stderr, "fopen(%s): %s\n", name, strerror(errno));
            goto pi_save_file_2d_leader_init_done;
        }

        rc = ftruncate(fd, wsiz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(%s): %s\n",
                    name, strerror(errno));
            close(fd);
            goto pi_save_file_2d_leader_init_done;
        }
        recvbuf = mmap(NULL, wsiz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
        if (recvbuf == MAP_FAILED) {
            fprintf(stderr, "mmap(%s): %s\n",
                    name, strerror(errno));
            recvbuf = NULL;
            close(fd);
            goto pi_save_file_2d_leader_init_done;
        }
pi_save_file_2d_leader_init_done:
        ;
    }

    /* Now all threads from job zero see the area mmaped by their leader
     * thread. It does not make sense on other jobs of course, since for
     * these, recvbuf is unused. */
    if (w->jrank == 0) {
        thread_broadcast(w, &recvbuf, 0);
    }

    /* Rather unfortunate, but error checking requires some checking. As
     * mentioned earlier, we don't feel concerned a lot by this, since
     * inour context I/Os are rare enough
     */
    int ok = recvbuf != NULL;
    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Bcast(&ok, 1, MPI_INT, 0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;
    if (!ok) return 0;

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        err = MPI_Gatherv(buf, mysize, MPI_BYTE,
                recvbuf, recvcounts, displs, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    // what a pain in the Xss. Because we don't exit with the mutex
    // locked, we can't do much... The main thread might end up calling
    // munmap before we complete the call to MPI_Gatherv.

    serialize_threads(w);

    free(recvcounts);
    free(displs);

    if (leader) {
        munmap(recvbuf, wsiz);
        rc = ftruncate(fd, siz);
        if (rc < 0) {
            fprintf(stderr, "ftruncate(): %s\n", strerror(errno));
            /* If only ftruncate failed, don't return an error */
        }
        close(fd);
    }

    return 1;
}

int pi_load_file(pi_wiring_ptr w, const char * name, void * buf, size_t mysize)
{
    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * sendcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements(w, mysize, displs, sendcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * sendbuf = NULL;
    int fd = -1;        // only used by leader
    int err;

    if (leader) {
        fd = open(name, O_RDONLY, 0666);
        DIE_ERRNO_DIAG(fd < 0, "fopen", name);
        sendbuf = mmap(NULL, wsiz, PROT_READ, MAP_SHARED, fd, 0);
        DIE_ERRNO_DIAG(sendbuf == MAP_FAILED, "mmap", name);
    }

    if (w->jrank == 0) {
        thread_broadcast(w, &sendbuf, 0);
    }

    /* For loading, in contrast to saving, we simply abort if reading
     * can't work. So no need to agree on an ok value. */

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        ASSERT_ALWAYS(mysize == (size_t) sendcounts[w->jrank]);
        err = MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE,
                buf, mysize, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    serialize_threads(w);

    free(sendcounts);
    free(displs);

    if (leader) {
        munmap(sendbuf, wsiz);
        close(fd);
    }
    return 1;
}

int pi_load_file_2d(parallelizing_info_ptr pi, int d, const char * name, void * buf, size_t mysize)
{
    pi_wiring_ptr w = pi->m;

    int * displs = (int *) malloc(w->njobs * sizeof(int));
    int * sendcounts = (int *) malloc(w->njobs * sizeof(int));
    size_t siz;
    siz = get_counts_and_displacements_2d(pi, d, mysize, displs, sendcounts);

    // the page size is always a power of two, so rounding to the next
    // multiple is easy.
    size_t wsiz = ((siz - 1) | (sysconf(_SC_PAGESIZE)-1)) + 1;
    int leader = w->jrank == 0 && w->trank == 0;
    void * sendbuf = NULL;
    int fd = -1;        // only used by leader
    int err;

    if (leader) {
        fd = open(name, O_RDONLY, 0666);
        DIE_ERRNO_DIAG(fd < 0, "fopen", name);
        sendbuf = mmap(NULL, wsiz, PROT_READ, MAP_SHARED, fd, 0);
        DIE_ERRNO_DIAG(sendbuf == MAP_FAILED, "mmap", name);
    }

    if (w->jrank == 0) {
        thread_broadcast(w, &sendbuf, 0);
    }

    /* For loading, in contrast to saving, we simply abort if reading
     * can't work. So no need to agree on an ok value. */

    SEVERAL_THREADS_PLAY_MPI_BEGIN(w) {
        ASSERT_ALWAYS(mysize == (size_t) sendcounts[w->jrank]);
        err = MPI_Scatterv(sendbuf, sendcounts, displs, MPI_BYTE,
                buf, mysize, MPI_BYTE,
                0, w->pals);
        ASSERT_ALWAYS(!err);
    }
    SEVERAL_THREADS_PLAY_MPI_END;

    serialize_threads(w);

    free(sendcounts);
    free(displs);

    if (leader) {
        munmap(sendbuf, wsiz);
        close(fd);
    }
    return 1;
}

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

int mpi_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    int ok = 0;
    int * pok = &ok;
    if (pi->m->trank == 0) {
        void * b_and = malloc(sz);
        void * b_or = malloc(sz);
        MPI_Allreduce(buffer, b_and, sz, MPI_BYTE, MPI_BAND, pi->m->pals);
        MPI_Allreduce(buffer, b_or, sz, MPI_BYTE, MPI_BOR, pi->m->pals);
        ok = memcmp(b_and, b_or, sz);
        free(b_and);
        free(b_or);
    }
    thread_broadcast(pi->m, (void**) &pok, 0);
    return *pok == 0;
}

int thread_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    char ** bufs;
    if (pi->m->trank == 0)
        bufs = malloc(pi->m->ncores * sizeof(char*));
    thread_broadcast(pi->m, (void**) &bufs, 0);
    bufs[pi->m->trank] = buffer;
    int ok = 0;
    int * pok = &ok;
    if (serialize_threads(pi->m)) {
        char * b_and = malloc(sz);
        char * b_or = malloc(sz);
        memset(b_and, 0, sz);
        memset(b_or, 0, sz);
        for(unsigned int i = 0 ; i < pi->m->ncores ; i++) {
            for(size_t k = 0 ; k < sz ; k++) {
                b_and[k] &= bufs[i][k];
                b_or[k]  |= bufs[i][k];
            }
        }
        ok = memcmp(b_and, b_or, sz);
        free(b_and);
        free(b_or);
    }
    thread_broadcast(pi->m, (void**) &pok, 0);
    return *pok == 0;
}

int global_data_eq(parallelizing_info_ptr pi, void *buffer, size_t sz)
{
    char ** bufs;
    if (pi->m->trank == 0)
        bufs = malloc(pi->m->ncores * sizeof(char*));
    thread_broadcast(pi->m, (void**) &bufs, 0);
    bufs[pi->m->trank] = buffer;
    int ok = 0;
    int * pok = &ok;
    if (serialize_threads(pi->m)) {
        char * b_and = malloc(sz);
        char * b_or = malloc(sz);
        memcpy(b_and, buffer, sz);
        memcpy(b_or, buffer, sz);
        for(unsigned int i = 0 ; i < pi->m->ncores ; i++) {
            for(size_t k = 0 ; k < sz ; k++) {
                b_and[k] &= bufs[i][k];
                b_or[k]  |= bufs[i][k];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, b_and, sz, MPI_BYTE, MPI_BAND, pi->m->pals);
        MPI_Allreduce(MPI_IN_PLACE, b_or,  sz, MPI_BYTE, MPI_BOR,  pi->m->pals);
        ok = memcmp(b_and, b_or, sz);
        /* We may use just any non-zero pointer */
        pok = ok == 0 ? (void*) pi : NULL;
        free(b_and);
        free(b_or);
    }
    thread_broadcast(pi->m, (void**) &pok, 0);
    ok = pok != NULL;
    return ok;
}

