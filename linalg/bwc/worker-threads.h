#ifndef CADO_WORKER_THREADS_H
#define CADO_WORKER_THREADS_H

/*
 * This API develops the concept of ``worker threads''. The control flow
 * of the main program is handled by the invoking (=controlling) thread
 * (which might or might not be the unique thread in the program, such a
 * distinction is not relevant here). Worker threads are used when some
 * clearly identified task has to be performed in parallel.
 * Schematically, the control flow goes like this.
 *
 *  (1)          (2)          (2')        (2)          (2')        (3)
 * ==+==(works)===#  (waits)   #===========#  (waits)   #===========+
 *   |                                                              |
 *   <  (waits)   +------------+  (waits)  +------------+  (waits)  >
 *   <  (waits)   +------------+  (waits)  +------------+  (waits)  >
 *
 * (1),(2),(2'),(3) are detailed below. M denotes the controlling thread,
 * T0 and T1 the two worker threads.
 *
 * (1) M does tg=worker_threads_init(2). T0 and T1 are created, and are
 * put on standby awaiting requests by M.
 *
 * (2) M calls worker_threads_do(tg, f, arg). This implies:
 *      T0 executes f(tg,0,arg)
 *      T1 executes f(tg,1,arg)
 *      M waits for T0 and T1 to complete.
 * (2') T0 and T1 have both completed. M resumes executions, while T0 and
 * T1 return to standby.
 *
 * Further occurrences of (2) are normal (otherwise pthread_create /
 * pthread_join would have done the job). The function f needs not be
 * constant across all worker_threads_do calls.
 *
 * (3) M calls worker_threads_clear(tg). T0 and T1 are terminated.
 */

/* Further notes.
 *
 * Sub-threads (T0, T1 above) may access some fiels from the tg structure while
 * executing. In particular, tg->mu is the mutex guarding the concurrency
 * spots of the control flow. It is safe to use it for protecting
 * critical parts in other locations. No thread keeps the mutex for long.
 */

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

struct worker_threads_group;

typedef void (*worker_func_t)(struct worker_threads_group *, int, void *);

struct worker_thread_info {
    int i;
    struct worker_threads_group * tg;
};

struct worker_threads_group {
    /* The mutex can be freely used by callees */
    pthread_mutex_t mu;

    /* Total number of threads */
    unsigned int n;

    /* Number of worker threads at work (not idle) */
    unsigned int working;

    /* Number of threads which have completed their current task */
    unsigned int done;

    /* practically everything else is meant to be private and/or useless
     * anyway to callees */
    pthread_cond_t cond_in;
    pthread_cond_t cond_out;

    struct worker_thread_info * thread_table;
    pthread_t * threads;

    worker_func_t f;
    void * arg;
};

extern void worker_threads_do(struct worker_threads_group * tg, worker_func_t f, void * arg);
extern struct worker_threads_group * worker_threads_init(unsigned int n);
extern void worker_threads_clear(struct worker_threads_group * tg); 

#ifdef __cplusplus
}
#endif

#endif	/* CADO_WORKER_THREADS_H */
