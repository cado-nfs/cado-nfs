#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <bits/types/struct_rusage.h>

#include <stdio.h>
#include <stdlib.h>     // exit
#include <stdint.h>
#include <inttypes.h>

#include <pthread.h>
#if defined(HAVE_GETRUSAGE)
#include <sys/resource.h> // IWYU pragma: keep
#endif

#include "timing.h"

/* We're not really checking the functionality here, in the sense
 * that timings are just for user reports, and even if they were
 * consistently printing 0, this would only be a low-priority
 * concern.
 */


/* some tests spawn multiple threads, here's the number */
#define NTHREADS_TEST 4

/* arrange so that this takes at least 20ms, so that scheduler gets a
 * chance to see it */
static unsigned long do_some_silly_thing()
{
    unsigned long p = 65521;
    volatile unsigned long a = 1009;
    for(int i = 0 ; i < 1000000 ; i++) {
        a = (a*a) % p;
    }
    return a;
}

/* this should take a wee bit more so that %.2f does not print a
 * seemingly zero value like .01 */
static unsigned long do_lots_of_silly_things()
{
    unsigned long v = 0;
    for(int i = 0; i < 10 ; i++) {
        v += do_some_silly_thing();
    }
    return v;
}

int test_microseconds ()
{
    uint64_t t = microseconds();
    do_some_silly_thing();
    printf("%s %" PRIu64"\n", __func__, microseconds()-t);
    return 1;
}

static void * test_microseconds_thread_subthread(void * arg)
{
    uint64_t t0 = microseconds_thread();
    do_some_silly_thing();
    uint64_t t1 = microseconds_thread();
    *(uint64_t*)arg = t1 - t0;
    return NULL;
}

static void generic_mutltithread_test(void*(*func)(void*), void ** args)
{
    pthread_t sub[NTHREADS_TEST];
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        pthread_create(sub + i, NULL, func, args[i]);
    }
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        pthread_join(sub[i], NULL);
    }
}

int test_microseconds_thread()
{
    uint64_t timers[NTHREADS_TEST]={0,};
    void * ptrs[NTHREADS_TEST];
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        ptrs[i] = timers + i;
    }
    generic_mutltithread_test(&test_microseconds_thread_subthread, ptrs);
    uint64_t s=0;
    for(int i = 0 ; i < NTHREADS_TEST ; s += timers[i++]);
    printf("%s %" PRIu64"\n", __func__, s);
    return 1;
}

int test_milliseconds ()
{
    unsigned long t = milliseconds();
    do_lots_of_silly_things();
    printf("%s %lu\n", __func__, milliseconds()-t);
    return 1;
}
int test_seconds ()
{
    double t = seconds();
    do_lots_of_silly_things();
    printf("%s %.2f\n", __func__, seconds()-t);
    return 1;
}

static void * test_seconds_thread_subthread(void * arg)
{
    double t = seconds_thread();
    do_some_silly_thing();
    *(double*)arg = seconds_thread()-t;
    return NULL;
}

int test_seconds_thread()
{
    double timers[NTHREADS_TEST]={0,};
    void * ptrs[NTHREADS_TEST];
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        ptrs[i] = timers + i;
    }
    generic_mutltithread_test(&test_seconds_thread_subthread, ptrs);
    double s=0;
    for(int i = 0 ; i < NTHREADS_TEST ; i++)
        s += timers[i];
    printf("%s %.2f\n", __func__, s);
    return 1;
}


int test_seconds_user_sys()
{
    double us[2];
    seconds_user_sys(us);
    return 1;
}

int test_wct_seconds()
{
    wct_seconds();
    return 1;
}

int test_print_timing_and_memory()
{
    print_timing_and_memory (stdout, seconds () - 17, wct_seconds() - 1000);
    return 1;
}

void * test_thread_seconds_user_sys_subthread(void * arg)
{
    thread_seconds_user_sys((double*)arg);
    return NULL;
}

int test_thread_seconds_user_sys()
{
    double timers[NTHREADS_TEST][2]={{0,}};
    void * ptrs[NTHREADS_TEST];
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        ptrs[i] = timers + i;
    }
    generic_mutltithread_test(&test_thread_seconds_user_sys_subthread, ptrs);
    double s=0,t=0;
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        s += timers[i][0];
        t += timers[i][1];
    }
    printf("%s %.2f %.2f\n", __func__, s, t);
    return 1;
}

void * test_timingstats_dict_subthread(void * arg)
{
    timingstats_dict_add_mythread((timingstats_dict_ptr) arg, "subthread");
    return NULL;
}

int test_timingstats_dict()
{
    timingstats_dict_t td;
    timingstats_dict_init(td);
#ifdef HAVE_GETRUSAGE
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    timingstats_dict_add(td, "Hello", res);
#endif
    timingstats_dict_add_myprocess(td, "main");
    pthread_t sub[NTHREADS_TEST];
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        pthread_create(sub + i, NULL, &test_timingstats_dict_subthread, td);
    }
    for(int i = 0 ; i < NTHREADS_TEST ; i++) {
        pthread_join(sub[i], NULL);
    }
    timingstats_dict_disp(td);
    timingstats_dict_clear(td);
    return 1;
}

int main()
{
    int status = 1;

    status = status && test_microseconds ();
    status = status && test_microseconds_thread();
    status = status && test_milliseconds ();
    status = status && test_seconds ();
    status = status && test_seconds_thread();
    status = status && test_seconds_user_sys();
    status = status && test_wct_seconds();
    status = status && test_print_timing_and_memory();
    status = status && test_thread_seconds_user_sys();
    status = status && test_timingstats_dict();

    return status ? EXIT_SUCCESS : EXIT_FAILURE;
}
