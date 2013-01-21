#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>		/* for cputime */
#ifdef HAVE_RESOURCE_H
#include <sys/resource.h>	/* for cputime */
#endif
#include <sys/time.h>	/* for gettimeofday */
#include "timing.h"
#include "memusage.h"

/* return total user time (all threads) */
uint64_t
microseconds (void)
{
#ifdef HAVE_GETRUSAGE
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    uint64_t r;
    r = (uint64_t) res->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec;
    return r;
#else
    return 0;
#endif
}

/* only consider user time of the current thread */
uint64_t
microseconds_thread (void)
{
#ifdef HAVE_GETRUSAGE
    struct rusage ru[1];
    uint64_t r;

#ifdef HAVE_RUSAGE_THREAD
    getrusage (RUSAGE_THREAD, ru);
#else
    getrusage (RUSAGE_SELF, ru);
#endif
    r = (uint64_t) ru->ru_utime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) ru->ru_utime.tv_usec;
    return r;
#else
    return 0;
#endif
}

/* cputime */
unsigned long 
milliseconds (void)
{
    return (unsigned long) (microseconds() / (uint64_t) 1000);
}

double
seconds (void)
{
    return (double) microseconds() / 1.0e6;
}

double
seconds_thread (void)
{
    return (double) microseconds_thread () / 1.0e6;
}

void
seconds_user_sys (double * res)
{
#ifdef HAVE_GETRUSAGE
    struct rusage ru[1];

    getrusage (RUSAGE_SELF, ru);
    res[0] = ru->ru_utime.tv_sec +  (double) ru->ru_utime.tv_usec / 1.0e6;
    res[1] = ru->ru_stime.tv_sec +  (double) ru->ru_stime.tv_usec / 1.0e6;
#else
    res[0] = res[1] = 0.;
#endif
}

/* returns the number of seconds since the Epoch (1970-01-01 00:00:00 +0000).
   Thus we have to call it twice and subtract to get the wall clock time of
   a given program. */
double
wct_seconds (void)
{
    struct timeval tv[1];
    gettimeofday (tv, NULL);
    return (double) tv->tv_sec + 1.0e-6 * tv->tv_usec;
}

void
print_timing_and_memory (double wct0)
{
  fprintf (stderr, "Total usage: time %1.0fs (cpu), %1.0fs (wct) ; "
           "memory %luM, peak %luM\n",
           seconds (), wct_seconds () - wct0,
           Memusage () >> 10, PeakMemusage () >> 10);
}

