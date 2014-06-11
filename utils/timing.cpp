#include "cado.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>		/* for cputime */
#ifdef HAVE_RESOURCE_H
#include <sys/resource.h>	/* for cputime */
#endif
#include <sys/time.h>	/* for gettimeofday */
#include "macros.h"
#include "timing.h"
#include "memusage.h"
#include "portability.h"

#ifdef HAVE_GETRUSAGE
/* I'm including some STL code for the timer info layer, but this could
 * equally well be done in C */
#include <map>
#include <string>
#endif

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
    res[0] = (double)ru->ru_utime.tv_sec + (double)ru->ru_utime.tv_usec/1.0e6;
    res[1] = (double)ru->ru_stime.tv_sec + (double)ru->ru_stime.tv_usec/1.0e6;
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
    return (double)tv->tv_sec + (double)tv->tv_usec*1.0e-6;
}

void
print_timing_and_memory (double wct0)
{
  fprintf (stderr, "Total usage: time %1.0fs (cpu), %1.0fs (wct) ; "
           "memory %luM, peak %luM\n",
           seconds (), wct_seconds () - wct0,
           Memusage () >> 10, PeakMemusage () >> 10);
}

/* We need some way to detect the time spent by threads. Unfortunately,
 * this is something not defined by POSIX.
 */

#if defined(HAVE_RUSAGE_THREAD)
void thread_seconds_user_sys(double * res)
{
    struct rusage ru[1];
    getrusage(RUSAGE_THREAD, ru);
    res[0] = (double)ru->ru_utime.tv_sec + (double)ru->ru_utime.tv_usec/1.0e6;
    res[1] = (double)ru->ru_stime.tv_sec + (double)ru->ru_stime.tv_usec/1.0e6;
}
#elif defined(__linux)
#include <unistd.h>
#include <sys/syscall.h>
#include <stdarg.h>
#include <string.h>

static inline pid_t gettid() { return syscall(SYS_gettid); }

/* On linux (well, using nptl and not linuxthreads, but we don't care
 * much), that's doable by parsing /proc/<tid>/stat
 */

/* statfields is fairly crude, for sure. ... has to be terminated by -1.
 * Arguments must come in groups of three, first the field index, then
 * its parsing format, then the adress of the pointer.
 */
static int statfields(pid_t t, ...)
{
    int nparsed = 0;
    va_list ap;
    va_start(ap, t);
    char * tmp;
    int rc = asprintf(&tmp,"/proc/%d/stat",t);
    if (rc < 0) return 0;
    char buf[1024];
    FILE * f = fopen(tmp,"r");
    char * s = fgets(buf, sizeof(buf), f);
    if (s == NULL) return 0;
    int j = va_arg(ap, int);
    for(int i = 0 ; j != -1 ; i++) {
        char * next;
        next = strchr(s, ' ');
        if (next == NULL)
            break;
        *next++='\0';
        if (i == j) {
            const char * fmt = va_arg(ap, const char *);
            void * ptr = va_arg(ap, void *);
            nparsed += sscanf(s, fmt, ptr);
            j = va_arg(ap, int);
        }
        s = next;
    }

    fclose(f);
    free(tmp);
    va_end(ap);

    return nparsed;
}

void thread_seconds_user_sys(double * res)
{
    unsigned long utime = 0;
    unsigned long stime = 0;
    statfields(gettid(), 13, "%lu", &utime, 14, "%lu", &stime, -1);
    res[0] = (double) utime / sysconf(_SC_CLK_TCK);
    res[1] = (double) stime / sysconf(_SC_CLK_TCK);
}

#else   /* Otherwise we'll do something stupid */
void thread_seconds_user_sys(double * res)
{
    seconds_user_sys(res); /* really stupid */
}
#endif

#ifdef HAVE_GETRUSAGE
typedef std::multimap<std::string, struct rusage> real_timingstats_dict_t;

void timingstats_dict_init(timingstats_dict_ptr p)
{
    *p = static_cast<void*>(new real_timingstats_dict_t());
}

void timingstats_dict_clear(timingstats_dict_ptr p)
{
    delete static_cast<real_timingstats_dict_t*>(*p);
}

void timingstats_dict_add(timingstats_dict_ptr p, const char * key, struct rusage * r)
{
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    real_timingstats_dict_t& s(*static_cast<real_timingstats_dict_t*>(*p));
    s.insert(std::make_pair(std::string(key), *r));
    pthread_mutex_unlock(&lock);
}

void timingstats_dict_add_mythread(timingstats_dict_ptr p, const char * key)
{
    struct rusage ru[1];
#ifdef HAVE_RUSAGE_THREAD
    getrusage (RUSAGE_THREAD, ru);
#else
    /* this will be plain bogus, but that's life */
    getrusage (RUSAGE_SELF, ru);
#endif
    timingstats_dict_add(p, key, ru);
}

void timingstats_dict_add_myprocess(timingstats_dict_ptr p, const char * key)
{
    struct rusage ru[1];
    getrusage (RUSAGE_SELF, ru);
    timingstats_dict_add(p, key, ru);
}

void timingstats_dict_disp(timingstats_dict_ptr p)
{
    real_timingstats_dict_t const& s(*static_cast<real_timingstats_dict_t*>(*p));
    /* multimap is sorted */
    typedef real_timingstats_dict_t::const_iterator it_t;
    for(it_t i = s.begin(), j ; i != s.end() ; i = j) {
        double tu = 0;
        double ts = 0;
        int n = 0;
        for(j = i ; j != s.end() && j->first == i->first ; j++, n++) {
            tu += (double)j->second.ru_utime.tv_sec
                + (double)j->second.ru_utime.tv_usec / 1.0e6;
            ts += (double)j->second.ru_stime.tv_sec
                + (double)j->second.ru_stime.tv_usec / 1.0e6;
        }
        printf("%s: %d process%s, total %.2fs+%.2fs on cpu\n",
                i->first.c_str(), n, n>1 ? "es" : "", tu, ts);
    }
}
#endif
