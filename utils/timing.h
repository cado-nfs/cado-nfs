#ifndef CADO_UTILS_TIMING_H_
#define CADO_UTILS_TIMING_H_

#include "macros.h"     /* MAYBE_UNUSED */
#include <stdint.h> /* for uint64_t */

#ifdef  HAVE_GETRUSAGE
#include <sys/resource.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
static inline uint64_t cputicks()
{
        uint64_t r;
        __asm__ __volatile__(
                "rdtsc\n\t"
                "shlq $32, %%rdx\n\t"
                "orq %%rdx, %%rax\n\t"
                : "=a"(r)
                :
                : "rdx");
        return r;
}
#endif

extern uint64_t microseconds (void);
extern uint64_t microseconds_thread (void);
extern unsigned long milliseconds (void);
extern unsigned long milliseconds_thread (void);
extern double seconds (void);
extern double seconds_thread (void);
extern void seconds_user_sys (double *);
extern double wct_seconds (void);
extern void print_timing_and_memory (FILE*, double);
extern void thread_seconds_user_sys(double *);

/* we provide an interface to collect the timings of the subprocesses.
 * Each subprocess is identified in the filter_io layer by some key, and
 * the caller has to call the proper function to display the tally of the
 * timings.
 *
 * This could be expanded to collect more versatile timing kinds.
 */

typedef void * timingstats_dict_t[1];
typedef void ** timingstats_dict_ptr;
#ifdef HAVE_GETRUSAGE   /* since this is only rusage-oriented... */
void timingstats_dict_init(timingstats_dict_ptr);
void timingstats_dict_clear(timingstats_dict_ptr);
/* display the tally */
void timingstats_dict_disp(timingstats_dict_ptr);
void timingstats_dict_add(timingstats_dict_ptr, const char * key, struct rusage * r);
void timingstats_dict_add_mythread(timingstats_dict_ptr, const char * key);
void timingstats_dict_add_myprocess(timingstats_dict_ptr, const char * key);
#else   /* provide all of these as no-ops */
static inline void timingstats_dict_init(timingstats_dict_ptr x MAYBE_UNUSED) {}
static inline void timingstats_dict_clear(timingstats_dict_ptr x MAYBE_UNUSED) {}
static inline void timingstats_dict_disp(timingstats_dict_ptr x MAYBE_UNUSED) {}
// sta inline ic void timingstats_dict_add(timingstats_dict_ptr x MAYBE_UNUSED, const char * key MAYBE_UNUSED, void * r MAYBE_UNUSED) {}
static inline void timingstats_dict_add_mythread(timingstats_dict_ptr x MAYBE_UNUSED, const char * key MAYBE_UNUSED) {}
static inline void timingstats_dict_add_myprocess(timingstats_dict_ptr x MAYBE_UNUSED, const char * key MAYBE_UNUSED) {}
#endif

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_TIMING_H_ */
