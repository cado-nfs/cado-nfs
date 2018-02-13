/*
 * rdtsc.h
 *
 *  Created on: 09.08.2016
 *      Author: kruppa
 */

#ifndef RDTSC_H_
#define RDTSC_H_

#include <stdint.h>
#include "macros.h"

/* Enable at most one of these three defines, preferably in the including
   compilation unit */
// #define USE_INTEL_PCM 1
// #define USE_PERF 1
// #define USE_JEVENTS 1

#include <stdlib.h>

uint64_t u32_to_64(uint32_t low, uint32_t high)
{ return ((uint64_t) high << 32) + (uint64_t) low; }

#ifdef USE_INTEL_PCM

#include <sched.h>
#include "cpucounters.h"
static PCM * m;
static CoreCounterState before_sstate, after_sstate;

#elif defined USE_PERF

#include "libperf.h"  /* standard libperf include */
static struct libperf_data* pd;
static uint64_t start_time, end_time;

#elif USE_JEVENTS

#ifdef __cplusplus
extern "C" {
#endif
#include "rdpmc.h"
#ifdef __cplusplus
}
#endif
struct rdpmc_ctx ctx;
#define IMMEDIATE_RDPMC 1
#ifdef IMMEDIATE_RDPMC
struct {unsigned lo, hi;} start_time;
static unsigned long end_time;
#else
static uint64_t start_time, end_time;
#endif

#else

#define START_COUNTER rdtscpl
#define END_COUNTER rdtscpl
static uint64_t start_time, end_time;

#endif

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void serialize()
{
    unsigned long id = 0;
    __asm__ volatile("CPUID\n\t"
		     : "+a" (id)
		     :: "%rbx", "%rcx", "%rdx"
	);
}

/* Use RDTSC to read CPU core clock cycles except that on contemporary CPUs,
   it actually reads some natural-time counter */
__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void rdtsc(uint32_t *low, uint32_t *high)
{
    __asm__ volatile("RDTSC\n\t"
		     : "=d" (*high), "=a" (*low)
	);
}

/* Use RDTSCP to serialize, then read TSC. Instructions following RDTSCP may
   start to execute before RDTSCP finishes. */
__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void rdtscp(uint32_t *low, uint32_t *high)
{
    __asm__ volatile("RDTSCP\n\t"
		     : "=d" (*high), "=a" (*low)
		     :: "%rcx"
	);
}

/* Read a performance measurement counter */
__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void
rdpmc(uint32_t *low, uint32_t *high, const unsigned int selector)
{
    __asm__ volatile("rdpmc" : "=a" (*low), "=d" (*high) : "c" (selector));
}

#else /* defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) */

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void serialize() {}

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void rdtsc(uint32_t *low, uint32_t *high)
{ *low = *high = 0; }

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void rdtscp(uint32_t *low, uint32_t *high)
{ *low = *high = 0; }

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void
rdpmc(uint32_t *low, uint32_t *high, const unsigned int selector MAYBE_UNUSED)
{ *low = *high = 0; }

#endif

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline uint64_t rdtscl()
{
    uint32_t high, low;
    rdtsc(&low, &high);
    return u32_to_64(low, high);
}

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline uint64_t rdtscpl()
{
    uint32_t high, low;
    rdtscp(&low, &high);
    return u32_to_64(low, high);
}

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline uint64_t rdpmcl(const unsigned int selector)
{
    uint32_t high, low;
    rdpmc(&low, &high, selector);
    return u32_to_64(low, high);
}

// rdpmcl_cycles uses a "fixed-function" performance counter to return
// the count of actual CPU core cycles executed by the current core.
// Core cycles are not accumulated while the processor is in the "HALT"
// state, which is used when the operating system has no task(s) to run
// on a processor core.
// Note that this counter continues to increment during system calls
// and task switches. As such, it may be unreliable for timing long
// functions where the CPU may serve an interrupt request or where
// the kernel may preempt execution and switch to another process.
// It is best used for timing short intervals which usually run
// uninterrupted, and where occurrences of interruption are easily
// detected by an abnormally large cycle count.

// The RDPMC instruction must be enabled for execution in user-space.
// This requires a total of three bits to be set in CR4 and MSRs of
// the CPU:
// Bit 1<<8 in CR4 must be set to 1. On Linux, this can be effected by
// executing as root:
//   echo 1 >> /sys/devices/cpu/rdpmc
// Bit 1<<33 must be set to 1 in the MSR_CORE_PERF_GLOBAL_CTRL
// (MSR address 0x38f). This enables the cycle counter so that it
// actually increment with each clock cycle; while this bit is 0,
// the counter value stays fixed.
// Bit 1<<5 must be set to 1 in the MSR_CORE_PERF_FIXED_CTR_CTRL
// (MSR address 0x38d) which causes the counter to be incremented
// on non-halted clock cycles that occur while the CPL is >0
// (user-mode). If bit 1<<4 is set to 1, then the counter will
// increment while CPL is 0 (kernel mode), e.g., during interrupts,
// etc.

// The only reliable way I found to enable all these bits is through
// the JEVENTS library which uses the PERF kernel interface which lets
// it know that the bits are supposed to stay on. Doing it manually
// with, e.g., a kernel module causes the kernel to clear the bits
// again next time it updates any of these CR/MSR.


__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline void
rdpmc_cycles(uint32_t *low, uint32_t *high)
{
    const unsigned c = (1U<<30) + 1; /* Second Fixed-function counter:
                                        clock cycles in non-HALT */
    rdpmc(low, high, c);
}

__attribute__((__unused__, __always_inline__)) ATTRIBUTE_ARTIFICIAL
static inline uint64_t rdpmcl_cycles()
{
    uint32_t low, high;
    rdpmc_cycles(&low, &high);
    return u32_to_64(low, high);
}

static inline void init_timing()
{
#ifdef USE_INTEL_PCM
    printf("# Using Intel PCM library\n");
    m = PCM::getInstance();
    if (m->program() != PCM::Success) {
	printf("Could not initialise PCM\n");
	exit(EXIT_FAILURE);
    }

    const bool have_smt = m->getSMT();
    if (have_smt) {
	printf("# CPU uses SMT\n");
    } else {
	printf("# CPU does not use SMT\n");
  }
#elif defined(USE_PERF)
    printf("# Using PERF library\n");
    pd = libperf_initialize(-1,-1); /* init lib */
    libperf_enablecounter(pd, LIBPERF_COUNT_HW_CPU_CYCLES);
    /* enable HW counter */
#elif defined(USE_JEVENTS)
    printf("# Using jevents library\n");
#ifdef TIMING_SERIALIZE
    printf("# Serializing with CPUID before timing start and with RDTSCP at "
	   "timing end\n");
#endif
    if (rdpmc_open(PERF_COUNT_HW_CPU_CYCLES, &ctx) < 0)
	exit(EXIT_FAILURE);
#else
    printf("# Using " CADO_STRINGIZE(START_COUNTER) " to start and "
	   CADO_STRINGIZE(END_COUNTER) " to end measurement\n");
#endif
}

static inline void clear_timing()
{
#ifdef USE_INTEL_PCM
    m->cleanup();
#elif defined(USE_PERF)
    libperf_close(pd);
    pd = NULL;
#elif defined(USE_JEVENTS)
    rdpmc_close (&ctx);
#else
    /* Timing with RDTSC[P] does not allocate anything */ 
#endif
}

static inline void start_timing()
{
#ifdef USE_INTEL_PCM
    const int cpu = sched_getcpu();
    before_sstate = getCoreCounterState(cpu);
#elif defined(USE_PERF)
    start_time = libperf_readcounter(pd, LIBPERF_COUNT_HW_CPU_CYCLES);
#elif defined(USE_JEVENTS)
#ifdef TIMING_SERIALIZE
    serialize();
#endif
#ifdef IMMEDIATE_RDPMC
    rdpmc_cycles(&start_time.lo, &start_time.hi);
#else
    start_time = rdpmc_read(&ctx);
#endif
#else
    start_time = START_COUNTER();
#endif
}

static inline void end_timing()
{
#ifdef USE_INTEL_PCM
    const int cpu = sched_getcpu();
    after_sstate = getCoreCounterState(cpu);
#elif defined(USE_PERF)
    end_time = libperf_readcounter(pd, LIBPERF_COUNT_HW_CPU_CYCLES);
#elif defined(USE_JEVENTS)
#ifdef TIMING_SERIALIZE
    rdtscpl();
#endif
#ifdef IMMEDIATE_RDPMC
    end_time = rdpmcl_cycles();
#else
    end_time = rdpmc_read(&ctx);
#endif
#else
    end_time = END_COUNTER();
#endif
}

static inline uint64_t get_diff_timing()
{
#ifdef USE_INTEL_PCM
    return getCycles(before_sstate,after_sstate);
#elif defined(USE_PERF)
    return end_time - start_time;
#elif defined(USE_JEVENTS)
#ifdef IMMEDIATE_RDPMC
    return end_time - u32_to_64(start_time.lo, start_time.hi);
#else
    return end_time - start_time;
#endif
#else
    return end_time - start_time;
#endif

}

#ifdef USE_JEVENTS
static uint64_t pmc0, pmc1, pmc2, pmc3;
#endif

void
readAllPmc()
{
#ifdef USE_JEVENTS
    pmc0 = rdpmcl(0);
    pmc1 = rdpmcl(1);
    pmc2 = rdpmcl(2);
    pmc3 = rdpmcl(3);
#endif
}

void
diffAllPmc()
{
#ifdef USE_JEVENTS
    pmc0 = rdpmcl(0) - pmc0;
    pmc1 = rdpmcl(1) - pmc1;
    pmc2 = rdpmcl(0) - pmc2;
    pmc3 = rdpmcl(1) - pmc3;
#endif
}

void
printAllPmc()
{
#ifdef USE_JEVENTS
    printf(" (pmc0: %lu, pmc1: %lu, pmc2: %lu, pmc3: %lu)", 
	   pmc0, pmc1, pmc2, pmc3);
#endif
}

#endif /* RDTSC_H_ */
