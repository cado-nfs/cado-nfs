#ifndef CADO_LAS_SMALLSIEVE_LOWLEVEL_HPP
#define CADO_LAS_SMALLSIEVE_LOWLEVEL_HPP

#include <cstddef>

#include "las-config.hpp"
#include "macros.h"

/* About row0_is_oddj: in order to check whether a j coordinate is even,
 * we need to take into account the bucket number, especially in case
 * buckets are as large as the sieve region. The row number corresponding
 * to a given i0 is i0/I, but we also need to add bucket_nr*bucket_size/I
 * to this, which is what this flag is for.  Sublat must also be taken
 * into account.
 */

/* This is copied from LOGNORM_FILL_COMMON_DEFS in las-norms.cpp ; from
 * logI, N, and LOG_BUCKET_REGION, define the integers i0, i1, j0, j1,
 * and I.
 */

#define SMALLSIEVE_COMMON_DEFS()                                         \
    const unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);\
    const unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);\
    const unsigned int regions_per_line = 1 << log_regions_per_line;           \
    const unsigned int region_rank_in_line = N & (regions_per_line - 1);       \
    const bool last_region_in_line MAYBE_UNUSED = region_rank_in_line == (regions_per_line - 1); \
    const unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    \
    const unsigned int j1 MAYBE_UNUSED = j0 + (1 << log_lines_per_region);    \
    const int I = 1 << logI;                                            \
    const int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          \
    const int i1 MAYBE_UNUSED = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     \
    /* those are (1,0,0) in the standard case */                        \
    const int sublatm MAYBE_UNUSED = sl.m ? sl.m : 1; \
    const unsigned int sublati0 MAYBE_UNUSED = sl.i0;       \
    const unsigned int sublatj0 MAYBE_UNUSED = sl.j0;       \
    const int row0_is_oddj MAYBE_UNUSED = (j0*sublatm + sublatj0) & 1;  \
    bool has_haxis = !j0;                                               \
    bool has_vaxis = region_rank_in_line == ((regions_per_line-1)/2);   \
    bool has_origin MAYBE_UNUSED = has_haxis && has_vaxis;              \
    do {} while (0)

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
/* {{{ preprocessor macro to generate assembly code for the small sieve */
/* {{{ Comments on the _OLDLOOP version (AF)
   0. The C code and the asm X86 code have the same algorithm.
   Read first the C code to understand easily the asm code.
   1. If there are less than 12 "T" in the line, the goal is to do
   only one jump (pipe-line breaks) and of course the instructions
   minimal number.
   2. 12 "T" seems the best in the critical loop. Before, gcc tries
   to optimize in a bad way the loop. For gcc generated code, the
   best is here a systematic code (12*2 instructions), like the C code.
   3. The asm X86 optimization uses addb logp,(pi+p_or_2p*[0,1,2])
   & three_p_or_2p = 3 * p_or_2p; the lea (add simulation) & real
   add interlace seems a bit interesting.
   So, the loop is smaller & faster (19 instructions versus 27 for
   gcc best X86 asm).
   4. Of course, the gain between the 2 versions is light, because
   the main problem is the access time of the L0 cache: read + write
   with sieve_increase(pi,logp,w), or *pi += logp in fact.
   }}} */
/* define sequences of operations rather systematically.
 *
 * ADD_NOINCR   *pi += logp.
 * ADD_INCR     *pi += logp; pi += incr
 * XADD_NOINCR  if (pi >= fence) goto exit; *pi += logp;
 * XADD_INCR    if (pi >= fence) goto exit; *pi += logp; pi += incr;
 */
/* {{{ *_ADD_NOINCR */
#define ONE_ADD_NOINCR(pi, logp)                                     	\
    "addb " logp ",(" pi ")\n"                /* pi[0] += logp */
#define TWO_ADDS_NOINCR(pi, incr, logp)                 		\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */
#define THREE_ADDS_NOINCR(pi, incr, logp)                 		\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */    	\
    ONE_ADD_NOINCR(pi "," incr ",2", logp) /* pi[2*incr] += logp */
/* }}} */
/* {{{ *_ADDS_INCR, including _PRECOMP variant */
#define ONE_ADD_INCR(pi, incr, logp)                             	\
    ONE_ADD_NOINCR(pi, logp)                                         	\
    "lea (" pi "," incr ",1)," pi "\n"        /* pi += incr */
#define TWO_ADDS_INCR(pi, incr, logp)                 			\
    ONE_ADD_NOINCR(pi, logp)               /* pi[0] += logp */       	\
    ONE_ADD_NOINCR(pi "," incr ",1", logp) /* pi[incr] += logp */    	\
    "lea (" pi "," incr ",2)," pi "\n"     /* pi += 2*incr */
/* The _PRECOMP variant uses a precomputed  bigincr = 3 * incr  */
#define THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                \
     THREE_ADDS_NOINCR(pi, incr, logp)                                  \
    "lea (" pi "," bigincr ",1)," pi "\n"     /* pi += 3*incr */
/* {{{ expand to define new lengths */
#define FOUR_ADDS_INCR(pi, incr, logp)                                  \
    TWO_ADDS_INCR(pi, incr, logp)                                       \
    TWO_ADDS_INCR(pi, incr, logp)
#define EIGHT_ADDS_INCR(pi, incr, logp)                                 \
    FOUR_ADDS_INCR(pi, incr, logp)                                      \
    FOUR_ADDS_INCR(pi, incr, logp)
#define TWELVE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)               \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)                    \
    THREE_ADDS_INCR_PRECOMP(pi, incr, bigincr, logp)
#define SIXTEEN_ADDS_INCR(pi, incr, logp)                               \
    EIGHT_ADDS_INCR(pi, incr, logp)                                     \
    EIGHT_ADDS_INCR(pi, incr, logp)
/* }}} */
/* }}} */
/* {{{ caution: *_XADDS_INCR */
#define ONE_XADD_NOINCR(pi, fence, logp, exit)               	\
    "cmp " fence ", " pi "\n"                 /* if (pi >= S1) break; */\
    "jae " exit "\n"                                                    \
    ONE_ADD_NOINCR(pi, logp)
#define ONE_XADD_INCR(pi, incr, fence, logp, exit)       		\
    ONE_XADD_NOINCR(pi, fence, logp, exit)                   		\
    "lea (" pi "," incr ",1)," pi "\n"        /* pi += incr */
#define TWO_XADDS_INCR(pi, incr, fence, logp, exit)             	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)                  	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)
/* no need for THREE_XADDs */
/* {{{ expand to define new lengths */
#define FOUR_XADDS_INCR(pi, incr, fence, logp, exit)     		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)          		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)
#define EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)    		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)         		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)
#define TWELVE_XADDS_INCR(pi, incr, fence, logp, exit) 			\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)        		\
    FOUR_XADDS_INCR(pi, incr, fence, logp, exit)
#define SIXTEEN_XADDS_INCR(pi, incr, fence, logp, exit)    		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)         		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)
/* }}} */
/* }}} */
/* {{{ We also have _NOLASTINCR variants. Those are used in the old
 * assembly loop. As a matter of fact, it is not a good idea to optimize
 * like this: we do need the end value of the pointer, and it's really
 * crucially important for the case there I>B. So using these functions
 * is part of the reason why the old code will not do in the longer term.
 */
#define TWO_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)       	\
    ONE_XADD_INCR(pi, incr, fence, logp, exit)                  	\
    ONE_XADD_NOINCR(pi, fence, logp, exit)
#define FOUR_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)   		\
    TWO_XADDS_INCR(pi, incr, fence, logp, exit)          		\
    TWO_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)
#define TWELVE_XADDS_NOLASTINCR(pi, incr, fence, logp, exit) 		\
    EIGHT_XADDS_INCR(pi, incr, fence, logp, exit)        		\
    FOUR_XADDS_NOLASTINCR(pi, incr, fence, logp, exit)
#define SMALLSIEVE_ASSEMBLY_OLD(_pi, _incr, _fence, _logp) do {             \
    unsigned char *dummy;                                              \
    size_t three_p_or_2p;                                               \
    __asm__ __volatile__ (                                              \
            "lea (%[incr],%[incr],2), %[three_p_or_2p]\n"     /* three_p_or_2p = p_or_2p * 3 */ \
            "lea (%[pi],%[three_p_or_2p],4), %[dummy]\n"     /* pi_end = pi + p_or_2p*12 */    \
            "cmp %[fence], %[dummy]\n"            /* if (pi_end > S1) no loop */    \
            "ja 1f\n"                                                  \
            ".balign 8\n"                                               \
            "0:\n"                                                      \
            TWELVE_ADDS_INCR_PRECOMP("%[pi]","%[incr]","%[three_p_or_2p]","%[logp]")               \
            "lea (%[pi],%[three_p_or_2p],4), %[dummy]\n"    /* if (pi+p_or_2p*12 > S1) break */\
            "cmp %[fence], %[dummy]\n"                                              \
            "jbe 0b\n"                                                  \
            "1:\n"                                                      \
            TWELVE_XADDS_INCR("%[pi]","%[incr]","%[fence]","%[logp]","2f")    \
            "2:\n"                                                      \
            : [dummy] "=&r"(dummy), [pi] "+r"(_pi), [three_p_or_2p]"=&r"(three_p_or_2p)            \
            : [incr]"r"(_incr), [_logp]"q"(logp), [fence]"r"(_fence) : "cc");                 \
} while (0)
/* }}} */
#if 0 /* {{{ unused GCC trick */
/* there is technically a way to drop the "volatile" qualifier, so that
 * gcc understands that we're touching memory, and which memory exactly.
 * Such a way can be found in the snippet below.
 * Alas, I don't know how it is possible to do so and at the same time
 * make sure that: 1 - the memory reference that is passed to the inline asm
 * is strictly of the form (%[register]), and 2 - we can actually get the
 * register back, or emit code that is similar to (%[register],%[displ])
 * from there.
 * And anyway, I don't really believe that it's important: this block is
 * really the expensive operation, what happens to the glue around it is
 * less important.
 */
#define SMALLSIEVE_ASSEMBLY_NEW(S0, pos, incr, fence, logp) do {        \
    /* see gcc 6.45.2.6 Clobbers */                                     \
    typedef unsigned char (*ss_area_t)[1<<16];                          \
    ss_area_t S0x = (ss_area_t) S0;                                     \
    asm (                                                               \
            "leaq (%[S0],%[pos]), %[pi]\n"                              \
            ONE_ADD_INCR("%[pi]","%[inc]", "%[logp]")            \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            ONE_XADD_INCR("%[pi]", "%[inc]", "%[ff]", "%[logp]", "0f") \
            "0:\n"                                                      \
            : [S0]"=m"(*S0x), [pi]"=&r"(pi)                             \
            : [pos]"r"(pos),                                            \
            [inc]"r"(incr), [logp]"q"(logp), [ff] "r" (fence));         \
} while (0)
#endif /* }}} */
/* {{{ asm volatile statements for all the variants */
/* {{{ unrolled statements, limited length */
#define SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)          \
            "0:\n"                                                      \
            : [pi]"+r"(_pi)                                             \
            : [inc]"r"(_incr), [logp]"q"(_logp), [ff] "r" (_fence)
#define SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_ADD_INCR("%[pi]","%[inc]","%[logp]")                    \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            TWO_ADDS_INCR("%[pi]","%[inc]","%[logp]")                   \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            FOUR_ADDS_INCR("%[pi]","%[inc]","%[logp]")                  \
            FOUR_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(_pi, _incr, _fence, _logp)      \
    asm volatile (                                                      \
            EIGHT_ADDS_INCR("%[pi]","%[inc]","%[logp]")                 \
            EIGHT_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(_pi, _incr, _fence, _logp)     \
    asm volatile (                                                      \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
/* this is the point where the complete unrolling starts being excessive */
#define SMALLSIEVE_ASSEMBLY_NEW_32_TO_64(_pi, _incr, _fence, _logp)     \
    asm volatile (                                                      \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
/* }}} */
/* {{{ loop8 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[pi], %[inc], 8), %%rdi\n"                          \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            EIGHT_ADDS_INCR("%[pi]","%[inc]","%[logp]")                 \
            "leaq (%%rdi, %[inc], 8), %%rdi\n"                          \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            EIGHT_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")   \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            : "rdi")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8(_pi, _incr, _fence, _logp)        \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8P0(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP8P1(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP8_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi, %[inc], 8)\n")
/* }}} */
/* {{{ loop12 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[inc], %[inc], 2), %%r10\n"                         \
            "leaq (%[pi], %%r10, 4), %%rdi\n"                           \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            TWELVE_ADDS_INCR_PRECOMP("%[pi]","%[inc]","%%r10","%[logp]")\
            "leaq (%%rdi, %%r10, 4), %%rdi\n"                           \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            TWELVE_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")  \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            : "rdi", "r10")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12(_pi, _incr, _fence, _logp)       \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P0(_pi, _incr, _fence, _logp)     \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P1(_pi, _incr, _fence, _logp)    \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi, %%r10, 4)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP12P2(_pi, _incr, _fence, _logp)    \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP12_COMMON(_pi, _incr, _fence, _logp,\
                "prefetch  (%%rdi, %%r10, 8)\n")
/* }}} */
/* {{{ loop16 */
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, _PP) \
    asm volatile (                                                      \
            "leaq (%[pi], %[inc2], 8), %%rdi\n"                         \
            "cmpq %[ff], %%rdi\n"                                       \
            "ja 1f\n"           /* if rdi > fence, no loop */           \
            ".balign 8\n"                                               \
            "2:\n"                                                      \
            _PP                                                         \
            SIXTEEN_ADDS_INCR("%[pi]","%[inc]","%[logp]")               \
            "leaq (%%rdi, %[inc2], 8), %%rdi\n"                         \
            "cmpq %[ff], %%rdi\n"                                       \
            "jbe 2b\n"                                                  \
            "1:\n"                                                      \
            SIXTEEN_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f") \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)      \
            , [inc2]"r"(2*(_incr)) : "rdi")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16(_pi, _incr, _fence, _logp)        \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16P0(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi)\n")
#define SMALLSIEVE_ASSEMBLY_NEW_LOOP16P1(_pi, _incr, _fence, _logp)      \
        SMALLSIEVE_ASSEMBLY_NEW_LOOP16_COMMON(_pi, _incr, _fence, _logp, \
                "prefetch  (%%rdi, %[inc2], 8)\n")
/*}}}*/
/* {{{ the routines below are unused. The might be, in case we insist on
 * using the same code for both even and odd lines (but then it would
 * probably be a loop anyway) */
#define SMALLSIEVE_ASSEMBLY_NEW_0_TO_1(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_0_TO_2(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")     \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp))
#define SMALLSIEVE_ASSEMBLY_NEW_0_TO_4(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            FOUR_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")    \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)) 
#define SMALLSIEVE_ASSEMBLY_NEW_1_TO_4(_pi, _incr, _fence, _logp)       \
    asm volatile (                                                      \
            ONE_ADD_INCR("%[pi]","%[inc]","%[logp]")                    \
            TWO_XADDS_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")     \
            ONE_XADD_INCR("%[pi]","%[inc]","%[ff]","%[logp]","0f")      \
            SMALLSIEVE_ASSEMBLY_TRAILER(_pi, _incr, _fence, _logp)) 
/*}}}*/
/*}}}*/
/* }}} */
#endif

/*{{{ function objects for all the routines that we have */
#define BEGIN_FOBJ(name_, compatibility_condition)       		\
    struct name_ {							\
        template<int b> struct is_compatible {                          \
            static_assert(compatibility_condition,                      \
                    # name_ " requires " # compatibility_condition);    \
            static const bool value = compatibility_condition;          \
        };                                                              \
        static constexpr const char * name = # name_;			\
        inline size_t operator()(					\
                unsigned char * S0,					\
                unsigned char * S1,					\
                size_t x0 MAYBE_UNUSED,					\
                size_t pos,						\
                size_t p_or_2p,						\
                unsigned char logp,					\
                where_am_I w MAYBE_UNUSED) const			\
        {								\
            unsigned char * pi = S0 + pos;                              \
            INTERMEDIARY_FOBJ()

#ifndef TRACE_K
#define INTERMEDIARY_FOBJ()       do {} while (0)
#else
#define INTERMEDIARY_FOBJ() do {					\
            if (trace_on_range_Nx(w->N, x0 + pos, x0 + S1 - S0)) {	\
                if ((trace_Nx.x - x0 - pos) % p_or_2p == 0) {           \
                    WHERE_AM_I_UPDATE(w, x, trace_Nx.x);		\
                    sieve_increase_logging(S0 + w->x - x0, logp, w);	\
                }							\
            }								\
} while (0)
#endif

#define END_FOBJ()							\
    return pi - S1;							\
}}
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
/* {{{ function objects for assembly code */
BEGIN_FOBJ(assembly_generic_loop8, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop8p0, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p0, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16p0, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16P0(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop8p1, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP8P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p1, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop16p1, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP16P1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_loop12p2, true);
SMALLSIEVE_ASSEMBLY_NEW_LOOP12P2(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly_generic_oldloop, true);
SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly0, b <= 0);
SMALLSIEVE_ASSEMBLY_NEW_0_TO_1(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly1, b == 1);
SMALLSIEVE_ASSEMBLY_NEW_1_TO_2(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly1x, b <= 1);
SMALLSIEVE_ASSEMBLY_NEW_0_TO_2(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly2, b == 2);
SMALLSIEVE_ASSEMBLY_NEW_2_TO_4(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly2x, b <= 2);
SMALLSIEVE_ASSEMBLY_NEW_0_TO_4(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly3, b == 3);
SMALLSIEVE_ASSEMBLY_NEW_4_TO_8(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly4, b == 4);
SMALLSIEVE_ASSEMBLY_NEW_8_TO_16(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly5, b == 5);
SMALLSIEVE_ASSEMBLY_NEW_16_TO_32(pi, p_or_2p, S1, logp);
END_FOBJ();
BEGIN_FOBJ(assembly6, b == 6);
SMALLSIEVE_ASSEMBLY_NEW_32_TO_64(pi, p_or_2p, S1, logp);
END_FOBJ();
/* }}} */
#endif
/* {{{ Function objects for C code, manually unrolled except for the
 * _loop versions. */
BEGIN_FOBJ(manual0, b <= 0);
if (pi < S1) { *pi += logp; pi += p_or_2p; }
END_FOBJ();
BEGIN_FOBJ(manual1, b == 1);
do {
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual2, b == 2);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual3, b == 3);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();
BEGIN_FOBJ(manual4, b == 4);
do {
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; pi += p_or_2p;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; if ((pi += p_or_2p) >= S1) break;
    *pi += logp; pi += p_or_2p;
} while (0);
END_FOBJ();

/* Gotcha: we're now preventively tracing the small sieve events outside
 * the asm code, which causes sieve_increase_logging_backend to be
 * called, and in turn test_divisible, which modifies traced norms (in
 * debug mode). So for consistency, we mustn't use sieve_increase here,
 * or the norm would be divided by p twice.
 */
BEGIN_FOBJ(manual_oldloop, true);
#define T do {                                                          \
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);                              \
    *pi += logp; /* sieve_increase (pi, logp, w); */ pi += p_or_2p;     \
} while(0)
while (UNLIKELY(pi + p_or_2p * 12 <= S1))
{ T; T; T; T; T; T; T; T; T; T; T; T; }
do {
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
    if (pi >= S1) break;
    T;
} while (0);
#undef T
END_FOBJ();
BEGIN_FOBJ(manual_oldloop_nounroll, true);
for ( ; pi < S1 ; pi += p_or_2p) {
    WHERE_AM_I_UPDATE(w, x, x0 + pi - S0);
    *pi += logp; // sieve_increase (pi, logp, w);
}
END_FOBJ();
/*}}}*/
/*}}}*/

/* {{{ our selection of best routines indexed by ceil(log2(I/p)) */
template<int bit> struct best_evenline;
template<int bit> struct best_oddline ;
#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && !defined(TRACK_CODE_PATH)
template<int bit> struct best_evenline {
    typedef int is_default;
    typedef manual_oldloop type;
};
template<int bit> struct best_oddline  {
    typedef int is_default;
    typedef manual_oldloop type;
};
template<> struct best_evenline<1> { typedef manual0 type; };
template<> struct best_evenline<2> { typedef assembly1 type; };
template<> struct best_evenline<3> { typedef assembly2 type; };
template<> struct best_evenline<4> { typedef assembly3 type; };
template<> struct best_evenline<5> { typedef assembly4 type; };
template<> struct best_evenline<6> { typedef assembly5 type; };

template<> struct best_oddline<1> { typedef assembly1 type; };
template<> struct best_oddline<2> { typedef assembly2 type; };
template<> struct best_oddline<3> { typedef assembly3 type; };
template<> struct best_oddline<4> { typedef assembly4 type; };
template<> struct best_oddline<5> { typedef assembly_generic_loop8 type; };
template<> struct best_oddline<6> { typedef assembly_generic_loop8 type; };

#if 1
template<> struct best_evenline<7> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<8> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<9> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<10> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<11> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<12> { typedef assembly_generic_loop12 type; };
template<> struct best_evenline<13> { typedef assembly_generic_loop16p0 type; };
template<> struct best_evenline<14> { typedef assembly_generic_loop16p1 type; };
template<> struct best_evenline<15> { typedef assembly_generic_loop16p0 type; };

template<> struct best_oddline<7> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<8> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<9> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<10> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<11> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<12> { typedef assembly_generic_loop12 type; };
template<> struct best_oddline<13> { typedef assembly_generic_loop16p0 type; };
template<> struct best_oddline<14> { typedef assembly_generic_loop16p1 type; };
template<> struct best_oddline<15> { typedef assembly_generic_loop16p0 type; };
typedef assembly_generic_oldloop default_smallsieve_inner_loop;
#endif
/* TODO: we could perhaps provide hints so that the generated code can
 * merge some rounds of the loop. Or maybe the compiler will do that ? */
#else
template<int bit> struct best_evenline { typedef manual_oldloop type; };
template<int bit> struct best_oddline  { typedef manual_oldloop type; };
typedef manual_oldloop default_smallsieve_inner_loop;
#endif
/* }}} */


#endif	/* CADO_LAS_SMALLSIEVE_LOWLEVEL_HPP */
