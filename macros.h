/* Common header file for the CADO project
 
Copyright 2007, 2008, 2009, 2010, 2011 Pierrick Gaudry, Alexander Kruppa,
                                       Emmanuel Thome, Paul Zimmermann

This file is part of the CADO project.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef CADO_MACROS_H
#define CADO_MACROS_H
// pragma no prototypes

/**********************************************************************/
/* Common asserting/debugging defines */
/* See README.macro_usage */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#ifdef __cplusplus
#include <stdexcept>
#endif
#include <gmp.h> /* for __GNU_MP_VERSION */

// scan-headers: stop here

#define ASSERT(x)	assert(x)

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#if defined(WANT_ASSERT_EXPENSIVE) || defined(STATIC_ANALYSIS)
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

#ifndef CPP_STRINGIFY
#define CPP_STRINGIFY0(x) #x
#define CPP_STRINGIFY(x) CPP_STRINGIFY0(x)
#endif
#ifndef CPP_PAD
#define CPP_PAD(x,y) x ## y
#endif

#define croak__(x,y)     						\
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",			\
                (x),__func__,__FILE__,__LINE__,(y))
#define croak_throw__(e, x)     					\
        throw e("code BUG() : condition " x            \
                " failed at " __FILE__ ":" CPP_STRINGIFY(__LINE__))

/* In C++ dtors which are not allowed to throw, use this variant instead.
 */
#define ASSERT_ALWAYS_NOTHROW(x)					\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)
// NOLINTBEGIN(readability-simplify-boolean-expr)
#ifdef __cplusplus
#define ASSERT_ALWAYS_OR_THROW(x, e)                                   \
    do {								\
        if (!(x)) 							\
            croak_throw__(e, #x);                                       \
    } while (0)
#define ASSERT_ALWAYS(x) ASSERT_ALWAYS_OR_THROW(x, std::runtime_error)
#else
#define ASSERT_ALWAYS(x) ASSERT_ALWAYS_NOTHROW(x)
#endif
// NOLINTEND(readability-simplify-boolean-expr)

/* never throw exceptions in that case, just exit */
#define FATAL_ERROR_CHECK(cond, msg)		        		\
    do {								\
      if (UNLIKELY((cond))) {                                           \
          croak__("Fatal error: ", msg);               			\
          abort();                                                      \
        }								\
    } while (0)

/* Note that string.h must be #included in order to use this macro */
#define DIE_ERRNO_DIAG(tst, fmt, ...) do {				\
    if (UNLIKELY(tst)) {				        	\
        fprintf(stderr, fmt ": %s\n", __VA_ARGS__, strerror(errno));    \
        abort();					        	\
    }							        	\
} while (0)

/* Note that string.h must be #included in order to use this macro */
#define WARN_ERRNO_DIAG(tst, fmt, ...) do {				\
    if (UNLIKELY(tst)) {				        	\
        fprintf(stderr, fmt ": %s\n", __VA_ARGS__, strerror(errno));    \
    }							        	\
} while (0)

/* This macro is used to guard against some trivial false positives
 * returned by static analyzer */
#if defined(__COVERITY__) || defined(STATIC_ANALYSIS)
#define ASSERT_FOR_STATIC_ANALYZER(x) do {                             \
    if (!(x)) {                                                        \
        abort();                                                       \
    }                                                                  \
} while (0)
#else
#define ASSERT_FOR_STATIC_ANALYZER(x)
#endif

#ifndef __cplusplus
/* c++: do not realloc at all, or if you really insist, use
 * checked_realloc from utils_cxx.hpp
 */
// NOLINTBEGIN(bugprone-macro-parentheses)
#define CHECKED_REALLOC(var, N, type)   do {				\
    if ((N) == 0) {							\
        free(var);							\
        (var) = NULL;							\
    } else {								\
        type * __p = (type *) realloc((var), (N) * sizeof(type));	\
        if (!__p && (var) != NULL)					\
            free((var));						\
        ASSERT_ALWAYS(__p != NULL);					\
        (var) = __p;							\
    }									\
} while (0)
// NOLINTEND(bugprone-macro-parentheses)
#endif


/*********************************************************************/
/* Helper macros */
/* See README.macro_usage */

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

/* Handy, and does not require libm */
#ifndef iceildiv
/* C99 defines division as truncating, i.e. rounding to 0. To get ceil in all
   cases, we have to distinguish by sign of arguments. Note that (afaik) the
   sign and truncation direction of the quotient with any negative operads was
   undefined in C89. */
/* iceildiv requires *unsigned* operands in spite of the name suggesting
   (signed) integer type. For negative operands, the result is wrong. We
   don't use the fairly common "(x + y - 1) / y" idiom because of issues
   with max values. Instead, we do (x-1)/y+1, with special adjustments
   for the x==0 case, which we want to do branch-less. */
// NOLINTBEGIN(readability-container-size-empty)
#define iceildiv(x,y) (((x)-((x)!=0))/(y)+((x)!=0))
/* siceildiv requires signed operands, or the compiler will throw warnings
   with -Wtype-limits */
#define siceildiv(x,y) ((x) == 0 ? 0 : ((x)<0) + ((y)<0) == 1 ? (x)/(y) : ((x)-1+2*((y)<0))/(y)+1)
// NOLINTEND(readability-container-size-empty)
#endif

/* Number of words holding B bits ; better naming sought. */
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))


#define LEXGE2(X,Y,A,B) ((X)>(A) || ((X) == (A) && (Y) >= (B)))
#define LEXGE3(X,Y,Z,A,B,C) ((X)>(A) || ((X) == (A) && LEXGE2((Y),(Z),(B),(C))))
#define LEXLE2(X,Y,A,B) LEXGE2((A),(B),(X),(Y))
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3((A),(B),(C),(X),(Y),(Z))

#ifndef GNUC_VERSION
#ifndef __GNUC__
#define GNUC_VERSION(X,Y,Z) 0
#else
#define GNUC_VERSION(X,Y,Z)     \
(__GNUC__ == (X) && __GNUC_MINOR__ == (Y) && __GNUC_PATCHLEVEL__ == (Z))
#endif
#endif

#ifndef GNUC_VERSION_ATLEAST
#ifndef __GNUC__
#define GNUC_VERSION_ATLEAST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,(X),(Y),(Z))
#endif
#endif

#ifndef GNUC_VERSION_ATMOST
#ifndef __GNUC__
#define GNUC_VERSION_ATMOST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,(X),(Y),(Z))
#endif
#endif

#ifndef INTEL_CC_VERSION
#ifndef __INTEL_COMPILER
#define INTEL_CC_VERSION(X, Y, Z) 0
#else
#define INTEL_CC_VERSION(X, Y, Z) \
    (__INTEL_COMPILER == (X) && __INTEL_COMPILER_UPDATE == (Y) && \
        __INTEL_COMPILER_BUILD_DATE == (Z))
#endif
#endif

#ifndef INTEL_CC_VERSION_ATLEAST
#ifndef __INTEL_COMPILER
#define INTEL_CC_VERSION_ATLEAST(X, Y, Z) 0
#else
#define INTEL_CC_VERSION_ATLEAST(X, Y, Z) \
    LEXGE3(__INTEL_COMPILER,    \
            __INTEL_COMPILER_UPDATE,    \
            __INTEL_COMPILER_BUILD_DATE,        \
            (X),(Y),(Z))
#endif
#endif

#ifndef INTEL_CC_VERSION_ATMOST
#ifndef __INTEL_COMPILER
#define INTEL_CC_VERSION_ATMOST(X, Y, Z) 0
#else
#define INTEL_CC_VERSION_ATMOST(X, Y, Z) \
    LEXGE3(__INTEL_COMPILER,    \
            __INTEL_COMPILER_UPDATE,    \
            __INTEL_COMPILER_BUILD_DATE,        \
            (X),(Y),(Z))
#endif
#endif


#ifndef GMP_VERSION_ATLEAST
#ifndef __GNU_MP_VERSION
#define GMP_VERSION_ATLEAST(X,Y,Z) 0
#else
#define GMP_VERSION_ATLEAST(X,Y,Z)     \
LEXGE3(__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL,(X),(Y),(Z))
#endif
#endif

#ifndef GMP_VERSION_ATMOST
#ifndef __GNU_MP_VERSION
#define GMP_VERSION_ATMOST(X,Y,Z) 0
#else
#define GMP_VERSION_ATMOST(X,Y,Z)     \
LEXLE3(__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL,(X),(Y),(Z))
#endif
#endif

/* Intel icc and Clang try to imitate all predefined macros present in
   FSF's GCC which can make it tricky to detect whether the compiler is
   genuine GCC. Thus we define GENUINE_GNUC by explicitly testing that
   the compiler is neither icc nor Clang. */
#ifndef GENUINE_GNUC
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
#define GENUINE_GNUC 1
#endif
#endif

#ifndef MPI_VERSION_ATLEAST
#define MPI_VERSION_ATLEAST(X,Y) LEXGE2(MPI_VERSION,MPI_SUBVERSION,(X),(Y))
#endif  /* MPI_VERSION_ATLEAST */

#ifndef MPI_VERSION_ATMOST
#define MPI_VERSION_ATMOST(X,Y) LEXLE2(MPI_VERSION,MPI_SUBVERSION,(X),(Y))
#endif  /* MPI_VERSION_ATMOST */

#ifndef MPI_VERSION_IS
#define MPI_VERSION_IS(X,Y) (((X) == MPI_VERSION) && ((Y) == MPI_SUBVERSION))
#endif  /* MPI_VERSION_IS */

#ifndef OMPI_VERSION_ATLEAST
#ifdef OPEN_MPI
#define OMPI_VERSION_ATLEAST(X,Y,Z)     \
    (LEXGE3(OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION,(X),(Y),(Z)))
#else
#define OMPI_VERSION_ATLEAST(X,Y,Z) 0
#endif
#endif  /* OMPI_VERSION_ATLEAST */

#ifndef OMPI_VERSION_ATMOST
#ifdef OPEN_MPI
#define OMPI_VERSION_ATMOST(X,Y,Z)     \
    (LEXLE3(OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION,(X),(Y),(Z)))
#else
#define OMPI_VERSION_ATMOST(X,Y,Z) 0
#endif
#endif  /* OMPI_VERSION_ATMOST */

#ifndef OMPI_VERSION_IS
#ifdef OPEN_MPI
#define OMPI_VERSION_IS(X,Y,Z)          \
    (((X) == OMPI_MAJOR_VERSION) &&     \
    ((Y) == OMPI_MINOR_VERSION) &&      \
    ((Z) == OMPI_RELEASE_VERSION))
#else
#define OMPI_VERSION_IS(X,Y,Z)          0
#endif
#endif  /* OMPI_VERSION_IS */


#ifndef MAYBE_UNUSED
#if GNUC_VERSION_ATLEAST(3,4,0)
/* according to
 * http://gcc.gnu.org/onlinedocs/gcc-3.1.1/gcc/Variable-Attributes.html#Variable%20Attributes
 * the 'unused' attribute already existed in 3.1.1 ; however the rules
 * for its usage remained quirky until 3.4.0, so we prefer to stick to
 * the more modern way of using the unused attribute, and recommend
 * setting the -Wno-unused flag for pre-3.4 versions of gcc
 */
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif
#endif
#if !(defined(GENUINE_GNUC) && GNUC_VERSION_ATMOST(11,9,9))
/* https://stackoverflow.com/questions/50646334/maybe-unused-on-member-variable-gcc-warns-incorrectly-that-attribute-is
 */
#define MAYBE_UNUSED_PRIVATE_DATA_MEMBER MAYBE_UNUSED
#else
#define MAYBE_UNUSED_PRIVATE_DATA_MEMBER
#endif

#if __STDC_VERSION__ >= 201112L
#define STATIC_ASSERT(COND,MSG) _Static_assert(COND, #MSG)
#else
/* Does not work in structs */
#define STATIC_ASSERT(COND,MSG)                                         \
    typedef char static_assertion_##MSG[(COND)?1:-1] MAYBE_UNUSED
#endif

#ifndef ATTRIBUTE_WARN_UNUSED_RESULT
/* https://gcc.gnu.org/onlinedocs/gcc-3.3.6/gcc/Function-Attributes.html#Function-Attributes
 * https://gcc.gnu.org/onlinedocs/gcc-3.4.6/gcc/Function-Attributes.html#Function-Attributes
 */
#if GNUC_VERSION_ATLEAST(3,4,0)
#define ATTRIBUTE_WARN_UNUSED_RESULT __attribute__ ((warn_unused_result))
#elif defined(__clang__)
#if __has_attribute(warn_unused_result)
#define ATTRIBUTE_WARN_UNUSED_RESULT __attribute__((warn_unused_result))
#else
#define ATTRIBUTE_WARN_UNUSED_RESULT
#endif
#else
#define ATTRIBUTE_WARN_UNUSED_RESULT
#endif
#endif

#ifndef ATTRIBUTE_NODISCARD
#if defined(__cplusplus)
#define ATTRIBUTE_NODISCARD [[nodiscard]]
#else
#define ATTRIBUTE_NODISCARD
#endif
#endif

#ifndef ATTRIBUTE_DEPRECATED
#if GNUC_VERSION_ATLEAST(3,1,1)
#define ATTRIBUTE_DEPRECATED __attribute__ ((deprecated))
#elif defined(__clang__)
#if __has_attribute(deprecated)
#define ATTRIBUTE_DEPRECATED __attribute__((deprecated))
#else
#define ATTRIBUTE_DEPRECATED
#endif
#else
#define ATTRIBUTE_DEPRECATED
#endif
#endif

#ifndef ATTRIBUTE_ARTIFICIAL
#if GNUC_VERSION_ATLEAST(4,3,0)
#define ATTRIBUTE_ARTIFICIAL __attribute__ ((__artificial__))
#elif defined(__clang__)
#if __has_attribute(artificial)
#define ATTRIBUTE_ARTIFICIAL __attribute__((artificial))
#else
#define ATTRIBUTE_ARTIFICIAL
#endif
#else
#define ATTRIBUTE_ARTIFICIAL
#endif
#endif

#ifndef ATTRIBUTE_NONNULL
#if GNUC_VERSION_ATLEAST(3,3,6)
#define ATTRIBUTE_NONNULL(which) __attribute__ ((__nonnull__ which))
#elif defined(__clang__)
#if __has_attribute(nonnull)
#define ATTRIBUTE_NONNULL(which) __attribute__((nonnull which))
#else
#define ATTRIBUTE_NONNULL(which)
#endif
#else
#define ATTRIBUTE_NONNULL(which)
#endif
#endif

#ifndef ATTRIBUTE_ALWAYS_INLINE
#if GNUC_VERSION_ATLEAST(3,1,0)
#define ATTRIBUTE_ALWAYS_INLINE __attribute__ ((__always_inline__))
#elif defined(__clang__)
#if __has_attribute(always_inline)
#define ATTRIBUTE_ALWAYS_INLINE __attribute__((always_inline))
#else
#define ATTRIBUTE_ALWAYS_INLINE
#endif
#else
#define ATTRIBUTE_ALWAYS_INLINE
#endif
#endif

#if defined(__GNUC__)

#ifndef NO_INLINE
#define NO_INLINE __attribute__ ((noinline))
#endif
#ifndef ATTR_EXPECT
#define ATTR_EXPECT(x,val)	__builtin_expect(x,val)
#endif
#ifndef ATTR_ALIGNED
#define ATTR_ALIGNED(x) __attribute__((aligned(x)))
#endif
#ifndef  HAVE_MINGW
#ifndef ATTR_PRINTF
#define ATTR_PRINTF(a,b) __attribute__((format(printf,a,b)))
#endif
#ifndef CONSTANT_P
#define CONSTANT_P(x) __builtin_constant_p(x)
#endif
#else
/* mingw's gcc is apparently unaware that the c99 format strings _may_ be
 * recognized by the win32 printf, for who asks nicely... */
#define ATTR_PRINTF(a,b) /**/
#endif  /* HAVE_MINGW */
/* Note that ATTRIBUTE is sort of a catch-all, but its use should be
 * discouraged, or at least limited to attributes which have been in gcc
 * versions for a very long time. For a newly introduced gcc version, it
 * is crucial to *NOT* use ATTRIBUTE() here, and instead define a macro
 * which is set depending on the GCC version (see further down in this
 * file).
 */
#ifndef ATTRIBUTE
#define ATTRIBUTE(x) __attribute__ (x)
#endif
#else
#ifndef NO_INLINE
#define NO_INLINE
#endif
#ifndef ATTR_EXPECT
#define	ATTR_EXPECT(x,val)	(x)
#endif
#ifndef ATTR_ALIGNED
#define ATTR_ALIGNED(x)
#endif
#ifndef ATTR_PRINTF
#define ATTR_PRINTF(a,b) /**/
#endif
#ifndef ATTRIBUTE
#define ATTRIBUTE(x)
#ifndef CONSTANT_P
#define CONSTANT_P(x) 0
#endif
#endif
#endif /* if defined(__GNUC__) */

/* These warnings are a nuisance, really. Not only do we have now to
 * add no_break() statements when we want switch cases fall through one
 * another, but on top of that, coverity wants the corresponding lines to
 * be preceded by a coverity[unterminated_case] comment...
 */
#if GNUC_VERSION_ATLEAST(7,0,0) && !defined(__ICC)
#define no_break() __attribute__ ((fallthrough))
#else
#define no_break()
#endif

/* as of version __INTEL_COMPILER==_ICC==1800, attribute assume_aligned
 * is not supported, even though the underlying gcc is 6.3...
 *
 * I'm flagging this as unsupported overall by icc. Maybe if someone
 * cares to check at some later point, we could have a finer grain test
 * case.
 */
#if GNUC_VERSION_ATLEAST(4,9,0) && !defined(__ICC)
#define ATTR_ASSUME_ALIGNED(x) __attribute__((assume_aligned(x)))
#else
#define ATTR_ASSUME_ALIGNED(x)
#endif

/* On 64 bit gcc (Ubuntu/Linaro 4.6.3-1ubuntu5) with -O3, the inline
   asm in ularith_div_2ul_ul_ul_r() is wrongly optimized (Alex Kruppa
   finds that gcc optimizes away the whole asm block and simply
   leaves a constant). */
/* On gcc 4.8.0, ..., 4.8.2, asm() blocks can be optimized away erroneously
   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58805 */
/* In both cases we force gcc not to omit any asm() blocks by declaring them
   volatile. */
#if defined(VOLATILE_IF_GCC_UBUNTU_BUG) || defined(VOLATILE_IF_GCC_58805_BUG)
#define __VOLATILE __volatile__
#else
#define __VOLATILE
#endif

#ifndef	LIKELY
#define LIKELY(x)	ATTR_EXPECT(x,1)
#endif
#ifndef	UNLIKELY
#define UNLIKELY(x)	ATTR_EXPECT(x,0)
#endif

#ifndef CADO_CONCATENATE
#define CADO_CONCATENATE_SUB(a,b) a##b // actually concatenate
#define CADO_CONCATENATE(a,b) CADO_CONCATENATE_SUB(a,b) // force expand
#define CADO_CONCATENATE3_SUB(a,b,c) a##b##c
#define CADO_CONCATENATE3(a,b,c) CADO_CONCATENATE3_SUB(a,b,c)
#define CADO_CONCATENATE4_SUB(a,b,c,d) a##b##c##d
#define CADO_CONCATENATE4(a,b,c,d) CADO_CONCATENATE4_SUB(a,b,c,d)
#endif

#ifndef CADO_STRINGIZE
#define CADO_STRINGIZE_(x) #x
#define CADO_STRINGIZE(x) CADO_STRINGIZE_(x)
#endif

#endif	/* CADO_MACROS_H */
