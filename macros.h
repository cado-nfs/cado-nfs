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

#ifndef CADO_MACROS_H_
#define CADO_MACROS_H_

/**********************************************************************/
/* Common asserting/debugging defines */
/* See README.macro_usage */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h> /* for __GNU_MP_VERSION */

#define ASSERT(x)	assert(x)

#define croak__(x,y) do {						\
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",			\
                (x),__func__,__FILE__,__LINE__,(y));			\
    } while (0)

#define ASSERT_ALWAYS(x)						\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)
#define FATAL_ERROR_CHECK(cond, msg)					\
    do {								\
      if (UNLIKELY((cond))) {                                           \
            croak__("Fatal error", msg);				\
          abort();                                                      \
        }								\
    } while (0)

#define DIE_ERRNO_DIAG(tst, func, arg) do {				\
    if (UNLIKELY(tst)) {				        	\
        fprintf(stderr, func "(%s): %s\n", arg, strerror(errno));       \
        exit(1);					        	\
    }							        	\
} while (0)

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
   (signed) integer type. For negative operands, the result is wrong. */
#define iceildiv(x,y) ((x) == 0 ? 0 : ((x)-1)/(y)+1)
/* siceildiv requires signed operands, or the compiler will throw warnings
   with -Wtype-limits */
#define siceildiv(x,y) ((x) == 0 ? 0 : ((x)<0) + ((y)<0) == 1 ? (x)/(y) : ((x)-1+2*((y)<0))/(y)+1)
#endif

#define LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && LEXGE2(Y,Z,B,C)))
#define LEXLE2(X,Y,A,B) LEXGE2(A,B,X,Y)
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3(A,B,C,X,Y,Z)

#ifndef GNUC_VERSION
#ifndef __GNUC__
#define GNUC_VERSION(X,Y,Z) 0
#else
#define GNUC_VERSION(X,Y,Z)     \
(__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z)
#endif
#endif

#ifndef GNUC_VERSION_ATLEAST
#ifndef __GNUC__
#define GNUC_VERSION_ATLEAST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif
#endif

#ifndef GNUC_VERSION_ATMOST
#ifndef __GNUC__
#define GNUC_VERSION_ATMOST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif
#endif

#ifndef GMP_VERSION_ATLEAST
#ifndef __GNU_MP_VERSION
#define GMP_VERSION_ATLEAST(X,Y,Z) 0
#else
#define GMP_VERSION_ATLEAST(X,Y,Z)     \
LEXGE3(__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL,X,Y,Z)
#endif
#endif

#ifndef GMP_VERSION_ATMOST
#ifndef __GNU_MP_VERSION
#define GMP_VERSION_ATMOST(X,Y,Z) 0
#else
#define GMP_VERSION_ATMOST(X,Y,Z)     \
LEXLE3(__GNU_MP_VERSION,__GNU_MP_VERSION_MINOR,__GNU_MP_VERSION_PATCHLEVEL,X,Y,Z)
#endif
#endif

#ifndef MPI_VERSION_ATLEAST
#define MPI_VERSION_ATLEAST(X,Y) LEXGE2(MPI_VERSION,MPI_SUBVERSION,X,Y)
#endif  /* MPI_VERSION_ATLEAST */

#ifndef MPI_VERSION_ATMOST
#define MPI_VERSION_ATMOST(X,Y) LEXLE2(MPI_VERSION,MPI_SUBVERSION,X,Y)
#endif  /* MPI_VERSION_ATMOST */

#ifndef MPI_VERSION_IS
#define MPI_VERSION_IS(X,Y) (((X) == MPI_VERSION) && ((Y) == MPI_SUBVERSION))
#endif  /* MPI_VERSION_IS */

#ifndef OMPI_VERSION_ATLEAST
#ifdef OPEN_MPI
#define OMPI_VERSION_ATLEAST(X,Y,Z)     \
    (LEXGE3(OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION,X,Y,Z))
#else
#define OMPI_VERSION_ATLEAST(X,Y,Z) 0
#endif
#endif  /* OMPI_VERSION_ATLEAST */

#ifndef OMPI_VERSION_ATMOST
#ifdef OPEN_MPI
#define OMPI_VERSION_ATMOST(X,Y,Z)     \
    (LEXLE3(OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION,X,Y,Z))
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

#ifndef ATTRIBUTE_DEPRECATED
#if GNUC_VERSION_ATLEAST(3,1,1)
#define ATTRIBUTE_DEPRECATED __attribute__ ((deprecated))
#else
#define ATTRIBUTE_DEPRECATED
#endif
#endif

#if defined(__GNUC__)

#ifndef NO_INLINE
#define NO_INLINE __attribute__ ((noinline))
#endif
#ifndef PACKED
#define PACKED __attribute__ ((packed))
#endif
#ifndef EXPECT
#define EXPECT(x,val)	__builtin_expect(x,val)
#endif
#ifndef  HAVE_MINGW
#ifndef ATTR_PRINTF
#define ATTR_PRINTF(a,b) __attribute__((format(printf,a,b)))
#endif
#else
/* mingw's gcc is apparently unaware that the c99 format strings _may_ be
 * recognized by the win32 printf, for who asks nicely... */
#define ATTR_PRINTF(a,b) /**/
#endif  /* HAVE_MINGW */
#ifndef ATTRIBUTE
#define ATTRIBUTE(x) __attribute__ (x)
#endif
#else
#ifndef NO_INLINE
#define NO_INLINE
#endif
#ifndef PACKED
#define PACKED
#endif
#ifndef EXPECT
#define	EXPECT(x,val)	(x)
#endif
#ifndef ATTR_PRINTF
#define ATTR_PRINTF(a,b) /**/
#endif
#ifndef ATTRIBUTE
#define ATTRIBUTE(x)
#endif
#endif /* if defined(__GNUC__) */

#if GNUC_VERSION_ATLEAST(7,0,0)
#define no_break() __attribute__ ((fallthrough))
#else
#define no_break()
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
#define LIKELY(x)	EXPECT(x,1)
#endif
#ifndef	UNLIKELY
#define UNLIKELY(x)	EXPECT(x,0)
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

#endif	/* CADO_MACROS_H_ */
