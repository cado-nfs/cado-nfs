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
/* unfortunately this fails miserably if x+y-1 overflows */
#define iceildiv(x,y)	(((x)+(y)-1)/(y))
#endif

#define LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && LEXGE2(Y,Z,B,C)))
#define LEXLE2(X,Y,A,B) LEXGE2(A,B,X,Y)
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3(A,B,C,X,Y,Z)

#ifndef GNUC_VERSION
#define GNUC_VERSION(X,Y,Z)     \
    defined(__GNUC__) &&        \
(__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z)
#endif

#ifndef GNUC_VERSION_ATLEAST
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif

#ifndef GNUC_VERSION_ATMOST
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
    defined(__GNUC__) &&        \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif

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
#ifndef ATTR_PRINTF
#define ATTR_PRINTF(a,b) __attribute__((format(printf,a,b)))
#endif
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
#endif

#ifndef	LIKELY
#define LIKELY(x)	EXPECT(x,1)
#endif
#ifndef	UNLIKELY
#define UNLIKELY(x)	EXPECT(x,0)
#endif

/* Portability macros */

/* MS VS and MinGW use the MS RTL (called MSVCRT for MinGW) which does not
   know the "%zu" format, they use "%Iu" instead. On MinGW, we use wrapper 
   functions that rewrite the %zu format accordingly, so the bulk of the
   code can continue to use C99 syntax.
   We do these renames only if stdio.h has been parsed before this file.
   Header files that need a certain include order are ugly, but that never
   stopped us and renaming printf() before parsing stdio.h would be "bad." 
   This way, when the renames are needed but don't happen, with any luck 
   gcc will complain about not understanding "%zu". */
#if defined(MINGW) && defined(_STDIO_H)
#define PRISIZ "Iu"
#define printf printf_subst_zu
#define fprintf fprintf_subst_zu
#define sprintf sprintf_subst_zu
#define snprintf snprintf_subst_zu
#endif

/* Handles portability cases which can be solved with a simple rename, 
   such as using the slower getc() instead of getc_unlocked() */
#ifndef HAVE_GETC_UNLOCKED
#define getc_unlocked getc
#endif
#ifndef HAVE_LRAND48
#define lrand48 rand
#endif

#endif	/* CADO_MACROS_H_ */
