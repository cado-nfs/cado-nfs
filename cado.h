/* Common header file for the CADO project

Copyright 2007-2015 Pierrick Gaudry, Alexander Kruppa, Emmanuel Thome, Paul Zimmermann

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

#ifndef CADO_H
#define CADO_H

// scan-headers: skip

/* The purpose of this header is to define some feature macros (and only
 * macros), which tweak the behaviour of include files. Our intent is to
 * define here in a unique place the required macros exposing the functions
 * we like to have.
 *
 * IT IS NECESSARY THAT THIS FILE APPEARS ONLY ON TOP OF THE COMPILATION
 * UNITS, AND BEFORE ANY OTHER HEADER. We promise to never include
 * another header file as a side-effect of this one (except
 * cado_config.h, which makes sense to include here as well).
 *
 * To make sure that the repository abides by the law above, see
 * scripts/check_compilation_units_policy.pl
 */

#if defined(__CYGWIN__) && defined(__STRICT_ANSI__)
/* The C library which comes with cygwin has no feature test macros. We do
   something ugly then.  */
#undef __STRICT_ANSI__
#endif

/* OpenBSD and FreeBSD expose *all* functions by default, and feature
 * macros are (apparently) used the other way around to restrict the
 * exposed interfaces.
 * FIXME: It's not entirely clear. Maybe it has been so in some version,
 * but it could also well be that I wholly misunderstood the thing.
 * OpenBSD 4.9 appears at least to grok _BSD_SOURCE as I expect (i.e.
 * _do_ expose BSD prototypes as an _addition_ to the rest).
 */
#ifdef __OpenBSD__
#define _BSD_SOURCE
#elif defined(__FreeBSD__)
/* XXX should check whether my former fear mentioned about turns out to
 * be true.
 */
#define _WITH_GETLINE /* See the COMPATIBILITY section of getline manpage
                         of FreeBSD. */
#else
#define _POSIX_C_SOURCE 200112L /* strtoumax */
/* POSIX: popen/pclose with -std=c99, -pedantic or -ansi (requires
 * _POSIX_C_SOURCE==2 ?) fileno */
#define _XOPEN_SOURCE   700     /* posix_memalign strndup */

/* This #define _FILE_OFFSET_BITS is useless, really: we never made the
 * effort to use off_t in cado-nfs code, and we always lived on the
 * assumption that size_t was certainly enough. Sadly untrue on 32-bit
 * machines: there, off_t + _FILE_OFFSET_BITS=64 are the only way
 * to deal with files above 2G. But at this point, it makes no sense to
 * amend the whole code to take this into account.
 */
#define _FILE_OFFSET_BITS 64
#define _BSD_SOURCE     /* M_LN2 gethostname strdup random */
#ifndef _ISOC99_SOURCE
#define _ISOC99_SOURCE 1  /* Sometimes there's link trickery which causes fscanf to be linked in *only* when this is defined */
#endif
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE 1 /* for glibc 2.20 and later. As per the man page
                             (feature_test_macros(7) ), the compilation
                             warning triggered by _BSD_SOURCE is silenced
                             if _DEFAULT_SOURCE is there too, so maybe we
                             can have both and not need a cmake-time check
                             for the libc version. */
#endif
#ifndef __cplusplus
#define _GNU_SOURCE         /* asprintf vasprintf */
#endif
#define _DARWIN_C_SOURCE    /* asprintf ; _ANSI_SOURCE must be undefined */
#define _NETBSD_SOURCE      /* asprintf vasprintf */
#endif

#ifdef __cplusplus
/* as per C99 standard:
 *      § 7.8.1 footnote 182            __STDC_FORMAT_MACROS
 *      § 7.18.2 footnote 217           __STDC_LIMIT_MACROS
 *      § 7.18.3 footnote 218           __STDC_LIMIT_MACROS
 *      § 7.18.4.1 footnote 220         __STDC_CONSTANT_MACROS
 *
 * Note though that C++11 items
 *      § 27.9.2.3                      __STDC_FORMAT_MACROS
 *      § 18.4.1.2                      the two other
 * explicitly say that despite what the C standard says, these macros
 * play no role in C++11. This might not be true for implementations
 * which do not conform completely to C++11.
 */
#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS    /* PRIu32 in lingen_mat_types.hpp */
#define __STDC_CONSTANT_MACROS
#endif

#if	defined(_AIX) && !defined(__cplusplus)
/* AIX seems to have weird rules. For example I don't understand exactly
 * how strndup is exposed. It seems that our choice of enabling feature
 * test macros above has the unfortunate effect of switching from the
 * default "expose-everything" mode for the libc headers to the more
 * restrictice "expose-only-what-has-been asked" mode. Putting
 * _ALL_SOURCE, which is obviously a catch-all, seems to maximize the
 * number of exported prototypes (it's clearly by lack of a better
 * solution). */
/* (for C++, _ALL_SOURCE is forcibly enabled by default) */
#define _ALL_SOURCE
#endif  /* _AIX */


/* This has the same effect as enforcing -Werror on the command line, so
 * as to trigger as many warnings as we can, and force ourselves to get
 * them fixed.
 * (useful to enable temporarily before releases, at least)
 *
 * note that clang and intel cc both define __GNUC__, so these flags
 * apply to them as well.
 */
// the problem is that for some given diagnostic, setting it as error
// once will apparently kill any opportunity to selectively disable it
// later on. So we now rely on the command line settings instead.
// #ifdef  __GNUC__
// #pragma GCC diagnostic error "-Wextra"
// #pragma GCC diagnostic error "-Wall"
// #endif

#if defined(__clang__)
#if __clang_major__ >= 17
/* See #30073. Sure, it's a bit of a pity, but I can't find a totally
 * satisfactory replacement for VLAs, which are used in somewhat critical
 * parts of the code.
 */
#pragma GCC diagnostic ignored "-Wvla-cxx-extension"
#endif
#endif

#include "cado_config.h"        // IWYU pragma: export

/* some stuff that we wish to enable only after we read our config-time
 * defines...
 */
#ifdef HAVE_MINGW
/* We define __USE_MINGW_ANSI_STDIO under MinGW to make MinGW provide */
/* C99-compliant printf() functions. I think that we should not have to
 * do this.
 */
#define    __USE_MINGW_ANSI_STDIO
#endif

/* It's really a bit of a pity, we don't even have our macros.h at this
 * point.
 */
#ifdef __cplusplus
#ifdef __GNUC__
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 9) || (__GNUC__ == 5 && __GNUC_MINOR__ <= 4)
#include <cstddef>	// problem with gcc 4.9.4 and gcc 5.4.0
#include <stdio.h>
#include <stdlib.h>
#endif
#endif
#endif

#endif  /* CADO_H */
