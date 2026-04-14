/* Portability header file for the CADO project
 
Copyright 2013 Pierrick Gaudry, Alexander Kruppa,
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

/* This header file defines macros (and perhaps static functions) to improve 
 * portability of the CADO code. They aim to provide a wrapper for some C99
 * and POSIX functionality for systems that lack those.
 */

#ifndef CADO_PORTABILITY_H
#define CADO_PORTABILITY_H

#ifndef CADO_CONFIG_H
#error cado_config.h must be included before portability.h
#endif

#include "macros.h"

// scan-headers: skip

#ifndef HAVE_STRDUP/*{{{*/
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline char * strdup(const char * const s)
{
    const size_t size = strlen(s) + 1;
    char * const r = (char *) malloc(size * sizeof(char));
    if (r != NULL)
        memcpy(r, s, size);
    return r;
}
#ifdef __cplusplus
}
#endif
#endif /* HAVE_STRDUP *//*}}}*/

#ifndef HAVE_STRNDUP/*{{{*/
/* strndup is posix-2008. providing a workalike is very easy */
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline char * strndup(const char * const a, const size_t n)
{
    const size_t l = strlen(a);
    const size_t size = (l < n ? l : n) + 1;
    char * const r = (char *) malloc(size * sizeof(char));
    if (r != NULL) {
        memcpy(r, a, size);
        r[size] = '\0';
    }
    return r;
}
#ifdef __cplusplus
}
#endif
#endif /* HAVE_STRNDUP *//*}}}*/

#ifndef HAVE_STRNLEN/*{{{*/
/* strnlen is posix-2008. providing a workalike is very easy */
#ifdef __cplusplus
extern "C" {
#endif
static inline size_t strnlen (const char *s, size_t maxlen)
{
  size_t ret = 0;
  for (; ret < maxlen && *s; s++)
    ret++;
  return ret;
}
#ifdef __cplusplus
}
#endif
#endif/*}}}*/

#ifndef HAVE_STRLCPY/*{{{*/
/* strlcpy is a bsd specificity, but we find it quite handy */

/*	$OpenBSD: strlcpy.c,v 1.11 2006/05/05 15:27:38 millert Exp $	*/

/*
 * Copyright (c) 1998 Todd C. Miller <Todd.Miller@courtesan.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <sys/types.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
/*
 * Copy src to string dst of size siz.  At most siz-1 characters
 * will be copied.  Always NUL terminates (unless siz == 0).
 * Returns strlen(src); if retval >= siz, truncation occurred.
 */
static inline size_t strlcpy(char *dst, const char *src, size_t size) ATTRIBUTE((__warn_unused_result__));
static inline size_t
strlcpy(char *dst, const char *src, size_t siz) ATTRIBUTE_WARN_UNUSED_RESULT;
static inline size_t
strlcpy(char *dst, const char *src, size_t siz)
{
	char *d = dst;
	const char *s = src;
	size_t n = siz;

	/* Copy as many bytes as will fit */
	if (n != 0) {
		while (--n != 0) {
			if ((*d++ = *s++) == '\0')
				break;
		}
	}

	/* Not enough room in dst, add NUL and traverse rest of src */
	if (n == 0) {
		if (siz != 0)
			*d = '\0';		/* NUL-terminate dst */
		while (*s++)
			;
	}

	return(s - src - 1);	/* count does not include NUL */
}
#ifdef __cplusplus
}
#endif
#endif/*}}}*/

#ifndef HAVE_STRLCAT/*{{{*/
/* strlcat is a bsd specificity, but we find it quite handy */

/*	$OpenBSD: strlcat.c,v 1.19 2019/01/25 00:19:25 millert Exp $	*/

/*
 * Copyright (c) 1998, 2015 Todd C. Miller <millert@openbsd.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <sys/types.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
/*
 * Appends src to string dst of size dsize (unlike strncat, dsize is the
 * full size of dst, not space left).  At most dsize-1 characters
 * will be copied.  Always NUL terminates (unless dsize <= strlen(dst)).
 * Returns strlen(src) + MIN(dsize, strlen(initial dst)).
 * If retval >= dsize, truncation occurred.
 */
static inline size_t
strlcat(char *dst, const char *src, size_t dsize) ATTRIBUTE_WARN_UNUSED_RESULT;
static inline size_t
strlcat(char *dst, const char *src, size_t dsize)
{
	const char *odst = dst;
	const char *osrc = src;
	size_t n = dsize;
	size_t dlen;

	/* Find the end of dst and adjust bytes left but don't go past end. */
	while (n-- != 0 && *dst != '\0')
		dst++;
	dlen = dst - odst;
	n = dsize - dlen;

	if (n-- == 0)
		return(dlen + strlen(src));
	while (*src != '\0') {
		if (n != 0) {
			*dst++ = *src;
			n--;
		}
		src++;
	}
	*dst = '\0';

	return(dlen + (src - osrc));	/* count does not include NUL */
}
#ifdef __cplusplus
}
#endif
#endif/*}}}*/

#ifndef HAVE_ASPRINTF/*{{{*/
/* Copied and improved from
 * http://mingw-users.1079350.n2.nabble.com/Query-regarding-offered-alternative-to-asprintf-td6329481.html
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline int 
vasprintf( char ** const sptr, const char *const fmt, va_list argv )
{
    int wanted = vsnprintf( *sptr = NULL, 0, fmt, argv );
    if (wanted<0)
        return -1;
    *sptr = (char *) malloc(1 + wanted);
    if (!*sptr)
        return -1;
    int rc = vsnprintf(*sptr, 1+wanted, fmt, argv );
    return rc;
}

static inline int 
asprintf( char ** const sptr, const char * const fmt, ... )
{
    int retval;
    va_list argv;
    va_start(argv, fmt);
    retval = vasprintf(sptr, fmt, argv);
    va_end(argv);
    return retval;
}
#ifdef __cplusplus
}
#endif
#endif  /* HAVE_ASPRINTF *//*}}}*/

#if defined(HAVE_SYSCONF)/*{{{*/
#include <unistd.h>
#endif/*}}}*/

static inline long pagesize ()/* {{{ */
{
#if defined(HAVE_SYSCONF)
  return sysconf (_SC_PAGESIZE);
#else
#error "Cannot determine page size"
#endif
}/*}}}*/


#endif /* ifndef CADO_PORTABILITY_H */
