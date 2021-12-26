// Copyright 2005 Caleb Epstein
// Copyright 2006 John Maddock
// Distributed under the Boost Software License, Version 1.0. (See accompany-
// ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

/* copied and modified from:
 *
 * http://boost.cvs.sourceforge.net/boost/boost/boost/detail/endian.hpp
 */

/*
 * Copyright (c) 1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

#ifndef CADO_ENDIAN_H_
#define CADO_ENDIAN_H_

// pragma no prototypes

#ifdef HAVE_GLIBC
// GNU libc offers the helpful header <endian.h> which defines
// __BYTE_ORDER
# include <endian.h>
# if (__BYTE_ORDER == __LITTLE_ENDIAN)
#  define CADO_LITTLE_ENDIAN
# elif (__BYTE_ORDER == __BIG_ENDIAN)
#  define CADO_BIG_ENDIAN
# elif (__BYTE_ORDER == __PDP_ENDIAN)
#  define CADO_PDP_ENDIAN
# else
#  error Unknown machine endianness detected.
# endif
# define CADO_BYTE_ORDER __BYTE_ORDER
/* There is no serious reason to think that _BIG_ENDIAN or _LITTLE_ENDIAN
 * being defined actually means that the machine is big (resp, little)
 * endian. Systems may like to unconditionally define these as constants
 * to correspond to some endianness, and define whichever user-exposed
 * constant they like (e.g. BYTE_ORDER) to either. After all, this is
 * more or less how it works for glibc above.
#elif defined(_BIG_ENDIAN)
# define CADO_BIG_ENDIAN
# define CADO_BYTE_ORDER 4321
#elif defined(_LITTLE_ENDIAN)
# define CADO_LITTLE_ENDIAN
# define CADO_BYTE_ORDER 1234
 */
#elif defined(__sparc) || defined(__sparc__) \
   || defined(_POWER) || defined(__powerpc__) \
   || defined(__ppc__) || defined(__hpux) \
   || defined(_MIPSEB) || defined(_POWER) \
   || defined(__AARCH64EB__) \
   || defined(__s390__)
# define CADO_BIG_ENDIAN
# define CADO_BYTE_ORDER 4321
#elif defined(__i386__) || defined(__alpha__) \
   || defined(__ia64) || defined(__ia64__) \
   || defined(_M_IX86) || defined(_M_IA64) \
   || defined(_M_ALPHA) || defined(__amd64) \
   || defined(__amd64__) || defined(_M_AMD64) \
   || defined(__x86_64) || defined(__x86_64__) \
   || defined(__AARCH64EL__) \
   || defined(_M_X64) \
   || defined(__LITTLE_ENDIAN__)
# define CADO_LITTLE_ENDIAN
# define CADO_BYTE_ORDER 1234
#else
# error The file cado-endian.h needs to be set up for your CPU type.
#endif
#endif
