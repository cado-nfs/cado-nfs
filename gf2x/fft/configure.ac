#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2013, 2014, 2015
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of either:
#   - If the archive contains a file named toom-gpl.c (not a trivial
#     placeholder), the GNU General Public License as published by the
#     Free Software Foundation; either version 3 of the License, or (at
#     your option) any later version.
#   - If the archive contains a file named toom-gpl.c which is a trivial
#     placeholder, the GNU Lesser General Public License as published by
#     the Free Software Foundation; either version 2.1 of the License, or
#     (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
#  
#  You should have received a copy of the GNU General Public License as
#  well as the GNU Lesser General Public License along with this program;
#  see the files COPYING and COPYING.LIB.  If not, write to the Free
#  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301, USA.

AC_INIT([gf2x-fft], [0.8])
: ${CFLAGS=-g -O4}
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CANONICAL_TARGET

AC_SUBST([gf2x_fft_lib_version], [0:8:0])

AM_INIT_AUTOMAKE

AC_LANG([C])

LT_INIT
LT_LANG([C])

AC_ARG_ENABLE([sse2],
              [AS_HELP_STRING([--enable-sse2],
                              [Turn on sse-2 code (default is yes)])])

AC_ARG_ENABLE([pclmul],
              [AS_HELP_STRING([--enable-pclmul],
                              [Turn on pclmul code (default is yes)])])

AC_PROG_CC
AC_COMPILE_WARNINGS
AC_PROG_CC_C99

GF2X_PROG_CC_FOR_BUILD
GF2X_PROG_EXEEXT_FOR_BUILD

# This macro is used for tuning
AM_PROG_CC_C_O


CHECK_SSE2_SUPPORT()
CHECK_PCLMUL_SUPPORT()

VERIFY_WORDSIZE()

AC_CONFIG_HEADERS([gf2xfft-config.h])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# vim: set sw=1:
