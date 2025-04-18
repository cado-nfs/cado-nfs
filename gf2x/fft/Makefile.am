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

ACLOCAL_AMFLAGS=-I config --install
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)

TESTS =
check_PROGRAMS=
check_LTLIBRARIES=

EXTRA_DIST=	\
	gf2x-fft-impl-utils.h

# ternary fft is included by default in gf2x.
check_LTLIBRARIES += libternary.la
libternary_la_SOURCES=gf2x-ternary-fft.c
libternary_la_CFLAGS= $(AM_CFLAGS) -DENGINE_TERNARY
libternary_la_LIBADD= $(top_builddir)/libgf2x.la
 
# The rest is optional

# we stress-test the libcantor implementation with the following
# extra compile-time flags

CCOMMON=-DENGINE_CANTOR
C64=$(CCOMMON) -DGF2X_CANTOR_BASE_FIELD_SIZE=64
C64nt=$(C64) -DWITHOUT_CANTOR_TRUNCATION
C64gm=$(C64) -DCANTOR_GM
C64gmt=$(C64gm) -DCANTOR_GM_TRUNCATE
C128=$(CCOMMON) -DGF2X_CANTOR_BASE_FIELD_SIZE=128
C128nt=$(C128) -DWITHOUT_CANTOR_TRUNCATION
C128gm=$(C128) -DCANTOR_GM
C128gmt=$(C128gm) -DCANTOR_GM_TRUNCATE

# while read x y ; do sed -e s/XXX/$x/g -e s/YYY/$y/g <<'EOF'
# check_LTLIBRARIES += libXXX.la
# libXXX_la_SOURCES=gf2x-cantor-fft.c
# libXXX_la_CFLAGS= $(AM_CFLAGS) $(YYY)
# check_PROGRAMS += gf2x_fft_check_XXX
# TESTS += gf2x_fft_check_XXX
# gf2x_fft_check_XXX_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_XXX_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE $(YYY)
# gf2x_fft_check_XXX_LDADD=../libgf2x.la
# 
# EOF
# done
# cantor64 C64
# cantor64nt C64nt
# cantor128 C128
# cantor128nt C128nt

# Gao-Mateer code is apparently buggy. No time to investigate.
# cantor64gm C64gm
# cantor64gmt C64gmt
# cantor128gm C128gm
# cantor128gmt C128gmt

check_LTLIBRARIES += libcantor64.la
libcantor64_la_SOURCES=gf2x-cantor-fft.c
libcantor64_la_CFLAGS= $(AM_CFLAGS) $(C64)
check_PROGRAMS += gf2x_fft_check_cantor64
TESTS += gf2x_fft_check_cantor64
gf2x_fft_check_cantor64_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor64_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE $(C64)
gf2x_fft_check_cantor64_LDADD=../libgf2x.la

check_LTLIBRARIES += libcantor64nt.la
libcantor64nt_la_SOURCES=gf2x-cantor-fft.c
libcantor64nt_la_CFLAGS= $(AM_CFLAGS) $(C64nt)
check_PROGRAMS += gf2x_fft_check_cantor64nt
TESTS += gf2x_fft_check_cantor64nt
gf2x_fft_check_cantor64nt_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor64nt_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE $(C64nt)
gf2x_fft_check_cantor64nt_LDADD=../libgf2x.la

check_LTLIBRARIES += libcantor128.la
libcantor128_la_SOURCES=gf2x-cantor-fft.c
libcantor128_la_CFLAGS= $(AM_CFLAGS) $(C128)
check_PROGRAMS += gf2x_fft_check_cantor128
TESTS += gf2x_fft_check_cantor128
gf2x_fft_check_cantor128_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor128_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE $(C128)
gf2x_fft_check_cantor128_LDADD=../libgf2x.la

check_LTLIBRARIES += libcantor128nt.la
libcantor128nt_la_SOURCES=gf2x-cantor-fft.c
libcantor128nt_la_CFLAGS= $(AM_CFLAGS) $(C128nt)
check_PROGRAMS += gf2x_fft_check_cantor128nt
TESTS += gf2x_fft_check_cantor128nt
gf2x_fft_check_cantor128nt_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor128nt_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE $(C128nt)
gf2x_fft_check_cantor128nt_LDADD=../libgf2x.la

# This last one is a duplicate of one of the other ones (presumably
# cantor128), and it exerts libgf2x-fft as it is shipped.
check_PROGRAMS += gf2x_fft_check
gf2x_fft_check_SOURCES=gf2x_fft_check.c
gf2x_fft_check_CFLAGS=$(AM_CFLAGS) -DGF2X_SOURCE_TREE
gf2x_fft_check_LDADD=../libgf2x.la
TESTS += gf2x_fft_check
