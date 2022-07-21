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

if HAVE_HIDDEN_VISIBILITY_ATTRIBUTE
VH_FFT_FLAG=-fvisibility=hidden -DGF2X_FFT_EXPORTED="__attribute__((visibility(\"default\")))"
endif

FFT_SOURCES=gf2x-fft.c
FFT_LIBS=
if ENABLE_FFT_INTERFACE
FFT_LIBS+=libgf2x-fft.la
FFT_SOURCES+=gf2x-fake-fft.c
FFT_SOURCES+=gf2x-ternary-fft.c
FFT_SOURCES+=gf2x-cantor-fft.c
gf2x_includedir=$(includedir)/gf2x
gf2x_include_HEADERS=	\
	gf2x-fft.h		\
	gf2x-fake-fft.h		\
	gf2x-cantor-fft.h	\
	gf2x-ternary-fft.h
lib_LTLIBRARIES=$(FFT_LIBS)

# The standard gf2x-fft library (installed by --enable-fft-interface) has
#
#  gf2x_ternary_fft
#  gf2x_fake_fft
#  gf2x_cantor_fft with -DCANTOR_BASE_FIELD_SIZE=128
#
libgf2x_fft_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DCANTOR_BASE_FIELD_SIZE=128
libgf2x_fft_la_SOURCES=$(FFT_SOURCES)
libgf2x_fft_la_LIBADD= $(top_builddir)/libgf2x.la
libgf2x_fft_la_LDFLAGS=-version-info @gf2x_fft_lib_version@ -no-undefined
endif

TESTS =
check_PROGRAMS=
check_LTLIBRARIES=

EXTRA_DIST=	\
	gf2x-fft-impl-utils.h		\
	gf2x-cantor-field-impl.h

# ternary fft is included by default in gf2x.
check_LTLIBRARIES += libternary.la
libternary_la_SOURCES=gf2x-ternary-fft.c
libternary_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) -DENGINE_TERNARY
libternary_la_LIBADD= $(top_builddir)/libgf2x.la
 
# The rest is optional

if ENABLE_FFT_INTERFACE

# With --enable-fft-interface, 
# we stress-test the libcantor implementation with the following
# extra features (compile-time flags) -- those lead to variants of
# libgf2x-fft.

CCOMMON=-DENGINE_CANTOR -DESCAPE_CANTOR_SYMBOLS
C64=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=64
C64nt=$(C64) -DWITHOUT_CANTOR_TRUNCATION
C64gm=$(C64) -DCANTOR_GM
C64gmt=$(C64gm) -DCANTOR_GM_TRUNCATE
C128=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=128
C128nt=$(C128) -DWITHOUT_CANTOR_TRUNCATION
C128gm=$(C128) -DCANTOR_GM
C128gmt=$(C128gm) -DCANTOR_GM_TRUNCATE

# while read x y ; do sed -e s/XXX/$x/g -e s/YYY/$y/g <<'EOF'
# check_LTLIBRARIES += libXXX.la
# libXXX_la_SOURCES=gf2x-cantor-fft.c
# libXXX_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(YYY)
# check_LTLIBRARIES += libgf2x-fft-XXX.la
# libgf2x_fft_XXX_la_SOURCES=$(FFT_SOURCES)
# libgf2x_fft_XXX_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(YYY)
# libgf2x_fft_XXX_la_LIBADD= $(top_builddir)/libgf2x.la
# check_PROGRAMS += gf2x_fft_check_XXX
# TESTS += gf2x_fft_check_XXX
# gf2x_fft_check_XXX_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_XXX_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(YYY)
# gf2x_fft_check_XXX_LDADD=libgf2x-fft-XXX.la ../libgf2x.la
# 
# EOF
# done
# cantor64 C64
# cantor64nt C64nt
# cantor64gm C64gm
# cantor64gmt C64gmt
# cantor128 C128
# cantor128nt C128nt
# cantor128gm C128gm
# cantor128gmt C128gmt

# Gao-Mateer code is apparently buggy. No time to investigate.
#
check_LTLIBRARIES += libcantor64.la
libcantor64_la_SOURCES=gf2x-cantor-fft.c
libcantor64_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C64)
check_LTLIBRARIES += libgf2x-fft-cantor64.la
libgf2x_fft_cantor64_la_SOURCES=$(FFT_SOURCES)
libgf2x_fft_cantor64_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C64)
libgf2x_fft_cantor64_la_LIBADD= $(top_builddir)/libgf2x.la
check_PROGRAMS += gf2x_fft_check_cantor64
TESTS += gf2x_fft_check_cantor64
gf2x_fft_check_cantor64_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor64_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C64)
gf2x_fft_check_cantor64_LDADD=libgf2x-fft-cantor64.la ../libgf2x.la

check_LTLIBRARIES += libcantor64nt.la
libcantor64nt_la_SOURCES=gf2x-cantor-fft.c
libcantor64nt_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C64nt)
check_LTLIBRARIES += libgf2x-fft-cantor64nt.la
libgf2x_fft_cantor64nt_la_SOURCES=$(FFT_SOURCES)
libgf2x_fft_cantor64nt_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C64nt)
libgf2x_fft_cantor64nt_la_LIBADD= $(top_builddir)/libgf2x.la
check_PROGRAMS += gf2x_fft_check_cantor64nt
TESTS += gf2x_fft_check_cantor64nt
gf2x_fft_check_cantor64nt_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor64nt_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C64nt)
gf2x_fft_check_cantor64nt_LDADD=libgf2x-fft-cantor64nt.la ../libgf2x.la

# check_LTLIBRARIES += libcantor64gm.la
# libcantor64gm_la_SOURCES=gf2x-cantor-fft.c
# libcantor64gm_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C64gm)
# check_LTLIBRARIES += libgf2x-fft-cantor64gm.la
# libgf2x_fft_cantor64gm_la_SOURCES=$(FFT_SOURCES)
# libgf2x_fft_cantor64gm_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C64gm)
# libgf2x_fft_cantor64gm_la_LIBADD= $(top_builddir)/libgf2x.la
# check_PROGRAMS += gf2x_fft_check_cantor64gm
# TESTS += gf2x_fft_check_cantor64gm
# gf2x_fft_check_cantor64gm_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_cantor64gm_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C64gm)
# gf2x_fft_check_cantor64gm_LDADD=libgf2x-fft-cantor64gm.la ../libgf2x.la
# 
# check_LTLIBRARIES += libcantor64gmt.la
# libcantor64gmt_la_SOURCES=gf2x-cantor-fft.c
# libcantor64gmt_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C64gmt)
# check_LTLIBRARIES += libgf2x-fft-cantor64gmt.la
# libgf2x_fft_cantor64gmt_la_SOURCES=$(FFT_SOURCES)
# libgf2x_fft_cantor64gmt_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C64gmt)
# libgf2x_fft_cantor64gmt_la_LIBADD= $(top_builddir)/libgf2x.la
# check_PROGRAMS += gf2x_fft_check_cantor64gmt
# TESTS += gf2x_fft_check_cantor64gmt
# gf2x_fft_check_cantor64gmt_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_cantor64gmt_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C64gmt)
# gf2x_fft_check_cantor64gmt_LDADD=libgf2x-fft-cantor64gmt.la ../libgf2x.la
# 
check_LTLIBRARIES += libcantor128.la
libcantor128_la_SOURCES=gf2x-cantor-fft.c
libcantor128_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C128)
check_LTLIBRARIES += libgf2x-fft-cantor128.la
libgf2x_fft_cantor128_la_SOURCES=$(FFT_SOURCES)
libgf2x_fft_cantor128_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C128)
libgf2x_fft_cantor128_la_LIBADD= $(top_builddir)/libgf2x.la
check_PROGRAMS += gf2x_fft_check_cantor128
TESTS += gf2x_fft_check_cantor128
gf2x_fft_check_cantor128_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor128_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C128)
gf2x_fft_check_cantor128_LDADD=libgf2x-fft-cantor128.la ../libgf2x.la

check_LTLIBRARIES += libcantor128nt.la
libcantor128nt_la_SOURCES=gf2x-cantor-fft.c
libcantor128nt_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C128nt)
check_LTLIBRARIES += libgf2x-fft-cantor128nt.la
libgf2x_fft_cantor128nt_la_SOURCES=$(FFT_SOURCES)
libgf2x_fft_cantor128nt_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C128nt)
libgf2x_fft_cantor128nt_la_LIBADD= $(top_builddir)/libgf2x.la
check_PROGRAMS += gf2x_fft_check_cantor128nt
TESTS += gf2x_fft_check_cantor128nt
gf2x_fft_check_cantor128nt_SOURCES=gf2x_fft_check.c
gf2x_fft_check_cantor128nt_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C128nt)
gf2x_fft_check_cantor128nt_LDADD=libgf2x-fft-cantor128nt.la ../libgf2x.la
# 
# check_LTLIBRARIES += libcantor128gm.la
# libcantor128gm_la_SOURCES=gf2x-cantor-fft.c
# libcantor128gm_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C128gm)
# check_LTLIBRARIES += libgf2x-fft-cantor128gm.la
# libgf2x_fft_cantor128gm_la_SOURCES=$(FFT_SOURCES)
# libgf2x_fft_cantor128gm_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C128gm)
# libgf2x_fft_cantor128gm_la_LIBADD= $(top_builddir)/libgf2x.la
# check_PROGRAMS += gf2x_fft_check_cantor128gm
# TESTS += gf2x_fft_check_cantor128gm
# gf2x_fft_check_cantor128gm_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_cantor128gm_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C128gm)
# gf2x_fft_check_cantor128gm_LDADD=libgf2x-fft-cantor128gm.la ../libgf2x.la
# 
# check_LTLIBRARIES += libcantor128gmt.la
# libcantor128gmt_la_SOURCES=gf2x-cantor-fft.c
# libcantor128gmt_la_CFLAGS= $(AM_CFLAGS) $(VH_FFT_FLAG) $(C128gmt)
# check_LTLIBRARIES += libgf2x-fft-cantor128gmt.la
# libgf2x_fft_cantor128gmt_la_SOURCES=$(FFT_SOURCES)
# libgf2x_fft_cantor128gmt_la_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) $(C128gmt)
# libgf2x_fft_cantor128gmt_la_LIBADD= $(top_builddir)/libgf2x.la
# check_PROGRAMS += gf2x_fft_check_cantor128gmt
# TESTS += gf2x_fft_check_cantor128gmt
# gf2x_fft_check_cantor128gmt_SOURCES=gf2x_fft_check.c
# gf2x_fft_check_cantor128gmt_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE $(C128gmt)
# gf2x_fft_check_cantor128gmt_LDADD=libgf2x-fft-cantor128gmt.la ../libgf2x.la
# 



# This last one is a duplicate of one of the other one (presumably
# cantor128), and it exerts libgf2x-fft as it is shipped.
check_PROGRAMS += gf2x_fft_check
gf2x_fft_check_SOURCES=gf2x_fft_check.c
gf2x_fft_check_CFLAGS=$(AM_CFLAGS) $(VH_FFT_FLAG) -DGF2X_SOURCE_TREE
gf2x_fft_check_LDADD=libgf2x-fft.la ../libgf2x.la
TESTS += gf2x_fft_check
endif
