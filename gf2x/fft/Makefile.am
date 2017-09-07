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
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir) -I$(srcdir)/mpfq

# We want to stress-test the libcantor_la implementation with the following
# extra features (compile-time flags)

CCOMMON=-DENGINE_CANTOR -DESCAPE_CANTOR_SYMBOLS
C64=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=64
C64nt=$(C64) -DWITHOUT_CANTOR_TRUNCATION
C64gm=$(C64) -DCANTOR_GM
C64gmt=$(C64gm) -DCANTOR_GM_TRUNCATE
C128=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=128
C128nt=$(C128) -DWITHOUT_CANTOR_TRUNCATION
C128gm=$(C128) -DCANTOR_GM
C128gmt=$(C128gm) -DCANTOR_GM_TRUNCATE

check_LTLIBRARIES=

EXTRA_DIST=	\
	mpfq/README.mpfq			\
	mpfq/i386/mpfq_2_128.c			\
	mpfq/i386/mpfq_2_128.h			\
	mpfq/i386/mpfq_2_64.c			\
	mpfq/i386/mpfq_2_64.h			\
	mpfq/mpfq.h				\
	mpfq/mpfq_gf2n_common.h			\
	mpfq/mpfq_name_K.h			\
	mpfq/x86_64/mpfq_2_128.c		\
	mpfq/x86_64/mpfq_2_128.h		\
	mpfq/x86_64/mpfq_2_64.c			\
	mpfq/x86_64/mpfq_2_64.h			\
	fft_adapter.hpp

# ternary fft is included by default in gf2x.
check_LTLIBRARIES += libternary.la
libternary_la_SOURCES=gf2x-ternary-fft.c
libternary_la_CFLAGS= $(AM_CFLAGS) -DENGINE_TERNARY
libternary_la_LIBADD= $(top_builddir)/libgf2x.la
 
# The rest is optional

if ENABLE_FFT_INTERFACE

check_LTLIBRARIES += libcantor64.la
libcantor64_la_SOURCES=gf2x-cantor-fft.c
libcantor64_la_CFLAGS= $(AM_CFLAGS) $(C64)

check_LTLIBRARIES += libcantor64nt.la
libcantor64nt_la_SOURCES=gf2x-cantor-fft.c
libcantor64nt_la_CFLAGS= $(AM_CFLAGS) $(C64nt)

check_LTLIBRARIES += libcantor64gm.la
libcantor64gm_la_SOURCES=gf2x-cantor-fft.c
libcantor64gm_la_CFLAGS= $(AM_CFLAGS) $(C64gm)

check_LTLIBRARIES += libcantor64gmt.la
libcantor64gmt_la_SOURCES=gf2x-cantor-fft.c
libcantor64gmt_la_CFLAGS= $(AM_CFLAGS) $(C64gmt)

check_LTLIBRARIES += libcantor128.la
libcantor128_la_SOURCES=gf2x-cantor-fft.c
libcantor128_la_CFLAGS= $(AM_CFLAGS) $(C128)

check_LTLIBRARIES += libcantor128nt.la
libcantor128nt_la_SOURCES=gf2x-cantor-fft.c
libcantor128nt_la_CFLAGS= $(AM_CFLAGS) $(C128nt)

check_LTLIBRARIES += libcantor128gm.la
libcantor128gm_la_SOURCES=gf2x-cantor-fft.c
libcantor128gm_la_CFLAGS= $(AM_CFLAGS) $(C128gm)

check_LTLIBRARIES += libcantor128gmt.la
libcantor128gmt_la_SOURCES=gf2x-cantor-fft.c
libcantor128gmt_la_CFLAGS= $(AM_CFLAGS) $(C128gmt)

endif
