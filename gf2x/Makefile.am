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

AUTOMAKE_OPTIONS=subdir-objects
ACLOCAL_AMFLAGS=-I config --install
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)

FFT_SOURCES=fft/gf2x-ternary-fft.c
if ENABLE_FFT_INTERFACE
AM_CPPFLAGS += -I$(srcdir)/mpfq
FFT_SOURCES+=fft/gf2x-fake-fft.c
FFT_SOURCES+=fft/gf2x-cantor-fft.c
AM_CPPFLAGS += -I$(top_srcdir)/fft/mpfq
gf2x_includedir=$(includedir)/gf2x
gf2x_include_HEADERS=	\
	fft/gf2x-fft.h			\
	fft/gf2x-fake-fft.h		\
	fft/gf2x-cantor-fft.h		\
	fft/gf2x-ternary-fft.h
endif


lib_LTLIBRARIES=libgf2x.la
if ENABLE_FFT_INTERFACE
# default which gets shipped
libgf2x_la_CFLAGS=$(AM_CFLAGS) -DCANTOR_BASE_FIELD_SIZE=128
endif
libgf2x_la_SOURCES=gf2x.c toom.c toom128.c toom-gpl.c $(FFT_SOURCES)




nobase_include_HEADERS=gf2x.h	\
		gf2x/gf2x-impl.h		\
		gf2x/gf2x-impl-export.h		\
		gf2x/gf2x-small.h

# All these are generated.
nobase_nodist_include_HEADERS=	\
		gf2x/gf2x-config.h		\
		gf2x/gf2x-config-export.h	\
		gf2x/gf2x-thresholds.h		\
		gf2x/gf2x_mul1.h		\
		gf2x/gf2x_mul2.h		\
		gf2x/gf2x_mul3.h		\
		gf2x/gf2x_mul4.h		\
		gf2x/gf2x_mul5.h		\
		gf2x/gf2x_mul6.h		\
	        gf2x/gf2x_mul7.h		\
		gf2x/gf2x_mul8.h		\
		gf2x/gf2x_mul9.h


libgf2x_la_LDFLAGS=-version-info @gf2x_lib_version@ -no-undefined

SUBDIRS=lowlevel src . fft tests $(MAYBE_APPS)

DIST_SUBDIRS=lowlevel src . fft tests apps

DISTCLEANFILES=gf2x/gf2x-thresholds.h	\
	       gf2x/gf2x_mul1.h		\
	       gf2x/gf2x_mul2.h		\
	       gf2x/gf2x_mul3.h		\
	       gf2x/gf2x_mul4.h		\
	       gf2x/gf2x_mul5.h		\
	       gf2x/gf2x_mul6.h		\
	       gf2x/gf2x_mul7.h		\
	       gf2x/gf2x_mul8.h		\
	       gf2x/gf2x_mul9.h

EXTRA_DIST=
EXTRA_DIST+=version.sh.in
EXTRA_DIST+=BUGS
EXTRA_DIST+=already_tuned
EXTRA_DIST+=config/configfsf.sub
EXTRA_DIST+=config/configfsf.guess
EXTRA_DIST+=toom-gpl-placeholder.c

# This is really a maintainer-only command. This can be used for creating
# a tarball with an LGPL-licensed gf2x. Note that when this is done,
# configure must be re-run !
untaint:
	cp $(srcdir)/toom-gpl-placeholder.c $(srcdir)/toom-gpl.c

dist-hook:
	-find $(distdir)/already_tuned/tuned \! -type d -a \! -name README | xargs rm
	-git log --boundary 29b13131cfaf95132f6a0602d3e0575694359d11..HEAD > $(distdir)/ChangeLog

# This one is simply taken from the generated Makefile. I agree it's a
# bit ugly.
dist-LGPL: distdir
	mv $(distdir)/ $(distdir)-LGPL/
	cp $(distdir)-LGPL/toom-gpl-placeholder.c $(distdir)-LGPL/toom-gpl.c
	tardir=$(distdir)-LGPL && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(distdir)-LGPL.tar.gz
	mv $(distdir)-LGPL/ $(distdir)/
	$(am__post_remove_distdir)

tune-lowlevel tune-toom tune-fft: all
	cd src; $(MAKE) $@

distclean-local:
	-find already_tuned/tuned -maxdepth 1 -name "gf2x*" | xargs -r rm -f
