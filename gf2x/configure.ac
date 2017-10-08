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

AC_INIT([gf2x], [1.2.1])

if test "x$CFLAGS" != x ; then
 has_CFLAGS_from_env=yes
else
 has_CFLAGS_from_env=no
fi

: ${CFLAGS=-g -O4}
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CANONICAL_TARGET

AC_LANG([C])

AC_FUNC_ALLOCA

AC_SUBST([gf2x_lib_version], [1:2:0])

AM_INIT_AUTOMAKE([1.13])

# disable for releases: AM_MAINTAINER_MODE([disable])
AM_MAINTAINER_MODE([disable])

LT_INIT()

AC_ARG_ENABLE([ntl-checks],
              [AS_HELP_STRING([--enable-ntl-checks],
                              [Turn on ntl checks])])

AC_ARG_ENABLE([hardware-specific-code],
              [AS_HELP_STRING([--enable-hardware-specific-code],
                              [Turn on all cpu-dependent optimized code (default is yes)])],
              [],[enable_hardware_specific_code=yes])

AC_ARG_ENABLE([sse2],
              [AS_HELP_STRING([--enable-sse2],
                              [Turn on sse-2 code (default is yes)])])

AC_ARG_ENABLE([sse3],
              [AS_HELP_STRING([--enable-sse3],
                              [Turn on sse-3 code (default is yes)])])

AC_ARG_ENABLE([ssse3],
              [AS_HELP_STRING([--enable-ssse3],
                              [Turn on ssse3 code (default is yes)])])

AC_ARG_ENABLE([sse41],
              [AS_HELP_STRING([--enable-sse41],
                              [Turn on sse-4.1 code (default is yes)])])

AC_ARG_ENABLE([pclmul],
              [AS_HELP_STRING([--enable-pclmul],
                              [Turn on pclmul code (default is yes)])])

if test "x$enable_ntl_checks" = xyes ; then
 MAYBE_APPS=apps
else
 MAYBE_APPS=
fi
AC_SUBST([MAYBE_APPS])

AC_ARG_ENABLE([fft-interface],
              [AS_HELP_STRING([--enable-fft-interface],
                [compile with the optional fft interface (default is no)])],
                [],[])
AM_CONDITIONAL([ENABLE_FFT_INTERFACE],[test "x$enable_fft_interface" = xyes])

AC_PROG_CC
AC_PROG_CXX
AC_COMPILE_WARNINGS
AC_PROG_CC_C99
AC_PROG_CPP

GF2X_PROG_CC_FOR_BUILD
GF2X_PROG_EXEEXT_FOR_BUILD

# This macro is used for tuning
AM_PROG_CC_C_O

# A fallback for generic builds. Otherwise unused.
AC_CHECK_SIZEOF([unsigned long])
gf2x_wordsize=`expr 8 \* $ac_cv_sizeof_unsigned_long`
AC_SUBST([GF2X_WORDSIZE], [$gf2x_wordsize])

AM_CONDITIONAL([WORDSIZE_IS_64],[test "x$gf2x_wordsize" = x64])


if test "x${enable_hardware_specific_code}" = xno ; then
 echo "Not using hardware-specific code due to --disable-hardware-specific-code flag"
else
 # It is necessary to make all tests. We do encounter in the wild
 # binutils (openbsd binutils 2.15, namely) which are buggy with ssse3,
 # and that isn't extremely quickly spotted by the later checks...
 CHECK_SSE2_SUPPORT()
 CHECK_SSE3_SUPPORT()
 CHECK_SSSE3_SUPPORT()
 CHECK_SSE41_SUPPORT()
 CHECK_PCLMUL_SUPPORT()
 # The undocumented environment variable DISABLE_MARCH_NATIVE_TEST can be
 # used to avoid any attempt to add -march=native flag to gcc. This is
 # unfortunately mandatory on some virtual machines for which the cpuid
 # capabilities are incorrectly reported (leading to i386 guests with
 # BMI2, or other similar absurdities).
 if test "x$DISABLE_MARCH_NATIVE_TEST" = x ; then
  if test "x$has_CFLAGS_from_env" = xno ; then
   CHECK_MARCH_NATIVE_SUPPORT()
   CHECK_MTUNE_NATIVE_SUPPORT()
  fi
 else
   echo "Not doing test for -march=native as per environment variable $DISABLE_MARCH_NATIVE_TEST"
  fi
fi

: ${ABI='default'}

if test x$ABI != xdefault ; then
 VERIFY_WORDSIZE([$ABI],[whether $CC and ABI=$ABI fit together])
fi

# Many values for target_cpu in the switch below are dead code, because
# for licensing reasons we are stuck with gmp's old config.guess from
# 2006. It's not that bad, since we check for the features we need by
# ourselves anyway.
if test x$hwdir = x ; then
 # Our _default_ setting for hwdir is to use 8*sizeof(ulong).
 hwdir=generic$gf2x_wordsize

 if test "x${enable_hardware_specific_code}" = xyes ; then
  # Yet there are several situations in which we select another directory.
  # This refined selection mechanism is impacted by the ABI= setting.
  case "$target_cpu" in
   # All 32-bit pentium patterns
   i?86|pentium*|athlon|prescott)
    if test "$gf2x_cv_cc_supports_sse2" = "no" ; then
     hwdir=generic$gf2x_wordsize
    else
     hwdir=x86_sse2
     # Maybe on some funny platform from hell, we could get i386 yet have a
     # 64-bit cpu ?
     if test x$ABI = x64 ; then
      hwdir=x86_64
     fi
    fi
    ;;
   core2|opteron|x86_64|nocona|k10)
    # Note: could set to something else if we have separate tuning files.
    if test "$gf2x_cv_cc_supports_pclmul" = "no" ; then
      hwdir=x86_64
    else
      hwdir=x86_64_pclmul
    fi
    if test x$ABI = x32 ; then
     if test "$gf2x_cv_cc_supports_sse2" = "no" ; then
      hwdir=generic32
     else
      hwdir=x86_sse2
     fi
    elif test x$ABI = xdefault ; then
     case x$gf2x_wordsize in
      x64|xdontknow) true ;;
      x32)
        AC_MSG_WARN([The default ABI for this compiler has $gf2x_wordsize-bit unsigned longs, not 64-bit as the CPU supports. Using $gf2x_wordsize-bit])
        if test "$gf2x_cv_cc_supports_sse2" = "no" ; then
         hwdir=generic$gf2x_wordsize
        else
         hwdir=x86_sse2
        fi;;
      *)
       AC_MSG_ERROR([The default ABI for this compiler has $gf2x_wordsize-bit unsigned longs, this is not supported])
       ;;
     esac
    fi
    ;;
   coreinhm) hwdir=x86_64;;
   coreiwsm) hwdir=x86_64_pclmul ;;
   coreisbr) hwdir=x86_64_pclmul ;;
   coreiibr) hwdir=x86_64_pclmul ;;
   coreihwl) hwdir=coreihwl ;;
   armv7l) hwdir=armv7l ;;
   neon) hwdir=armv7l ;;
   # TODO: add bobcat, bulldozer.
  esac
 fi
else
 AC_MSG_NOTICE([Using supplied hwdir value $hwdir])
fi
case "$hwdir" in
 x86_64*|coreihwl)
    if test "$gf2x_cv_cc_supports_sse2" = "no" ; then
      # this is only if we've been given an explicit --disable-sse2
      hwdir=x86_64_nosse
     fi;;
esac

AC_MSG_NOTICE([using ABI="$ABI"])
AC_MSG_NOTICE([      CC="$CC"])
AC_MSG_NOTICE([      CFLAGS="$CFLAGS"])
AC_MSG_NOTICE([      CPPFLAGS="$CPPFLAGS"])
AC_MSG_NOTICE([      hwdir="$hwdir"])

for f in gf2x-thresholds.h      \
        gf2x_mul1.h gf2x_mul2.h gf2x_mul3.h gf2x_mul4.h \
        gf2x_mul5.h gf2x_mul6.h gf2x_mul7.h gf2x_mul8.h \
        gf2x_mul9.h ; do
 who=tuned
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=$hwdir; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=generic$gf2x_wordsize; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=generic; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then AC_MSG_ERROR([$f not found]); fi
 AC_CONFIG_LINKS([gf2x/$f:already_tuned/$who/$f])
 if test "$f" = "gf2x-thresholds.h" ; then
    tuned_nbits=[`sed -n 's/^#define[ 	][ 	]*GF2X_WORDSIZE[ 	][ 	]*\([0-9][0-9]*\).*$/\1/p' $srcdir/already_tuned/$who/$f`]
 fi
done

if test x$tuned_nbits = x ; then
 tuned_nbits=$gf2x_wordsize
else
 if test x$ABI != xdefault ; then
  if test x$tuned_nbits != x$ABI ; then
   AC_MSG_ERROR([already_tuned/$hwdir/ assumes $tuned_nbits-bits unsigned longs, which conflicts with ABI=$ABI])
  fi
 else
  # At this point we haven't run the standard check.
  VERIFY_WORDSIZE([$tuned_nbits],[whether already_tuned/$hwdir/ is right assuming $tuned_nbits-bits unsigned longs])
 fi
fi

# This is used in src/Makefile.am
AM_CONDITIONAL([GF2X_32BIT_SOURCES],[test "x$tuned_nbits" = x32])
AM_CONDITIONAL([GF2X_64BIT_SOURCES],[test "x$tuned_nbits" = x64])
AM_CONDITIONAL([GF2X_SSE2_AVAILABLE],[expr "x$gf2x_cv_cc_supports_sse2" : "xyes.*"])
AM_CONDITIONAL([GF2X_PCLMUL_AVAILABLE],[expr "x$gf2x_cv_cc_supports_pclmul" : "xyes"])
AM_CONDITIONAL([GPL_CODE_PRESENT],[grep -q "released under the GPL" "$srcdir/toom-gpl.c"])
if grep -q "released under the GPL" "$srcdir/toom-gpl.c" ; then
 echo "File $srcdir/toom-gpl.c is present, using a GPL-licensed gf2x"
 AC_DEFINE([GPL_CODE_PRESENT],1,[Have GPL code])
else
 echo "File $srcdir/toom-gpl.c is a placeholder, using an LGPL-licensed gf2x"
 AC_DEFINE([GPL_CODE_PRESENT],0,[Have LGPL code (no GPL taint)])
fi

AC_CONFIG_HEADERS([gf2x/gf2x-config.h gf2x/gf2x-config-export.h])

#if test -d "$srcdir/gf2x-fft" ; then
# AC_MSG_NOTICE([Enabling optional gf2x-fft])
# AC_CONFIG_SUBDIRS([gf2x-fft])
# MAYBE_GF2X_FFT=gf2x-fft
#else
# AC_MSG_NOTICE([optional gf2x-fft not found])
# MAYBE_GF2X_FFT=
#fi
#AC_SUBST([MAYBE_GF2X_FFT])

AC_CONFIG_FILES([version.sh Makefile tests/Makefile lowlevel/Makefile
                 src/Makefile fft/Makefile apps/Makefile])

AC_OUTPUT

# vim: set sw=1:
