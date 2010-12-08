#  This file is part of the gf2x library.
# 
#  Copyright 2007, 2008, 2009, 2010
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
# 
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2.1 of the License, or (at
#  your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License for more details.
#  
#  You should have received a copy of the GNU Lesser General Public
#  License along with CADO-NFS; see the file COPYING.  If not, write to
#  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
#  Boston, MA 02110-1301, USA.

AC_DEFUN([WORDSIZE_CODE],[
/* We check wraparound rather than zero, because that's the only thing
   the norm guarantees (C99) -- UINT_MAX isn't committed to being a power
   of two */
#include <stdio.h>
int main() {
    unsigned long x = 1UL;
    unsigned long y;
    FILE * f = fopen("conftest.out","w");
    int i = 1;
    for( ; ; i++) {
        y = x << 1;
        if (y < x) {
            break;
        }
        x = y;
    }
    fprintf(f,"%d\n",i);
    fclose(f);
    return 0;
}
])

AC_DEFUN([RUNTIME_ULONG_BITS],[
    if test x$gf2x_cv_ulongbits = x ; then
    AC_CACHE_CHECK([the number of bits in an unsigned long],
        [gf2x_cv_ulongbits],[
        AC_RUN_IFELSE([WORDSIZE_CODE()],[
            detected=`cat conftest.out | tr -d -c 0-9`
            if test x$detected = x ; then
                AC_MSG_ERROR([test program failed])
            else
                gf2x_cv_ulongbits=$detected
            fi
        ],[
            AC_MSG_FAILURE([cannot compile/run test program])
        ],[
            AC_MSG_NOTICE([check skipped because of cross-compiling])
            gf2x_cv_ulongbits=dontknow
        ])
    ])
    fi
])

AC_DEFUN([VERIFY_WORDSIZE],[
    RUNTIME_ULONG_BITS()
    AC_MSG_CHECKING([$2])
    case x$gf2x_cv_ulongbits in
        xdontknow) AC_MSG_NOTICE([cannot tell (cross-compiling)]);;
        x$1) AC_MSG_RESULT([yes]);;
        *)   AC_MSG_ERROR([no, $gf2x_cv_ulongbits-bit. Please provide appropriate \$CC variable]);;
    esac
])

AC_DEFUN([SSE2_EXAMPLE],[
#include <emmintrin.h>
__v2di x;
int main() {}
])

# Check whether we need some flag such as -msse2 in order to enable sse-2
# support
AC_DEFUN([CHECK_SSE2_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile sse-2 code], [gf2x_cv_cc_supports_sse2],[
  gf2x_cv_cc_supports_sse2=no
  if test "x${enable_sse2}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   AC_COMPILE_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cc_supports_sse2=yes
   ],[
    CFLAGS="$ac_save_CFLAGS -msse2"
    AC_COMPILE_IFELSE([SSE2_EXAMPLE()],[
     gf2x_cv_cc_supports_sse2="requires -msse2"
    ],[
     gf2x_cv_cc_supports_sse2=no
    ])
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse2" = "requires -msse2" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -msse2 is also needed by the preprocessor],
   [gf2x_cv_cpp_requires_msse2_flag],[
   AC_PREPROC_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cpp_requires_msse2_flag=no
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS -msse2"
    AC_PREPROC_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cpp_requires_msse2_flag=yes
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse sse-2 !])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse2" = "requires -msse2" ;then
  CFLAGS="$CFLAGS -msse2"
 fi
 if test "$gf2x_cv_cpp_requires_msse2_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -msse2"
 fi
 if test "$gf2x_cv_cc_supports_sse2" != "no" ;then
  AC_DEFINE([HAVE_SSE2_SUPPORT],[1],[Define if sse-2 code as present in the source tree is supported by the compiler])
 fi
])# CHECK_SSE2_SUPPORT



AC_DEFUN([PCLMUL_EXAMPLE],[
#include <wmmintrin.h>
#include <assert.h>
int main() {
assert(sizeof(unsigned long) == 8); /* assume 64-bit */
typedef union { __v2di s; unsigned long x[[2]]; } __v2di_proxy;
__v2di xx, yy;
__v2di_proxy zz;
xx = (__v2di) { 23, 0 };
yy = (__v2di) { 47, 0 };
zz.s = _mm_clmulepi64_si128(xx, yy, 0);
return zz.x[[0]] - 61;
}
])

# Check whether we need some flag such as -mpclmul in order to enable pclmulqdq
# support
AC_DEFUN([CHECK_PCLMUL_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile pclmulqdq and if it is supported by the hardware], [gf2x_cv_cc_supports_pclmul],[
  gf2x_cv_cc_supports_pclmul=no
  if test "x${enable_pclmul}" = xno ; then
   echo $ECHO_N " disabled, "
  else
   AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cc_supports_pclmul=yes
   ],[
    CFLAGS="$ac_save_CFLAGS -mpclmul"
    AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
     gf2x_cv_cc_supports_pclmul="requires -mpclmul"
    ],[
     gf2x_cv_cc_supports_pclmul=no
    ])
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_pclmul" = "requires -mpclmul" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -mpclmul is also needed by the preprocessor],
   [gf2x_cv_cpp_requires_mpclmul_flag],[
   AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpclmul_flag=no
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS -mpclmul"
    AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_requires_mpclmul_flag=yes
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse pclmul !])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_pclmul" = "requires -mpclmul" ;then
  CFLAGS="$CFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cpp_requires_mpclmul_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cc_supports_pclmul" != "no" ;then
  AC_DEFINE([HAVE_PCLMUL_SUPPORT],[1],[Define if pclmul as present in the source tree is supported by the compiler and hardware])
 fi
])# CHECK_PCLMUL_SUPPORT




AC_DEFUN([AC_COMPILE_WARNINGS], [
AC_MSG_CHECKING([warning verbosity option])
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_CXX])

  AC_ARG_WITH([compile-warnings],
              AS_HELP_STRING([--without-compile-warnings],
                             [Disable warning verbosity]),
              [ac_compile_warnings_on="$withval"],
              [ac_compile_warnings_on=""])

  if test x"$ac_compile_warnings_on" = xno
  then
    ac_compile_warnings_msg=no
  else
    if test -n "$CXX"
    then
      if test "$GXX" = "yes"
      then
        ac_compile_warnings_opt='-Wall -W'
      fi
      CXXFLAGS="$ac_compile_warnings_opt $CXXFLAGS"
      ac_compile_warnings_msg="$ac_compile_warnings_opt for C++"
    fi

  if test -n "$CC"
  then
    if test "$GCC" = "yes"
    then
      ac_compile_warnings_opt='-Wall -W'
    fi
    CFLAGS="$ac_compile_warnings_opt $CFLAGS"
    ac_compile_warnings_msg="$ac_compile_warnings_msg $ac_compile_warnings_opt for C"
  fi
  fi
  AC_MSG_RESULT([$ac_compile_warnings_msg])
  unset ac_compile_warnings_msg
  unset ac_compile_warnings_opt
])


