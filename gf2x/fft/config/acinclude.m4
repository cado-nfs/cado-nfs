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

AC_DEFUN([WORDSIZE_CODE],[AC_LANG_SOURCE([
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
])])

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
    # A fallback for generic builds. Otherwise unused.
    AC_CHECK_SIZEOF([unsigned long])
    gf2x_wordsize=`expr 8 \* $ac_cv_sizeof_unsigned_long`
    # AC_SUBST([GF2X_WORDSIZE_VALUE], [$gf2x_wordsize])
    AC_DEFINE_UNQUOTED([GF2X_WORDSIZE], [$gf2x_wordsize])
    RUNTIME_ULONG_BITS()
    AC_MSG_CHECKING([whether $CC indeed has $gf2x_wordsize-bit unsigned longs])
    case x$gf2x_cv_ulongbits in
        xdontknow) AC_MSG_NOTICE([cannot tell (cross-compiling)]);;
        x$gf2x_wordsize) AC_MSG_RESULT([yes]);;
        *)   AC_MSG_ERROR([no, $gf2x_cv_ulongbits-bit. Please provide appropriate \$CC variable]);;
    esac
    AH_TEMPLATE([GF2X_WORDSIZE],[Define this to the number of bits in an unsigned long])
])

AC_DEFUN([SSE2_EXAMPLE],[AC_LANG_SOURCE([
#include <emmintrin.h>
int main() {
    __m128i x = _mm_setr_epi64((__m64) (uint64_t) 42, (__m64) (uint64_t) 17);
    __m128d g = _mm_set_pd(42.0, 17.0);
    __m128i shift = _mm_setr_epi64((__m64) (uint64_t) 42, (__m64) (uint64_t) 17);
    x = _mm_srl_epi64(x, shift);
    return 0;
}
])])

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
 AM_CONDITIONAL([GF2X_SSE2_AVAILABLE],[test "x$gf2x_cv_cc_supports_sse2" != xno])
])# CHECK_SSE2_SUPPORT



AC_DEFUN([PCLMUL_EXAMPLE],[AC_LANG_SOURCE([
#include <wmmintrin.h>
#include <assert.h>
int main() {
    assert(sizeof(unsigned long) == 8); /* assume 64-bit */
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
    /* _m128i from 2 int64_t's */
#define _gf2x_mm_setr_epi64(lo, hi)                                     \
    _mm_setr_epi64(                                                     \
            _gf2x_mm_cvtsi64_m64((int64_t) (lo)),                       \
            _gf2x_mm_cvtsi64_m64((int64_t) (hi))                        \
            )
    /* _m128i from 1 int64_t's */
#define _gf2x_mm_set1_epi64(u) _mm_set1_epi64( _gf2x_mm_cvtsi64_m64((int64_t) (u)))
    __m128i a = _gf2x_mm_set1_epi64(17);
    __m128i b = _gf2x_mm_set1_epi64(42);
    union { __m128i s; unsigned long x[[2]]; } proxy;
    proxy.s = _mm_clmulepi64_si128(a, b, 0);
    return proxy.x[[0]] - 650;
}
])])

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
   ],[
   echo $ECHO_N " cross-compiling, "
   gf2x_cv_cc_supports_pclmul=no
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
 AM_CONDITIONAL([GF2XFFT_PCLMUL_SOURCES],[test "$gf2x_cv_cc_supports_pclmul" != "no"])
])# CHECK_PCLMUL_SUPPORT




AC_DEFUN([AC_COMPILE_WARNINGS], [
AC_MSG_CHECKING([warning verbosity option])
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH([compile-warnings],
              AS_HELP_STRING([--without-compile-warnings],
                             [Disable warning verbosity]),
              [ac_compile_warnings_on="$withval"],
              [ac_compile_warnings_on=""])

  if test x"$ac_compile_warnings_on" = xno
  then
    ac_compile_warnings_msg=no
  else

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



AC_DEFUN([GF2X_PROG_CC_TRY_MORE_CFLAGS],[
    AC_REQUIRE([AC_PROG_CC])
    for i in "march" "mtune" ; do
        flag="-$i=native"
        # Do not add another march/mtune flag if one is already present !
        if ! echo $CFLAGS | grep -q $i ; then
        GF2X_PROG_CC_FOR_BUILD_WORKS([$CC_FOR_BUILD $CFLAGS $flag],[
            CFLAGS="$CFLAGS $flag"
        ])
        fi
    done
])



dnl -- taken from gmp-4.2.1, LGPL v2.1+ --
dnl -- renamed GMP_ to GF2X_ --
dnl
dnl
dnl
dnl  GF2X_PROG_CC_FOR_BUILD
dnl  ---------------------
dnl  Establish CC_FOR_BUILD, a C compiler for the build system.
dnl
dnl  If CC_FOR_BUILD is set then it's expected to work, likewise the old
dnl  style HOST_CC, otherwise some likely candidates are tried, the same as
dnl  configfsf.guess.

AC_DEFUN([GF2X_PROG_CC_FOR_BUILD],
[AC_REQUIRE([AC_PROG_CC])
if test -n "$CC_FOR_BUILD"; then
  GF2X_PROG_CC_FOR_BUILD_WORKS($CC_FOR_BUILD,,
    [AC_MSG_ERROR([Specified CC_FOR_BUILD doesn't seem to work])])
elif test -n "$HOST_CC"; then
  GF2X_PROG_CC_FOR_BUILD_WORKS($HOST_CC,
    [CC_FOR_BUILD=$HOST_CC],
    [AC_MSG_ERROR([Specified HOST_CC doesn't seem to work])])
else
  for i in "$CC" "$CC $CFLAGS $CPPFLAGS" cc gcc c89 c99; do
    GF2X_PROG_CC_FOR_BUILD_WORKS($i,
      [CC_FOR_BUILD=$i
       break])
  done
  if test -z "$CC_FOR_BUILD"; then
    AC_MSG_ERROR([Cannot find a build system compiler])
  fi
fi
    
AC_ARG_VAR(CC_FOR_BUILD,[build system C compiler])
AC_SUBST(CC_FOR_BUILD)
GF2X_PROG_CC_TRY_MORE_CFLAGS()
])

dnl  GF2X_PROG_CC_FOR_BUILD_WORKS(cc/cflags[,[action-if-good][,action-if-bad]])
dnl  -------------------------------------------------------------------------
dnl  See if the given cc/cflags works on the build system.
dnl
dnl  It seems easiest to just use the default compiler output, rather than
dnl  figuring out the .exe or whatever at this stage.

AC_DEFUN([GF2X_PROG_CC_FOR_BUILD_WORKS],
[AC_MSG_CHECKING([build system compiler $1])
# remove anything that might look like compiler output to our "||" expression
rm -f conftest* a.out b.out a.exe a_out.exe
cat >conftest.c <<EOF
int
main ()
{
  exit(0);
}
EOF
gf2x_compile="$1 conftest.c"
cc_for_build_works=no
if AC_TRY_EVAL(gf2x_compile); then
  if (./a.out || ./b.out || ./a.exe || ./a_out.exe || ./conftest) >&AC_FD_CC 2>&1; then
    cc_for_build_works=yes
  fi
fi
rm -f conftest* a.out b.out a.exe a_out.exe
AC_MSG_RESULT($cc_for_build_works)
if test "$cc_for_build_works" = yes; then
  ifelse([$2],,:,[$2])
else
  ifelse([$3],,:,[$3])
fi
])

dnl  GF2X_PROG_EXEEXT_FOR_BUILD
dnl  -------------------------
dnl  Determine EXEEXT_FOR_BUILD, the build system executable suffix.
dnl
dnl  The idea is to find what "-o conftest$foo" will make it possible to run
dnl  the program with ./conftest.  On Unix-like systems this is of course
dnl  nothing, for DOS it's ".exe", or for a strange RISC OS foreign file
dnl  system cross compile it can be ",ff8" apparently.  Not sure if the
dnl  latter actually applies to a build-system executable, maybe it doesn't,
dnl  but it won't hurt to try.

AC_DEFUN([GF2X_PROG_EXEEXT_FOR_BUILD],
[AC_REQUIRE([GF2X_PROG_CC_FOR_BUILD])
AC_CACHE_CHECK([for build system executable suffix],
               gf2x_cv_prog_exeext_for_build,
[cat >conftest.c <<EOF
int
main ()
{
  exit (0);
}
EOF
for i in .exe ,ff8 ""; do
  gf2x_compile="$CC_FOR_BUILD conftest.c -o conftest$i"
  if AC_TRY_EVAL(gf2x_compile); then
    if (./conftest) 2>&AC_FD_CC; then
      gf2x_cv_prog_exeext_for_build=$i
      break
    fi
  fi
done
rm -f conftest*
if test "${gf2x_cv_prog_exeext_for_build+set}" != set; then
  AC_MSG_ERROR([Cannot determine executable suffix])
fi
])
AC_SUBST(EXEEXT_FOR_BUILD,$gf2x_cv_prog_exeext_for_build)
])
