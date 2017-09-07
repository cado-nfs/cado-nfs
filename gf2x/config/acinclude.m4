#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2012, 2013
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
            # see bug #15631 and autoconf manual about tr.
            # detected=`cat conftest.out | tr -d -c 0-9`
            detected=`cat conftest.out`
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

# pepper this check with more sse-2 only statements, or we might be
# fooled by some early athlon64 cpus supporting extended 3dnow, which
# includes a subset of sse-2, but do not support the full sse-2 insn set.
AC_DEFUN([SSE2_EXAMPLE],[AC_LANG_SOURCE([
#include <emmintrin.h>
#include <stdint.h>
#include <emmintrin.h>

int main(int argc, char * argv[[]]) {
    volatile int a0 = 17;
    volatile int a1 = 42;
    __m128i foo = _mm_setr_epi32(argc, argc + 1, argc + 2, argc + 3);
    __m128i bar = _mm_setr_epi32(argc + 3, argc + 2, argc + 1, argc);
    __m128i x = _mm_setr_epi32(a1, 0, a0, 0);
    __m128d g = _mm_set_pd((double) a1, (double) a0);
    x = _mm_srl_epi64(x, _mm_setr_epi32(2,0,0,0));
    foo = _mm_mullo_epi16(foo, bar);
    foo = _mm_slli_epi64(foo, 1);
    foo = _mm_xor_si128(bar, _mm_unpacklo_epi32(foo, bar));
    foo = _mm_srli_epi64(foo, 1);
    foo = _mm_mullo_epi16(foo, bar);
    foo = _mm_shuffle_epi32(foo, 78);
    foo = _mm_xor_si128(bar, _mm_unpacklo_epi32(foo, bar));
    foo = _mm_srli_si128(foo, 1);
    foo = _mm_xor_si128(foo, x);

    return _mm_extract_epi16(foo, 0) & (argc - 1);
}
])])

AC_DEFUN([SSE3_EXAMPLE],[AC_LANG_SOURCE([
/* This source file is our test case for sse-3 support. */
#include <stdint.h>
#include <string.h>
#include <pmmintrin.h>

int main()
{
    volatile double a0 = 12.34;
    volatile double a1 = 56.78;
    __m128d x = _mm_setr_pd(a0, 34.12);
    __m128d y = _mm_setr_pd(78.56, a1);
    double a[[2]], b[[2]] = { 78.56 + 56.78, 12.34 + 34.12 };

    y = _mm_hadd_pd(y, x);
    memcpy(a, &y, 16);
    return (a[[0]] != b[[0]] || a[[1]] != b[[1]]);
}
])])

AC_DEFUN([SSSE3_EXAMPLE],[AC_LANG_SOURCE([
/* This source file is our test case for ssse3 support. */
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>

int main()
{
    volatile uint32_t a0 = 0x03020100;
    volatile uint32_t a1 = 0x1F1E1D1C;
    __m128i x = _mm_setr_epi32(a0, 0x07060504, 0x0B0A0908, 0x0F0E0D0C);
    __m128i y = _mm_setr_epi32(0x13121110, 0x17161514, 0x1B1A1918, a1);
    uint64_t a[[2]], b[[2]] = { 0x0C0B0A0908070605, 0x14131211100F0E0D };
    y = _mm_alignr_epi8(y, x, 0x5);
    memcpy (a, &y, 16);
    return(a[[0]] != b[[0]] || a[[1]] != b[[1]]);
}
])])

AC_DEFUN([SSE41_EXAMPLE],[AC_LANG_SOURCE([
#include <stdint.h>
#include <stdlib.h>
#include <smmintrin.h>

int main() {
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    volatile int a0 = 17;
    volatile int a1 = 42;
    __m128i x = _mm_setr_epi32(a1, 0, a0, 0);
    // x = 0 0x2a 0 0x11
    __m128i y = _mm_setr_epi32(42, 0, 17, 0);
    // y = 0 0x2a 0 0x11
    __m128i ma = _mm_max_epi32(x, y);
    __m128i mi = _mm_min_epi32(x, y);
    __m128i z = _mm_xor_si128(mi, ma);
    int ok0 = _mm_testz_si128(z, z);
    __m128i c = _mm_cmpeq_epi64(x, y);
    int ok1 = _mm_extract_epi32(c, 0) && _mm_extract_epi32(c, 1);
    return (ok0 && ok1) ? EXIT_SUCCESS : EXIT_FAILURE;
}
])])
AC_DEFUN([PCLMUL_EXAMPLE],[AC_LANG_SOURCE([
#include <stdint.h>
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
    volatile int a0 = 17;
    volatile int a1 = 42;
    __m128i a = _gf2x_mm_set1_epi64(a0);
    __m128i b = _gf2x_mm_set1_epi64(a1);
    union { __m128i s; unsigned long x[[2]]; } proxy;
    proxy.s = _mm_clmulepi64_si128(a, b, 0);
    return proxy.x[[0]] - 650;
}
])])

# Check whether we need some flag such as -msse2 in order to enable sse-2
# support
AC_DEFUN([CHECK_SSE2_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile and run sse-2 code], [gf2x_cv_cc_supports_sse2],[
  gf2x_cv_cc_supports_sse2=no
  if test "x${enable_sse2}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   CFLAGS="$ac_save_CFLAGS -msse2"
   AC_RUN_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cc_supports_sse2=yes
   ],[
    CFLAGS="$ac_save_CFLAGS"
    AC_RUN_IFELSE([SSE2_EXAMPLE()],[
     gf2x_cv_cc_supports_sse2="yes, but without -msse2"
    ],[
     gf2x_cv_cc_supports_sse2=no
    ])
   ],[
   echo $ECHO_N "cross-compiling, "
   if test "x${enable_sse2}" = xyes ; then
    echo $ECHO_N "explicitly enabled, "
    gf2x_cv_cc_supports_sse2=yes
   else
    gf2x_cv_cc_supports_sse2=no
   fi
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse2" = "yes" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -msse2 is also understood by the preprocessor],
   [gf2x_cv_cpp_understands_msse2_flag],[
   CPPFLAGS="$ac_save_CPPFLAGS -msse2"
   AC_PREPROC_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse2_flag=yes
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_PREPROC_IFELSE([SSE2_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse2_flag=no
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse sse-2!])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse2" = "yes" ;then
  CFLAGS="$CFLAGS -msse2"
 fi
 if test "$gf2x_cv_cpp_understands_msse2_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -msse2"
 fi
 if test "$gf2x_cv_cc_supports_sse2" != "no" ;then
  AC_DEFINE([GF2X_HAVE_SSE2_SUPPORT],[1],[Define if sse-2 code as present in the source tree is supported by the compiler])
 fi
])# CHECK_SSE2_SUPPORT

AC_DEFUN([CHECK_SSE3_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile and run sse-3 code], [gf2x_cv_cc_supports_sse3],[
  gf2x_cv_cc_supports_sse3=no
  if test "x${gf2x_cv_cc_supports_sse2}" = xno ; then
   echo $ECHO_N "skipped, "
  elif test "x${enable_sse3}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   CFLAGS="$ac_save_CFLAGS -msse3"
   AC_RUN_IFELSE([SSE3_EXAMPLE()],[
    gf2x_cv_cc_supports_sse3=yes
   ],[
    CFLAGS="$ac_save_CFLAGS"
    AC_RUN_IFELSE([SSE3_EXAMPLE()],[
     gf2x_cv_cc_supports_sse3="yes, but without -msse3"
    ],[
     gf2x_cv_cc_supports_sse3=no
    ])
   ],[
   echo $ECHO_N "cross-compiling, "
   if test "x${enable_sse3}" = xyes ; then
    echo $ECHO_N "explicitly enabled, "
    gf2x_cv_cc_supports_sse3=yes
   else
    gf2x_cv_cc_supports_sse3=no
   fi
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse3" = "yes" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -msse3 is also understood by the preprocessor],
   [gf2x_cv_cpp_understands_msse3_flag],[
   CPPFLAGS="$ac_save_CPPFLAGS -msse3"
   AC_PREPROC_IFELSE([SSE3_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse3_flag=yes
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_PREPROC_IFELSE([SSE3_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse3_flag=no
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse sse-3!])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse3" = "yes" ;then
  CFLAGS="$CFLAGS -msse3"
 fi
 if test "$gf2x_cv_cpp_understands_msse3_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -msse3"
 fi
 if test "$gf2x_cv_cc_supports_sse3" != "no" ;then
  AC_DEFINE([GF2X_HAVE_SSE3_SUPPORT],[1],[Define if sse-3 code as present in the source tree is supported by the compiler])
 fi
])# CHECK_SSE3_SUPPORT

AC_DEFUN([CHECK_SSSE3_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile and run ssse3 code], [gf2x_cv_cc_supports_ssse3],[
  gf2x_cv_cc_supports_ssse3=no
  if test "x${gf2x_cv_cc_supports_sse3}" = xno ; then
   echo $ECHO_N "skipped, "
  elif test "x${enable_ssse3}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   CFLAGS="$ac_save_CFLAGS -mssse3"
   AC_RUN_IFELSE([SSSE3_EXAMPLE()],[
    gf2x_cv_cc_supports_ssse3=yes
   ],[
    CFLAGS="$ac_save_CFLAGS"
    AC_RUN_IFELSE([SSSE3_EXAMPLE()],[
     gf2x_cv_cc_supports_ssse3="yes, but without -mssse3"
    ],[
     gf2x_cv_cc_supports_ssse3=no
    ])
   ],[
   echo $ECHO_N "cross-compiling, "
   if test "x${enable_ssse3}" = xyes ; then
    echo $ECHO_N "explicitly enabled, "
    gf2x_cv_cc_supports_ssse3=yes
   else
    gf2x_cv_cc_supports_ssse3=no
   fi
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_ssse3" = "yes" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -mssse3 is also understood by the preprocessor],
   [gf2x_cv_cpp_understands_mssse3_flag],[
   CPPFLAGS="$ac_save_CPPFLAGS -mssse3"
   AC_PREPROC_IFELSE([SSSE3_EXAMPLE()],[
    gf2x_cv_cpp_understands_mssse3_flag=yes
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_PREPROC_IFELSE([SSSE3_EXAMPLE()],[
    gf2x_cv_cpp_understands_mssse3_flag=no
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse ssse3!])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_ssse3" = "yes" ;then
  CFLAGS="$CFLAGS -mssse3"
 fi
 if test "$gf2x_cv_cpp_understands_mssse3_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -mssse3"
 fi
 if test "$gf2x_cv_cc_supports_ssse3" != "no" ;then
  AC_DEFINE([GF2X_HAVE_SSSE3_SUPPORT],[1],[Define if ssse3 code as present in the source tree is supported by the compiler])
 fi
])# CHECK_SSSE3_SUPPORT


AC_DEFUN([CHECK_SSE41_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile and run sse-4.1 code], [gf2x_cv_cc_supports_sse41],[
  gf2x_cv_cc_supports_sse41=no
  if test "x${gf2x_cv_cc_supports_ssse3}" = xno ; then
   echo $ECHO_N "skipped, "
  elif test "x${enable_sse41}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   CFLAGS="$ac_save_CFLAGS -msse4.1"
   AC_RUN_IFELSE([SSE41_EXAMPLE()],[
    gf2x_cv_cc_supports_sse41=yes
   ],[
    CFLAGS="$ac_save_CFLAGS"
    AC_RUN_IFELSE([SSE41_EXAMPLE()],[
     gf2x_cv_cc_supports_sse41="yes, but without -msse4.1"
    ],[
     gf2x_cv_cc_supports_sse41=no
    ])
   ],[
   echo $ECHO_N "cross-compiling, "
   if test "x${enable_sse41}" = xyes ; then
    echo $ECHO_N "explicitly enabled, "
    gf2x_cv_cc_supports_sse41=yes
   else
    gf2x_cv_cc_supports_sse41=no
   fi
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse41" = "yes" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -msse4.1 is also understood by the preprocessor],
   [gf2x_cv_cpp_understands_msse41_flag],[
   CPPFLAGS="$ac_save_CPPFLAGS -msse4.1"
   AC_PREPROC_IFELSE([SSE41_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse41_flag=yes
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_PREPROC_IFELSE([SSE41_EXAMPLE()],[
    gf2x_cv_cpp_understands_msse41_flag=no
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse sse-4.1!])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_sse41" = "yes" ;then
  CFLAGS="$CFLAGS -msse4.1"
 fi
 if test "$gf2x_cv_cpp_understands_msse41_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -msse4.1"
 fi
 if test "$gf2x_cv_cc_supports_sse41" != "no" ;then
  AC_DEFINE([GF2X_HAVE_SSE41_SUPPORT],[1],[Define if sse-4.1 code as present in the source tree is supported by the compiler])
 fi
])# CHECK_SSE41_SUPPORT


AC_DEFUN([CHECK_PCLMUL_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC can compile pclmulqdq and if it is supported by the hardware], [gf2x_cv_cc_supports_pclmul],[
  gf2x_cv_cc_supports_pclmul=no
  if test "x${gf2x_cv_cc_supports_sse41}" = xno ; then
   echo $ECHO_N "skipped, "
  elif test "x${enable_pclmul}" = xno ; then
   echo $ECHO_N "explicitly disabled, "
  else
   CFLAGS="$ac_save_CFLAGS -mpclmul"
   AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cc_supports_pclmul=yes
   ],[
    CFLAGS="$ac_save_CFLAGS"
    AC_RUN_IFELSE([PCLMUL_EXAMPLE()],[
     gf2x_cv_cc_supports_pclmul="yes, but without -mpclmul"
    ],[
     gf2x_cv_cc_supports_pclmul=no
    ])
   ],[
   echo $ECHO_N "cross-compiling, "
   if test "x${enable_pclmul}" = xyes ; then
    echo $ECHO_N "explicitly enabled, "
    gf2x_cv_cc_supports_pclmul=yes
   else
    gf2x_cv_cc_supports_pclmul=no
   fi
   ])
  fi
 ])
 ac_save_CPPFLAGS=$CPPFLAGS
 if test "$gf2x_cv_cc_supports_pclmul" = "yes" ;then
  # Tweaking CFLAGS is often not enough.
  AC_CACHE_CHECK([whether -mpclmul is also understood by the preprocessor],
   [gf2x_cv_cpp_understands_mpclmul_flag],[
   CPPFLAGS="$ac_save_CPPFLAGS -mpclmul"
   AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_understands_mpclmul_flag=yes
   ],[
    CPPFLAGS="$ac_save_CPPFLAGS"
    AC_PREPROC_IFELSE([PCLMUL_EXAMPLE()],[
    gf2x_cv_cpp_understands_mpclmul_flag=no
    ],[
     AC_MSG_ERROR([Sorry, the preprocessor can't parse pclmul!])
    ])
   ])
  ])
 fi
 CFLAGS=$ac_save_CFLAGS
 CPPFLAGS=$ac_save_CPPFLAGS
 if test "$gf2x_cv_cc_supports_pclmul" = "yes" ;then
  CFLAGS="$CFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cpp_understands_mpclmul_flag" = "yes" ; then
  CPPFLAGS="$CPPFLAGS -mpclmul"
 fi
 if test "$gf2x_cv_cc_supports_pclmul" != "no" ;then
  AC_DEFINE([GF2X_HAVE_PCLMUL_SUPPORT],[1],[Define if pclmul code as present in the source tree is supported by the compiler])
 fi
])# CHECK_PCLMUL_SUPPORT

AC_DEFUN([HELLO_WORLD_EXAMPLE],[AC_LANG_SOURCE([
#include <stdio.h>

int main() {
    printf("hello\n");
    return 0;
}
])])


AC_DEFUN([CHECK_MARCH_NATIVE_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 special_double_setting="yes, via -march=x86-64 -march=native"
 AC_CACHE_CHECK([whether $CC understands -march=native], [gf2x_cv_cc_supports_march_native],[
  gf2x_cv_cc_supports_march_native=no
  CFLAGS="$ac_save_CFLAGS -march=native"
  AC_COMPILE_IFELSE(
      [HELLO_WORLD_EXAMPLE()],
      [
      gf2x_cv_cc_supports_march_native=yes
      ],
      [
      CFLAGS="$ac_save_CFLAGS -march=x86-64 -march=native"
      AC_COMPILE_IFELSE(
          [HELLO_WORLD_EXAMPLE()],
          [
          gf2x_cv_cc_supports_march_native="$special_double_setting"
          ],
          [AC_MSG_RESULT(no)])
      ]
  )
  CFLAGS=$ac_save_CFLAGS
  ])
  if test "$gf2x_cv_cc_supports_march_native" = "$special_double_setting" ;then
    CFLAGS="$CFLAGS -march=x86-64 -march=native"
  elif test "$gf2x_cv_cc_supports_march_native" = "yes" ;then
    CFLAGS="$CFLAGS -march=native"
  fi
])# CHECK_MARCH_NATIVE_SUPPORT

AC_DEFUN([CHECK_MTUNE_NATIVE_SUPPORT],[
 ac_save_CFLAGS=$CFLAGS
 AC_CACHE_CHECK([whether $CC understands -mtune=native], [gf2x_cv_cc_supports_mtune_native],[
  gf2x_cv_cc_supports_mtune_native=no
  CFLAGS="$ac_save_CFLAGS -mtune=native"
  AC_COMPILE_IFELSE(
      [HELLO_WORLD_EXAMPLE()],
      [
      gf2x_cv_cc_supports_mtune_native=yes
      ])
  CFLAGS=$ac_save_CFLAGS
  ])
  if test "$gf2x_cv_cc_supports_mtune_native" = "yes" ;then
    CFLAGS="$CFLAGS -mtune=native"
  fi
])# CHECK_MTUNE_NATIVE_SUPPORT


# It is necessary to make all tests. We do encounter in the wild binutils
# (openbsd binutils 2.15, namely) which are buggy with ssse3, and that
# isn't extremely quickly spotted by the later checks...
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
      CXXFLAGS="$CXXFLAGS $ac_compile_warnings_opt"
      ac_compile_warnings_msg="$ac_compile_warnings_opt for C++"
    fi

  if test -n "$CC"
  then
    if test "$GCC" = "yes"
    then
      ac_compile_warnings_opt='-Wall -W'
    fi
    CFLAGS="$CFLAGS $ac_compile_warnings_opt"
    ac_compile_warnings_msg="$ac_compile_warnings_msg $ac_compile_warnings_opt for C"
  fi
  fi
  AC_MSG_RESULT([$ac_compile_warnings_msg])
  unset ac_compile_warnings_msg
  unset ac_compile_warnings_opt
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

