changequote(<<, >>)dnl

dnl Usage for the FEATURE_CHECK code generation macro:
dnl       SSE2    # SSE2_EXAMPLE
dnl       sse-2   # human-readable feature name
dnl       -msse2  # for compiler flag
dnl       sse2    # for shell variable gf2x_cv_cc_supports_sse2
dnl       msse2   # for shell variable gf2x_cv_cpp_needs_msse2
dnl       parent  # name of parent feature. checks gf2x_cv_cc_supports_parent
dnl
dnl


m4_define(<<FEATURE_CHECK>>,<<
  # Hello world, inserting code for >><<$1>><<
  # First check the compiler support
    AC_CACHE_CHECK([whether $CC can compile and run >><<$2>><< code], [gf2x_cv_cc_supports_>><<$4>><<],[
      gf2x_cv_cc_supports_>><<$4>><<=no
      if >>m4_if(<<$6>>,<<>>,<<>>,<<test "x${gf2x_cv_cc_supports_>><<$6>><<}" = xno ; then
        echo $ECHO_N "skipped>><<<<,>>>><< "
      elif >>)<<
      test "x${enable_>><<$4>><<}" = xno ; then
        echo $ECHO_N "explicitly disabled, "
      else
        AC_RUN_IFELSE([>><<$1>><<_EXAMPLE()],[
          gf2x_cv_cc_supports_>><<$4>><<=yes
        ],[
          case "x$CFLAGS" in #((
           *-march=*)
             echo $ECHO_N "(not testing >><<$3>><< because -march is already present), "
             gf2x_cv_cc_supports_>><<$4>><<="no";;
           *)
             ac_save_CFLAGS=$CFLAGS
             CFLAGS="$ac_save_CFLAGS >><<$3>><<"
             AC_RUN_IFELSE([>><<$1>><<_EXAMPLE()],[
              gf2x_cv_cc_supports_>><<$4>><<="yes, with >><<$3>><<"
             ],[
              CFLAGS="$ac_save_CFLAGS"
              gf2x_cv_cc_supports_>><<$4>><<=no
             ])
             ;;
          esac
        ],[
          echo $ECHO_N "cross-compiling, "
          if test "x${enable_>><<$4>><<}" = xyes ; then
            echo $ECHO_N "explicitly enabled, "
            gf2x_cv_cc_supports_>><<$4>><<=yes
          fi
        ])
      fi
    ])
  # Post-interpret our findings.
    CFLAGS=$ac_save_CFLAGS
    if test "$gf2x_cv_cc_supports_>><<$4>><<" = "yes, with >><<$3>><<" ;then
      CFLAGS="$CFLAGS >><<$3>><<"
    fi
  # Publish this as a config.h macro
    case "$gf2x_cv_cc_supports_>><<$4>><<" in #((
      yes*)
        AC_DEFINE([GF2X_HAVE_>><<$1>><<_SUPPORT],[1],[Define if >><<$2>><< code as present in the source tree is supported by the compiler])
        ;;
      *) : ;;
    esac
  # Next, see if we need to do something with the preprocessor flags.
    ac_save_CPPFLAGS=$CPPFLAGS
    if test "$gf2x_cv_cc_supports_>><<$4>><<" = "yes, with >><<$3>><<" ;then
      # Tweaking CFLAGS is often not enough.
      AC_CACHE_CHECK([whether >><<$3>><< is also needed by the preprocessor],
        [gf2x_cv_cpp_needs_>><<$5>><<_flag],[
          CPPFLAGS="$ac_save_CPPFLAGS"
          AC_PREPROC_IFELSE([>><<$1>><<_EXAMPLE()],[
            gf2x_cv_cpp_needs_>><<$5>><<_flag=no
          ],[
            CPPFLAGS="$ac_save_CPPFLAGS >><<$3>><<"
            AC_PREPROC_IFELSE([>><<$1>><<_EXAMPLE()],[
              gf2x_cv_cpp_needs_>><<$5>><<_flag=yes
            ],[
              AC_MSG_ERROR([Sorry, the preprocessor can't parse >><<$2>><<!])
            ])
        ])
      ])
    fi
  # Post-interpret our findings.
    CPPFLAGS=$ac_save_CPPFLAGS
    if test "$gf2x_cv_cpp_needs_>><<$5>><<_flag" = "yes" ; then
      CPPFLAGS="$CPPFLAGS >><<$3>><<"
    fi
>>)
changequote([, ])dnl

AC_DEFUN([CHECK_SSE2_SUPPORT],
 [FEATURE_CHECK(SSE2,sse-2,-msse2,sse2,msse2)])
AC_DEFUN([CHECK_SSE3_SUPPORT],
 [FEATURE_CHECK(SSE3,sse-3,-msse3,sse3,msse3,sse2)])
AC_DEFUN([CHECK_SSSE3_SUPPORT],
 [FEATURE_CHECK(SSSE3,ssse-3,-mssse3,ssse3,mssse3,sse3)])
AC_DEFUN([CHECK_SSE41_SUPPORT],
 [FEATURE_CHECK(SSE41,sse-4.1,-msse4.1,sse41,msse41,ssse3)])
AC_DEFUN([CHECK_PCLMUL_SUPPORT],
 [FEATURE_CHECK(PCLMUL,pclmul,-mpclmul,pclmul,mpclmul,sse41)])
# Only checked by apps/
AC_DEFUN([CHECK_BMI2_SUPPORT],
 [FEATURE_CHECK(BMI2,bmi2,-mbmi2,bmi2,mbmi2,sse41)])

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


AC_DEFUN([BMI2_EXAMPLE],[AC_LANG_SOURCE([
#include <x86intrin.h> /* for _pdep_u64 */
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char * argv[[]]) {
    volatile uint64_t a = UINT64_C(0x3456789012391873);
    volatile uint64_t b = UINT64_C(0xe92183172937191e);
    a = _pdep_u64 (a, UINT64_C(0x5555555555555555));
    b = _pdep_u64 (b, UINT64_C(0xaaaaaaaaaaaaaaaa));
    volatile uint64_t c = a ^ b;
    return c == UINT64_C(0x9860f6b03c217ad) ? EXIT_SUCCESS : EXIT_FAILURE;
}
])])

# vim: fdm=indent sw=1 sta et:
