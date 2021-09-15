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
#include <stdlib.h>
int
main ()
{
  exit(0);
}
EOF
gf2x_compile="$1 conftest.c"
cc_for_build_works=no
if AC_TRY_EVAL(gf2x_compile); then
  if (./a.out || ./b.out || ./a.exe || ./a_out.exe || ./conftest) >&AS_MESSAGE_LOG_FD 2>&1; then
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
#include <stdlib.h>
int
main ()
{
  exit (0);
}
EOF
for i in .exe ,ff8 ""; do
  gf2x_compile="$CC_FOR_BUILD conftest.c -o conftest$i"
  if AC_TRY_EVAL(gf2x_compile); then
    if (./conftest) 2>&AS_MESSAGE_LOG_FD; then
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


AC_DEFUN([GF2X_CHECK_VISIBILITY_HIDDEN],
[
AC_CACHE_CHECK([for __attribute__((visibility("hidden")))],
    gf2x_cv_hidden_visibility_attribute, [
    cat > conftest.c <<EOF
int __attribute__ ((visibility ("hidden"))) foo (void) { return 1; }
int __attribute__ ((visibility ("default"))) bar (void) { return 1; }
int baz (void) { return 1; }
EOF
    gf2x_cv_hidden_visibility_attribute=no
    if AC_TRY_COMMAND(${CC-cc} -fvisibility=hidden -Werror -S conftest.c -o conftest.s 1>&AS_MESSAGE_LOG_FD);
    then
        if (grep '\.hidden.*foo' conftest.s && grep '\.hidden.*baz' conftest.s && ! grep '\.hidden.*bar' conftest.s) >/dev/null;
        then
            gf2x_cv_hidden_visibility_attribute=yes
        fi
    fi
    rm -f conftest.*
    ])
if test $gf2x_cv_hidden_visibility_attribute = yes;
then
    AC_DEFINE(HAVE_HIDDEN_VISIBILITY_ATTRIBUTE, 1,
          [Define if __attribute__((visibility("hidden"))) is supported.])
fi
AM_CONDITIONAL([HAVE_HIDDEN_VISIBILITY_ATTRIBUTE],[test "x$gf2x_cv_hidden_visibility_attribute" = xyes])
])

AC_DEFUN([GF2X_TEST_CLOCK_NOT_CONSTANT_CODE],[AC_LANG_SOURCE([
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
int main()
{
    clock_t c0 = clock(), c1;
    long count = 0;
    for( ; count < 1000000 ; count++) {
        if ((c1 = clock()) != c0)
            return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
])])
    
AC_DEFUN([GF2X_TEST_CLOCK_NOT_CONSTANT],
[
    AC_CACHE_CHECK([that clock() returns non-constant values],
    gf2x_cv_clock_is_not_constant, [
        AC_RUN_IFELSE([GF2X_TEST_CLOCK_NOT_CONSTANT_CODE()],[
            gf2x_cv_clock_is_not_constant=yes
        ],[
            gf2x_cv_clock_is_not_constant=no
        ],[
            AC_MSG_NOTICE([check skipped because of cross-compiling])
            gf2x_cv_clock_is_not_constant=dontknow
        ])
    ])
    if test $gf2x_cv_clock_is_not_constant = yes;
    then
        AC_DEFINE(HAVE_NONCONSTANT_CLOCK, 1,
              [Define if clock() returns non-constant values.])
    fi
    AM_CONDITIONAL([HAVE_NONCONSTANT_CLOCK],[test "x$gf2x_cv_clock_is_not_constant" = xyes])
])
