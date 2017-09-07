#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2013
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


# Be aware that ``make check'' from here cannot ensure that the
# upper-level library is up-to-date.

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/fft -I$(top_srcdir)/fft/mpfq -I$(top_builddir)
AM_LDFLAGS = -no-install

check_PROGRAMS=check-addmul
check_addmul_SOURCES=check-addmul.c
check_addmul_LDADD=$(top_builddir)/libgf2x.la

SYNTHESIZED_TESTS =
TESTS = check-addmul
TEST_EXTENSIONS =

EXTRA_DIST=test-tools.h check-and-bench-common.h

############# mul checks

check_PROGRAMS += check-mul
check_mul_SOURCES=check-mul.c test-tools.c
check_mul_LDADD=$(top_builddir)/libgf2x.la
SYNTHESIZED_TESTS += \
        check_1_1_1_819ea234.mul                \
        check_2_2_1_928e7e94.mul                \
        check_3_3_1_14aac2f1.mul                \
        check_4_4_1_fbc8042f.mul                \
        check_5_5_1_049e6e28.mul                \
        check_6_6_1_f7fa48ac.mul                \
        check_7_7_1_39d1a8b9.mul                \
        check_8_8_1_11334230.mul                \
        check_9_9_1_a499e414.mul                \
        check_20_20_1_d40a21e3.mul              \
        check_100_20_1_6dc5ac13.mul             \
        check_100_100_1_5094c182.mul            \
        check_1000_20_1_6ce07fed.mul            \
        check_1000_100_1_08825281.mul           \
        check_1000_1000_1_c28bc2e3.mul

# Do more tests only on 64-bit machines.
if WORDSIZE_IS_64
SYNTHESIZED_TESTS += \
        check_5437_20_1_ad84f617.mul            \
        check_5437_100_1_8eb064ef.mul           \
        check_5437_1000_1_c546a34c.mul          \
        check_5437_5437_1_f09c5936.mul          \
        check_10000_20_1_4d3b9b2d.mul           \
        check_10000_100_1_f2e7cc20.mul          \
        check_10000_1000_1_3b8778b1.mul         \
        check_10000_5437_1_b337c8eb.mul         \
        check_10000_10000_1_2706ca53.mul        \
        check_832040_1346269_1_8513d1fb.mul     \
        check_785672_1271244_1_842e0d19.mul
endif

TEST_EXTENSIONS += .mul
# We'll run the *.mul tests above as:
#       ./do-check-mul.sh -m check_10000_5437_64751d23.mul
MUL_LOG_COMPILER = $(srcdir)/do-check-mul.sh
AM_MUL_LOG_FLAGS = -b $(builddir) -m

############# fft interface checks

general_fftcheck_sources=check.c        \
                      check-and-bench-common.c          \
                      test-tools.c
# This must be exactly the same list of configurations as in fft/Makefile.am
# the flags WITHOUT_CANTOR_TRUNCATION, CANTOR_GM, and CANTOR_GM_TRUNCATE
# are irrelevant here, but we lazily copy the config block from
# fft/Makefile.am
CCOMMON=-DENGINE_CANTOR -DESCAPE_CANTOR_SYMBOLS
C64=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=64
C64nt=$(C64) -DWITHOUT_CANTOR_TRUNCATION
C64gm=$(C64) -DCANTOR_GM
C64gmt=$(C64gm) -DCANTOR_GM_TRUNCATE
C128=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=128
C128nt=$(C128) -DWITHOUT_CANTOR_TRUNCATION
C128gm=$(C128) -DCANTOR_GM
C128gmt=$(C128gm) -DCANTOR_GM_TRUNCATE
TERNARY=-DENGINE_TERNARY
TERNARY_NS=-DENGINE_TERNARY -DARTIFICIAL_NON_SPLIT_VERSION

# The list of tests is auto-generated with perl, and edited afterwards.
# The edits are:
#  - link check_ternaryns with libternary (since this is the same
#  library, really).
#  - add bench_ternary (introduced 7b9a9b45)
# -- begin pattern --
# check_PROGRAMS+=check_ARG0
# check_ARG0_SOURCES=$(general_fftcheck_sources)
# check_ARG0_CFLAGS=$(AM_CFLAGS) ARG1
# check_ARG0_LDADD=$(top_builddir)/fft/libARG0.la
# SYNTHESIZED_TESTS+=   \
#     check_ARG0_1000_1_c28bc2e3.fft            \
#     check_ARG0_1000_2_77db95cc.fft            \
#     check_ARG0_1000_3_6dc06b2f.fft
# if WORDSIZE_IS_64
# SYNTHESIZED_TESTS+=   \
#     check_ARG0_10000_1_2706ca53.fft           \
#     check_ARG0_10000_2_3f32a516.fft           \
#     check_ARG0_10000_3_fd7bba07.fft           \
#     check_ARG0_100000_1_fca70ced.fft          \
#     check_ARG0_100000_2_cd1293a5.fft          \
#     check_ARG0_100000_3_91e7d307.fft
# endif
# -- end pattern --
# -- begin arguments --
#  ternary $(TERNARY)
#  ternaryns $(TERNARY_NS)
# if ENABLE_FFT_INTERFACE
#  cantor64 $(C64)
#  cantor64nt $(C64nt)
#  cantor64gm $(C64gm)
#  cantor64gmt $(C64gmt)
#  cantor128 $(C128)
#  cantor128nt $(C128nt)
#  cantor128gm $(C128gm)
#  cantor128gmt $(C128gmt)
# endif
# -- end arguments --
# -- begin generated code --
check_PROGRAMS+=check_ternary
check_ternary_SOURCES=$(general_fftcheck_sources)
check_ternary_CFLAGS=$(AM_CFLAGS) $(TERNARY)
check_ternary_LDADD=$(top_builddir)/fft/libternary.la
check_PROGRAMS+=bench_ternary
bench_ternary_SOURCES=bench.c check-and-bench-common.c test-tools.c
bench_ternary_CFLAGS=$(AM_CFLAGS) $(TERNARY)
bench_ternary_LDADD=$(top_builddir)/fft/libternary.la
SYNTHESIZED_TESTS+=   \
    check_ternary_1000_1_c28bc2e3.fft            \
    check_ternary_1000_2_77db95cc.fft            \
    check_ternary_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_ternary_10000_1_2706ca53.fft           \
    check_ternary_10000_2_3f32a516.fft           \
    check_ternary_10000_3_fd7bba07.fft           \
    check_ternary_100000_1_fca70ced.fft          \
    check_ternary_100000_2_cd1293a5.fft          \
    check_ternary_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_ternaryns
check_ternaryns_SOURCES=$(general_fftcheck_sources)
check_ternaryns_CFLAGS=$(AM_CFLAGS) $(TERNARY_NS)
check_ternaryns_LDADD=$(top_builddir)/fft/libternary.la
SYNTHESIZED_TESTS+=   \
    check_ternaryns_1000_1_c28bc2e3.fft            \
    check_ternaryns_1000_2_77db95cc.fft            \
    check_ternaryns_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_ternaryns_10000_1_2706ca53.fft           \
    check_ternaryns_10000_2_3f32a516.fft           \
    check_ternaryns_10000_3_fd7bba07.fft           \
    check_ternaryns_100000_1_fca70ced.fft          \
    check_ternaryns_100000_2_cd1293a5.fft          \
    check_ternaryns_100000_3_91e7d307.fft
endif
if ENABLE_FFT_INTERFACE
check_PROGRAMS+=check_cantor64
check_cantor64_SOURCES=$(general_fftcheck_sources)
check_cantor64_CFLAGS=$(AM_CFLAGS) $(C64)
check_cantor64_LDADD=$(top_builddir)/fft/libcantor64.la
SYNTHESIZED_TESTS+=   \
    check_cantor64_1000_1_c28bc2e3.fft            \
    check_cantor64_1000_2_77db95cc.fft            \
    check_cantor64_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor64_10000_1_2706ca53.fft           \
    check_cantor64_10000_2_3f32a516.fft           \
    check_cantor64_10000_3_fd7bba07.fft           \
    check_cantor64_100000_1_fca70ced.fft          \
    check_cantor64_100000_2_cd1293a5.fft          \
    check_cantor64_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor64nt
check_cantor64nt_SOURCES=$(general_fftcheck_sources)
check_cantor64nt_CFLAGS=$(AM_CFLAGS) $(C64nt)
check_cantor64nt_LDADD=$(top_builddir)/fft/libcantor64nt.la
SYNTHESIZED_TESTS+=   \
    check_cantor64nt_1000_1_c28bc2e3.fft            \
    check_cantor64nt_1000_2_77db95cc.fft            \
    check_cantor64nt_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor64nt_10000_1_2706ca53.fft           \
    check_cantor64nt_10000_2_3f32a516.fft           \
    check_cantor64nt_10000_3_fd7bba07.fft           \
    check_cantor64nt_100000_1_fca70ced.fft          \
    check_cantor64nt_100000_2_cd1293a5.fft          \
    check_cantor64nt_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor64gm
check_cantor64gm_SOURCES=$(general_fftcheck_sources)
check_cantor64gm_CFLAGS=$(AM_CFLAGS) $(C64gm)
check_cantor64gm_LDADD=$(top_builddir)/fft/libcantor64gm.la
SYNTHESIZED_TESTS+=   \
    check_cantor64gm_1000_1_c28bc2e3.fft            \
    check_cantor64gm_1000_2_77db95cc.fft            \
    check_cantor64gm_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor64gm_10000_1_2706ca53.fft           \
    check_cantor64gm_10000_2_3f32a516.fft           \
    check_cantor64gm_10000_3_fd7bba07.fft           \
    check_cantor64gm_100000_1_fca70ced.fft          \
    check_cantor64gm_100000_2_cd1293a5.fft          \
    check_cantor64gm_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor64gmt
check_cantor64gmt_SOURCES=$(general_fftcheck_sources)
check_cantor64gmt_CFLAGS=$(AM_CFLAGS) $(C64gmt)
check_cantor64gmt_LDADD=$(top_builddir)/fft/libcantor64gmt.la
SYNTHESIZED_TESTS+=   \
    check_cantor64gmt_1000_1_c28bc2e3.fft            \
    check_cantor64gmt_1000_2_77db95cc.fft            \
    check_cantor64gmt_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor64gmt_10000_1_2706ca53.fft           \
    check_cantor64gmt_10000_2_3f32a516.fft           \
    check_cantor64gmt_10000_3_fd7bba07.fft           \
    check_cantor64gmt_100000_1_fca70ced.fft          \
    check_cantor64gmt_100000_2_cd1293a5.fft          \
    check_cantor64gmt_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor128
check_cantor128_SOURCES=$(general_fftcheck_sources)
check_cantor128_CFLAGS=$(AM_CFLAGS) $(C128)
check_cantor128_LDADD=$(top_builddir)/fft/libcantor128.la
SYNTHESIZED_TESTS+=   \
    check_cantor128_1000_1_c28bc2e3.fft            \
    check_cantor128_1000_2_77db95cc.fft            \
    check_cantor128_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor128_10000_1_2706ca53.fft           \
    check_cantor128_10000_2_3f32a516.fft           \
    check_cantor128_10000_3_fd7bba07.fft           \
    check_cantor128_100000_1_fca70ced.fft          \
    check_cantor128_100000_2_cd1293a5.fft          \
    check_cantor128_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor128nt
check_cantor128nt_SOURCES=$(general_fftcheck_sources)
check_cantor128nt_CFLAGS=$(AM_CFLAGS) $(C128nt)
check_cantor128nt_LDADD=$(top_builddir)/fft/libcantor128nt.la
SYNTHESIZED_TESTS+=   \
    check_cantor128nt_1000_1_c28bc2e3.fft            \
    check_cantor128nt_1000_2_77db95cc.fft            \
    check_cantor128nt_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor128nt_10000_1_2706ca53.fft           \
    check_cantor128nt_10000_2_3f32a516.fft           \
    check_cantor128nt_10000_3_fd7bba07.fft           \
    check_cantor128nt_100000_1_fca70ced.fft          \
    check_cantor128nt_100000_2_cd1293a5.fft          \
    check_cantor128nt_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor128gm
check_cantor128gm_SOURCES=$(general_fftcheck_sources)
check_cantor128gm_CFLAGS=$(AM_CFLAGS) $(C128gm)
check_cantor128gm_LDADD=$(top_builddir)/fft/libcantor128gm.la
SYNTHESIZED_TESTS+=   \
    check_cantor128gm_1000_1_c28bc2e3.fft            \
    check_cantor128gm_1000_2_77db95cc.fft            \
    check_cantor128gm_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor128gm_10000_1_2706ca53.fft           \
    check_cantor128gm_10000_2_3f32a516.fft           \
    check_cantor128gm_10000_3_fd7bba07.fft           \
    check_cantor128gm_100000_1_fca70ced.fft          \
    check_cantor128gm_100000_2_cd1293a5.fft          \
    check_cantor128gm_100000_3_91e7d307.fft
endif
check_PROGRAMS+=check_cantor128gmt
check_cantor128gmt_SOURCES=$(general_fftcheck_sources)
check_cantor128gmt_CFLAGS=$(AM_CFLAGS) $(C128gmt)
check_cantor128gmt_LDADD=$(top_builddir)/fft/libcantor128gmt.la
SYNTHESIZED_TESTS+=   \
    check_cantor128gmt_1000_1_c28bc2e3.fft            \
    check_cantor128gmt_1000_2_77db95cc.fft            \
    check_cantor128gmt_1000_3_6dc06b2f.fft
if WORDSIZE_IS_64
SYNTHESIZED_TESTS+=   \
    check_cantor128gmt_10000_1_2706ca53.fft           \
    check_cantor128gmt_10000_2_3f32a516.fft           \
    check_cantor128gmt_10000_3_fd7bba07.fft           \
    check_cantor128gmt_100000_1_fca70ced.fft          \
    check_cantor128gmt_100000_2_cd1293a5.fft          \
    check_cantor128gmt_100000_3_91e7d307.fft
endif
endif
# -- end generated code --

TEST_EXTENSIONS += .fft
# We'll run the *.fft tests above as:
#       ./do-check-fft.sh -m check_cantor64_100_1_00000000.fft
FFT_LOG_COMPILER = $(srcdir)/do-check-fft.sh
AM_FFT_LOG_FLAGS=-b $(builddir) -m




$(SYNTHESIZED_TESTS) : ; @touch $@
TESTS += $(SYNTHESIZED_TESTS)


EXTRA_DIST+=do-check-mul.sh do-check-fft.sh
