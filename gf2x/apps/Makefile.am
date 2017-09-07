# put here the location of the NTL installation (for example /usr/local)
NTL=

# add -mbmi2 to use the pdep_u64 instruction
AM_CPPFLAGS = -I$(top_srcdir) -T$(top_builddir) -fopenmp -std=c++11 -I$(NTL)/include

check_PROGRAMS=factor
factor_SOURCES=factor.cpp halfgcd.cpp
noinst_HEADERS=halfgcd.hpp
factor_LDFLAGS=-fopenmp
factor_LDADD=../libgf2x.la $(NTL)/lib/libntl.a -lgmp

LOG_COMPILER=$(SHELL) 2>/dev/null
TESTS=check-factor-000.sh check-factor-001.sh	\
      check-factor-002.sh check-factor-003.sh

EXTRA_DIST=$(TESTS)
