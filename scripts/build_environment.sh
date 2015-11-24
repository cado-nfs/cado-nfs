#!/usr/bin/env bash

# This file is sourced from call_cmake.sh, and sets all the shell
# variables which are useful to cmake. Whe call_cmake.sh sources this
# files, it sets $# to zero.
#
# cado-nfs.py also *runs* this file, setting $# to 1, and $1 to "--show".
# The output is then parsed by python in order to guess variables such as
# the build tree.

# This is called only from within the source tree.

########################################################################

# find the absolute path of the source tree, and deduce the default
# location of the build directory (by default,
# $cado_source_tree/build/`hostname`, but customizing it is easy).


# This is readlink -f on many unices. Alas, not on mac.
if [ "`uname -s`" = Darwin ] ; then
    # use code from https://github.com/mkropat/sh-realpath, MIT-licensed.
    source "`dirname $0`/realpath.sh"
    readlink_f() { realpath "$@" ; }
else
    readlink_f() { readlink -f "$@" ; }
fi


up_path="$(readlink_f `dirname $0`/..)"
if [ "$(readlink_f $up_path/scripts/`basename $0`)" != "$(readlink_f $0)" ] ; then
    echo "Error: `basename $0` must be inside the cado-nfs source tree"
    exit 1
fi
pwdP="pwd -P"
if ! $pwdP >/dev/null 2>&1 ; then
    pwdP=pwd
fi
called_from="`pwd`"
absolute_path_of_source="$(readlink_f $up_path)"
relative_path_of_cwd="${called_from##$absolute_path_of_source}"

########################################################################
# Set some default variables which can be overridden from local.sh

# The default behaviour
# build_tree="/tmp/cado-build-`id -un`"
: ${build_tree:="${up_path}/build/`hostname`"}

# Joe gets confused by " in prev line; this line fixes it

# By default, we also avoid /usr/local ; of course, it may be overridden.
# Note that cmake does not really have phony targets, so one must avoid
# basenames which correspond to targets !
: ${PREFIX:="$absolute_path_of_source/installed"}

# XXX XXX XXX LOOK LOOK LOOK: here you've got an entry point for customizing.
# The source directory may contain a hint script with several useful
# preferences. Its job is to put environment variables. Quite notably,
# $build_tree is amongst these.
if [ -f "${up_path}/local.sh" ] ; then
    . "${up_path}/local.sh"
fi

# If no CFLAGS have been set yet, set something sensible: get optimization by
# default, as well as asserts.  If you want to disable this, use either
# local.sh or the environment to set an environment variable CFLAGS to be
# something non-empty, like CFLAGS=-g or CFLAGS=-O0 ; Simply having CFLAGS=
# <nothing> won't do, because bash makes no distinction between null and unset
# here.
: ${CFLAGS:=-O2}
: ${CXXFLAGS:=-O2}
: ${PTHREADS:=1}

########################################################################
# Arrange so that relevant stuff is passed to cmake -- the other end of
# the magic is in CMakeLists.txt. The two lists must agree.
# (no, it's not as simple. As long as the cmake checks care about
# *environment variables*, we are here at the right place for setting
# them. However, we might also be interested in having cmake export test
# results to scripts. This is done by cmake substitutions, but the
# corresponding names need not match the ones below).

export CC
export CXX
export CFLAGS
export CXXFLAGS
export FLAGS_SIZE
export LDFLAGS
export PREFIX
export CADO_DIST_ARCHIVE_NAME
export MPI
export GF2X_CONFIGURE_EXTRA
export GMP
export GMP_INCDIR
export GMP_LIBDIR
export MPIR
export MPIR_INCDIR
export MPIR_LIBDIR
export PTHREADS
export CURL
export CURL_INCDIR
export CURL_LIBDIR
export HWLOC
export HWLOC_INCDIR
export HWLOC_LIBDIR
export NUMA
export NUMA_INCDIR
export NUMA_LIBDIR
export GF2X_CONFIGURE_EXTRA_FLAGS
export CMAKE_DUMP_VARIABLES
export ENABLE_SHARED
export NO_PYTHON_CHECK
export NO_SSE
export NO_INLINE_ASSEMBLY
export NO_GAS_ASSEMBLY
export CHECKS_EXPENSIVE
export BWC_GF2_MATMUL_BACKENDS
export BWC_GFP_MATMUL_BACKENDS
export BWC_GF2_ARITHMETIC_BACKENDS
export BWC_GFP_ARITHMETIC_BACKENDS
export BWC_EXTRA_BACKENDS

########################################################################
# make show
if [ "$1" = "--show" ] ; then
    echo "build_tree=\"$build_tree\""
    echo "up_path=\"$up_path\""
    echo "called_from=\"$called_from\""
    echo "absolute_path_of_source=\"$absolute_path_of_source\""
    echo "relative_path_of_cwd=\"$relative_path_of_cwd\""
    echo "CC=\"$CC\""
    echo "CXX=\"$CXX\""
    echo "CFLAGS=\"$CFLAGS\""
    echo "CXXFLAGS=\"$CXXFLAGS\""
    echo "FLAGS_SIZE=\"$FLAGS_SIZE\""
    echo "LDFLAGS=\"$CFLAGS\""
    echo "PREFIX=\"$PREFIX\""
    echo "CADO_DIST_ARCHIVE_NAME=\"$CADO_DIST_ARCHIVE_NAME\""
    echo "MPI=\"$MPI\""
    echo "GF2X_CONFIGURE_EXTRA=\"$GF2X_CONFIGURE_EXTRA\""
    echo "PTHREADS=\"$PTHREADS\""
    echo "CMAKE_GENERATOR=\"$CMAKE_GENERATOR\""
    echo "ENABLE_SHARED=\"$ENABLE_SHARED\""
    echo "CHECKS_EXPENSIVE=\"$CHECKS_EXPENSIVE\""
    echo "BWC_GF2_MATMUL_BACKENDS=\"$BWC_GF2_MATMUL_BACKENDS\""
    echo "BWC_GFP_MATMUL_BACKENDS=\"$BWC_GFP_MATMUL_BACKENDS\""
    echo "BWC_GF2_ARITHMETIC_BACKENDS=\"$BWC_GF2_ARITHMETIC_BACKENDS\""
    echo "BWC_GFP_ARITHMETIC_BACKENDS=\"$BWC_GFP_ARITHMETIC_BACKENDS\""
    echo "BWC_EXTRA_BACKENDS=\"$BWC_EXTRA_BACKENDS\""
    exit 0
fi
