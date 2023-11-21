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
    source `dirname "$0"`/realpath.sh
    readlink_f() { realpath "$@" ; }
else
    readlink_f() { readlink -f "$@" ; }
fi


up_path=$(readlink_f `dirname "$0"`/..)
path_inside_tree=$(readlink_f "$up_path/scripts/`basename $0`")
path_myself=$(readlink_f "$0")

if [ "$path_inside_tree" != "$path_myself" ] ; then
    echo "Error: `basename $0` must be inside the cado-nfs source tree"
    exit 1
fi

case "$path_myself" in
    *\ *) echo "You are trying to compile cado-nfs from a directory with spaces in its full path. This is not supported. (Mostly because the GNU autotools don't support it either.)" >&2
        exit 1
        ;;
esac

pwdP="pwd -P"
if ! $pwdP >/dev/null 2>&1 ; then
    pwdP=pwd
fi
#
# When "make" is run from the subdirectory $CADO/linalg/bwc, we set:
#
# called_from to be the absolute path of $CADO/linalg/bwc
# absolute_path_of_source to be the absolute path of $CADO
# relative_path_of_cwd to be linalg/bwc
#
called_from=$(readlink_f "`pwd`")
absolute_path_of_source=$(readlink_f "$up_path")
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
if [ "$LOCAL_SH" ] ; then
    if [ -f "$LOCAL_SH" ] ; then
        . "$LOCAL_SH"
    else
        echo "Weird, LOCAL_SH=$LOCAL_SH points to nonexisting file" >&2
    fi
elif [ -f "${up_path}/local.sh" ] ; then
    . "${up_path}/local.sh"
fi

if [ "$MPI" ] && [ "$MPI" != 0 ] && ! [[ "$build_tree" =~ \.mpi ]] ; then
    build_tree="$build_tree".mpi
fi

if [ "$force_build_tree" ] ; then
        build_tree="$force_build_tree"
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

passed_env=(
CC
CXX
CFLAGS
CXXFLAGS
FLAGS_SIZE
LDFLAGS
PREFIX
CADO_DIST_ARCHIVE_NAME
MPI
GF2X_CONFIGURE_EXTRA
GMP
GMP_INCDIR
GMP_LIBDIR
MPIR
MPIR_INCDIR
MPIR_LIBDIR
PTHREADS
HWLOC
HWLOC_INCDIR
HWLOC_LIBDIR
GMPECM
GMPECM_INCDIR
GMPECM_LIBDIR
JEVENTS
JEVENTS_INCDIR
JEVENTS_LIBDIR
NUMA
NUMA_INCDIR
NUMA_LIBDIR
GF2X_CONFIGURE_EXTRA_FLAGS
CMAKE_DUMP_VARIABLES
ENABLE_SHARED
NO_PYTHON_CHECK
NO_SSE
NO_INLINE_ASSEMBLY
NO_GAS_ASSEMBLY
CHECKS_EXPENSIVE
BWC_GF2_MATMUL_BACKENDS
BWC_GFP_MATMUL_BACKENDS
BWC_GF2_ARITHMETIC_BACKENDS
BWC_GFP_ARITHMETIC_BACKENDS
BWC_EXTRA_BACKENDS
TIMEOUT_SCALE
)

for e in "${passed_env[@]}" ; do export $e ; done

########################################################################
# make show
if [ "$1" = "--show" ] ; then
    env_me=(build_tree up_path called_from absolute_path_of_source relative_path_of_cwd)
    for e in "${env_me[@]}" ; do
        echo "${e}=\"${!e}\""
    done
    for e in "${passed_env[@]}" ; do
        echo "${e}=\"${!e}\""
    done
    exit 0
fi
