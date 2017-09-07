#!/usr/bin/env bash

# For debug, uncomment:
# set -x

########################################################################
# This script is responsible of handing over the build process, in a
# proper out of source build directory. It takes care of calling cmake
# first if needed, and then cd's into the proper sub-directory of the
# build tree, and runs make there. The intent is that this script is
# called from a Makefile within the source tree.
# In particular, the following tasks are done,
#  - check if the calling path is correct?
#  - if exists, parse the file ${up_path}/local.sh
#  - check if cmake is installed, if not install it.
#  - "cmake" to generate Makefile
#  - "make"
#
# Type "make ?" for more options.
########################################################################

: ${MAKE=make}
export MAKE

if echo "${MAKEFLAGS}" | grep -q "jobserver-fds=0," ; then
    echo "# You are calling the top-level cado makefile with file descriptor 0 closed.">&2
    echo "# This is unsupported (at least for a parallel build), because in that">&2
    echo "# case GNU Make opens and uses a pipe on file descriptor 0, and we">&2
    echo "# suspect that cmake closes it right away, causing the compilation to">&2
    echo "# fail.">&2
    echo "#">&2
    echo "# Simple fix: make -j \$number_of_cpus < /dev/null">&2
    echo "#">&2
    exit 42
fi

args=("$@")
if [ "$1" = "show" ] ; then
    if [ "$#" != 1 ] ; then
        echo "Argument 'show' must be alone on the command line of $0" >&2
        fi
    set -- --show
else
    set --
fi
source "$(dirname $0)/build_environment.sh"
if [ "$1" ] ; then
    # we've done our deeds, finish.
    exit 0
fi
set -e
set "${args[@]}"

########################################################################
# "make ?" or "make help" (when the Makefile does not exist)
info=""
function make_usage {
    echo "-------------------------------------------------------------- "
    echo "[Options] (see $0 for details)"
    echo "-------------------------------------------------------------- "
    echo " $info \"make ?\"       -- this help"
    echo " $info \"make show\"    -- only show env variables"
    echo " $info \"make cmake\"   -- only run cmake to generate Makefile"
    echo " $info \"make\"         -- run cmake first and then make"
    echo " $info \"make tidy\"    -- delete folder $build_tree (dangerous)"
    echo " $info  Any other options will be passed to the actual make followed."
    echo "-------------------------------------------------------------- "
    exit 0
}
if [ "$1" == "?" ] ; then
    make_usage
fi
if [ "$1" == "help" ] && [ ! -f "$build_tree/Makefile" ] ; then
    make_usage
fi

########################################################################
# make tidy (warn, this delete the whole build folder)
wn="[warning]"
if [ "$1" == "tidy" ] ; then
    echo "$wn this deletes the whole folder $build_tree. Do you want to continue? (n/Y)"
    if [ -e "`tty`" ] ; then
        read TIDY_BUILD
    else
        echo "$wn no input terminal, assuming no"
        TIDY_BUILD=n
    fi
    if [ "$TIDY_BUILD" == "Y" ]; then
        echo "$wn wiping out $build_tree"
        rm -rf "$build_tree"
    else
        echo "$wn no action and quit now"
    fi
    exit 0
fi

########################################################################
# Make sure we have cmake, by the way !
:  ${cmake_path:="`which cmake 2>/dev/null`"}
cmake_companion_install_location="$absolute_path_of_source/cmake-installed"
if [ "$?" != "0" ] || ! [ -x "$cmake_path" ] ; then
    echo "CMake not found" >&2
    cmake_path=
# Recall that (some versions of) bash do not want quoting for regex patterns.
elif [[ "`"$cmake_path" --version`" =~ ^cmake\ version\ [012] ]] && ! ( [[ "`"$cmake_path" --version`" =~ ^cmake\ version\ 2.(9|8.11|8.12) ]]) ; then
    echo "CMake found, but not with version 2.8.11 or newer" >&2
    cmake_path=
fi
if ! [ "$cmake_path" ] ; then
    cmake_path="$cmake_companion_install_location/bin/cmake"
    if [ -x "$cmake_path" ] ; then
        echo "Using custom cmake in $cmake_companion_install_location" >&2
    else
        echo "I am about to download and compile a compatible version of Cmake."
        echo "Do you want to continue ? (y/n)"
        if [ -e "`tty`" ] ; then
            read INSTALL_CMAKE
        else
            echo "No input terminal, assuming yes"
            INSTALL_CMAKE=y
        fi
        if [ ! "$INSTALL_CMAKE" = "y" ]; then
            echo "Please install a compatible version of Cmake."
            exit 1
        fi
        echo "Need to get cmake first -- this takes long !"
        cd "$up_path"
        if ! scripts/install-cmake.sh "$cmake_companion_install_location" ; then
            echo "cmake install Failed, sorry" >&2
            exit 1
        fi
        cd "$called_from"
    fi
fi


########################################################################
# handle "make clean"
if [ "$1" == "clean" ] && [ ! -f "$build_tree/Makefile" ] ; then
    echo "There is no $build_tree/Makefile. Nothing to clean."
    exit 0
fi

########################################################################
# call cmake (if Makefile does not exist)
if [ "$1" = "cmake" ] || [ ! -f "$build_tree/Makefile" ] ; then
    mkdir -p "$build_tree"
    absolute_path_of_build_tree="`cd "$build_tree" ; $pwdP`"
    if [ ! "x$CMAKE_GENERATOR" == "x" ] ; then
      CMAKE_GENERATOR_OPT="-G$CMAKE_GENERATOR"
    else
      unset CMAKE_GENERATOR_OPT
    fi
    (cd "$absolute_path_of_build_tree" ; "$cmake_path" "$CMAKE_GENERATOR_OPT" $CMAKE_EXTRA_ARGS "$absolute_path_of_source")
fi

if [ "$1" = "cmake" ] ; then
    exit 0
fi

########################################################################
# Now cd into the target directory, and build everything required.
# Note that it's useful to kill MAKELEVEL, or otherwise we'll get scores
# and scores of ``Entering directory'' messages (sure, there's the
# --no-print-directory option -- but it's not the right cure here).
# env | grep -i make
unset MAKELEVEL
absolute_path_of_build_tree="`cd "$build_tree" ; $pwdP`"
(cd "$absolute_path_of_build_tree$relative_path_of_cwd" ; ${MAKE} "$@")
