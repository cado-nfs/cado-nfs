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
#  - check if cmake is installed
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
    if [ "$TIDY_BUILD" == "Y" ] ; then
        echo "$wn wiping out $build_tree"
        rm -rf "$build_tree"
    else
        echo "$wn no action and quit now"
    fi
    exit 0
fi

########################################################################
# Make sure we have cmake, by the way !

set +e

if [ "$cmake_path" ] ; then
    if ! [ -x "$cmake_path" ] ; then
        echo "cmake_path=$cmake_path points to non-existing (or non-executable) file" >&2
        exit 1
    fi
else
    # just try with the PATH.
    : ${cmake_path="`type -p cmake 2>/dev/null`"}
fi

cmake_version=

if [ "$cmake_path" ] ; then
cmake_version=$("$cmake_path" --version || :)
fi

if ! [ "$cmake_version" ] ; then
    echo "CMake not found" >&2
    cmake_path=
# Recall that (some versions of) bash do not want quoting for regex patterns.
elif [[ "$cmake_version" =~ ^cmake\ version\ [012] ]] ; then
    echo "CMake found, but not with version 3.5 or newer" >&2
    cmake_path=
elif [[ "$cmake_version" =~ ^cmake\ version\ 3\.[01234]\. ]] ; then
    echo "CMake found, but not with version 3.5 or newer" >&2
    cmake_path=
fi

set -e

cmake_companion_install_location="$absolute_path_of_source/cmake-installed"

if ! [ "$cmake_path" ] ; then
    if [ "$CI_JOB_NAME" ] ; then
        echo "Refusing to auto-install cmake in CI builds" >&2
        exit 1
    fi
    cmake_path="$cmake_companion_install_location/bin/cmake"
    if [ -x "$cmake_path" ] ; then
        echo "Using custom cmake in $cmake_companion_install_location" >&2
    else
        echo "No cmake binary was found on your system."
        echo
        echo "Most probably, you want to rely on your system distribution to"
        echo "provide cmake in some way. Luckily, the package is usually called"
        echo "cmake. Version **at least** 3.5 is necessary."
        echo
        echo "You may try to have cado-nfs download and install some version of"
        echo "cmake for you. THIS IS A PRIORI A VERY BAD IDEA, and we advise"
        echo "you to prefer the system-level option above (because you will"
        echo "most likely get a more recent cmake version)."
        echo
        echo "If, despite this word of warning, you want to auto-install cmake,"
        echo "you can do so with scripts/install-cmake.sh \"$cmake_companion_install_location\""
        echo "(note that compiling cmake in itself takes a bit of time)"
        echo
        exit 1

        echo "Do you want to continue ? (y/n)"
        if [ -e "`tty`" ] ; then
            read INSTALL_CMAKE
        else
            echo "No input terminal, assuming yes"
            INSTALL_CMAKE=y
        fi
        if [ ! "$INSTALL_CMAKE" = "y" ] ; then
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

if [ "$cmake_path" ] ; then
    # we want the absolute path to cmake_path, since cmake gets called in
    # a directory which is not the cwd, and the current logic will fail
    # if cmake_path is a relative path (or if cmake_path=cmake and a
    # relative path is in the PATH)
    absolute_cmake_path="$(type -p "$cmake_path" || which "$cmake_path" 2>/dev/null)"
    if [ "$absolute_cmake_path" ] ; then
        absolute_cmake_path="$(realpath $absolute_cmake_path)"
    fi
    if [ $? != 0 ] ; then
        echo "Cannot find absolute path of $cmake_path" >&2
    else
        cmake_path="$absolute_cmake_path"
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
    cmake_gen=()
    if [ "$CMAKE_GENERATOR" ] ; then
        cmake_gen+=("-G$CMAKE_GENERATOR")
    fi
    if [ "$(bash -c 'echo ${CC}')" ] ; then
        cmake_overrides+=(-DCMAKE_C_COMPILER="$CC")
    fi
    if [ "$(bash -c 'echo ${CXX}')" ] ; then
        cmake_overrides+=(-DCMAKE_CXX_COMPILER="$CXX")
    fi
    if [ "$(bash -c 'echo ${MAKE}')" ] ; then
        cmake_overrides+=(-DCMAKE_MAKE_PROGRAM="$MAKE")
    fi
    (cd "$absolute_path_of_build_tree" ; "$cmake_path" "${cmake_gen[@]}" $CMAKE_EXTRA_ARGS "${cmake_overrides[@]}" "$absolute_path_of_source")
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

callit_args=("$@")
callit() {
(cd "$absolute_path_of_build_tree$relative_path_of_cwd" ; ${MAKE} "${callit_args[@]}")
}
if [ "$1" = "check" ] && ! [ "$ctest_filter" = "no" ] ; then
    set -o pipefail
    # the ctest_filter groks -nc, -q, -v
    # fd 3 is another stdout that is guaranteed to escape the ctest
    # filter. We may want to use it.
    (callit | "$absolute_path_of_source/scripts/filter-ctest.pl" $ctest_filter) 3>&1
else
    callit
fi
