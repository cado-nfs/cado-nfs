# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent, because we're read by 00-dockerfile.sh. Note
# though that 00-dockerfile.sh sets HUSH_STDOUT, so that ECHO_E is
# actually a no-op (i.e. ECHO_E or major_message are ok to use safely,
# just don't use plain echo)

if [ "$DISPLAY_CONFIG" ] ; then
    display_config() { major_message "$@" ; }
else
    display_config() { : ; }
fi

if ! [ "$BUILD_NAME" ] && [ "$1" ] ; then
    BUILD_NAME="$1"
    display_config "Setting BUILD_NAME=\"$1\""
fi

#### Set BUILD_NAME to either CI_BUILD_NAME or GITHUB_JOB
if [ "$BUILD_NAME" ] ; then
    :
elif [ "$CI_BUILD_NAME" ] ; then
    BUILD_NAME="$CI_BUILD_NAME"
elif [ "$GITHUB_JOB" ] ; then
    BUILD_NAME="`echo $GITHUB_JOB | tr - ' '`"
else
    # this triggers a failure down the line
    $ECHO_E "${CSI_RED}This set of scripts really really expect that BUILD_NAME is set to something!${CSI_RESET}"
fi
display_config "BUILD_NAME=$BUILD_NAME"

case "$BUILD_NAME" in
    *"coverage tests"*)
    : ${CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    : ${CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    coverage=1
    display_config coverage reports enabled
    ;;
esac
case "$BUILD_NAME" in
    *"with gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    gcc=1
    display_config compiler: using gcc
    ;;
esac
case "$BUILD_NAME" in
    *"with 32"[-\ ]"bit gcc"*)
        : ${CC=gcc}
        : ${CXX=g++}
        : ${CFLAGS="$CFLAGS -m32"}
        : ${CXXFLAGS="$CXXFLAGS -m32"}
        GMP="/usr/local/gmp-6.2.1.abi32"
        export GMP
        gcc32=1
        display_config 32-bit build
        ;;
    *) display_config default abi build
        ;;
esac
case "$BUILD_NAME" in
    *"shared libs"*)
    ENABLE_SHARED=1
    shared_libs=1
    display_config shared libraries enabled
    ;;
esac
case "$BUILD_NAME" in
    *"with clang"*)
    : ${CC=clang}
    : ${CXX=clang++}
    clang=1
    # We want to recognize "clangNN" or "clangdev" as monikers for
    # specific versions of clang.
    case "$BUILD_NAME" in
        *"with clangdev"*) clang=dev;;
        *"with clang12"*) clang=12;;
        *"with clang13"*) clang=13;;
        *"with clang14"*) clang=14;;
        *"with clang15"*) clang=15;;
        *"with clang16"*) clang=16;;
        *"with clang17"*) clang=17;;
    esac
    display_config compiler: using clang-$clang
    ;;
esac
case "$BUILD_NAME" in
    *"with icc"*)
    : ${CC=icc}
    : ${CXX=icpc}
    icc=1
    display_config compiler: using icc
    ;;
esac
case "$BUILD_NAME" in
    *"checks"*)
        checks=1
    ;;
esac
case "$BUILD_NAME" in
    *"coverity"*)
        coverity=1
        display_config producing static analysis data for coverity
    ;;
esac
case "$BUILD_NAME" in
    *"using cmake directly"*)
        using_cmake_directly=1
        # use build_tree in this case, which matches the variable that
        # call_cmake.sh uses, by the way.
        source_tree="$PWD"
        export source_tree
        if [ "$BASH_VERSION" ] ; then
            if ! [ "$build_tree" ] ; then
                build_tree="/tmp/$BUILD_NAME"
                # spaces in dir names don't work, mostly because of libtool
                # (look at gf2x/fft/libgf2x-fft.la)
                # This substitution is bash-only, but this should be fine to 
                # have in a conditional that non-bash skips over
                build_tree="${build_tree// /_}"
                export build_tree
            fi
            if ! [ -d "$build_tree" ] ; then
                mkdir "$build_tree"
            fi
        else
            # just a safeguard
            build_tree=/no/build_tree/set/because/we/require/bash/for/that
        fi
        display_config using cmake directly
    ;;
esac
case "$BUILD_NAME" in
    *"expensive checks"*|*"checks expensive"*)
        export CHECKS_EXPENSIVE=1
        display_config doing expensive checks
    ;;
esac

