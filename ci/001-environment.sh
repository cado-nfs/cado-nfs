# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent, because we're read by 00-dockerfile.sh. Note
# though that 00-dockerfile.sh sets HUSH_STDOUT, so that ECHO_E is
# actually a no-op (i.e. ECHO_E or major_message are ok to use safely,
# just don't use plain echo)

export CLICOLOR_FORCE=1

# Note that our set of scripts reacts on BUILD_NAME, and the logic for
# this is here. BUILD_NAME must follow the regexp below.
#
# (note that we **CANNOT** do regexp matching in this script because
# we're /bin/sh, not bash !)
#
# (container for )?(coverage tests on )?((expensive )?checks)?( on (osx|alpine|(debian|fedora|centos|freebsd)[0-9]*) system)?( with ((32[- ]bit )?gcc|clang|icc))?
#
# with one exception which is "merge coverage tests"
#
# a downside is that in the gitlab page, this makes many pipeline steps
# with similar names. We could change to abbreviated names (e.g.
# "coverage tests on " would be V, "container for " would be L, "checks"
# and "expensive checks" would be C and XC, and so on. But it would be
# really cryptic. to have, e.g. LC/debian10+gcc ; wouldn't it ?

if type -p hostname > /dev/null 2>&1 ; then
    HOSTNAME=$(hostname)
else
    HOSTNAME="[[placeholder]]"
fi
case "$HOSTNAME" in
    # some of our very slow machines have so little ram that clearly, we
    # must not tax them too much.
    genepi|calva|pine64) export NCPUS_FAKE=1;;
    *) : ;;
esac

if ! [ "$BUILD_NAME" ] && [ "$1" ] ; then
    BUILD_NAME="$1"
    $ECHO_E "${CSI_BLUE}Setting BUILD_NAME=\"$1\"${CSI_RESET}"
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
major_message "BUILD_NAME=$BUILD_NAME"

#### set COMMIT_SHORT_SHA to CI_COMMIT_SHORT_SHA or GITHUB_SHA
if [ "$CI_COMMIT_SHORT_SHA" ] ; then
    COMMIT_SHORT_SHA="$CI_COMMIT_SHORT_SHA"
elif [ "$GITHUB_SHA" ] ; then
    COMMIT_SHORT_SHA="$GITHUB_SHA"
elif [ -d .git ] && type -p git > /dev/null 2>&1 ; then
    COMMIT_SHORT_SHA="$(git rev-parse --short HEAD)"
    $ECHO_E "${CSI_BLUE}Setting COMMIT_SHORT_SHA=\"$COMMIT_SHORT_SHA\"${CSI_RESET}"
fi
major_message "COMMIT_SHORT_SHA=$COMMIT_SHORT_SHA"

#### set JOB_ID to either CI_JOB_ID or GITHUB_RUN_ID
if [ "$CI_JOB_ID" ] ; then
    JOB_ID="$CI_JOB_ID"
elif [ "$GITHUB_RUN_ID" ] ; then
    JOB_ID="$GITHUB_RUN_ID"
else
    JOB_ID=0
    $ECHO_E "${CSI_BLUE}Setting JOB_ID=\"$JOB_ID\"${CSI_RESET}"
fi
major_message "JOB_ID=$JOB_ID"

### set REPOSITORY to either $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME or GITHUB_REPOSITORY

if [ "$CI_PROJECT_NAMESPACE" ] && [ "$CI_PROJECT_NAME" ] ; then
    REPOSITORY="$CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME"
elif [ "$GITHUB_REPOSITORY" ] ; then
    REPOSITORY="$GITHUB_REPOSITORY"
else
    # no default
    REPOSITORY=
fi
major_message "REPOSITORY=$REPOSITORY"
    
case "$BUILD_NAME" in
    *"coverage tests"*)
    : ${CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    : ${CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    coverage=1
    major_message coverage reports enabled
    ;;
esac
case "$BUILD_NAME" in
    *"with gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    gcc=1
    major_message compiler: using gcc
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
        major_message 32-bit build
        ;;
    *) major_message default abi build
        ;;
esac
case "$BUILD_NAME" in
    *"shared libs"*)
    ENABLE_SHARED=1
    shared_libs=1
    major_message shared libraries enabled
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
    major_message compiler: using clang-$clang
    ;;
esac
case "$BUILD_NAME" in
    *"with icc"*)
    : ${CC=icc}
    : ${CXX=icpc}
    icc=1
    major_message compiler: using icc
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
        major_message producing static analysis data for coverity
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
        major_message using cmake directly
    ;;
esac
case "$BUILD_NAME" in
    *"expensive checks"*)
        export CHECKS_EXPENSIVE=1
        major message doing expensive checks
    ;;
esac

if [ -x /opt/homebrew/bin/brew ] ; then
    eval `/opt/homebrew/bin/brew shellenv`
fi

MAKE=make
if type -p gmake > /dev/null 2>&1 ; then
    MAKE=gmake
fi

export CC CXX
export CFLAGS CXXFLAGS
export MAKE
export ENABLE_SHARED
