# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent, because we're read by 00-dockerfile.sh. Note
# though that 00-dockerfile.sh sets HUSH_STDOUT, so that ECHO_E is
# actually a no-op (i.e. ECHO_E is ok to use safely, just don't use plain
# echo)

export CLICOLOR_FORCE=1


# Note that our set of scripts reacts on CI_JOB_NAME, and the logic for
# this is here. CI_JOB_NAME must follow the regexp below.
#
# (note that we **CANNOT** do regexp matching in this script because
# we're /bin/sh, not bash !)
#
# (container for )?(coverage tests on )?((expensive )?checks)?( on (osx|alpine|(debian|fedora|centos|freebsd)[0-9]*) system)?( with ((32-bit )?gcc|clang|icc))?
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

if ! [ "$CI_JOB_NAME" ] && [ "$1" ] ; then
    CI_JOB_NAME="$1"
    $ECHO_E "${CSI_BLUE}Setting CI_JOB_NAME=\"$1\"${CSI_RESET}"
fi

if ! [ "$CI_JOB_NAME" ] ; then
    $ECHO_E "${CSI_RED}This set of scripts really really expect that CI_JOB_NAME is set to something!${CSI_RESET}"
fi

if ! [ "$CI_COMMIT_SHORT_SHA" ] && [ -d .git ] && type -p git > /dev/null 2>&1 ; then
    CI_COMMIT_SHORT_SHA="$(git rev-parse --short HEAD)"
    $ECHO_E "${CSI_BLUE}Setting CI_COMMIT_SHORT_SHA=\"$CI_COMMIT_SHORT_SHA\"${CSI_RESET}"
fi
    
if ! [ "$CI_JOB_ID" ] ; then
    CI_JOB_ID=0
    $ECHO_E "${CSI_BLUE}Setting CI_JOB_ID=\"$CI_JOB_ID\"${CSI_RESET}"
fi
    
case "$CI_JOB_NAME" in
    *"coverage tests"*)
    : ${CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    : ${CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    coverage=1
    ;;
esac
case "$CI_JOB_NAME" in
    *"with gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    gcc=1
    ;;
esac
case "$CI_JOB_NAME" in
    *"with 32-bit gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    : ${CFLAGS="$CFLAGS -m32"}
    : ${CXXFLAGS="$CXXFLAGS -m32"}
    GMP="/usr/local/gmp-6.2.1.abi32"
    export GMP
    gcc32=1
    ;;
esac
case "$CI_JOB_NAME" in
    *"and architecture "*)
    set ${CI_JOB_NAME#*and architecture }
    : ${CFLAGS="$CFLAGS -march=$1"}
    : ${CXXFLAGS="$CXXFLAGS -march=$1"}
    ;;
esac
case "$CI_JOB_NAME" in
    *"with clang"*)
    : ${CC=clang}
    : ${CXX=clang++}
    clang=1
    ;;
esac
case "$CI_JOB_NAME" in
    *"with icc"*)
    : ${CC=icc}
    : ${CXX=icpc}
    icc=1
    ;;
esac
case "$CI_JOB_NAME" in
    *"checks"*) checks=1 ;;
esac
case "$CI_JOB_NAME" in
    *"coverity"*) coverity=1 ;;
esac

case "$CI_JOB_NAME" in
    *"expensive checks"*) export CHECKS_EXPENSIVE=1 ;;
esac

source_tree="$PWD"
export source_tree

case "$CI_JOB_NAME" in
    *"out-of-source"*) out_of_source=1 ;;
esac
case "$CI_JOB_NAME" in
    *"make-dist"*"tarball"*) do_make_dist=1 ;;
esac
case "$CI_JOB_NAME" in
    *"LGPL code"*|*"LGPL tarball"*) remove_gpl_sources=1 ;;
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
