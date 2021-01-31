# This file must be sourced.
#
# We're /bin/sh, not bash.
#
# We **must** be silent !

export CLICOLOR_FORCE=1

if type -p hostname > /dev/null 2>&1 ; then
    HOSTNAME=$(hostname)
else
    HOSTNAME="[[placeholder]]"
fi
case "$HOSTNAME" in
    # docker runners (gitlab ones) are apparently called "runner-*"
    runner-*) ;;
    docker-script-*) ;;
    cado-*) ;;
    raclette|fondue|tartiflette|berthoud) ;;
    plafrim|fcatrel|fnancy|catrel-*|miriel*|mistral*|bora*) ;;
    gcc*) ;;
    poire*) ;;
    macintosh*home) ;;
    fedora*|debian*|ubuntu*|centos*|freebsd*|openbsd*|netbsd*) ;;
    # some of our very slow machines have so little ram that clearly, we
    # must not tax them too much.
    genepi|calva|pine64) export NCPUS_FAKE=1;;
esac

case "$CI_BUILD_NAME" in
    *"coverage tests"*)
    : ${CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    : ${CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    coverage=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"with gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    gcc=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"with 32-bit gcc"*)
    : ${CC=gcc}
    : ${CXX=g++}
    : ${CFLAGS="$CFLAGS -m32"}
    : ${CXXFLAGS="$CXXFLAGS -m32"}
    GMP="/usr/local/gmp-6.1.2.abi32"
    export GMP
    gcc32=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"with clang"*)
    : ${CC=clang}
    : ${CXX=clang++}
    clang=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"with icc"*)
    : ${CC=icc}
    : ${CXX=icpc}
    icc=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"checks"*)
        checks=1
    ;;
esac
case "$CI_BUILD_NAME" in
    *"expensive checks"*)
        export CHECKS_EXPENSIVE=1
    ;;
esac

if [ -x /opt/homebrew/bin/brew ] ; then
    eval `/opt/homebrew/bin/brew shellenv`
fi

export CC CXX
export CFLAGS CXXFLAGS
