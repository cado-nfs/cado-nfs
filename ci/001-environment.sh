# This file must be sourced.
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
    genepi|calva) ;;
    plafrim|fcatrel|fnancy|catrel-*|miriel*|mistral*|bora*) ;;
    gcc*) ;;
    poire*) ;;
    macintosh*home) ;;
    fedora*|debian*|ubuntu*|centos*|freebsd*|openbsd*|netbsd*) ;;
    pine64) export NCPUS_FAKE=1;;
    *)
    if ! [ "$DOCKER_SCRIPT" ] ; then
        echo "${CSI_RED}Hostname is $HOSTNAME ; that is unexpected${CSI_RESET}" >&2
    fi
    # Make this a non-fatal event now.
    # exit 1
    ;;
esac

# prefer if's to case switches, as it is possible that several such cases
# are matches.

if [[ $CI_BUILD_NAME =~ "coverage tests" ]] ; then
    : ${CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
    : ${CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"}
fi

if [[ $CI_BUILD_NAME =~ "with clang" ]] ; then
    : ${CC=clang}
    : ${CXX=clang++}
fi

if [[ $CI_BUILD_NAME =~ "with icc" ]] ; then
    : ${CC=icc}
    : ${CXX=icpc}
fi

if [[ $CI_BUILD_NAME =~ "expensive checks" ]] ; then
    export CHECKS_EXPENSIVE=1
fi

export CC CXX
export CFLAGS CXXFLAGS
