#!/usr/bin/env bash

set -e

: ${wdir:?missing}

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        BUILD_DIR="$1"
        shift
    elif [ "$1" = "-poly" ] ; then
        shift
        POLY="$1"
        shift
    elif [ "$1" = "-lpbs" ] ; then
        shift
        LPBS="$1"
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

if ! [ "$BUILD_DIR" ] ; then
    BUILD_DIR=${PROJECT_BINARY_DIR?missing}
fi
: ${POLY:?missing}
: ${LPBS:?missing}

RENUMBER="${wdir}/renumber.gz"

${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER} \
                           -lpbs "$LPBS" 

${BUILD_DIR}/misc/debug_renumber -poly ${POLY} -renumber ${RENUMBER} -check -quiet
