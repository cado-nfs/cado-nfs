#!/usr/bin/env bash

set -e

: ${WORKDIR?missing}

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

: ${BUILD_DIR?missing}
: ${POLY?missing}
: ${LPBS?missing}

RENUMBER="${WORKDIR}/renumber.gz"

${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER} \
                           -lpbs "$LPBS" 

${BUILD_DIR}/misc/debug_renumber -poly ${POLY} -renumber ${RENUMBER} -check -quiet
