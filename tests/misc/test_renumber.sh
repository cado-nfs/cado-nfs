#!/usr/bin/env bash

set -e

: ${WORKDIR?missing}

FREERELS_EXTRA=()
LC=1

CHECK_MULTI=

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
    elif [ "$1" = "-lcideals" ] ; then
        LC=1
        shift
    elif [ "$1" = "-nolcideals" ] ; then
        LC=
        shift
    elif [ "$1" = "-check-multi" ] ; then
        CHECK_MULTI=1
        shift
    elif [ "$1" = "-renumber_format" ] ; then
        FREERELS_EXTRA+=("$1" "$2")
        shift
        shift
    elif [ "$1" = "-renumber_cache_bits" ] ; then
        FREERELS_EXTRA+=("$1" "$2")
        shift
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

if [ "$LC" ] ; then FREERELS_EXTRA+=(-lcideals) ; fi

: ${BUILD_DIR?missing}
: ${POLY?missing}
: ${LPBS?missing}

RENUMBER="${WORKDIR}/renumber.gz"

${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER} \
                           -lpbs "$LPBS" "${FREERELS_EXTRA[@]}"

if [ "$CHECK_MULTI" ] ; then
    SHA1BIN=sha1sum
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then
        echo "Could not find a SHA-1 checksumming binary !" >&2
        exit 1
    fi
    ${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER}.multi \
                               -lpbs "$LPBS" "${FREERELS_EXTRA[@]}"
    hash_reference=`gzip -dc "$RENUMBER" | $SHA1BIN`
    hash_multi=`cat "${RENUMBER}.multi"/*[0-9] | $SHA1BIN`
    if [ "$hash_reference" != "$hash_multi" ] ; then
        echo "multi-I/O and single-end I/O do not give consistent results" >&2
        exit 1
    fi
fi

${BUILD_DIR}/misc/debug_renumber -poly ${POLY} -renumber ${RENUMBER} -check -quiet
