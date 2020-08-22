#!/usr/bin/env bash

# This does a quick check that dup2 works without crashing, based on a
# poly file and some sample relations out of las (or dup1).

set -e

: ${WORKDIR?missing}

LC=(-lcideals)

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
    elif [ "$1" = "-rels" ] ; then
        shift
        RELS="$1"
        shift
    elif [ "$1" = "-lcideals" ] ; then
        LC=(-lcideals)
        shift
    elif [ "$1" = "-nolcideals" ] ; then
        LC=()
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

: ${BUILD_DIR?missing}
: ${POLY?missing}
: ${LPBS?missing}
: ${RELS?missing}

RENUMBER="${WORKDIR}/renumber.gz"
WORK_RELS="${WORKDIR}/rels.gz"

cp "$RELS" "$WORK_RELS"

common=(-poly "$POLY" -renumber "${RENUMBER}")

"${BUILD_DIR}/sieve/freerel" "${common[@]}" -out /dev/null \
                           -lpbs "$LPBS" "${LC[@]}"
# bail out early if debug_renumber sees an inconsistency.
"${BUILD_DIR}/misc/debug_renumber" "${common[@]}" -check -quiet
"${BUILD_DIR}/filter/dup2" "${common[@]}"       \
                    -nrels $(zcat "$WORK_RELS" | wc -l) "${WORK_RELS}"
