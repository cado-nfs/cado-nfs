#!/usr/bin/env bash

# This does a quick check that dup2 works without crashing, based on a
# poly file and some sample relations out of las (or dup1).

set -e

: ${WORKDIR?missing}

DL=()

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        BUILD_DIR="$1"
        shift
    elif [ "$1" = "-poly" ] ; then
        shift
        POLY="$1"
        shift
    elif [ "$1" = "-S" ] ; then
        shift
        REFERENCE_SHA1="$1"
        shift
    elif [ "$1" = "-R" ] ; then
        shift
        REFERENCE_REVISION="$1"
        shift
    elif [ "$1" = "-lpbs" ] ; then
        shift
        LPBS="$1"
        shift
    elif [ "$1" = "-rels" ] ; then
        shift
        RELS="$1"
        shift
    elif [ "$1" = "-dl" ] ; then
        DL=(-dl)
        shift
    elif [ "$1" = "-nodl" ] ; then
        DL=()
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

"${BUILD_DIR}/sieve/freerel" "${common[@]}" \
                           -lpbs "$LPBS"
# bail out early if debug_renumber sees an inconsistency.
"${BUILD_DIR}/misc/debug_renumber" "${common[@]}" -check -quiet
"${BUILD_DIR}/filter/dup2" "${common[@]}"       \
                    -nrels $(gzip -dc "$WORK_RELS" | wc -l) "${DL[@]}" "${WORK_RELS}"

if [ "$REFERENCE_SHA1" ] ; then
    SHA1BIN=sha1sum
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then
        echo "Could not find a SHA-1 checksumming binary !" >&2
        exit 1
    fi

    # This is taken from sievetest.sh
    export LC_ALL=C
    export LANG=C
    export LANGUAGE=C

    sort_rels() {
        read -s -r -d '' perl_code <<- 'EOF'
            /^[^#]/ or next;
            chomp($_);
            my ($ab,@sides) = split(":", $_, 3);
            for (@sides) {
                $_ = join(",", sort({ hex($a) <=> hex($b) } split(",",$_)));
            }
            print join(":", ($ab, @sides)), "\n";
            
EOF
        perl -ne "$perl_code" "$@" | sort -n
    }

    if [[ "$WORK_RELS" =~ \.gz$ ]] ; then
        SHA1=`gzip -dc "${WORK_RELS}" | grep "^[^#]" | sort_rels | ${SHA1BIN}` || exit 1
    else
        SHA1=`grep "^[^#]" "${WORK_RELS}" | sort_rels | ${SHA1BIN}` || exit 1
    fi
    SHA1="${SHA1%% *}"
    echo "$0: Got SHA1 of ${SHA1}"
    echo "$0: expected ${REFERENCE_SHA1}"
    if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
      if [ -n "${REFERENCE_REVISION}" ] ; then
        REFMSG=", as created by Git revision ${REFERENCE_REVISION}"
      fi
      if [ "$CADO_DEBUG" ] ; then
          REFMSG=". Files remain in ${WORKDIR}"
      else
          REFMSG=". Set CADO_DEBUG=1 to examine log output"
      fi
      echo "$0: Got SHA1(sort(${WORK_RELS}))=${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}"
      gzip -dc ${WORK_RELS}
      exit 1
    fi
fi
