#!/usr/bin/env bash

set -e
set -x
# Create a fake sequence

# Note that if we arrive here, we are 64-bit only, since the GFP backends
# for bwc are explicitly disabled on i386 (for now -- most probably
# forever too).

wordsize=64

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

bindir="$1"
shift

dotest() {
    REFERENCE_SHA1="$1"
    shift

    : ${TMPDIR:=/tmp}
    TMPDIR=`mktemp -d $TMPDIR/lingen-test.XXXXXXXXXX`

    m="$1"; shift
    n="$1"; shift
    length="$1"; shift
    seed="$1"; shift

    # this is just copying the argument array here. Could do more if we
    # considered having influential parameters here.
    args=()
    mt_args=()
    for x in "$@" ; do
        case "$x" in
            *) args+=("$x");;
        esac
    done

    F="$TMPDIR/base"
    "`dirname $0`"/perlrandom.pl $((m*n*length/8)) $seed > $F
    G="$TMPDIR/seq.bin"
    cat $F $F $F > $G
    rm -f $F

    $bindir/linalg/bwc/lingen m=$m n=$n prime=2 --lingen-input-file $G --lingen-output-file $G.gen "${args[@]}"
    [ -f "$G.gen" ]
    SHA1=$($SHA1BIN < $G.gen)
    SHA1="${SHA1%% *}"

    if [ "$REFERENCE_SHA1" ] ; then
        reftab=(`echo $REFERENCE_SHA1 | tr ',' ' '`)
        for ref in "${reftab[@]}" ; do
            if [ "${SHA1}" = "$ref" ] ; then
                match="$ref"
            fi
        done
        if ! [ "$match" ] ; then
            echo "$0: Got SHA1 of ${SHA1}; matches none of ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}" >&2
            exit 1
        fi
    else
        echo "========= $SHA1 ========"
    fi
    rm -f $G.gen

    rm -rf "$TMPDIR"
}

dotest "$@"

