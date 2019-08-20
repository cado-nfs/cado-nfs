#!/usr/bin/env bash

set -e
if [ "$CADO_DEBUG" ] ; then set -x ; fi
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

while [ $# -gt 0 ] ; do
    if [[ "$1" =~ ^(seed|m|n|sequence_length|expect_sha1_F)=[0-9a-f]+$ ]] ; then
        eval "$1"
        shift
        continue
    elif [[ "$1" =~ ^wdir=(.+)$ ]] ; then
        wdir="${BASH_REMATCH[1]}"
        if ! [ -d "$wdir" ] ; then
            echo "wdir $wdir does not exist" >&2
            exit 1
        fi
        shift
        continue
    elif [[ "$1" =~ ^bindir=(.+)$ ]] ; then
        bindir="${BASH_REMATCH[1]}"
        if ! [ -d "$bindir" ] ; then
            echo "bindir $bindir does not exist" >&2
            exit 1
        fi
        shift
        continue
    elif [ "$1" = -- ]  ; then
        shift
        break
    else
        echo "argument $1 not understood" >&2
        exit 1
    fi
done

tail_args=("$@")

for v in m n sequence_length seed wdir bindir expect_sha1_F ; do
    if ! [ "${!v}" ] ; then
        echo "\$$v must be provided" >&2
        exit 1
    fi
done

if ! [ -d "$bindir" ] ; then
    echo "bindir $bindir does not exist" >&2
    exit 1
fi

if ! [ -d "$wdir" ] ; then
    echo "wdir $wdir does not exist" >&2
    exit 1
fi

TMPDIR="$wdir"
REFERENCE_SHA1="$expect_sha1_F"
length="$sequence_length"

dotest() {
    F="$TMPDIR/base"
    "`dirname $0`"/perlrandom.pl $((m*n*length/8)) $seed > $F
    G="$TMPDIR/seq.bin"
    cat $F $F $F > $G
    rm -f $F

    $bindir/linalg/bwc/lingen m=$m n=$n prime=2 --lingen-input-file $G --lingen-output-file $G.gen "${tail_args[@]}"
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
}

dotest "$@"
