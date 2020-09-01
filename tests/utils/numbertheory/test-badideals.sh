#!/usr/bin/env bash

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 --binary <binary> [--wdir <directory>] <test sample name>" >&2
    for x in "$@" ; do echo "$x" ; done
    exit 1
}

tests=()

while [ $# -gt 0 ] ; do
    if [ "$1" = "--wdir" ] ; then
        shift
        wdir=$1
        shift
    elif [ "$1" = "--binary" ] ; then
        shift
        binary=$1
        shift
    else
        tests+=("$1")
        shift
    fi
done

if ! [ "$binary" ] ; then
    usage
fi

: ${TMPDIR=/tmp}
: ${wdir=$(mktemp -d  $TMPDIR/cado-nfs.badideals.XXXXXXXX)}

check_there() {
    u="$1"
    cd "$wdir"
    "$binary" $(cat "$u.input") > "$u.stdout"
    for e in stdout badideals badidealinfo ; do
        if ! [ -f "$u.$e" ] ; then continue; fi
        if ! [ -f "$u.expect_$e" ] ; then usage "$t.$e: missing file" ; fi
        grep -v '^#' "$u.$e" > "$u.$e.expunged" || :
        grep -v '^#' "$u.expect_$e" > "$u.expect_$e.expunged" || :
        if ! diff -q "$u.expect_$e.expunged" "$u.$e.expunged" ; then
            echo "Difference for test $t"
            diff "$u.expect_$e.expunged" "$u.$e.expunged"
            exit 1
        fi
    done
}

for t in "${tests[@]}" ; do
    for e in input ; do
        if [ -f "$t.$e" ] ; then
            cp -f "$t.$e" "$wdir"
        else
            usage "$t.$e: missing file"
        fi
    done
    u=$(basename "$t")
    if [ -f "$t.poly" ] ; then
        cp -f "$t.poly" "$wdir/"
    fi
    for e in stdout badideals badidealinfo ; do
        if [ -f "$t.$e" ] ; then
            cp -f "$t.$e" "$wdir/$u.expect_$e"
        fi
    done
    check_there "$u"
done


