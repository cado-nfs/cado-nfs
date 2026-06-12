#!/usr/bin/env bash

# This test is meant to trigger an almost certain failure, but with some
# randomness. It works as a companion to the other torture test that is
# frozen.
#
: "${wdir?missing}"

bindir="${PROJECT_BINARY_DIR?missing}/filter"

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

if type -p shuf > /dev/null 2>&1 ; then
    fshuf() { shuf ; }
elif ! sort --version 2>&1 | grep -q -e -R ; then
    fshuf() { sort -R ; }
else
    echo "no shuffle function found"
    exit 1
fi

$(dirname "$0")/create_purge_test_file.py 192 3 320 32 32 > $wdir/relations_base

# every once in a while, it does occur that all tests pass, even with a
# somewhat buggy purge. No big deal.
ntests=64
if [ "$VALGRIND" ] ; then
    ntests=16
fi
for i in `seq 1 $ntests` ; do
    fshuf < $wdir/relations_base > $wdir/relations.$i
    if ! "$bindir/purge" -keep 0 -col-min-index 192 $wdir/relations.$i > $wdir/out-err.$i 2>&1 ; then
        cat $wdir/out-err.$i
        echo "randomization #$i failed" >&2
        exit 1
    fi
    echo "randomization #$i ok"
done
