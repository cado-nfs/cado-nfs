#!/usr/bin/env bash

FB="$1"
LAS="$2"
DUPSUP="$3"
SRCDIR="$4"

TMPFILE=`mktemp ${TMPDIR-/tmp}/cadotest.dupsup.XXXXXXXXXX`

cmd="${LAS} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -lim0 2300000 -lim1 1200000 -lpb0 28 -lpb1 28 -mfb0 56 -mfb1 56 \
    -lambda0 2.1 -lambda1 2.1 \
    -I 12 -fb1 ${FB} -q0 250000000 -nq 1 \
    -ncurves0 10 -ncurves1 10"
echo "$cmd -dup"
ndup_ref=`$cmd -dup | grep "DUPE " | wc -l`
# remove leading spaces (for openbsd 5.3)
let ndup_ref=ndup_ref
if [ "$ndup_ref" != "12" ]; then
    echo "Wrong number of duplicates in las -dup"
    echo "Expected 12, got $ndup_ref"
    exit 1
fi

cmd="${LAS} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -lim0 2300000 -lim1 1200000 -lpb0 28 -lpb1 28 -mfb0 56 -mfb1 56 \
    -lambda0 2.1 -lambda1 2.1 \
    -I 12 -fb1 ${FB} -q0 250000000 -nq 1 \
    -ncurves0 10 -ncurves1 10 -out $TMPFILE"
echo $cmd
$cmd

cmd="${DUPSUP} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -lim0 2300000 -lim1 1200000 -lpb0 28 -lpb1 28 -mfb0 56 -mfb1 56 \
    -I 12 -ncurves0 10 -ncurves1 10 -skew 1.0 \
    -lambda0 2.1 -lambda1 2.1 \
    $TMPFILE"
echo $cmd
ndup=`$cmd | grep "DUPE " | wc -l`
# remove leading spaces (for openbsd 5.3)
let ndup=ndup
if [ "$ndup" != "12" ]; then
    echo "Wrong number of duplicates in dupsup"
    echo "Expected 12, got $ndup"
    exit 1
fi

rm $TMPFILE

exit 0
