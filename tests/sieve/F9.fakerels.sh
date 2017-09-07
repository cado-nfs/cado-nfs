#!/usr/bin/env bash

FB="$1"
LAS="$2"
SRCDIR="$3"

BINDIR="`dirname "$2"`"
FREEREL=${BINDIR}/freerel
FAKERELS=${BINDIR}/fake_rels

TMPSAMPLE=`mktemp ${TMPDIR-/tmp}/cadotest.fakerel.sample.XXXXXXXXXX`
TMPRENUMBER=`mktemp ${TMPDIR-/tmp}/cadotest.fakerel.renumber.XXXXXXXXXX`

cmd="${LAS} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -lim0 100000 -lim1 100000 -lpb0 23 -lpb1 23 -mfb0 46 -mfb1 46 \
    -I 13 -fb1 ${FB} -sqside 0 \
    -q0 1000000 -q1 1001000 -random-sample 4 -dup"
echo $cmd
$cmd > $TMPSAMPLE

cmd="${FREEREL} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -renumber $TMPRENUMBER -out /dev/null -lpb0 23 -lpb1 23 -pmax 1"
echo $cmd
$cmd

cmd="${FAKERELS} -poly ${SRCDIR}/parameters/polynomials/F9.poly \
    -lpb0 23 -lpb1 23 -q0 1000000 -q1 1001000 -sqside 0 \
    -renumber $TMPRENUMBER -sample $TMPSAMPLE"
echo $cmd
nfake=`$cmd | grep -c "^[^#]"`

# remove leading spaces (for openbsd 5.3)
let nfake=nfake
nfake_exp=111
if [ "$nfake" != "$nfake_exp" ]; then
    echo "Wrong number of fake relations ($nfake, expected $nfake_exp)"
    exit 1
fi

rm $TMPSAMPLE
rm $TMPRENUMBER

exit 0
