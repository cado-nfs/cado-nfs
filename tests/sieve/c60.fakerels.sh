#!/usr/bin/env bash

# See cado-nfs-discuss/2020-September/001310.html

BINDIR="$1"
SRCDIR="$2"

FREEREL=${BINDIR}/sieve/freerel
FAKERELS=${BINDIR}/sieve/fake_rels

TMPRENUMBER="${WORKDIR:?missing}/cadotest.fakerel.renumber"

cmd="${FREEREL} -poly ${SRCDIR}/parameters/polynomials/c60.poly \
    -renumber $TMPRENUMBER -out /dev/null -lpb0 18 -lpb1 19 -pmax 1"
echo $cmd
$cmd

# The following was segfaulting
cmd="${FAKERELS} -poly ${SRCDIR}/parameters/polynomials/c60.poly \
    -lpb0 18 -lpb1 19 -q0 61961 -q1 62961 -sqside 1 -shrink-factor 1 \
    -renumber $TMPRENUMBER -sample ${SRCDIR}/tests/sieve/c60.sample"
echo $cmd
$cmd
