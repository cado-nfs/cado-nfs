#!/usr/bin/env bash

# See https://sympa.inria.fr/sympa/arc/cado-nfs/2020-09/msg00020.html

BINDIR="$1"
SRCDIR="$2"

FREEREL=${BINDIR}/sieve/freerel
FAKERELS=${BINDIR}/sieve/fake_rels

TMPRENUMBER="${WORKDIR:?missing}/cadotest.fakerel.renumber"

cmd="${FREEREL} -poly ${SRCDIR}/parameters/polynomials/c60.poly \
    -renumber $TMPRENUMBER -lpb0 18 -lpb1 19 -pmax 1"
echo $cmd
$cmd

# The following was segfaulting
cmd="${FAKERELS} -poly ${SRCDIR}/parameters/polynomials/c60.poly \
    -lpb0 18 -lpb1 19 -q0 61961 -q1 62961 -sqside 1 -shrink-factor 1 \
    -renumber $TMPRENUMBER -sample ${SRCDIR}/tests/sieve/c60.sample"
echo $cmd
$cmd
