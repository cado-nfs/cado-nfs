#!/usr/bin/env bash

# See https://sympa.inria.fr/sympa/arc/cado-nfs/2020-09/msg00020.html

: ${CADO_NFS_BINARY_DIR:?missing}
: ${CADO_NFS_SOURCE_DIR:?missing}

FREEREL=${CADO_NFS_BINARY_DIR}/sieve/freerel
FAKERELS=${CADO_NFS_BINARY_DIR}/sieve/fake_rels

TMPRENUMBER="${wdir:?missing}/cadotest.fakerel.renumber"

cmd="${FREEREL} -poly ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/c60.poly \
    -renumber $TMPRENUMBER -lpb0 18 -lpb1 19 -pmax 1"
echo $cmd
$cmd

# The following was segfaulting
cmd="${FAKERELS} -poly ${CADO_NFS_SOURCE_DIR}/parameters/polynomials/c60.poly \
    -lpb0 18 -lpb1 19 -q0 61961 -q1 62961 -sqside 1 -shrink-factor 1 \
    -renumber $TMPRENUMBER -sample ${CADO_NFS_SOURCE_DIR}/tests/sieve/c60.sample"
echo $cmd
$cmd
