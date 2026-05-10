#!/usr/bin/env bash

usage() {
    cat <<EOF
Usage: $0 <N> <ncols_dense> <nrows_base> [<ccsize0> <ccsize1> ...]

This script creates a set of dummy relations for testing filter/purge.

A dense part consisting of ncols_dense columns is created, with indices
between 0 and N-1. nrows_base such "plain rows" are created. Then we also
create rows, modeled in the same way, that form connected components of
sizes ccsize0, ccsize1, etc. The resulting set of relations is then shuffled.
EOF
}

oneline() {
    N=$1
    nj=$2
    js=($(for i in `seq 1 $nj`; do echo $((RANDOM%N)) ; done))
    ab=$(printf "%x,%x" $RANDOM $RANDOM)
    sep=':'
    for c in "${js[@]}"; do
        ab="${ab}${sep}`printf %x $c`"
        sep=,
    done
    echo $ab
}

onecc() {
    N=$1
    nj=$2
    base=$3
    ccsize=$4
    let i=$((base))
    for k in `seq 1 $((ccsize-1))` ; do
        let j=$((i+1))
        echo `oneline $N $nj`,`printf %x,%x $i $j`
        i=$j
    done
    echo `oneline $N $nj`,`printf %x,%x $i $base`
}

buildmatrix() {
    N=$1; shift
    nj=$1; shift
    nr=$1; shift
    if ! [ "$N$nj$nr" ] ; then usage; exit 1 ; fi
    for i in `seq 1 $nr` ; do oneline $N $nj ; done
    b=$N
    while [ $# -gt 0 ] ; do
        onecc $N $nj $b $1
        let b+=$1
        shift
    done
}

if type -p shuf > /dev/null 2>&1 ; then
    fshuf() { shuf ; }
elif ! sort --version 2>&1 | grep -q -e -R ; then
    fshuf() { sort -R ; }
else
    echo "no shuffle function found"
    exit 1
fi
buildmatrix "$@" | fshuf
