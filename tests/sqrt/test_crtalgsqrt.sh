#!/usr/bin/env bash

BINARY=$1
DEP=$2
RATDEP=$3
POLY=$4

set -- $DEP $RATDEP $POLY
# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    set -- "$@" mpi="$mpi"
fi

run=("${mpirun[@]}" $BINARY "$@")
echo "${run[@]}"
"${run[@]}" 2>&1 | egrep '(44371162641954939938390944368|40462797324737803355716975649)'
