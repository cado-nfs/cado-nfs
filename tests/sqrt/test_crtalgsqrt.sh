#!/usr/bin/env bash

DEP=$1
RATDEP=$2
POLY=$3

: ${CADO_NFS_BINARY_DIR:?missing}

set -- $DEP $RATDEP $POLY
# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    set -- "$@" mpi="$mpi"
fi

run=("${mpirun[@]}" "${CADO_NFS_BINARY_DIR}/sqrt/crtalgsqrt" "$@")
echo "${run[@]}"
"${run[@]}" 2>&1 | egrep '(44371162641954939938390944368|40462797324737803355716975649)'
