#!/usr/bin/env bash

set -e

if ! [ "$WDIR" ] ; then
    echo "Want \$WDIR" >&2
    exit 1
fi

echo "Using mpi=$mpi; $mpirun; $mpi_extra_args"

args=("$@")

set -- "${args[@]}" tuning_schedule_filename="$WDIR/ts.txt"

"`dirname $0`"/test-plingen.sh "$@" --tune
cpdir="$WDIR/cp"
mkdir "$cpdir"
set -- "$@" tuning_quiet=1 checkpoint_directory="$cpdir"
"`dirname $0`"/test-plingen.sh "$@"
# osx does not have -r option to xargs.
ls -l "$cpdir"/pi*aux
find "$cpdir" -name 'pi.0.*' | xargs rm -f || :
nfiles=$(ls -rt "$cpdir"/pi*aux | wc -l)
ls -t "$cpdir"/pi*aux | head -n $((2*nfiles/3)) | xargs rm -vf || :
ls -l "$cpdir"/pi*aux
"`dirname $0`"/test-plingen.sh "$@"
