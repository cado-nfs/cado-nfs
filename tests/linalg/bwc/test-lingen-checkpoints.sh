#!/usr/bin/env bash

set -e

# fetch the parameters that are needed to call guess_mpi_configs
bindir="$1"
for x in "$@" ; do
    case "$x" in
        mpi_magic=*) mpi_magic="${x#mpi_magic=}";;
    esac
done
mpi=
mpirun=()
. "`dirname $0`/guess_mpi_configs.sh"
set_mpi_derived_variables

if ! [ "$WDIR" ] ; then
    echo "Want \$WDIR" >&2
    exit 1
fi

args=("$@")

set -- "${args[@]}" tuning_schedule_filename="$WDIR/ts.txt"

"`dirname $0`"/test-plingen.sh "$@" --tune
cpdir="$WDIR/cp"
mkdir "$cpdir"
"`dirname $0`"/test-plingen.sh "$@" checkpoint_directory="$cpdir"
# osx does not have -r option to xargs.
ls -l "$cpdir"/pi*aux
find "$cpdir" -name 'pi.0.*' | xargs rm -f || :
nfiles=$(ls -rt "$cpdir"/pi*aux | wc -l)
ls -t "$cpdir"/pi*aux | head -n $((2*nfiles/3)) | xargs rm -vf || :
ls -l "$cpdir"/pi*aux
"`dirname $0`"/test-plingen.sh "$@" checkpoint_directory="$cpdir"
