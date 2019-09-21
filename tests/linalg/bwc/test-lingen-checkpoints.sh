#!/bin/bash

set -e

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
find "$cpdir" -name 'pi.[012].*' | xargs rm -f || :
nfiles=$(ls -rt "$cpdir/pi*" | wc -l)
ls -rt "$cpdir/pi*" | head -n $((nfiles/3)) | xargs rm -f || :
"`dirname $0`"/test-plingen.sh "$@" checkpoint_directory="$cpdir"
