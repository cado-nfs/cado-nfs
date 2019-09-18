#!/bin/bash

set -e

if ! [ "$WDIR" ] ; then
    echo "Want \$WDIR" >&2
    exit 1
fi

args="$@"

"`dirname $0`"/test-plingen.sh "$@" --tune
cpdir="$WDIR/cp"
mkdir "$cpdir"
"`dirname $0`"/test-plingen.sh "$@" checkpoint_directory="$cpdir"
find "$cpdir" -name 'pi.[012].*' | xargs -r rm -f
nfiles=$(ls -rt "$cpdir/pi*" | wc -l)
ls -rt "$cpdir/pi*" | head -n $((nfiles/3)) | xargs -r rm -f
"`dirname $0`"/test-plingen.sh "$@" checkpoint_directory="$cpdir"
