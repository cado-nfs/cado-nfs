#!/usr/bin/env bash

set -e
set -o pipefail

# The goal is to minimize the context, so that silly little changes to
# the ci/ tree don't necessarily trigger a full rebuild of the
# containers.

IMAGE="$1"

if [ "$tmp" ] ; then
    external=1
    tmp=$(mktemp -d $tmp/context.XXXXXXXXXX)
else
    external=
    tmp=$(mktemp -d /tmp/XXXXXXXXXX)
    trap "rm -rf $tmp" EXIT
fi

mkdir $tmp/context

needed_files() {
    cat <<EOF
000-functions.sh
001-environment.sh
005-build-environment.sh
00-prepare-docker.sh
utilities/ncpus.sh
EOF
    find phony-packages -type f
}

(cd "$(dirname $0)" ; needed_files | xargs tar cf $tmp/context.tar.gz)


(cd "$tmp"/context ; tar xf $tmp/context.tar.gz)

(cd "$(dirname $0)" ; ci/00-dockerfile.sh) > "$tmp/context/Dockerfile"

if [ "$GITHUB_OUTPUT" ] ; then
    echo "context=$tmp/context" >> $GITHUB_OUTPUT
fi

if ! [ "$external" ] ; then
    docker build --pull -t $IMAGE --cache-from $IMAGE:latest $tmp/context
fi
