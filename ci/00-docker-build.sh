#!/bin/sh

# This wrapper is run on the hosts that build the containers, i.e.  from
# the "docker" docker image, using /bin/sh. We do like to have
# "pipefail", although it's not obvious that it will work with /bin/sh.
# As a matter of fact, currently it does work, so we're happy.

set -e
set -o pipefail

# The goal is to minimize the context, so that silly little changes to
# the ci/ tree don't necessarily trigger a full rebuild of the
# containers.

IMAGE="$1"

tmp=$(mktemp -d /tmp/XXXXXXXXXX)
trap "rm -rf $tmp" EXIT

mkdir $tmp/context

needed_files() {
    cat <<EOF
000-functions.sh
001-environment.sh
00-prepare-docker.sh
utilities/ncpus.sh
EOF
}

(cd "$(dirname $0)" ; needed_files | xargs tar cf $tmp/context.tar.gz)


(cd "$tmp"/context ; tar xf $tmp/context.tar.gz)

ci/00-dockerfile.sh > "$tmp/context/Dockerfile"

docker build --pull -t $IMAGE --cache-from $IMAGE:latest $tmp/context
