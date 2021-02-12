#!/bin/bash
#
# This script is a debugging helper that can be used to mimick what is
# done in the build jobs triggered by .gitlab-ci.yml
#
# On a machine with docker installed, this script can be invoked with a
# command line argument that copies the name of a given test in
# .gitlab-ci.yml ; note that the goal is just to give you a shell in the
# given environment, not to run the tests!
#
#
# Examples:
#
#  ./ci/debug.sh "checks with gcc"
#  ./ci/debug.sh "checks with clang"
#  ./ci/debug.sh "coverage tests with gcc"
#  ./ci/debug.sh "checks on fedora31 system with clang"
#  ./ci/debug.sh "fedora31 system with clang"
# (the latter example is just a shorthand of the preceding one).
#
# The example below is special, because it needs a **LOT** of disk space
# (about 20G).
#
#  ./ci/debug.sh "checks with icc"

# This leaves you in a container, within directory /host, which is
# directly mapped to the directory from which you're calling this
# script. Your edits are directly reflected in the directory on the host,
# and you can even do git commits. Beyond that, the script is pretty
# frugal and doesn't claim to do much.
#
# The (pristine) disk image of your container takes the name that you
# provide in argument, with "debug_" as a prefix.
#
# You can sudo from within the container. However, the modifications that
# you do to the disk space of the container itself are lost when your
# shell exits (but remember that /host is mapped to your original file
# system anyway).
#
# While the container still runs, if you want to open a second connection
# to it, you might find the following commands useful:
#       docker container list
#       docker  exec -ti $CONTAINER_ID bash -i
# (where $CONTAINER_ID is obtained from the first command)
#
#
# Caveat: absolute symlinks, or symlinks to outside the filesystem
# hierarchy under the path from which this script is called, cannot work.
# For this reason, the scripts force the build tree to a temporary
# location in the container's /tmp directory
#
#
# Because the disk image of your container takes a predictible name, a
# subsequent call to the same script should bring you back to exactly the
# same system !
#
# The downside of this is that the images need to be purged every so
# often. A command for this might be:
#
# docker image prune -f --filter "until=`date -Is --date='1 week ago'`"
#
#
# NOTE: to debug the freebsd tests, this wrapper script cannot be used
# since it has not been adapted to this case. The following command line
# is potentially a good start (gives a bash shell on a fresh tree in sync
# with current git HEAD).
#    DOCKER_SCRIPT=1 CI_BUILD_NAME='checks on freebsd13 with gcc' ci/50-libvirt-wrap-tests.sh freebsd:13.0
# but there's no funny volume mounting and so on. Getting to the shell
# prompt takes roughly two minutes in this case.

set -e
export DOCKER_SCRIPT=1
export CI_BUILD_NAME="$1"
tmp=$(mktemp -d /tmp/XXXXXXXXXXXXXX)
trap "rm -rf $tmp" EXIT
cat > $tmp/prepare.sh <<EOF
# This installs packages just as our CI jobs do. However, because of the
# DOCKER_SCRIPT environment variable, some bonus packages are installed
# as well (vim, gdb).
/host/ci/00-prepare-docker.sh
. /host/ci/000-functions.sh
. /host/ci/001-environment.sh
. /host/ci/999-debug.sh
EOF
ci/00-dockerfile.sh > $tmp/Dockerfile
imagename="$1"
shift
: ${imagename=docker-image-$RANDOM}
imagename="${imagename// /_}"
imagename="${imagename//:/_}"
imagename="${imagename//-/_}"
imagename="debug_$imagename"
docker build -t "$imagename" -f $tmp/Dockerfile ci
echo "# NOTE: docker image is $imagename"
echo "# NOTE: this image contains a few extra debug tools"
DARGS=()
if ! [ "$NO_REMOVE" ] ; then
    DARGS+=(--rm)
else
    echo "# NOTE: the container will remain up on script exit"
fi
docker run "${DARGS[@]}" -ti --hostname docker-script-$RANDOM --volume $PWD:/host "$imagename" /host/ci/999-debug.sh "$@"
