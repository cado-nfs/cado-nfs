#!/bin/bash
# This script is a debugging helper that can be used to mimick what is
# done in the build jobs triggered by .gitlab-ci.yml
#
# On a machine with docker installed, this script can be invoked with the
# following command lines:
#
# Example with the "gcc" image that is used by some of the .gitlab-ci.yml jobs:
#     docker run --rm -ti --hostname docker-script-$RANDOM --volume $PWD:/host gcc bash /host/.docker-script.sh
#
# Example with another base image, here "fedora". The script below has
# provision for fedora as well.
#     docker run --rm -ti --hostname docker-script-$RANDOM --volume $PWD:/host fedora bash /host/.docker-script.sh
#
# Third example, on an image with intel icc / icpc installed (you need
# about 20G of disk space for this)
#     docker run --rm -ti --hostname docker-script-$RANDOM --volume $PWD:/host intel/oneapi-hpckit bash /host/.docker-script.sh
#
# This leaves you in a container, within directory /host, which is
# directly mapped to the directory from which you're calling this
# script. Your edits are directly reflected in the directory on the host,
# and you can even do git commits. Beyond that, the script is pretty
# frugal and doesn't claim to do much.
#
# All traces of the container are removed when the script exit (as per --rm)
#
# While the container still runs, if you want to open a second connection
# to it, you might find the following commands useful:
#       docker container list
#       docker  exec -ti $CONTAINER_ID bash -i
# (where $CONTAINER_ID is obtained from the first command)


# 
# Caveat: absolute symlinks, or symlinks to outside the filesystem
# hierarchy under the path from which this script is called, cannot work.
# For this reason, the scripts force the build tree to a temporary
# location in the container's /tmp directory
#
# Note that on purpose, this script reinstalls the needed dependencies.
# This is because the gitlab-ci.yml jobs do so as well, and we want to be
# able to be in sync with what they do. Of course, if the docker
# preparation time annoys you, it is also possibel to
# base work on an image that has these dependencies pre-installed, still
# in pretty much the same way that this script works. For this, you may
# use one of the Dockerfile examples under ci/, e.g.
#
#       docker build -t cado-nfs-debug -f ci/Dockerfile.debian ci
#       docker run --rm -ti --hostname docker-script-$RANDOM --volume $PWD:/host cado-nfs-debug bash /host/.docker-script.sh
#

export DOCKER_SCRIPT=1

# This installs packages just as our CI jobs do. However, because of the
# DOCKER_SCRIPT environment variable, some bonus packages are installed
# as well (vim, gdb).
/host/ci/00-prepare-docker.sh

. /host/ci/000-functions.sh
. /host/ci/001-environment.sh
. /host/ci/999-debug.sh
