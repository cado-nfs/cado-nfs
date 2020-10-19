#!/bin/bash
# This script is a debugging helper that can be used to mimick what is
# done in the build jobs triggered by .gitlab-ci.yml
#
# On a machine with docker installed, this script can be invoked with the
# following command lines:
#
# Example with the "gcc" image that is used by the .gitlab-ci.yml jobs:
#     docker run --rm -ti  --volume $PWD:/host gcc bash /host/.docker-script.sh
#
# Example with another base image, here "fedora". The script below has
# provision for fedora as well.
#     docker run --rm -ti  --volume $PWD:/host fedora bash /host/.docker-script.sh
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
# In particular, you might find it useful to do things such as:
#       export force_build_tree=/tmp/b
#       make cmake
#       [...]
#
# Note that on purpose, this script reinstalls the needed dependencies.
# This is because the gitlab-ci.yml jobs do so as well, and we want to be
# able to be in sync with what they do. Of course, it is also possibel to
# base work on an image that has these dependencies pre-installed, but
# that would not suit the desired purpose here.

if [ -f /etc/debian_version ] ; then
    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen
    env DEBIAN_FRONTEND=noninteractive apt-get -y update
    env DEBIAN_FRONTEND=noninteractive apt-get -y install bc locales cmake libhwloc-dev libgmp-dev python3 git vim gdb g++
elif [ -f /etc/fedora-release ] ; then
    dnf -y install bc cmake hwloc-devel gmp-devel python git vim gdb gcc g++ make
fi
cd /host
uid=$(stat --printf '%u' .)
gid=$(stat --printf '%g' .)
groupadd -g $gid hostgroup
useradd -m -g $gid -u $uid hostuser
if ! [ -f /etc/fedora-release ] ; then
    cat <<EOF
# NOTE: You might need to type:
#       exec bash -i
# first, so as to get a proper shell
#
EOF
fi
if ! [ -f build ] && ! [ -e build ] ; then
    cat <<EOF
# NOTE: ./build appears to be a symlink that goes outside the
# current container.
# You might want to set the build directory by hand to work around this,
# e.g. with:
#       export force_build_tree=/tmp/b
#
EOF
fi
exec su hostuser
