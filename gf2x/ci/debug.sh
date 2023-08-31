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
# NOTE: to debug the freebsd tests, this wrapper script tries to mimick
# what we do in gitlab-ci as well. But it doesn't work in exactly the
# same way. So, you might want to use this script, and you might choose
# not to. The possible things to try are:
#
#  ./ci/debug.sh "checks on freebsd13.0 system with gcc"
#     --> This will fire up a freebsd13.0 virtual machine, and reach a
#     shell within an sshfs-mounted copy of your current directory. As a
#     consequence, please note that symlink traversal won't work, and
#     sshfs might incur some delays.
#
#  DOCKER_SCRIPT=1 CI_JOB_NAME='checks on freebsd13 with gcc' ci/50-libvirt-wrap-tests.sh freebsd:13.0
#     --> This second option is more faithful to what we do with
#     gitlab-ci, and gives a bash shell on a fresh tree in sync with
#     current git HEAD. But there's no funny volume mounting and so on.
#
# In either case, getting to the shell prompt takes from 2 to 10 minutes
# depending on a variety of factors, and on the progress of caching in
# the tanker script (at the moment: none, but I'm thinking of it).

set -e
export DOCKER_SCRIPT=1
export CI_JOB_NAME="$1"
tmp=$(mktemp -d /tmp/XXXXXXXXXXXXXX)
trap "rm -rf $tmp" EXIT

DARGS=()
if ! [ "$NO_REMOVE" ] ; then
    DARGS+=(--rm)
else
    echo "# NOTE: the container will remain up on script exit"
fi

if [[ "$CI_JOB_NAME" =~ freebsd([0-9]+\.[0-9]+) ]] ; then
    # use tanker script instead
    IMAGE_NAME=freebsd:"${BASH_REMATCH[1]}"
    if [[ "$CI_JOB_NAME" =~ '32-bit freebsd' ]] ; then
        IMAGE_NAME+="?arch=i386"
    fi
    . "$(dirname $0)/002-tanker.bash"
    # create base image. As in the current gitlab-ci case, we have no
    # caching, which is a pity.
    myimage=gf2x-"$(id -n -u)"-"$IMAGE_NAME"
    tanker vm build -R -t $myimage "$IMAGE_NAME" @host "$(dirname $0)/" env "${exports[@]}" ./00-prepare-docker.sh

    cat > $tmp/getaddress.xsl <<'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<stylesheet version="1.0" xmlns="http://www.w3.org/1999/XSL/Transform">
<output method="text" /><template match="/"><for-each select="//network/ip">
<value-of select="@address" /><text>
</text></for-each></template></stylesheet>
EOF
    libvirt_host=$(xsltproc $tmp/getaddress.xsl <(virsh -c qemu:///system net-dumpxml default))
    random=$(uuidgen  | sha256sum | cut -c1-10)
    # contrary to what I found in some places, user access to fues mounts
    # is not a question of group membership. /dev/fuse is 0666 anyway.
    # However the sysctl matters!
    # pw groupmod operator -m user --
    force_build_tree=/tmp/b
    exports+=(force_build_tree=$force_build_tree)
    echo "# NOTE: gf2x build tree has just been set to $force_build_tree"
    commands=(
        @guest root@ env ASSUME_ALWAYS_YES=yes pkg install fusefs-sshfs \;
                     kldload fusefs \;
                     sysctl vfs.usermount=1 --
        @guest user@ mkdir /tmp/$random \;
                     sshfs -o idmap=user,StrictHostkeyChecking=no $(id -u -n)@$libvirt_host:$PWD /tmp/$random \;
                     cd /tmp/$random \;
                     env "${exports[@]}" bash
    )
    tanker vm run "${DARGS[@]}" -t $myimage "${commands[@]}"
else
    # remove CI_JOB_NAME from the args!
    shift
    : ${imagename=docker-image-$RANDOM}
    # run with bash instead of sh. the "docker" docker image has a
    # /bin/sh shell that groks "set -o pipefail", which appears in
    # ci/00-docker-build.sh ; such is not the case of /bin/sh on debian,
    # at least.
    ci/00-docker-build.sh "$imagename" "$@"
    echo "# NOTE: docker image is $imagename"
    echo "# NOTE: this image contains a few extra debug tools"
    # CI_JOB_NAME is passed to the script via 00-dockerfile.sh
    docker run "${DARGS[@]}" -ti --hostname docker-script-$RANDOM --volume $PWD:/host "$imagename" /host/ci/999-debug.sh "$@"
fi
