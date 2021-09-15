#!/bin/bash

# Usage:
# ./ci/50-libvirt-wrap-tests.sh [IMAGE NAME] [[SCRIPT]]

IMAGE_NAME="$1"
shift

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"
. "$(dirname $0)/002-tanker.bash"
. "$(dirname $0)/003-trap-add.bash"

# if ! tanker image list | grep -q gf2x-"$IMAGE_NAME" ; then
    enter_section configuration "Creating base image gf2x-$IMAGE_NAME"
    tanker vm build -R -t gf2x-"$(id -n -u)"-"$IMAGE_NAME" "$IMAGE_NAME" @host "$(dirname $0)/" env "${exports[@]}" ./00-prepare-docker.sh
    leave_section
# fi

tty=()
context=./
if [ $# = 0 ] ; then
    $ECHO_E "${CSI_BLUE}Starting an interactive shell. Transporting only git archive${CSI_RESET}"
    tty+=(-t)
    tmp=$(mktemp -d /tmp/XXXXXXXXXXXXXXX)
    git archive --format=tar.gz HEAD > $tmp/git.tar.gz
    context=$tmp/git.tar.gz
    trap_add "rm -rf $tmp" EXIT
    set -- bash
fi
tanker vm run --rm "${tty[@]}" gf2x-"$(id -n -u)"-"$IMAGE_NAME" @host user@ "$context" env "${exports[@]}" "$@"
