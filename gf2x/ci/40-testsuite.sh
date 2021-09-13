#!/usr/bin/env bash

CI_PATH=$(cd "$(dirname $0)" > /dev/null ; echo $PWD)

. "${CI_PATH}/000-functions.sh"
. "${CI_PATH}/001-environment.sh"
. "${CI_PATH}/003-trap-add.bash"
. "${CI_PATH}/004-disksize-watchdog.bash"

export TMPDIR=$(mktemp -d /tmp/XXXXXXXXXXXX)
trap_add "rm -rf $TMPDIR" EXIT

# Note that we now **SOURCE** all these scripts. This makes it possible
# to play tricks with what the source and build directories are.

export sourced_from_testsuite=1
. "${CI_PATH}"/01-conf.sh
. "${CI_PATH}"/02-build.sh
. "${CI_PATH}"/03-check.sh
