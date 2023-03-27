#!/usr/bin/env bash

env
ls -ld "$GITHUB_WORKSPACE"
id -a

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

enter_section configuration Configuring
if [ "$using_cmake_directly" ] ; then
    (cd "$build_tree" ; cmake "$source_tree")
else
    "${MAKE}" cmake
fi
leave_section

