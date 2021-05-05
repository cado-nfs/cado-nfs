#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
enter_section build2 "Building test depedencies"
if [ "$using_cmake_directly" ] ; then
    (cd "$build_tree" ; "${MAKE}" -j$NCPUS all_test_dependencies)
else
    "${MAKE}" -j$NCPUS all_test_dependencies
fi
leave_section

