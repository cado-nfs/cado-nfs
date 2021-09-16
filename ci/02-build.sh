#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
enter_section build Building
if [ "$using_cmake_directly" ] ; then
    SOURCEDIR="$PWD"
    (cd "$build_tree" ; "${MAKE}" -j$NCPUS)
else
    "${MAKE}" -j$NCPUS
fi
leave_section

enter_section build2 "Building test depedencies"
if [ "$using_cmake_directly" ] ; then
    SOURCEDIR="$PWD"
    (cd "$build_tree" ; "${MAKE}" -j$NCPUS all_test_dependencies)
else
    "${MAKE}" -j$NCPUS all_test_dependencies
fi
leave_section

