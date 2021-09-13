#!/usr/bin/env bash

# This is actually __ONLY__ called from ci/40-testsuite.sh now

if ! [ "$sourced_from_testsuite" ] ; then
    . "${CI_PATH:-$(dirname $0)}/000-functions.sh"
    . "${CI_PATH:-$(dirname $0)}/001-environment.sh"
fi

NCPUS=`"${CI_PATH:-$(dirname $0)}/utilities/ncpus.sh"`
enter_section build Building
(cd "$build_tree" ; "${MAKE}" -j$NCPUS)
leave_section
