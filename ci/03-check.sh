#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
export NCPUS
if [ "$CHECKS_EXPENSIVE" ] ; then
    enter_section "xtest" "Running expensive tests"
else
    enter_section "test" "Running tests"
fi
env OMP_DYNAMIC=true STATS_PARSING_ERRORS_ARE_FATAL=1 make check ARGS="-j$NCPUS"
leave_section
