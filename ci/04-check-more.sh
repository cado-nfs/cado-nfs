#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
enter_section "test" "Running tests"
env OMP_DYNAMIC=true STATS_PARSING_ERRORS_ARE_FATAL=1 make check ARGS="-j$NCPUS"
leave_section
