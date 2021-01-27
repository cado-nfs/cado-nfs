#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
enter_section build2 "Building test depedencies"
make -j$NCPUS all_test_dependencies
leave_section

