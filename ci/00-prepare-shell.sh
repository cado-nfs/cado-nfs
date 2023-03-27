#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

echo "Enter CI script for $REPOSITORY ; $BUILD_NAME"

enter_section preparation "System preparation (shell).  Checking required and optional software"

check_mandatory_tools bc cmake
if [ "$checks" ] ; then
check_mandatory_tools xsltproc
fi

# These are directories where we want to search for include files, but
# without any assurance that these directories exist.
tested_include_paths=(/usr/include /opt/homebrew/include)

# check_mandatory_file "${tested_include_paths[@]/%//gmp.h}"
check_optional_file "${tested_include_paths[@]/%//hwloc.h}"
check_mandatory_nonzero_output_shell find "${tested_include_paths[@]}" -name gmp.h

leave_section
