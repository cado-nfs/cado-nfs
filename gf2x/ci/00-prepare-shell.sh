#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

echo "Enter CI script for $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME, stage $CI_JOB_STAGE ; $CI_JOB_NAME"

enter_section preparation "System preparation (shell).  Checking required and optional software"

check_mandatory_tools autoconf autoreconf

# These are directories where we want to search for include files, but
# without any assurance that these directories exist.
tested_include_paths=(/usr/include /opt/homebrew/include)

# check_mandatory_file "${tested_include_paths[@]/%//gmp.h}"
check_mandatory_nonzero_output_shell find "${tested_include_paths[@]}" -name gmp.h

leave_section
