#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

echo "Enter CI script for $CI_PROJECT_NAMESPACE/$CI_PROJECT_NAME, stage $CI_JOB_STAGE ; $CI_BUILD_NAME"

enter_section preparation "System preparation (shell)"
id -a
export
leave_section

enter_section install_packages "Checking required and optional software"

check_mandatory_tools bc cmake
if [ "$checks" ] ; then
check_mandatory_tools xsltproc
fi
# check_mandatory_files /usr/include/gmp.h
check_optional_files /usr/include/hwloc.h
check_mandatory_nonzero_output_shell 'find /usr/include -name gmp.h'

