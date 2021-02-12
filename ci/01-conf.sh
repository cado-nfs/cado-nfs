#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

enter_section configuration Configuring
"${MAKE}" cmake
leave_section

