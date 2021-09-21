#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"
. "$(dirname $0)/003-trap-add.bash"
. "$(dirname $0)/004-disksize-watchdog.bash"

"$(dirname $0)"/01-conf.sh
"$(dirname $0)"/02-build.sh
"$(dirname $0)"/03-check.sh
