#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

"$(dirname $0)"/01-conf.sh
"$(dirname $0)"/02-build1.sh
"$(dirname $0)"/02-build2.sh
"$(dirname $0)"/03-check.sh
