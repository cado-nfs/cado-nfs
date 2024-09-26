#!/usr/bin/env bash

# We need this silly script inside the test definition just in order to
# catch potential precommands...

"${TEST_PRECOMMAND[@]}" "$@"
