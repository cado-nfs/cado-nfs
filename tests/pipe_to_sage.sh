#!/usr/bin/env bash

set -e
: ${wdir:?missing}
"$@" > "$wdir/test.sage"
"$SAGE" "$wdir/test.sage"
