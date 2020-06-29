#!/usr/bin/env bash

set -e

: ${WORKDIR?missing}

SOURCE_TEST_DIR="`dirname "$0"`"
BUILD_DIR="$1"

# Test MNFS with 5 polys (interesting primes are 2 and 3)
POLY="${SOURCE_TEST_DIR}/test_renumber.data/mnfs5.poly"
RENUMBER="${WORKDIR}/mnfs5.renumber.gz"

${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER} -out /dev/null \
                           -lpbs 11,10,10,10,10 -lcideals

${BUILD_DIR}/misc/debug_renumber -poly ${POLY} -renumber ${RENUMBER} -check -quiet
