#!/usr/bin/env bash

SOURCE_TEST_DIR="`dirname "$0"`"
BUILD_DIR="$1"
WORKDIR=`mktemp -d ${TMPDIR-/tmp}/cadotest.XXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${WORKDIR}"

# Test MNFS with 5 polys (interesting primes are 2 and 3)
POLY="${SOURCE_TEST_DIR}/test_renumber.data/mnfs5.poly"
BADIDEALS="${SOURCE_TEST_DIR}/test_renumber.data/mnfs5.badideals"
RENUMBER="${WORKDIR}/mnfs5.renumber.gz"

${BUILD_DIR}/sieve/freerel -poly ${POLY} -renumber ${RENUMBER} -out /dev/null \
                           -lpbs 11,10,10,10,10 -badideals ${BADIDEALS} \
                           -lcideals
if [ ! $? ] ; then
  echo "$0: freerel binary failed. Files remain in ${WORKDIR}"
  exit 1
fi
${BUILD_DIR}/misc/debug_renumber -poly ${POLY} -renumber ${RENUMBER} -check
if [ ! $? ] ; then
  echo "$0: debug_renumber binary failed. Files remain in ${WORKDIR}"
  exit 1
fi

rm -rf "${WORKDIR}"
