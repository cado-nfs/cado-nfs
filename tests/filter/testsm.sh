#!/usr/bin/env bash

# The main functions to compute SM are tested in utils/test_sm_utils.
# This test is here to check that the multi-threaded and mono-threaded versions
# give the same results and that the -nsm option works correctly.

SM="$1"
SOURCE_TEST_DIR="`dirname "$0"`"
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${TMPDIR}"

poly="${SOURCE_TEST_DIR}/testsm.p59.poly"
purged="${SOURCE_TEST_DIR}/testsm.p59.purged.gz"
id="${SOURCE_TEST_DIR}/testsm.p59.index.gz"
go="2926718140519"


args="-poly ${poly} -purged ${purged} -index ${id} -gorder ${go} -nsm 0,4"
args2="-poly ${poly} -purged ${purged} -index ${id} -gorder ${go} -nsm 0,2"

#with -nsm 0,4 (ie nsm = deg F-1) and -t 1
if ! "${SM}" ${args} -out "${TMPDIR}/sm.4.1" -t 1 ; then
  echo "$0: sm binary failed with -t 1 and without -nsm. Files remain in ${TMPDIR}"
  exit 1
fi
#with -nsm 0,4 (ie nsm = deg F-1) and -t 1
if ! "${SM}" ${args} -out "${TMPDIR}/sm.4.2" -t 2 ; then
  echo "$0: sm binary failed with -t 2 and without -nsm. Files remain in ${TMPDIR}"
  exit 1
fi

#with -nsm 0,2 and -t 1
if ! "${SM}" ${args2} -out "${TMPDIR}/sm.2.1" -t 1 ; then
  echo "$0: sm binary failed with -t 1 and -nsm1 2. Files remain in ${TMPDIR}"
  exit 1
fi
#with -nsm 0,2 and -t 2
if ! "${SM}" ${args2} -out "${TMPDIR}/sm.2.2" -t 2 ; then
  echo "$0: sm binary failed with -t 2 and -nsm1 2. Files remain in ${TMPDIR}"
  exit 1
fi


if ! diff -b "${TMPDIR}/sm.4.1" "${TMPDIR}/sm.4.2" > /dev/null ; then
  echo "$0: Mono-threaded and multi-threaded versions do not match (without -nsm). Files remain in ${TMPDIR}"
  exit 1
fi

if !  diff -b "${TMPDIR}/sm.2.1" "${TMPDIR}/sm.2.2" > /dev/null ; then
  echo "$0: Mono-threaded and multi-threaded versions do not match (with -nsm1 2). Files remain in ${TMPDIR}"
  exit 1
fi

tail -n +2 "${TMPDIR}/sm.4.1" | cut -d ' ' -f1-2 > "${TMPDIR}/sm.4.1.short"
tail -n +2 "${TMPDIR}/sm.2.1" > "${TMPDIR}/sm.2.1.short"

if ! diff -b "${TMPDIR}/sm.4.1.short" "${TMPDIR}/sm.2.1.short" > /dev/null ; then
  echo "$0: First two SMs computed without -nsm do not match SMs computed with -nsm 0,2). Files remain in ${TMPDIR}"
  exit 1
fi

rm -f "${TMPDIR}/sm.4.1" "${TMPDIR}/sm.4.2"
rm -f "${TMPDIR}/sm.2.1" "${TMPDIR}/sm.2.2"
rm -f "${TMPDIR}/sm.4.1.short" "${TMPDIR}/sm.2.1.short"
rmdir "${TMPDIR}"
