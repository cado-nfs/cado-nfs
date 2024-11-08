#!/usr/bin/env bash

# The main functions to compute SM are tested in utils/test_sm_utils.
# This test is here to check that the multi-threaded and mono-threaded versions
# give the same results and that the -nsms option works correctly.

build_tree="$PROJECT_BINARY_DIR"

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        build_tree="$1"
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

SM="$build_tree/filter/sm"
SOURCE_TEST_DIR="`dirname "$0"`"

: ${wdir:?missing}

poly="${SOURCE_TEST_DIR}/testsm.p59.poly"
purged="${SOURCE_TEST_DIR}/testsm.p59.purged.gz"
id="${SOURCE_TEST_DIR}/testsm.p59.index.gz"
go="2926718140519"


args="-poly ${poly} -purged ${purged} -index ${id} -ell ${go} -nsms 0,4"
args2="-poly ${poly} -purged ${purged} -index ${id} -ell ${go} -nsms 0,2"

#with -nsms 0,4 (ie nsms = deg F-1)
if ! "${SM}" ${args} -out "${wdir}/sm.4.1" ; then
  echo "$0: sm binary failed without -nsms. Files remain in ${wdir}"
  exit 1
fi

#with -nsms 0,2 and -t 1
if ! "${SM}" ${args2} -out "${wdir}/sm.2.1" ; then
  echo "$0: sm binary failed with -nsm1 2. Files remain in ${wdir}"
  exit 1
fi


tail -n +2 "${wdir}/sm.4.1" | cut -d ' ' -f1-2 > "${wdir}/sm.4.1.short"
tail -n +2 "${wdir}/sm.2.1" > "${wdir}/sm.2.1.short"

if ! diff -b "${wdir}/sm.4.1.short" "${wdir}/sm.2.1.short" > /dev/null ; then
  echo "$0: First two SMs computed without -nsms do not match SMs computed with -nsms 0,2). Files remain in ${wdir}"
  exit 1
fi
