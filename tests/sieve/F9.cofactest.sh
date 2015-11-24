#!/usr/bin/env bash

# This file checks las with a -file-cofact command-line option.

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
CHECKSUM_FILE="$4"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 4

REFERENCE_SHA1="e04ae591795ff30af703a2235a4f9fd09d17b380"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="bfc0fd6210b9642e82343b4a28a6808ed5fb8b2f" # 431 relations

rlim=2300000
alim=1200000
lpbr=26
lpba=26
maxbits=10
mfbr=52
mfba=52
rlambda=2.1
alambda=2.2
I=12
q0=1200000
q1=1200200

export rlim alim lpbr lpba maxbits mfbr mfba rlambda alambda I q0 q1
"${SOURCE_TEST_DIR}"/sievetest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "${CHECKSUM_FILE}" -v -v "$@" || exit 1

exit 0
