#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls batchtest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
CHECKSUM_FILE="$4"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 4

# the checksum should not change in batch mode, which should find the maximum
# number of relations
REFERENCE_SHA1="ada1f8268b249109dd3e4d3bab5abf0513e817ac" # 272 relations

rlim=2300000
alim=1200000
lpbr=26
lpba=24
maxbits=10
mfbr=52
mfba=72
rlambda=2.2
alambda=2.9
I=12
q0=1200000
q1=1200200

export rlim alim lpbr lpba maxbits mfbr mfba rlambda alambda I q0 q1
"${SOURCE_TEST_DIR}"/batchtest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "${CHECKSUM_FILE}" "$@" || exit 1

exit 0
