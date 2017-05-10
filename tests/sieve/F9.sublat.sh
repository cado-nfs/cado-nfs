#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
CHECKSUM_FILE="$4"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 4

REFERENCE_SHA1="090c7a9172ba1e70a459df8b0d7e024c4289d093"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="604a709c968e05ec09423dba1a977e100f038a90"

lim0=2300000
lim1=1200000
lpb0=26
lpb1=26
maxbits=10
mfb0=52
mfb1=52
lambda0=2.1
lambda1=2.2
I=10
q0=1200000
q1=1200100

export lim0 lim1 lpb0 lpb1 maxbits mfb0 mfb1 lambda0 lambda1 I q0 q1
"${SOURCE_TEST_DIR}"/sievetest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "${CHECKSUM_FILE}" -v -v --adjust-strategy 0 -sublat 3 "$@" || exit 1

exit 0
