#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
CHECKSUM_FILE="$4"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 4

REFERENCE_SHA1="c12c10808914fb445dd32e6801ffb8f320db36e6"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="0bd0c87d93c709279cc5c3dd471292b28ae04198" # 430 relations
# Previous revisions and number of relations found (most recent first):
# 9253f54766c5438fe3198053f9edd811ee254a17 (421 relations)

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
