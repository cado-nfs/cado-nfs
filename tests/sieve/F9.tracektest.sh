#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 3

REFERENCE_SHA1="0ac3255855fe39fe729b771147a36dd72fc60d1b"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="9253f54766c5438fe3198053f9edd811ee254a17"

# Set default values for variables not set by caller
: ${rlim:=2300000}
: ${alim:=1200000}
: ${lpbr:=26}
: ${lpba:=26}
: ${maxbits:=10}
: ${mfbr:=52}
: ${mfba:=52}
: ${rlambda:=2.1}
: ${alambda:=2.2}
: ${I:=12}
: ${q0:=1200007}
: ${rho:=554209}

export rlim alim lpbr lpba maxbits mfbr mfba rlambda alambda I q0 rho
"${SOURCE_TEST_DIR}"/sievetest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "" "$@" || exit 1

exit 0
