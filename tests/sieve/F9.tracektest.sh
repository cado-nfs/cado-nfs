#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 3

REFERENCE_SHA1="d84c58cae72d27c85778bbce7c84e2eb671b336e"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="b97cad99b07072a6ed47a0d12b03eb768604ab1e"
# REFERENCE_SHA1="0ac3255855fe39fe729b771147a36dd72fc60d1b"
# The Git revision that created the REFERENCE_SHA1 hash
# REFERENCE_REVISION="9253f54766c5438fe3198053f9edd811ee254a17"

# Set default values for variables not set by caller
: ${lim0:=2300000}
: ${lim1:=1200000}
: ${lpb0:=26}
: ${lpb1:=26}
: ${maxbits:=10}
: ${mfb0:=52}
: ${mfb1:=52}
: ${lambda0:=2.1}
: ${lambda1:=2.2}
: ${I:=12}
: ${q0:=1200007}
: ${rho:=554209}

export lim0 lim1 lpb0 lpb1 maxbits mfb0 mfb1 lambda0 lambda1 I q0 rho
"${SOURCE_TEST_DIR}"/sievetest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "" --adjust-strategy 0 "$@" || exit 1

exit 0
