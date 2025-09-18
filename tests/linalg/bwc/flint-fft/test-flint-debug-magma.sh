#!/usr/bin/env bash

: ${wdir=/tmp}

TMPDIR="${wdir}"
export TMPDIR

set -e

magma="$1"
shift
magmascript="$1"
shift
binary="$1"
shift
binary_args=("$@")

"$binary" "${binary_args[@]}" > $TMPDIR/binary_output.m
cd $TMPDIR
$magma -b binary_output.m $magmascript < /dev/null | tee magma_output.txt
if grep -q 'Assertion failed' magma_output.txt ; then
    exit 1
fi
if grep -qi 'error' magma_output.txt ; then
    exit 1
fi

