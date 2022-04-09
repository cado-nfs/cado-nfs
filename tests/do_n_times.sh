#!/usr/bin/env bash
n="$1"
shift
set -e
for i in `seq 1 $n` ; do
    "$@"
done
