#!/usr/bin/env bash

echo "programs: $all_programs" >&2
echo "sizes: $all_sizes" >&2

for binary in $all_programs ; do
    for size in $all_sizes ; do
        "$tests_location/$binary" --test $size > /tmp/toto
        if `magma < debug/check_cantor.m | grep true`; then
            echo "ok $binary $size"
        else
            echo "failed $binary $size"
            exit 1
        fi
    done
done

