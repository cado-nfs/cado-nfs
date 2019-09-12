#!/usr/bin/env bash

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

seed=0
m=64
n=64
random_stem=30000
sequence_length=200

while [ $# -gt 0 ] ; do
    if [[ "$1" =~ ^(seed|m|n|random_stem|sequence_length|expect_crc_[a-zA-Z]+)=[0-9a-f]+$ ]] ; then
        eval "$1"
        shift
        continue
    elif [[ "$1" =~ ^(expect_crc_[a-zA-Z]+)=([0-9a-f]+(,[0-9a-f]+)*)$ ]] ; then
        eval "`echo "${BASH_REMATCH[1]}=(${BASH_REMATCH[2]})" | tr , ' '`"
        shift
        continue
    elif [[ "$1" =~ ^wdir=(.+)$ ]] ; then
        wdir="${BASH_REMATCH[1]}"
        if ! [ -d "$wdir" ] ; then
            echo "wdir $wdir does not exist" >&2
            exit 1
        fi
        shift
        continue
    elif [[ "$1" =~ ^bindir=(.+)$ ]] ; then
        bindir="${BASH_REMATCH[1]}"
        if ! [ -d "$bindir" ] ; then
            echo "bindir $bindir does not exist" >&2
            exit 1
        fi
        shift
        continue
    else
        echo "argument $1 not understood" >&2
        exit 1
    fi
done

if ! [ "${#expect_crc_pi}" ] || ! [ "${#expect_crc_F}" ] ; then
    echo "Please set expect_crc_pi and expect_crc_F on the command line" >&2
    exit 1
fi

if ! [ -d "$bindir" ] ; then
    echo "bindir $bindir does not exist" >&2
    exit 1
fi

if ! [ -d "$wdir" ] ; then
    echo "wdir $wdir does not exist" >&2
    exit 1
fi

"$(dirname $0)/perlrandom.pl" $random_stem $seed > $wdir/X
A_length=$((m*n*sequence_length/8))
nn=0
while [ $nn -lt $A_length ] ; do
    cat $wdir/X
    let nn+=$random_stem
done | dd ibs=1 count=$A_length > $wdir/seq


# we used to use process substitution here to have lingen write to >(tee
# $wdir/output); our goal being to both catch failure of the lingen
# command, and display+save its output. Alas, there's a race condition,
# as the subprocess may not have completed when we move on to the check
# that comes afterwards.
# 
# https://stackoverflow.com/questions/4489139/bash-process-substitution-and-syncing
# 
# among the possible solutions, I quite like the pipefail option.
set -o pipefail

if ! "$bindir/lingen" m=$m n=$n lingen-input-file=$wdir/seq lingen-output-file=$wdir/seq.gen 2>&1 | tee $wdir/output ; then
    echo "Error running lingen >&2"
    exit 1
fi

rc=0
gotpi=
gotF=
for c in "${expect_crc_F[@]}" ; do
    if grep -q "crc(F)=$c" $wdir/output ; then
        gotF=1
        break
    fi
done
for c in "${expect_crc_pi[@]}" ; do
    if grep -q "crc(pi)=$c" $wdir/output ; then
        gotpi=1
        break
    fi
done
if ! [ "$gotpi" ] ; then
    echo "No accepted crc(pi) found in $wdir/output (expected: ${expect_crc_pi[*]})" >&2
    rc=1
fi
if ! [ "$gotF" ] ; then
    echo "No accepted crc(F) found in $wdir/output (expected: ${expect_crc_F[*]})" >&2
    rc=1
fi

exit $rc
