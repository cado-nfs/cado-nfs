#!/usr/bin/env bash

base_args=("$@")

set -e

pre_doubledash=()
post_doubledash=()
while [ $# -gt 0 ] ; do
    x="$1"
    if [ "$x" = -- ] ; then
        break
    fi
    shift
    if [[ $x =~ ^wdir=(.*) ]] ; then
        wdir="${BASH_REMATCH[1]}"
    elif [[ $x =~ ^interval=(.*) ]] ; then
        interval="${BASH_REMATCH[1]}"
        continue
    elif [[ $x =~ ^seed=(.*) ]] ; then
        seed="${BASH_REMATCH[1]}"
    fi
    pre_doubledash+=("$x")
done
post_doubledash=("$@")
: ${wdir:?missing}
: ${seed:?missing}
: ${interval:?missing}

args1=(
    "${pre_doubledash[@]}"
    interval=$interval
    script_steps=wipecheck,matrix,bwc.pl/prep/secure
    "${post_doubledash[@]}"
    check_stops=$interval
)
args2=(
    "${pre_doubledash[@]}"
    interval=$((4*interval)) start=$interval
    script_steps=keepdir,bwc.pl/secure
    "${post_doubledash[@]}"
    check_stops=$((interval/2)),$interval,$((2*interval))
)
args3=(
    "${pre_doubledash[@]}"
    interval=$((4*interval))
    script_steps=keepdir,bwc.pl/secure
    "${post_doubledash[@]}"
    check_stops=$((interval/2)),$interval,$((2*interval))
)
"`dirname $0`"/bwc-ptrace.sh "${args1[@]}"
rm -f "$wdir"/C[rdt]*
"`dirname $0`"/bwc-ptrace.sh "${args2[@]}"
mkdir "$wdir/saved_check"
mv "$wdir"/C[rvdt]* "$wdir/saved_check"
"`dirname $0`"/bwc-ptrace.sh "${args3[@]}"

failed=
for f in `(cd "$wdir" ; ls C[rvdt]* ; cd "$wdir/saved_check" ; ls C[rvdt]*) | sort -u` ; do
    if ! diff -q "$wdir/$f" "$wdir/saved_check/$f" ; then
        echo "Files $wdir/$f and $wdir/saved_check/$f differ" >&2
        sha1sum "$wdir/$f" "$wdir/saved_check/$f"
        failed=1
    fi
done
if ! [ "$failed" ] ; then
    echo "Check files are consistent, good"
else
    exit 1
fi
