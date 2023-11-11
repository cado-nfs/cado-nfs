#!/usr/bin/env bash

# bwc.pl and linalg/bwc/secure both have a notion of default check stops.
# This test makes sure that they agree on that. It's not terribly
# important.

base_args=("$@")

set -e

tmp=$(mktemp -d /tmp/cado-nfs.XXXXXXXXXXXXXX)

if ! [ "$CADO_DEBUG" ] ; then
    trap "rm -rf $tmp" EXIT
fi

pre_doubledash=()
post_doubledash=()
restricted_for_direct_call=()
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
    if [[ $x =~ ^(seed|m|n|mn|matrix|wdir|nullspace)=(.*) ]] ; then
        restricted_for_direct_call+=("$x")
    fi
    pre_doubledash+=("$x")
done
post_doubledash=("$@")
: ${wdir:?missing}
: ${seed:?missing}
: ${interval:?missing}

# these arguments will be fed to bwc.pl, which will add its own
# check_stops parameter
args1=(
    "${pre_doubledash[@]}"
    interval=$interval
    script_steps=wipecheck,matrix,bwc.pl/prep/secure
    "${post_doubledash[@]}"
)
"`dirname $0`"/bwc-ptrace.sh "${args1[@]}" | tee "$tmp/run.log"
mkdir "$wdir/saved_check"
mv "$wdir"/C[rvdt]* "$wdir/saved_check"
cmd=($(perl -ne 'm{(\S+/linalg/bwc/secure.*)$} && do { $_=$1; s{check_stops=\S+}{}; print "$_\n"; };' "$tmp/run.log"))
[ "${cmd[*]}" ]
"${cmd[@]}"


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
