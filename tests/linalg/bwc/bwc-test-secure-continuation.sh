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

CSI_RED="[01;31m"
CSI_BLUE="[01;34m"
CSI_RESET="[00;39m[m"

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

echo -e "${CSI_BLUE}RUNNING linalg/bwc/secure FOR THE 1st TIME${CSI_RESET}"
"`dirname $0`"/bwc-ptrace.sh "${args1[@]}"
# rm -f "$wdir"/C[d]*
echo -e "${CSI_BLUE}RUNNING linalg/bwc/secure FOR THE 2nd TIME${CSI_RESET}"
"`dirname $0`"/bwc-ptrace.sh "${args2[@]}"
mkdir "$wdir/saved_check"
mv "$wdir"/C[rvdt]* "$wdir/saved_check"
echo -e "${CSI_BLUE}RUNNING linalg/bwc/secure FOR THE 3rd TIME${CSI_RESET}"
"`dirname $0`"/bwc-ptrace.sh "${args3[@]}"

echo -e "${CSI_BLUE}COMPARING THE DIFFERENT CHECK FILES${CSI_RESET}"

failed=
passed=0
for f in `(cd "$wdir" ; ls C[rvdt]* ; cd "$wdir/saved_check" ; ls C[rvdt]*) | sort -u` ; do
    if [ -f "$wdir/$f" ] && [ -f "$wdir/saved_check/$f" ] ; then
        if ! diff -q "$wdir/$f" "$wdir/saved_check/$f" ; then
            echo -e "${CSI_RED}Files $wdir/$f and $wdir/saved_check/$f differ${CSI_RESET}" >&2
            sha1sum "$wdir/$f" "$wdir/saved_check/$f"
            failed=1
        else
            echo "Files $wdir/$f and $wdir/saved_check/$f are consistent"
            let passed+=1
        fi
    fi
done
if [ "$passed" -gt 0 ] && ! [ "$failed" ] ; then
    echo -e "${CSI_BLUE}Check files are consistent ($passed checks done), good${CSI_RESET}"
elif [ "$passed" -eq 0 ] && ! [ "$failed" ] ; then
    echo -e "${CSI_RED}No checks done. This is weird${CSI_RESET}"
    exit 1
else
    exit 1
fi
