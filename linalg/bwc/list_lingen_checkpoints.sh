#!/usr/bin/env bash

# This script can be used to list all possible lists of file stems that
# can be passed to lingen_verify_checkpoints. The directory in $1 is
# scanned for lingen checkpoints, and each output line is a triple (or
# pair) of polynomials that are eligible for checkpoint verification.

cd "$1"
ls E*aux | while read x ; do
    read _ level t0 aux < <(echo $x  | tr . ' ')
    pi0=`ls pi.$level.$t0.*aux 2>/dev/null`
    pi1=`ls pi.$((level+1)).$t0.*aux 2>/dev/null`
    read uu0 uu1 uu2 th uu3 < <(echo $pi1 | tr . ' ')
    pi2=`ls pi.$((level+1)).$th.*aux 2>/dev/null`
    if [ "$pi0" ] ; then
        echo ${x%%.aux} ${pi0%%.aux}
    fi
    if [ "$pi2" ] ; then
        echo  ${pi1%%.aux} ${pi2%%.aux} ${pi0%%.aux}
    fi
done

