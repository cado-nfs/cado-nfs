#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then set -x ; fi

#merge the strategies for each r0 and each r1 to compute the best
#strategies for each pair of cofactors (r0, r1).

#At the end, we will optain in the directory 'res_matrix' some files
#where each file is a set of optimal strategies for a pair (r0,r1).

lim0=$1 # say 2^27
mfb0=$2 # say 110
lim1=$3 # say 2^29
mfb1=$4 # say 110

#the directory where the strategies related on each value of r were
#stored.
in=$5 # say 'res_precompt_st'
#the directory where the files results will be stored. One file will
#be generated for each pair of bit size (r0,r1).
out=$6 # say 'res_matrix'

let "r0=0"
while [ $r0 -le $mfb0 ] ; do
let "r1=0"
while [ $r1 -le $mfb1 ] ; do
$CADO_NFS/sieve/strategies/gst -gst -mfb0 $mfb0 -lim0 $lim0 -mfb1 $mfb1 -lim1 $lim1 -r0 $r0 -r1 $r1 -in $in -out $out
let "r1 = r1 + 1"
done
let "r0 = r0 + 1"
done
