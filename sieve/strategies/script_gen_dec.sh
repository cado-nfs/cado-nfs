#!/usr/bin/env bash

#to precompute all decompostions of cofactors!

if [ "$CADO_DEBUG" ] ; then set -x ; fi

lim=$1
x=$2
mfb=$3
while [ $x -le $mfb ]
do
    $CADO_NFS/sieve/strategies/gst -gdc -lim0 $lim -mfb0 $x -out decomp_tmp/decomp_${lim}_${x}
    x=$[$x+1]
done
