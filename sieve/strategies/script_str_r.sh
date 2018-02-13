#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then set -x ; fi

#precompute for each size of cofactors and each side (RAT or ALG), the
#convex hull of strategies!


#the sieving region
lim0=$1 # say 2^27
lpb0=$2 # say 37
r0=$3  #say 0
mfb0=$4 # say 110
ncurves=$5 # say 10

#directory which store our decomposition files!
decomp=$6 # say '/localdisk/trichard/results/decomp_cofactor/decomp_tmp/'
#file which contains our factoring methods (PM1, PP1, EC)!
methods=$7 # say 'All_methods'
#directory where the output file will be stored, because this binary
#generates a file for each r.
out=$8 # say 'res_precompt_st/'

while [ $r0 -le $mfb0 ]; do
pathname_decomp=$decomp"/decomp_"$lim0"_"$r0
$CADO_NFS/sieve/strategies/gst -gst_r -lim0 $lim0 -lpb0 $lpb0 -r0 $r0 -ncurves $ncurves -in $methods -decomp $pathname_decomp -out $out &
let "r0 = r0 + 1"
done

