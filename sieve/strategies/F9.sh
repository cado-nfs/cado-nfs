#!/usr/bin/env bash

# this file generates a strategy file for F9 (fbb=21, lpb=26)

# set CADO_NFS to the CADO-NFS repository
# set CADO_BIN to the directory where binaries are

if [ "$CADO_DEBUG" ] ; then set -x ; fi

fbb=21
lpb=26

#select our factoring methods
$CADO_BIN/sieve/strategies/gfm -m PM1 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 50 -cmin 10 -cmax 50 -cstep 10 -out pm1
$CADO_BIN/sieve/strategies/gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 50 -cmin 10 -cmax 50 -cstep 10 -out pp1

$CADO_BIN/sieve/strategies/gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 100 -b1max 500 -b1step 50 -cmin 10 -cmax 100 -cstep 10 -out ecm-m12

$CADO_BIN/sieve/strategies/gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 100 -b1max 500 -b1step 50 -cmin 10 -cmax 100 -cstep 10 -out ecm-m16

$CADO_BIN/sieve/strategies/gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 50 -cmin 10 -cmax 100 -cstep 10 -out ecm-b12

$CADO_BIN/sieve/strategies/gfm -fch -fch_in pm1 -fch_out pm1_ch
$CADO_BIN/sieve/strategies/gfm -fch -fch_in pp1 -fch_out pp1_ch
$CADO_BIN/sieve/strategies/gfm -fch -fch_in ecm-b12 -fch_out ecm-b12_ch
$CADO_BIN/sieve/strategies/gfm -fch -fch_in ecm-m12 -fch_out ecm-m12_ch
$CADO_BIN/sieve/strategies/gfm -fch -fch_in ecm-m16 -fch_out ecm-m16_ch

$CADO_BIN/sieve/strategies/benchfm -in pm1_ch -p -t -lb 21 -out pm1_pt  &
$CADO_BIN/sieve/strategies/benchfm -in pp1_ch -p -t -lb 21 -out pp1_pt  &
$CADO_BIN/sieve/strategies/benchfm -in ecm-b12_ch -p -t -lb 21 -out ecm-b12_pt
$CADO_BIN/sieve/strategies/benchfm -in ecm-m12_ch -p -t -lb 21 -out ecm-m12_pt & 
$CADO_BIN/sieve/strategies/benchfm -in ecm-m16_ch -p -t -lb 21 -out ecm-m16_pt  &

cat pm1_pt pp1_pt ecm-b12_pt ecm-m12_pt ecm-m16_pt > All_methods

mkdir decomp
./script_gen_dec.sh 2300000 43 52&
./script_gen_dec.sh 1200000 41 52&

#create strategies for one length!
mkdir res_precompt_st
./script_str_r.sh 2300000 26 0 20 10 decomp All_methods res_precompt_st
./script_str_r.sh 2300000 26 21 40 10 decomp All_methods res_precompt_st
./script_str_r.sh 2300000 26 41 52 10 decomp All_methods res_precompt_st
./script_str_r.sh 1200000 26 0  20 10 decomp All_methods res_precompt_st
./script_str_r.sh 1200000 26 21 40 10 decomp All_methods res_precompt_st
./script_str_r.sh 1200000 26 41 52 10 decomp All_methods res_precompt_st

/bin/rm -fr decomp

mkdir res_matrix
./script_str.sh 2300000 52 1200000 52 res_precompt_st res_matrix

/bin/rm -fr res_precompt_st

$CADO_BIN/sieve/las -poly F9.poly -fb1 F9.roots -I 12 -lim0 2300000 -lpb0 26 -mfb0 52 -lambda0 2.1 -lim1 1200000 -lpb1 26 -mfb1 52 -lambda1 2.2 -q0 1200000 -q1 1200200 -v -v -unsievethresh 1  -stats-cofact cofactors.stats
# Total cpu time 2.14s [norm 0.01+0.0, sieving 0.9 (0.7 + 0.0 + 0.2), factor 1.2 (1.2 + 0.0)]
# Total elapsed time 2.16s, per special-q 0.27059s, per relation 0.00503423s
# Total 430 reports [0.00497s/r, 53.8r/sq]

$CADO_BIN/sieve/strategies/finalst -st res_matrix -dist cofactors.stats -t 2.14 -mfb0 52 -mfb1 52 -out final_st
/bin/rm -fr res_matrix
# Y = 422.847013 relations, T = 2.161445 s., yt = 0.0051116476 s/rel.

# Total cpu time 2.16s [norm 0.00+0.0, sieving 1.0 (0.7 + 0.1 + 0.2), factor 1.2 (1.2 + 0.0)]
# Total elapsed time 2.19s, per special-q 0.273291s, per relation 0.00508448s
# Total 430 reports [0.00502s/r, 53.8r/sq]
