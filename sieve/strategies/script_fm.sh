#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then set -x ; fi

fbb=$1
lpb=$2

#$fbb-$lpb

./gfm -m PM1 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PM1 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PM1 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PM1 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PM1 -lb $fbb -ub $lpb -b1min 600 -b1max 1000 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PM1 -lb $fbb -ub $lpb -b1min 600 -b1max 1000 -b1step 20 -cmin 50 -cmax 200 -cstep 10&

./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M12 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 50 -cmax 300 -cstep 10

./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-M16 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &

./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 300 -cstep 10 &
./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m ECM-B12 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 50 -cmax 300 -cstep 10 

./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-27 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &

./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 50 -b1max 300 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 300 -b1max 600 -b1step 20 -cmin 50 -cmax 200 -cstep 10 &
./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 10 -cmax 50 -cstep 5 &
./gfm -m PP1-65 -lb $fbb -ub $lpb -b1min 600 -b1max 1500 -b1step 20 -cmin 50 -cmax 200 -cstep 10


