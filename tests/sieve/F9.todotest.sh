#!/usr/bin/env bash

FB="$1"
LAS="$2"
SRCDIR="$3"
TODO="$4"

cmd="${LAS} -poly ${SRCDIR}/parameters/polynomials/F9.poly -lim0 2300000 -lim1 1200000 -lpb0 26 -lpb1 26 -mfb0 52 -mfb1 52 -lambda0 2.0 -lambda1 2.0 -I 12 -fb1 ${FB} -todo ${TODO} --adjust-strategy 0"
echo $cmd
$cmd | grep -e "Total 420 reports" || exit 1

exit 0
