#!/usr/bin/env bash

FB="$1"
LAS="$2"
SRCDIR="$3"
HINT="$4"

${LAS} -poly "${SRCDIR}/parameters/polynomials/F9.poly" \
    -lim0 2300000 -lim1 1200000 -lpb0 26 -lpb1 26 -mfb0 52 -mfb1 52 \
    -I 12 -fb1 ${FB} -q0 235086167414177 -rho 46506392392611 \
    -descent-hint ${HINT} -allow-largesq -ncurves0 10 -ncurves1 10 || exit 1

exit 0
