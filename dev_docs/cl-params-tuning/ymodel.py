#!/usr/bin/env python3
"""Calibrated yield model: relations per special-q."""
import math
from dickman import psmooth
from alpha import alpha
from qsbench import discriminant
C0, C1 = 1.43, 1.17      # X_eff = ln(M sqrt|D|) + 0.4 + C0 + C1*alpha  (fit on d80, d100 samples)
def xeff(absD, I, al=0.0):
    return (I - 1) * math.log(2) + 0.5 * math.log(absD) + 0.4 + C0 + C1 * al
def yield_sq(absD, I, n, lim, lpb, al=0.0, mfb=None):
    X = xeff(absD, I, al)
    P = psmooth(X, math.log(lim), lpb * math.log(2), (mfb or 2 * lpb) * math.log(2))
    return 2.0 ** (I + n - 1) * sum(P), P
if __name__ == "__main__":
    import json
    T = json.load(open("timings.json"))
    for key, v in T.items():
        d, I, n, lim, lpb, mfb, nc = map(int, key.split(","))
        D = discriminant(d); al = alpha(D)
        y, P = yield_sq(abs(D), I, n, lim, lpb, al)
        nrel = v["rels_per_sq"] * 2 * v["n2"]
        print(f"d{d} I={I} n={n} lim={lim} lpb={lpb} alpha={al:+.2f}: measured {v['rels_per_sq']:.3f}/sq (~{nrel:.0f} rels) model {y:.3f} ratio {v['rels_per_sq']/y:.2f}")
