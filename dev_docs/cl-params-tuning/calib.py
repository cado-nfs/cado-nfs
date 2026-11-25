#!/usr/bin/env python3
"""Compare measured yields (per special-q, from subsampling) with the Dickman-based model."""
import sys, json, math, collections
import numpy as np
from study import parse_rel
from dickman import psmooth

def yields(samplejson, lims, lpbs):
    sj = json.load(open(samplejson))
    cnt = collections.Counter(); types = collections.defaultdict(lambda: [0, 0, 0])
    with open(sj["file"]) as f:
        for line in f:
            ps = sorted(parse_rel(line), reverse=True)[:3]
            for lim in lims:
                big = [p for p in ps if p >= lim]
                if len(big) > 2: continue
                for lpb in lpbs:
                    if big and big[0] >= (1 << lpb): continue
                    cnt[(lim, lpb)] += 1; types[(lim, lpb)][len(big)] += 1
    return sj, {k: v / sj["nq"] for k, v in cnt.items()}, types

def model(sj, lim, lpb, delta=0.0):
    D = abs(int(sj["D"])); I, n = sj["I"], sj["nfac"]
    X = math.log(2.0 ** (I - 1)) + 0.5 * math.log(D / 2) - delta
    P = psmooth(X, math.log(lim), lpb * math.log(2))
    npos = 2.0 ** (I + n - 1)
    return npos * sum(P), P

if __name__ == "__main__":
    for fn in sys.argv[1:]:
        sj0 = json.load(open(fn))
        lims = [l for l in (50000, 100000, 200000, 400000, 800000, 1600000, 3200000) if l <= sj0["lim"]]
        lpbs = [l for l in range(20, 34, 2) if l <= sj0["lpb"]] + ([sj0["lpb"]] if sj0["lpb"] % 2 else [])
        sj, y, types = yields(fn, lims, lpbs)
        print(fn.split("/")[-1])
        for lim in lims:
            for lpb in sorted(set(lpbs)):
                if (lim, lpb) not in y: continue
                m, P = model(sj, lim, lpb)
                t = types[(lim, lpb)]; tot = sum(t)
                # delta that would match
                lo, hi = -10.0, 20.0
                for _ in range(40):
                    mid = (lo + hi) / 2
                    if model(sj, lim, lpb, mid)[0] < y[(lim, lpb)]: lo = mid
                    else: hi = mid
                print(f"  lim={lim:>8} lpb={lpb}: y={y[(lim,lpb)]:9.3f} model={m:9.3f} ratio={y[(lim,lpb)]/m:5.2f} delta_fit={lo:5.2f} "
                      f"types meas 0/1/2LP = {t[0]/tot:.3f}/{t[1]/tot:.3f}/{t[2]/tot:.3f} model = {P[0]/sum(P):.3f}/{P[1]/sum(P):.3f}/{P[2]/sum(P):.3f}")
