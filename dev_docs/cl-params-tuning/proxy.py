#!/usr/bin/env python3
"""Build a proxy relation sample for a larger discriminant from a real sample at a smaller one.
The large-prime type (0/1/2 LP) is drawn from the smoothness model at the target size; the
factor-base part of a template of that type is padded with extra factor-base primes (drawn from
the pooled factor-base primes of the templates, restricted to < lim) until the size deficit
dX (in nats) is filled.
   make_proxy(template_file, out, dX, lim, lpb, Xt, nout, qmax)"""
import math, random, sys, json
from study import parse_rel
from dickman import psmooth

def make_proxy(tfile, out, dX, lim, lpb, Xt, nout, qlist=(), seed=1, mfb=None):
    rng = random.Random(seed)
    L = 1 << lpb
    mfb = mfb or 2 * lpb
    bytype = {0: [], 1: [], 2: []}
    pool = []
    qset = set(qlist)
    for line in open(tfile):
        ps = parse_rel(line)
        big = [p for p in ps if p >= lim]
        if len(big) > 2 or (big and max(big) >= L) or (big and math.prod(big) >> mfb):
            continue
        bytype[len(big)].append(ps)
        pool.extend(p for p in ps if p < lim and p not in qset)
    P = psmooth(Xt, math.log(lim), lpb * math.log(2), mfb * math.log(2))
    w = [P[0], P[1], P[2]]
    for k in range(3):
        if not bytype[k]:
            w[k] = 0.0
    meas = [len(bytype[k]) for k in range(3)]
    with open(out, "w") as f:
        for i in range(nout):
            k = rng.choices((0, 1, 2), weights=w)[0]
            ps = list(rng.choice(bytype[k]))
            s = 0.0
            while True:
                p = rng.choice(pool)
                lp = math.log(p)
                if s + lp <= dX:
                    ps.append(p); s += lp
                else:
                    # stochastic rounding so that E[sum] = dX
                    if rng.random() < (dX - s) / lp:
                        ps.append(p)
                    break
            ps.sort()
            f.write(f"{i+1:x},1:" + ",".join(f"{p:x}" for p in ps) + ":\n")
    tot = sum(meas)
    return dict(template_types=[m / tot for m in meas], target_types=[x / sum(w) for x in w], ntemplates=tot,
                avg_fb_primes_template=len(pool) / tot)
