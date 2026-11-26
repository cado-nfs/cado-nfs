#!/usr/bin/env python3
"""Re-evaluate the stored filtering simulations (final/*.jsonl) with given cost assumptions and print
the optimum per size.  reeval.py [alpha=0] [dup=..] [nbits=704] [sscale=1] [lascale=1]
sscale / lascale: factors on the sieving / linear algebra times, e.g. those printed by calibrate.py."""
import json, math, sys, collections, glob
from opt import la_cost
from ymodel import yield_sq
import timing2
GEOM = {80: (13, 13), 90: (13, 14), 100: (13, 16), 110: (13, 16), 120: (15, 16), 130: (15, 16), 140: (15, 16),
        150: (15, 16), 160: (15, 16)}
FILES = {80: "final/d80.I13.n13.a0.jsonl", 90: "final/d90.I13.n14.a0.jsonl", 100: "final/d100.I13.n16.a0.jsonl",
         110: "final/d110.I13.n16.a0.jsonl", 120: "final/d120.I13.n18.a0.jsonl", 130: "final/d130.I13.n18.a0.jsonl",
         140: "final/d140.I13.n18.a0.jsonl", 150: "final/d150.I13.n18.a0.jsonl", 160: "final/d160.I13.n18.a0.jsonl"}
DUP = {80: 0.07, 90: 0.07, 100: 0.08, 110: 0.08, 120: 0.10, 130: 0.10, 140: 0.10, 150: 0.10, 160: 0.10}
# per-position sieve time (ns) = a + b*sqrt(lim/1e6)
TPOS = {(13, 13): (0.55, 10.5), (13, 14): (1.23, 6.85), (13, 16): (1.93, 4.17), (13, 18): (1.53, 3.1), (15, 16): (1.2, 2.15)}
def t_sq(d, I, n, lim):
    T = timing2.load()
    best = None
    for key, v in T.items():
        kd, kI, kn, klim = map(int, key.split(",")[:4])
        if (kI, kn, klim) == (I, n, lim) and abs(kd - d) <= 20 and not (kI == 16):
            if best is None or abs(kd - d) < best[0]:
                best = (abs(kd - d), v["core_s_per_sq"])
    if best:
        return best[1]
    a, b = TPOS[(I, n)]
    return (a + b * math.sqrt(lim / 1e6)) * 1e-9 * 2.0 ** (I + n - 1)
def evaluate(d, alpha=0.0, nbits=704, lascale=1.0, sscale=1.0, dup=None, neq=2, beta=1.95, geom=None):
    I, n = geom or GEOM[d]
    dup = DUP[d] if dup is None else dup
    Dt = 10 ** (d - 0.5)
    rows = []
    for l in open(FILES[d]):
        r = json.loads(l)
        if "N" not in r: continue
        lim, lpb, R = r["lim"], r["lpb"], r["R"]
        y = yield_sq(Dt, I, n, lim, lpb, alpha)[0]
        ts = t_sq(d, I, n, lim)
        N, W = r["N_scaled"], r["W_scaled"]
        if r["shrink"] > 1 and N > 8000: W = max(W, 170.0 * N)
        raw = R / (1 - dup)
        sv = raw / y * ts * sscale
        la, P = la_cost(N, W, nbits=nbits, neq=neq, beta=beta); la *= lascale
        rows.append(dict(d=d, I=I, nfac=n, lim=lim, lpb=lpb, R=R, raw=raw, nsq=raw / y, y=y, t_sq=ts, N=N, w=W / N, P=P,
                         sieve=sv, la=la, total=sv + la))
    return rows
if __name__ == "__main__":
    kw = dict(a.split("=") for a in sys.argv[1:])
    kw = {k: (int(v) if k in ("nbits", "neq") else float(v)) for k, v in kw.items()}
    allbest = []
    for d in sorted(FILES):
        rows = evaluate(d, **kw)
        g = collections.defaultdict(list)
        for r in rows: g[(r["lim"], r["lpb"])].append(r)
        print(f"d{d} geometry I={rows[0]['I']} nfac={rows[0]['nfac']}")
        bests = []
        for key in sorted(g):
            r = min(g[key], key=lambda r: r["total"]); bests.append(r)
            print(f"   lim={key[0]:>9} lpb={key[1]} y={r['y']:.3g} t_sq={r['t_sq']:.3g} uniq={r['R']:.3g} raw={r['raw']:.3g} nsq={r['nsq']:.3g} N={r['N']:.0f} P={r['P']} "
                  f"sieve={r['sieve']/86400:.4g} la={r['la']/86400:.4g} total={r['total']/86400:.4g} cd")
        allbest.append(min(bests, key=lambda r: r["total"]))
    json.dump(allbest, open("allbest.json", "w"), indent=1)
    print("\nsummary (core-days):")
    for r in allbest:
        print(f" d{r['d']}: I={r['I']} n={r['nfac']} lim={r['lim']} lpb={r['lpb']} raw rels={r['raw']:.3g} nsq={r['nsq']:.3g} N={r['N']:.0f} sieve={r['sieve']/86400:.4g} la={r['la']/86400:.4g} total={r['total']/86400:.4g}")
