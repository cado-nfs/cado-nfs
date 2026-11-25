#!/usr/bin/env python3
"""Final per-size optimisation (typical discriminant, alpha given).
   final.py d I n template.json lims=.. lpbs=.. [alpha=0] [nbits=704] [maxsim=2.5e6] [step=1.25]
Sieve cost: yield model (ymodel) x measured per-position sieve time; filtering: proxy relations
(template padded to the target size) -> purge/merge-dl; linear algebra: la_cost."""
import os, sys, json, math
from proxy import make_proxy
from simfilter import simulate
from opt import la_cost, pi_approx
from ymodel import yield_sq, xeff
import timing2

# per-position sieve time (ns) = a + b*sqrt(lim/1e6), fitted per number of special-q factors
TPOS = {13: (0.55, 10.5), 14: (1.23, 6.85), 15: (1.6, 5.4), 16: (1.93, 4.17), 17: (1.7, 3.6), 18: (1.53, 3.1)}
def t_sq(I, n, lim, d=None):
    T = timing2.load()
    if d is not None:
        for key, v in T.items():
            kd, kI, kn, klim = map(int, key.split(",")[:4])
            if (kd, kI, kn, klim) == (d, I, n, lim):
                return v["core_s_per_sq"]
    a, b = TPOS[n]
    return (a + b * math.sqrt(lim / 1e6)) * 1e-9 * 2.0 ** (I + n - 1)

def typical_D(d):
    return 10 ** (d - 0.5)

def run(d, I, n, tpl, lims, lpbs, alpha=0.0, nbits=704, maxsim=2.5e6, step=1.25, neq=2, out=None, beta=1.95,
        quiet=False):
    sj = json.load(open(tpl))
    D0 = abs(int(sj["D"])); Dt = typical_D(d)
    dX = 0.5 * math.log(Dt / D0) + (I - sj["I"]) * math.log(2)
    out = out or f"final/d{d}.I{I}.n{n}.a{alpha:g}.jsonl"
    os.makedirs("final", exist_ok=True)
    done = {}
    if os.path.exists(out):
        for l in open(out):
            r = json.loads(l); done[(r["lim"], r["lpb"], r["k"])] = r
    res = []
    with open(out, "a") as fo:
        for lim in lims:
            ts = t_sq(I, n, lim, d)
            for lpb in lpbs:
                if lim > sj["lim"] or lpb > sj["lpb"]:
                    print(f"template too small for lim={lim} lpb={lpb}"); continue
                y, P = yield_sq(Dt, I, n, lim, lpb, alpha)
                Xt = xeff(Dt, I, alpha)
                proxyf = f"samples/tmp.proxy.{os.getpid()}"
                made = False
                R0 = 3.0 * pi_approx(lim) / 2
                best = None; rising = 0; nsucc = 0
                for k in range(45):
                    R = R0 * step ** k
                    if (lim, lpb, k) in done:
                        r = done[(lim, lpb, k)]
                    else:
                        if not made:
                            info = make_proxy(sj["file"], proxyf, max(dX, 0.0), lim, lpb, Xt, 80000,
                                              qlist=range(sj["qmin"], sj["qmax"]))
                            made = True
                        shrink = 1
                        while R / shrink > maxsim:
                            shrink *= 2
                        r = simulate(proxyf, R, lpb, 30, shrink=shrink, tag=f"f{d}")
                        r.update(lim=lim, lpb=lpb, k=k, d=d, I=I, n=n, alpha=alpha, y=y, t_sq=ts)
                        for kk in ("purge_tail", "merge_tail"):
                            r.pop(kk, None)
                        fo.write(json.dumps(r) + "\n"); fo.flush()
                    if "N" not in r:
                        continue
                    nsucc += 1
                    N, W = r["N_scaled"], r["W_scaled"]
                    if r["shrink"] > 1 and N > 8000:
                        W = max(W, 170.0 * N)     # small simulated matrices stop merging early
                    sv = R / y * ts; la, Pn = la_cost(N, W, nbits=nbits, neq=neq, beta=beta)
                    tot = sv + la
                    rr = dict(lim=lim, lpb=lpb, R=R, nsq=R / y, y=y, t_sq=ts, N=N, w=W / N, P=Pn, sieve=sv, la=la,
                              total=tot, purged=r["p_rows_scaled"], shrink=r["shrink"])
                    res.append(rr)
                    if not quiet:
                        print(f"d{d} lim={lim} lpb={lpb} R={R:.3g} nsq={R/y:.0f} y={y:.3g} N={N:.0f} w={W/N:.0f} P={Pn} "
                              f"sieve={sv/86400:.3g}cd la={la/86400:.3g}cd total={tot/86400:.4g}cd", flush=True)
                    if best is None or tot < best:
                        best = tot; rising = 0
                    else:
                        rising += 1
                    if rising >= 2 and nsucc >= 4:
                        break
                if made:
                    os.unlink(proxyf)
    return res

def summarize(res):
    import collections
    g = collections.defaultdict(list)
    for r in res:
        g[(r["lim"], r["lpb"])].append(r)
    for key in sorted(g):
        r = min(g[key], key=lambda r: r["total"])
        Rs = sorted(x["R"] for x in g[key])
        edge = " (low edge)" if r["R"] == Rs[0] else (" (HIGH edge)" if r["R"] == Rs[-1] else "")
        print(f"  lim={key[0]:>9} lpb={key[1]} y={r['y']:.3g} t_sq={r['t_sq']:.3g} R={r['R']:.3g}{edge} nsq={r['nsq']:.0f} purged={r['purged']:.0f} "
              f"N={r['N']:.0f} w={r['w']:.0f} P={r['P']} sieve={r['sieve']/86400:.3g} la={r['la']/86400:.3g} total={r['total']/86400:.4g} core-days")

if __name__ == "__main__":
    d, I, n, tpl = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    kw = dict(a.split("=") for a in sys.argv[5:])
    lims = [int(float(x)) for x in kw.pop("lims").split(",")]
    lpbs = [int(x) for x in kw.pop("lpbs").split(",")]
    kw = {k: (int(v) if k in ("nbits", "neq") else float(v)) for k, v in kw.items()}
    res = run(d, I, n, tpl, lims, lpbs, **kw)
    print(f"== d{d} I={I} n={n} {kw}")
    summarize(res)
