#!/usr/bin/env python3
"""Slope timing of siqs: core-seconds per special-q with all 4 cores busy (2 procs x 2 threads),
   startup excluded.  Results cached in timings.json keyed by 'd,I,nfac,lim,lpb,mfb,ncurves'."""
import os, sys, json
from concurrent.futures import ThreadPoolExecutor
from qsbench import *
HERE = os.path.dirname(os.path.abspath(__file__))
TJ = HERE + "/timings.json"


def load():
    return json.load(open(TJ)) if os.path.exists(TJ) else {}


def tsq(d, I, nfac, lim, lpb, mfb=None, ncurves=20, target=20.0, k=0, force=False):
    mfb = mfb or 2 * lpb
    key = f"{d},{I},{nfac},{lim},{lpb},{mfb},{ncurves}"
    T = load()
    if key in T and not force:
        return T[key]
    D = discriminant(d, k)
    qmin, qmax, pool = choose_qfac(D, I, nfac)

    def run(n, off):
        def job(i):
            return run_siqs(D, I, lim, lpb, mfb, ncurves, nfac, qmin, qmax, nq=n, idx0=off + i * n, threads=2)
        with ThreadPoolExecutor(2) as ex:
            rs = list(ex.map(job, range(2)))
        return (sum(r["elapsed"] for r in rs) * 2, sum(r["cpu"] for r in rs), sum(r["rels"] for r in rs),
                max(r["mem"] for r in rs))
    # pilot
    e, c, r, m = run(2, 0)
    est = max(e / 4 / 2, 1e-3)          # rough core-s per sq per process pair
    n1 = max(2, int(target / 4 / est))
    n2 = 3 * n1
    e1, c1, r1, m1 = run(n1, 100)
    e2, c2, r2, m2 = run(n2, 100)
    res = dict(core_s_per_sq=(e2 - e1) / (2 * (n2 - n1)), cpu_s_per_sq=(c2 - c1) / (2 * (n2 - n1)),
               rels_per_sq=r2 / (2 * n2), mem=max(m1, m2), n1=n1, n2=n2, startup=(e1 - (e2 - e1) / 2) / 4)
    T = load(); T[key] = res
    json.dump(T, open(TJ, "w"), indent=1)
    return res


if __name__ == "__main__":
    d, I, nfac, lpb = map(int, sys.argv[1:5])
    for lim in sys.argv[5].split(","):
        r = tsq(d, I, nfac, int(lim), lpb)
        print(d, I, nfac, lim, lpb, {k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()}, flush=True)
