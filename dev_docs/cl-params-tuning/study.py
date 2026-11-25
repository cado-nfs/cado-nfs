#!/usr/bin/env python3
"""Per-size study helpers.
  sample(d, I, nfac, lim, lpb, nq, ...): sieve the first nq special-q's with 2 procs x 2 threads,
      store unique relations in samples/d{d}.I{I}.n{nfac}.lim{lim}.lpb{lpb}.rels + .json (timing)
  subsample(file, lim, lpb, mfb): relations that would also be found with smaller bounds
"""
import os, re, subprocess, sys, json, math, time
from concurrent.futures import ThreadPoolExecutor
from qsbench import *

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = HERE + "/samples"
os.makedirs(SAMPLES, exist_ok=True)


def sample(d, I, nfac, lim, lpb, nq, mfb=None, ncurves=20, nproc=2, thr=2, k=0, tag="", idx0=0,
           qmin=None, qmax=None):
    D = discriminant(d, k)
    if qmin is None:
        qmin, qmax, pool = choose_qfac(D, I, nfac)
    mfb = mfb or 2 * lpb
    base = f"{SAMPLES}/d{d}.I{I}.n{nfac}.lim{lim}.lpb{lpb}.mfb{mfb}{tag}"
    if os.path.exists(base + ".json"):
        return json.load(open(base + ".json"))
    per = nq // nproc

    def job(i):
        return run_siqs(D, I, lim, lpb, mfb, ncurves, nfac, qmin, qmax, nq=per, idx0=idx0 + i * per,
                        threads=thr, keep=f"{base}.part{i}")
    t0 = time.time()
    with ThreadPoolExecutor(nproc) as ex:
        rs = list(ex.map(job, range(nproc)))
    wct = time.time() - t0
    for r in rs:
        if r["rc"]:
            raise RuntimeError(r["err"])
    seen = set()
    nraw = 0
    with open(base + ".rels", "w") as out:
        for i in range(nproc):
            with open(f"{base}.part{i}") as f:
                for line in f:
                    if line[0] == "#":
                        continue
                    nraw += 1
                    key = line.split(":", 1)[0]
                    if key in seen:
                        continue
                    seen.add(key)
                    out.write(line)
            os.unlink(f"{base}.part{i}")
    res = dict(d=d, D=str(D), I=I, nfac=nfac, lim=lim, lpb=lpb, mfb=mfb, ncurves=ncurves, nq=per * nproc,
               qmin=qmin, qmax=qmax, raw=nraw, uniq=len(seen), cpu=sum(r["cpu"] for r in rs),
               elapsed=sum(r["elapsed"] for r in rs), thr=thr, nproc=nproc, wct=wct,
               mem=max(r["mem"] for r in rs), file=base + ".rels")
    res["core_s_per_sq"] = res["elapsed"] * thr / res["nq"]   # core occupancy
    res["cpu_s_per_sq"] = res["cpu"] / res["nq"]
    json.dump(res, open(base + ".json", "w"), indent=1)
    return res


def parse_rel(line):
    ps = line.rstrip("\n").split(":")[1].split(",")
    return [int(x, 16) for x in ps if x]


def subsample(fn, lim, lpb, mfb=None, out=None, qmax=0):
    """relations of fn that have at most mfb bits of primes >= lim, each < 2^lpb, and at most
    mfb/lpb such primes. Returns count (and writes them to out)."""
    mfb = mfb or 2 * lpb
    L = 1 << lpb
    nlp_max = round(mfb / lpb)
    n = 0
    o = open(out, "w") if out else None
    with open(fn) as f:
        for line in f:
            ps = parse_rel(line)
            big = [p for p in ps if p >= lim]
            if len(big) > nlp_max or (big and max(big) >= L):
                continue
            if big and math.prod(big) >> mfb:
                continue
            n += 1
            if o:
                o.write(line)
    if o:
        o.close()
    return n


def timing(d, I, nfac, lim, lpb, nq=8, mfb=None, ncurves=20, k=0, idx0=0, qmin=None, qmax=None):
    """core-seconds per special-q (2 procs x 2 threads, all 4 cores busy)"""
    D = discriminant(d, k)
    if qmin is None:
        qmin, qmax, pool = choose_qfac(D, I, nfac)
    mfb = mfb or 2 * lpb
    per = max(1, nq // 2)

    def job(i):
        return run_siqs(D, I, lim, lpb, mfb, ncurves, nfac, qmin, qmax, nq=per, idx0=idx0 + i * per, threads=2)
    with ThreadPoolExecutor(2) as ex:
        rs = list(ex.map(job, range(2)))
    return dict(core_s_per_sq=sum(r["elapsed"] for r in rs) * 2 / (2 * per),
                cpu_s_per_sq=sum(r["cpu"] for r in rs) / (2 * per),
                rels_per_sq=sum(r["rels"] for r in rs) / (2 * per), mem=max(r["mem"] for r in rs))


if __name__ == "__main__":
    kw = dict(a.split("=") for a in sys.argv[2:])
    kw = {k: int(v) for k, v in kw.items()}
    if sys.argv[1] == "sample":
        r = sample(**kw)
        print({k: r[k] for k in ("d", "I", "nfac", "lim", "lpb", "nq", "raw", "uniq", "cpu", "elapsed",
                                 "core_s_per_sq", "cpu_s_per_sq", "mem", "qmin", "qmax")})
    elif sys.argv[1] == "timing":
        print(timing(**kw))
