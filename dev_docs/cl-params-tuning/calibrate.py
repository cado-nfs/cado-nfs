#!/usr/bin/env python3
"""Calibrate the timing estimates on another machine (e.g. a grdix or grvingt node).

Loads the node the way a production run would (one 2-thread process per two
cores) and measures
  * the d120 sieve: core-seconds per special-q, start-up excluded
    (reference machine: 5.36);
  * the 704-bit Krylov iteration on a random 40000 x 40000 matrix:
    core-nanoseconds per matrix coefficient (reference machine: 12.2 with
    this matrix; real merged matrices run at 10.4 there).
The two ratios are the factors to apply to the "sieving" and "lin. alg."
columns of README.md.

    export CADO_BUILD=/path/to/cado-nfs/build/<host>
    python3 calibrate.py [--cores N] [--secs 40] [--skip-la] [--skip-sieve]

Memory: about 2.8 GB per sieving process (so 1.4 GB per core).
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from qsbench import discriminant, choose_qfac, run_siqs
import labench2

HERE = os.path.dirname(os.path.abspath(__file__))
REF_SIEVE = 5.36      # core-s per special-q, d120 file, reference machine
REF_LA = 12.2         # core-ns per coefficient per iteration, 704 bits, this random matrix
D120 = dict(d=120, I=15, nfac=16, lim=3200000, lpb=29)


def sieve(nproc, secs):
    p = D120
    D = discriminant(p["d"])
    qmin, qmax, _ = choose_qfac(D, p["I"], p["nfac"])

    def run(n, off):
        def job(i):
            return run_siqs(D, p["I"], p["lim"], p["lpb"], 2 * p["lpb"], 20, p["nfac"], qmin, qmax,
                            nq=n, idx0=off + i * n, threads=2)
        with ThreadPoolExecutor(nproc) as ex:
            rs = list(ex.map(job, range(nproc)))
        bad = [r for r in rs if r["rc"]]
        if bad:
            sys.exit("siqs failed:\n" + bad[0]["err"])
        return sum(r["elapsed"] for r in rs) * 2, sum(r["cpu"] for r in rs), sum(r["rels"] for r in rs)

    # the factor base file is created by the first call; do it before timing
    run_siqs(D, p["I"], p["lim"], p["lpb"], 2 * p["lpb"], 20, p["nfac"], qmin, qmax, nq=1, threads=2)
    n1 = max(2, int(secs / REF_SIEVE / 2))
    n2 = 3 * n1
    e1, c1, _ = run(n1, 1000)
    e2, c2, r2 = run(n2, 1000)
    core = (e2 - e1) / (nproc * (n2 - n1))
    cpu = (c2 - c1) / (nproc * (n2 - n1))
    return dict(core_s_per_sq=core, cpu_s_per_sq=cpu, rels_per_sq=r2 / (nproc * n2), nproc=nproc,
                special_q=nproc * (n1 + n2))


def linalg(nproc, secs, N=40000, w=170):
    mat = f"{HERE}/la/calib{N}.bin"
    os.makedirs(os.path.dirname(mat), exist_ok=True)
    if not os.path.exists(mat):
        subprocess.run([sys.executable, f"{HERE}/randmat.py", str(N), str(w), mat], check=True,
                       stdout=subprocess.DEVNULL)
    import numpy as np
    nnz = int(np.fromfile(mat[:-4] + ".rw.bin", dtype=np.uint32).sum())
    # two runs of different lengths: the difference removes the start-up (balancing, ...)
    def iters(t):
        rs = labench2.bench(mat, 704, t, nproc, 2)
        ok = [r for r in rs if r.get("iters_done")]
        if not ok:
            sys.exit("bwc failed:\n" + str(rs[0].get("tail")))
        return ok
    a = iters(secs / 2)
    b = iters(secs / 2 + secs)
    d_it = sorted(r["iters_done"] for r in b)[len(b) // 2] - sorted(r["iters_done"] for r in a)[len(a) // 2]
    if d_it <= 0:
        sys.exit("krylov made no progress; increase --la-secs")
    return dict(core_ns_per_coeff=2 * secs / d_it / nnz * 1e9,
                thread_ns_item=sorted(r["ns_item"] for r in b)[len(b) // 2], nproc=len(b), N=N, nnz=nnz)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cores", type=int, default=os.cpu_count())
    ap.add_argument("--secs", type=float, default=60.0, help="target duration of the short sieving run")
    ap.add_argument("--la-secs", type=float, default=60.0, help="krylov is timed over about this long")
    ap.add_argument("--skip-la", action="store_true")
    ap.add_argument("--skip-sieve", action="store_true")
    a = ap.parse_args()
    nproc = max(1, a.cores // 2)
    out = {"host": os.uname().nodename, "cores_used": 2 * nproc}
    if not a.skip_sieve:
        t0 = time.time()
        s = sieve(nproc, a.secs)
        out["sieve"] = s
        print(f"sieve d120: {s['core_s_per_sq']:.2f} core-s per special-q ({s['cpu_s_per_sq']:.2f} cpu), "
              f"{s['rels_per_sq']:.1f} relations per special-q, {nproc} processes x 2 threads, "
              f"{time.time() - t0:.0f} s  ->  factor {s['core_s_per_sq'] / REF_SIEVE:.2f} on the sieving column")
    if not a.skip_la:
        l = linalg(nproc, a.la_secs)
        out["linalg"] = l
        print(f"krylov 704 bits, N={l['N']}: {l['core_ns_per_coeff']:.1f} core-ns per coefficient "
              f"(thread time {l['thread_ns_item']:.1f} ns), {l['nproc']} processes x 2 threads"
              f"  ->  factor {l['core_ns_per_coeff'] / REF_LA:.2f} on the lin. alg. column")
    json.dump(out, open(f"{HERE}/calibration.{out['host']}.json", "w"), indent=1)
