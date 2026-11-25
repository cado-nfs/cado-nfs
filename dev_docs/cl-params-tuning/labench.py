#!/usr/bin/env python3
"""Benchmark one determinant computation (bwc.pl :determinant) on a matrix.
   labench.py matrix.bin nbits [m] [n] [thr]"""
import os, re, subprocess, sys, json, shutil, time
from qsbench import is_prime
# cado-nfs build tree: $CADO_BUILD, or the first directory under <source>/build/
import glob as _glob
_here = os.path.dirname(os.path.abspath(__file__))
CADO = os.environ.get("CADO_BUILD") or (_glob.glob(os.path.join(_here, "..", "..", "build", "*", "")) + [""])[0].rstrip("/")
HERE = os.path.dirname(os.path.abspath(__file__))


def next_prime(n):
    n += 1 + (n % 2 == 0) * 0
    if n % 2 == 0:
        n += 1
    while not is_prime(n):
        n += 2
    return n


def bench(matrix, nbits, m=2, n=1, thr=4, keep=False, interval=None):
    p = next_prime(2 ** (nbits - 1))
    wdir = f"{HERE}/la/{os.path.basename(matrix)}.{nbits}.{m}.{n}.{os.getpid()}"
    os.makedirs(os.path.dirname(wdir), exist_ok=True)
    cmd = [f"{CADO}/linalg/bwc/bwc.pl", ":determinant", f"lingen_thr={thr}", f"thr={thr}",
           f"m={m}", f"n={n}", "nullspace=right", f"matrix={matrix}", f"prime={p}", f"wdir={wdir}"]
    if interval:
        cmd.append(f"interval={interval}")
    t0 = time.time()
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    res = {"wct": time.time() - t0, "rc": r.returncode, "nbits": nbits, "m": m, "n": n, "thr": thr}
    txt = r.stdout
    for step in ("prep", "secure", "krylov", "lingen_pz", "acollect", "det_from_lingen"):
        for t in ("wct", "cpu"):
            mm = re.search(rf"Timings for {step}: \({t}\) ([\d.]+)", txt)
            if mm:
                res[f"{step}_{t}"] = float(mm.group(1))
    mm = re.search(r"Target iteration is (\d+)", txt)
    if mm:
        res["iters"] = int(mm.group(1))
    ms = re.findall(r"CPU0: ([\d.]+) s/iter \[[\d.]+\.\.[\d.]+\], ([\d.]+) ns", txt)
    if ms:
        res["s_iter"] = float(ms[-1][0]); res["ns_item"] = float(ms[-1][1])
    mm = re.search(r"determinant modulo p: (\S+)", txt)
    res["det"] = mm.group(1)[:20] if mm else None
    if r.returncode != 0 or not mm:
        res["tail"] = txt[-1500:]
    if not keep:
        shutil.rmtree(wdir, ignore_errors=True)
    return res


if __name__ == "__main__":
    a = sys.argv
    r = bench(a[1], int(a[2]), int(a[3]) if len(a) > 3 else 2, int(a[4]) if len(a) > 4 else 1,
              int(a[5]) if len(a) > 5 else 4)
    print(json.dumps(r))
