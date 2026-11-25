#!/usr/bin/env python3
"""Filtering simulation: synthetic relations -> purge -> merge-dl.
   simulate(sample, R, lpb, T, shrink=1, density=170, keep=1) -> dict"""
import os, re, subprocess, sys, json, glob, shutil, time
# cado-nfs build tree: $CADO_BUILD, or the first directory under <source>/build/
import glob as _glob
_here = os.path.dirname(os.path.abspath(__file__))
CADO = os.environ.get("CADO_BUILD") or (_glob.glob(os.path.join(_here, "..", "..", "build", "*", "")) + [""])[0].rstrip("/")
HERE = os.path.dirname(os.path.abspath(__file__))
SIMDIR = os.environ.get("SIMDIR", HERE + "/sim")


def purge_merge(files, wdir, density=170, keep=1, threads=4, merge=True):
    res = {}
    purged = f"{wdir}/purged"
    t0 = time.time()
    r = subprocess.run([f"{CADO}/filter/purge", "-out", purged, "-keep", str(keep), "-t", str(threads),
                        "-required_excess", "0.0"] + files, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    res["t_purge"] = time.time() - t0
    txt = r.stdout
    m = re.search(r"Sing\. rem\.: begin with: nrows=(\d+) ncols=(\d+) excess=(-?\d+)", txt)
    if m:
        res["in_rows"], res["in_cols"], res["in_excess"] = map(int, m.groups())
    ms = re.findall(r"Sing\. rem\.:\s+iter \d+: nrows=(\d+) ncols=(\d+) excess=(-?\d+)", txt)
    # first singleton removal pass ends at the first "No more singletons"
    first = txt.split("No more singletons")[0]
    ms = re.findall(r"Sing\. rem\.:\s+iter \d+: nrows=(\d+) ncols=(\d+) excess=(-?\d+)", first)
    if ms:
        res["sr_rows"], res["sr_cols"], res["sr_excess"] = map(int, ms[-1])
    m = re.search(r"Final values:\s*\nnrows=(\d+) ncols=(\d+) excess=(-?\d+)", txt)
    if m:
        res["p_rows"], res["p_cols"], res["p_excess"] = map(int, m.groups())
    m = re.search(r"weight=(\d+) weight", txt)
    res["purge_rc"] = r.returncode
    if r.returncode != 0 or "p_rows" not in res:
        res["purge_tail"] = txt[-600:]
        return res
    if not merge:
        return res
    t0 = time.time()
    r = subprocess.run([f"{CADO}/filter/merge-dl", "-mat", purged, "-out", f"{wdir}/history",
                        "-skip", "0", "-target_density", str(density), "-t", str(threads)],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    res["t_merge"] = time.time() - t0
    m = re.search(r"Final matrix has N=(\d+) nc=(\d+) \(\d+\) W=(\d+)", r.stdout)
    if m:
        res["N"], res["nc"], res["W"] = map(int, m.groups())
    else:
        res["merge_tail"] = r.stdout[-600:]
    return res


def simulate(sample, R, lpb, T, shrink=1, density=170, keep=1, eps=0.05, seed=1, tag="x", merge=True,
             keepdir=False):
    wdir = f"{SIMDIR}/{tag}.{os.getpid()}"
    os.makedirs(wdir, exist_ok=True)
    nf = 4
    subprocess.run([f"{HERE}/fakeqs", sample, str(int(R / shrink)), str(lpb), str(T), str(eps),
                    str(shrink), str(seed), f"{wdir}/rels", str(nf)], check=True,
                   stderr=subprocess.DEVNULL)
    files = [f"{wdir}/rels.{i}" for i in range(nf)]
    res = purge_merge(files, wdir, density, keep, merge=merge)
    res.update(R=R, lpb=lpb, T=T, shrink=shrink, density=density)
    for k in ("in_rows", "in_cols", "sr_rows", "sr_cols", "sr_excess", "p_rows", "p_cols", "N", "W"):
        if k in res and k != "sr_excess":
            res[k + "_scaled"] = res[k] * shrink
    if not keepdir:
        shutil.rmtree(wdir)
    else:
        res["wdir"] = wdir
    return res


if __name__ == "__main__":
    sample, lpb, T, shrink = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    for R in sys.argv[5:]:
        r = simulate(sample, float(R), lpb, T, shrink)
        print(json.dumps(r), flush=True)


def make_matrix(wdir, nsq=2):
    """replay-dl on a kept simulation directory -> wdir/sparse.bin (+ sparse.1.bin)"""
    r = subprocess.run([f"{CADO}/filter/replay-dl", "-purged", f"{wdir}/purged", "-ideals", f"{wdir}/ideal",
                        "-his", f"{wdir}/history", "-index", f"{wdir}/index", "-out", f"{wdir}/sparse.bin",
                        "-nsquare_matrices", str(nsq)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode:
        raise RuntimeError(r.stdout[-500:])
    return f"{wdir}/sparse.bin"
