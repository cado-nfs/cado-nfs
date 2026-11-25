#!/usr/bin/env python3
"""Sampling harness for the cado-nfs siqs siever (class group computations).

Runs siqs on a random sample of special-q's and reports yield / timing.
"""
import os, re, subprocess, sys, json, math, random, hashlib, gzip

# cado-nfs build tree: $CADO_BUILD, or the first directory under <source>/build/
import glob as _glob
_here = os.path.dirname(os.path.abspath(__file__))
CADO = os.environ.get("CADO_BUILD") or (_glob.glob(os.path.join(_here, "..", "..", "build", "*", "")) + [""])[0].rstrip("/")
WDIR = os.environ.get("QSBENCH_WDIR", os.path.dirname(os.path.abspath(__file__)) + "/work")
os.makedirs(WDIR, exist_ok=True)


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def discriminant(digits, k=0):
    """k-th test discriminant with the given number of digits: -p, p = 3 mod 4
    prime, p pseudo-random of exactly `digits` digits (deterministic)."""
    r = int.from_bytes(hashlib.shake_256(f"cado-cl-d{digits}-{k}".encode()).digest(128), "big")
    lo, hi = 10 ** (digits - 1), 10 ** digits
    # leading part uniform in log scale would bias; take uniform in [2,8]*10^(d-1)
    p = 2 * lo + r % (6 * lo)
    p -= p % 4
    p += 3
    while not is_prime(p):
        p += 4
    return -p


def polyfile(D):
    fn = f"{WDIR}/d{len(str(-D))}_{-D % 10**8}.poly"
    if not os.path.exists(fn):
        with open(fn, "w") as f:
            f.write(f"n: {D}\nskew: {math.isqrt(-D)}\npoly0: {-D},0,1\n")
    return fn


def fbfile(D, lim, I):
    poly = polyfile(D)
    fn = poly[:-5] + f".lim{lim}.I{I}.roots.gz"
    if not os.path.exists(fn):
        subprocess.run([f"{CADO}/sieve/makefb", "-poly", poly, "-lim", str(lim),
                        "-maxbits", str(I), "-out", fn, "-t", "2"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return fn


def legendre(a, p):
    return pow(a % p, (p - 1) // 2, p)


def qfac_pool(D, qmin, qmax):
    return [p for p in range(qmin | 1, qmax, 2) if is_prime(p) and legendre(D, p) == 1]


def choose_qfac(D, I, nfac, npool=None, spread=None):
    """Choose [qfac_min, qfac_max) so that the product of nfac primes of the
    pool is close to the optimum sqrt(|D|/2)/2^(I-1), with a pool of about
    npool split primes (default: enough for many special-q's)."""
    target = math.sqrt(-D / 2) / 2 ** (I - 1)
    g = target ** (1.0 / nfac)
    if npool is None:
        npool = nfac + max(12, nfac)
    # take the npool split primes closest (in log scale) to g, arranged so
    # that the first special-q's (smallest primes) are slightly below target
    lo = int(g)
    primes = []
    p = lo | 1
    below, above = [], []
    q = p
    while len(above) < npool:
        if is_prime(q) and legendre(D, q) == 1:
            above.append(q)
        q += 2
    q = p - 2
    while len(below) < npool and q > 2:
        if is_prime(q) and legendre(D, q) == 1:
            below.append(q)
        q -= 2
    # the special-q's are enumerated starting from the smallest primes; we
    # want the nfac smallest to have geometric mean slightly below g.
    nb = min(len(below), nfac // 2 + 1)
    pool = sorted(below[:nb]) + above[:npool - nb]
    return pool[0], pool[-1] + 1, pool


def run_siqs(D, I, lim, lpb, mfb, ncurves, nfac, qmin, qmax, nq=20, idx0=0, idx1=None,
             seed=1, threads=1, extra=(), keep=None, timeout=None, lam=None):
    poly = polyfile(D)
    fb = fbfile(D, lim, I)
    cmd = [f"{CADO}/sieve/siqs", "-poly", poly, "-I", str(I), "-lim0", str(lim),
           "-lpb0", str(lpb), "-mfb0", str(mfb), "-ncurves0", str(ncurves),
           "-fb0", fb, "-qfac-min", str(qmin), "-qfac-max", str(qmax),
           "-qfac-nfac", str(nfac), "-galois", "_y", "-t", str(threads),
           "-qidx0", str(idx0)]
    if lam is not None:
        cmd += ["-lambda0", str(lam)]
    if idx1 is None:
        cmd += ["-nq", str(nq)]
    else:
        cmd += ["-qidx1", str(idx1), "-random-sample", str(nq), "-seed", str(seed)]
    cmd += list(extra)
    out = keep or f"{WDIR}/tmp.{os.getpid()}.out"
    cmd += ["-out", out]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                       timeout=timeout)
    res = {"cmd": " ".join(cmd), "rc": r.returncode}
    if r.returncode != 0:
        res["err"] = r.stdout[-2000:]
        return res
    txt = open(out).read()
    if keep is None:
        os.unlink(out)
    m = re.search(r"# Total cpu time ([\d.]+)s", txt)
    res["cpu"] = float(m.group(1))
    m = re.search(r"# Total elapsed time ([\d.]+)s", txt)
    res["elapsed"] = float(m.group(1))
    m = re.search(r"# Total (\d+) reports \[([\d.e+-]+)s/r, ([\d.]+)r/sq\]", txt)
    res["rels"] = int(m.group(1))
    m = re.search(r"# Discarded (\d+) special-q's out of (\d+) pushed", txt)
    res["nq"] = int(m.group(2)) - int(m.group(1))
    m = re.search(r"# PeakMemusage \(MB\) = (\d+)", txt)
    res["mem"] = int(m.group(1)) if m else None
    for key in ("after_sieve", "enter_cofactoring", "smooth"):
        ms = re.findall(r"# survivors %s: (\d+)" % key, txt)
        if ms:
            res[key] = int(ms[-1])
    m = re.search(r"# Total cpu time [\d.]+s \[norm ([\d.]+)\+([\d.]+), sieving ([\d.]+) "
                  r"\(([\d.]+)\+([\d.]+) \+ ([\d.]+)\), factor ([\d.]+) \(([\d.]+)\+([\d.]+) \+ ([\d.]+)\), rest ([\d.]+)\]", txt)
    if m:
        res["t_sieve"] = float(m.group(3))
        res["t_factor"] = float(m.group(7))
    res["rels_per_sq"] = res["rels"] / max(1, res["nq"])
    res["s_per_sq"] = res["cpu"] / max(1, res["nq"])
    res["s_per_rel"] = res["cpu"] / max(1, res["rels"])
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("digits", type=int)
    ap.add_argument("-I", type=int, default=14)
    ap.add_argument("-lim", type=int, required=True)
    ap.add_argument("-lpb", type=int, required=True)
    ap.add_argument("-mfb", type=int)
    ap.add_argument("-ncurves", type=int, default=10)
    ap.add_argument("-nfac", type=int, required=True)
    ap.add_argument("-qmin", type=int)
    ap.add_argument("-qmax", type=int)
    ap.add_argument("-nq", type=int, default=10)
    ap.add_argument("-idx1", type=int)
    ap.add_argument("-k", type=int, default=0)
    ap.add_argument("-t", type=int, default=1)
    ap.add_argument("-keep")
    ap.add_argument("extra", nargs="*")
    a = ap.parse_args()
    D = discriminant(a.digits, a.k)
    if a.qmin is None:
        a.qmin, a.qmax, pool = choose_qfac(D, a.I, a.nfac)
    r = run_siqs(D, a.I, a.lim, a.lpb, a.mfb or 2 * a.lpb, a.ncurves, a.nfac, a.qmin, a.qmax,
                 nq=a.nq, idx1=a.idx1, threads=a.t, keep=a.keep, extra=sum([("-" + e).split() for e in a.extra if e], []))
    r["D"] = D
    r["qmin"], r["qmax"] = a.qmin, a.qmax
    print(json.dumps(r, indent=1))
