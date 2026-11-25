#!/usr/bin/env python3
"""Krylov speed under load: nproc concurrent bwc determinant runs (thr threads each), killed after
`secs`; reports seconds/iteration and thread-ns per matrix coefficient.
   labench2.py matrix.bin nbits [secs=90] [nproc=2] [thr=2]"""
import os, re, subprocess, sys, json, shutil, time, signal
from labench import next_prime, CADO, HERE
def bench(matrix, nbits, secs=90, nproc=2, thr=2):
    procs = []
    p = 2 ** (nbits - 1)
    for i in range(nproc):
        p = next_prime(p)
        wdir = f"{HERE}/la/{os.path.basename(matrix)}.{nbits}.b{i}.{os.getpid()}"
        os.makedirs(os.path.dirname(wdir), exist_ok=True)
        cmd = [f"{CADO}/linalg/bwc/bwc.pl", ":determinant", f"lingen_thr={thr}", f"thr={thr}", "m=2", "n=1",
               "nullspace=right", f"matrix={matrix}", f"prime={p}", f"wdir={wdir}"]
        out = open(wdir + ".log", "w")
        procs.append((subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, start_new_session=True), wdir, out))
    time.sleep(secs)
    res = []
    for pr, wdir, out in procs:
        try:
            os.killpg(pr.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        pr.wait(); out.close()
        txt = open(wdir + ".log").read()
        ms = re.findall(r"N=(\d+) ; CPU0: ([\d.]+) s/iter \[[\d.]+\.\.[\d.]+\], ([\d.]+) ns", txt)
        ws = re.findall(r"krylov: N=\d+ ; ETA \(N=(\d+)\): .* \[([\d.]+) s/iter\]", txt)
        if ms:
            res.append(dict(iters_done=int(ms[-1][0]), cpu_s_iter=float(ms[-1][1]), ns_item=float(ms[-1][2]),
                            wct_s_iter=float(ws[-1][1]) if ws else None, target=int(ws[-1][0]) if ws else None))
        else:
            res.append(dict(tail=txt[-400:]))
        shutil.rmtree(wdir, ignore_errors=True); os.unlink(wdir + ".log")
    return res
if __name__ == "__main__":
    a = sys.argv
    r = bench(a[1], int(a[2]), *(int(x) for x in a[3:]))
    print(json.dumps(r))
