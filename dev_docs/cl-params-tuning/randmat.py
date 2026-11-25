#!/usr/bin/env python3
"""Random sparse matrix in cado bwc format (with coefficients), structure similar to a merged
relation matrix: N x N, row weight ~w, column j chosen with density ~ 1/(j+c0).
   randmat.py N w out.bin"""
import sys, numpy as np
N, w, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(1)
c0 = N / 200
rw = np.zeros(N, dtype=np.uint32); cw = np.zeros(N, dtype=np.uint32)
with open(out, "wb") as f:
    for i in range(N):
        u = rng.random(int(w * 1.15))
        cols = np.unique((c0 * ((1 + N / c0) ** u - 1)).astype(np.int64).clip(0, N - 1))
        cols[0] = i if i not in cols else cols[0]   # ensure non-singular-ish diagonal presence
        cols = np.unique(cols)
        e = rng.choice(np.array([-1, 1, -1, 1, 2, -2], dtype=np.int32), size=len(cols))
        row = np.empty(1 + 2 * len(cols), dtype=np.int32)
        row[0] = len(cols); row[1::2] = cols; row[2::2] = e
        row.tofile(f)
        rw[i] = len(cols); np.add.at(cw, cols, 1)
rw.tofile(out[:-4] + ".rw.bin"); cw.tofile(out[:-4] + ".cw.bin")
print("N", N, "avg w", rw.mean())
