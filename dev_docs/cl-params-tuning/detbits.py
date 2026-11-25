#!/usr/bin/env python3
"""log2|det| of a cado sparse matrix with coefficients (uint32 n; n x (uint32 col, int32 coeff))"""
import sys, numpy as np
def load(fn):
    a = np.fromfile(fn, dtype=np.int32)
    rows = []; i = 0
    while i < len(a):
        n = a[i]; rows.append(a[i + 1:i + 1 + 2 * n].reshape(n, 2)); i += 1 + 2 * n
    return rows
def detbits(fn):
    rows = load(fn)
    N = len(rows)
    nc = max(int(r[:, 0].max()) for r in rows if len(r)) + 1
    M = np.zeros((N, max(N, nc)))
    W = 0
    for i, r in enumerate(rows):
        np.add.at(M[i], r[:, 0], r[:, 1]); W += len(r)
    M = M[:, :N]
    s, l = np.linalg.slogdet(M)
    had = sum(0.5 * np.log2((r[:, 1].astype(float) ** 2).sum()) for r in rows)
    return N, nc, W, l / np.log(2), had, s
if __name__ == "__main__":
    for fn in sys.argv[1:]:
        N, nc, W, b, had, s = detbits(fn)
        print(f"{fn.split('/')[-1]}: N={N} ncols={nc} w={W/N:.1f} log2|det|={b:.0f} = {b/N:.3f} bits/row; Hadamard {had/N:.3f} bits/row sign={s}")
