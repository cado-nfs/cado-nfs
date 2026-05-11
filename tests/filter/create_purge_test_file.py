#!/usr/bin/env python3

import sys
from random import randrange, shuffle


def usage():
    print("""
Usage: $0 <N> <ncols_dense> <nrows_base> [<ccsize0> <ccsize1> ...]

This script creates a set of dummy relations for testing filter/purge.

A dense part consisting of ncols_dense columns is created, with indices
between 0 and N-1. nrows_base such "plain rows" are created. Then we also
create rows, modeled in the same way, that form connected components of
sizes ccsize0, ccsize1, etc. The resulting set of relations is then shuffled.
""".strip())


def oneline(N, nj):
    ab = f"{randrange(1024):x},{randrange(1024):x}"
    sep = ':'
    for c in range(nj):
        ab += sep + f"{randrange(N):x}"
        sep = ','
    return ab


def onecc(N, nj, base, ccsize):
    i = base
    C = []
    for k in range(ccsize - 1):
        j = i + 1
        C.append(f"{oneline(N, nj)},{i:x},{j:x}")
        i = j
    C.append(f"{oneline(N, nj)},{i:x},{base:x}")
    return C


def buildmatrix(N, nj, nr, *ccs):
    if not N or not nj or not nr:
        usage
        raise RuntimeError("missing argument")

    L = [oneline(N, nj) for i in range(nr)]
    b = N
    for c in ccs:
        L += onecc(N, nj, b, c)
        b += c
    shuffle(L)
    for r in L:
        print(r)


buildmatrix(*[int(c) for c in sys.argv[1:]])
