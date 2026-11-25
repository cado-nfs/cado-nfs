#!/usr/bin/env python3
"""Linear algebra cost model for class group computations (determinants by CRT)."""
import math


def la_cost(N, W, nbits=704, beta=1.95, neq=2, nmat=2):
    """core-seconds (this machine, all cores busy) for nmat determinants by CRT.
    krylov: 1.5*N iterations (m=2,n=1), c core-s per matrix coefficient per iteration."""
    c = {64: 2.8e-9, 256: 5.5e-9, 704: 10.4e-9, 960: 14.5e-9}[nbits]
    P = math.ceil((beta * N + 1) / nbits) + neq
    lingen = 3.6 * (nbits / 704) * (N / 4200) * (math.log(N) / math.log(4200)) ** 2
    return nmat * P * (1.5 * N * W * c + lingen), P


def pi_approx(x):
    l = math.log(x)
    return x / l * (1 + 1 / l + 2 / l**2 + 6 / l**3)
