#!/usr/bin/env python3
"""Smoothness model with up to k large primes (Dickman rho based)."""
import math, numpy as np

_H = 1e-3
_UMAX = 40.0
def _rho_table():
    n = int(_UMAX / _H) + 1
    r = np.ones(n)
    k1 = int(1 / _H)
    # u rho(u) = int_{u-1}^{u} rho(t) dt  ->  solve incrementally with trapezoid
    S = np.zeros(n)  # cumulative integral
    for i in range(1, k1 + 1):
        S[i] = S[i - 1] + _H
    for i in range(k1 + 1, n):
        u = i * _H
        # S[i] = S[i-1] + H/2 (r[i-1] + r[i]);  u r[i] = S[i] - S[i-k1]
        # => r[i] (u - H/2) = S[i-1] + H/2 r[i-1] - S[i-k1]
        r[i] = (S[i - 1] + _H / 2 * r[i - 1] - S[i - k1]) / (u - _H / 2)
        S[i] = S[i - 1] + _H / 2 * (r[i - 1] + r[i])
    return r
_RHO = _rho_table()
def rho(u):
    u = np.asarray(u, dtype=float)
    x = np.clip(u, 0, _UMAX - 2 * _H) / _H
    i = x.astype(int); f = x - i
    v = _RHO[i] * (1 - f) + _RHO[i + 1] * f
    return np.where(u <= 1, 1.0, v)

def psmooth(X, b1, b2, mfb=None, nlp=2, npts=80):
    """probability that an integer of natural-log size X is e^b1-smooth except for at most nlp
    primes in (e^b1, e^b2] whose product is < e^mfb. Returns (P0, P1, P2)."""
    if mfb is None:
        mfb = nlp * b2
    P0 = float(rho(X / b1))
    if b2 <= b1 or nlp == 0:
        return P0, 0.0, 0.0
    t = np.linspace(b1, b2, npts + 1)
    tm = (t[:-1] + t[1:]) / 2; dt = t[1] - t[0]
    P1 = float(np.sum(rho((X - tm) / b1) / tm) * dt)
    P2 = 0.0
    if nlp >= 2:
        T1, T2 = np.meshgrid(tm, tm)
        ok = (T1 + T2) <= mfb
        P2 = 0.5 * float(np.sum(np.where(ok, rho((X - T1 - T2) / b1) / (T1 * T2), 0.0)) * dt * dt)
    return P0, P1, P2
if __name__ == "__main__":
    for u in (1, 2, 3, 5, 10):
        print(u, float(rho(u)))
    print(psmooth(math.log(2.0**180), math.log(8e5), 24 * math.log(2)))
