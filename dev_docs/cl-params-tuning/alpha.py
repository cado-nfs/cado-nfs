#!/usr/bin/env python3
"""alpha(D): average log-size correction of the values X^2-D compared to random integers
(negative = more small prime factors than random = higher yield)."""
import math, sys
from qsbench import discriminant, legendre
def primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\0\0"
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
PR = primes(100000)
def alpha(D, B=100000):
    a = 0.0
    for p in PR:
        if p > B: break
        if p == 2:
            # average 2-adic valuation of X^2-D over X (D odd, D=1 mod 4 here)
            m = D % 8
            v = {1: 2.0, 5: 1.0}.get(m, 0.5 if D % 4 == 3 else 1.0)   # E[v2]: D=1 mod 8: 1/2*(3+...)=2 ; 5 mod 8: 1/2*2=1
            a += (1.0 - v) * math.log(2)
            continue
        if D % p == 0:
            chi = 0
        else:
            chi = 1 if legendre(D, p) == 1 else -1
        # expected valuation: random integer 1/(p-1); X^2-D: (1+chi) * p/(p^2-1)*... ~ (1+chi)/(p-1) for chi=1
        ev = (2.0 / (p - 1)) if chi == 1 else (1.0 / p if chi == 0 else 0.0)
        a += (1.0 / (p - 1) - ev) * math.log(p)
    return a
if __name__ == "__main__":
    for d in (60, 80, 90, 100, 110, 120, 130, 140, 150, 160):
        D = discriminant(d)
        sp = [p for p in PR[:12] if p > 2 and legendre(D, p) == 1]
        print(f"d{d}: D mod 8 = {D % 8}, alpha = {alpha(D):+.2f}, small split primes: {sp}")
    import random
    random.seed(1)
    al = sorted(alpha(discriminant(100, k)) for k in range(1, 201))
    print("distribution over 200 random 100-digit -p (p=3 mod 4): min %.2f 10%% %.2f median %.2f 90%% %.2f max %.2f mean %.2f" % (al[0], al[20], al[100], al[180], al[-1], sum(al)/len(al)))
