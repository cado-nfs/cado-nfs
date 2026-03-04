#!/usr/bin/env python3

import enum
import math
import random

# Compatibility matrix between Algorithm and Computation
#         NFS    QS
# FACT    ✅     ✅
# DLP     ✅     ❌
# CL      ❌     ✅


class Algorithm(enum.Enum):
    """
    Enum for type of algorithm that can be used in cado-nfs.

    Two possible values:
        - NFS, for the Number Field Sieve;
        - QS, for the Quadratic Sieve.

    >>> Algorithm("NFS")
    <Algorithm.NFS: 'NFS'>

    >>> Algorithm("QS")
    <Algorithm.QS: 'QS'>

    >>> Algorithm("donotexist")
    Traceback (most recent call last):
        ...
    ValueError: 'donotexist' is not a valid Algorithm

    >>> Algorithm("NFS").name
    'NFS'

    >>> Algorithm("NFS") == "NFS", "NFS" == Algorithm("NFS")
    (True, True)

    >>> Algorithm("QS") == "QS", "QS" == Algorithm("QS")
    (True, True)

    >>> Algorithm("QS") == Algorithm("NFS")
    False

    >>> Algorithm("QS") == 'donotexist', Algorithm("QS") != 'donotexist'
    (False, True)

    >>> Algorithm("QS") == "NFS", "QS" == Algorithm("NFS")
    (False, False)
    """
    NFS = "NFS"
    QS = "QS"

    @classmethod
    def _missing_(cls, name):
        name = name.upper()
        for member in cls:
            if member.name == name:
                return member
        return None

    def __eq__(self, other):
        try:
            return self.name == type(self)(other).name
        except ValueError:
            return NotImplemented


class Computation(enum.Enum):
    """
    Enum for type of computation that can be performed in cado-nfs.

    Three possible values:
        - FACT, for integer factorization;
        - DLP, for finite field discrete logarithm;
        - CL, for class groups computation.

    >>> Computation("FACT")
    <Computation.FACT: ('Factorization', 'factor', 'c')>

    >>> Computation("DLP")
    <Computation.DLP: ('Discrete logarithm', 'dlp', 'p')>

    >>> Computation("donotexist")
    Traceback (most recent call last):
        ...
    ValueError: 'donotexist' is not a valid Computation

    >>> Computation("CL").name
    'CL'

    >>> Computation("FACT") == "FACT", "FACT" == Computation("FACT")
    (True, True)

    >>> Computation("DLP") == "DLP", "DLP" == Computation("DLP")
    (True, True)

    >>> Computation("CL") == Computation("FACT")
    False

    >>> Computation("CL") == "FACT", "FACT" == Computation("DLP")
    (False, False)

    >>> Computation("DLP") == 'donotexist', Computation("DLP") != 'donotexist'
    (False, True)
    """
    FACT = "Factorization", "factor", "c"
    DLP = "Discrete logarithm", "dlp", "p"
    CL = "Class group", "cl", "d"

    def __init__(self, desc, shortname, letter):
        self.desc = desc
        self.param_dirname = self.shortname = shortname
        self.letter = letter

    @classmethod
    def _missing_(cls, name):
        name = name.upper()
        for member in cls:
            if member.name == name:
                return member
        return None

    def __eq__(self, other):
        try:
            return self.name == type(self)(other).name
        except ValueError:
            return NotImplemented


def xgcd(a, b):
    """
    Compute the extended GCD of the two integers a and b.

    Return a tuple (d, u, v) such that d is the GCD of a and b, and u and v are
    such that u*a + v*b == d

    >>> xgcd(17, 42)
    (1, 5, -2)

    >>> xgcd(42, 17)
    (1, -2, 5)

    >>> xgcd(56, 44)
    (4, 4, -5)

    >>> xgcd(-56, 44)
    (4, -4, -5)

    >>> xgcd(56, -44)
    (4, 4, 5)

    >>> xgcd(-56, -44)
    (4, -4, 5)

    >>> xgcd(4, 8)
    (4, 1, 0)

    >>> xgcd(0, 17)
    (17, 0, 1)

    >>> xgcd(42, 0)
    (42, 1, 0)

    >>> xgcd(0, 0)
    (0, 0, 0)

    >>> xgcd(505753929367110, 391303222658950)
    (10, 8244650487681, -10656095168522)

    >>> import random
    >>> a, b = random.randrange(10**100), random.randrange(10**120)
    >>> d, u, v = xgcd(a, b)
    >>> d == u*a+b*v, d == math.gcd(a, b), a, b  # doctest: +ELLIPSIS
    (True, True, ...)
    """
    if a == 0 and b == 0:
        return (0, 0, 0)
    u0, v0, r0 = 1, 0, a
    u1, v1, r1 = 0, 1, b
    while r1 != 0:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        u0, u1 = u1, u0 - q * u1
        v0, v1 = v1, v0 - q * v1
    if r0 < 0:
        u0, v0, r0 = -u0, -v0, -r0
    return r0, u0, v0


def CRT(r0, m0, r1, m1, signed=False):
    """
    Return a tuple (r, m) such that
      - r = ri mod mi for i=0,1
      - m = m0*m1
    If signed is True, return r in [-m/2,m/2[, else return r in [0,m[.

    Assumes that m0 and m1 are coprime, a ValueError is raised otherwise.

    Compute CRT of 2 mod 3 and 1 mod 5
    >>> CRT(2, 3, 1, 5)
    (11, 15)
    >>> CRT(2, 3, 1, 5, signed=True)
    (-4, 15)

    Compute CRT of 13 mod 100 and 20 mod 301
    >>> CRT(13, 100, 20, 301)
    (28013, 30100)
    >>> CRT(13, 100, 20, 301, signed=True)
    (-2087, 30100)

    Compute CRT of 1 mod 5 and 2 mod 7
    >>> CRT(1, 5, 2, 7)
    (16, 35)
    >>> CRT(1, 5, 2, 7, signed=True)
    (16, 35)

    CRT of 4 mod 8 and 6 mod 12 should fail
    >>> CRT(4, 8, 6, 12)
    Traceback (most recent call last):
      ...
    ValueError: 8 and 12 are not coprime, gcd is 4
    """
    d, u0, u1 = xgcd(m0, m1)
    if d != 1:
        raise ValueError("%d and %d are not coprime, gcd is %d" % (m0, m1, d))
    m = m0*m1
    r = (r0*m1*u1 + r1*m0*u0) % m
    if signed:
        r = r if 2*r <= m else r-m  # return in [-m/2,m/2[
    return (r, m)


def is_prime(n, niter=10):
    """
    Miller-Rabin primality test.

    Return False if n is known to be composite, return True if n is probably
    prime.

    >>> [is_prime(m) for m in range(10)]
    [False, False, True, True, False, True, False, True, False, False]

    >>> is_prime(-3)
    False

    >>> len([_ for n in range(100) if is_prime(n)])
    25

    >>> len([_ for n in range(10000) if is_prime(n)])
    1229

    test all factors of the fixed bases for the 64-bit case
    >>> all(is_prime(p) for p in [2, 3, 5, 13, 19, 73, 193, 407521, 299210837])
    True

    64-bit primes generated with SageMath
    >>> all(is_prime(p) for p in [13456180279374568867, 4536562633469233429,
    ...     4147834561926919979, 7813982791570798891, 17366322772979549569])
    True

    64-bit composite generated with SageMath
    >>> any(is_prime(c) for c in [2463229809936173032, 3059833457856409124,
    ...     11375980981929811056, 11769466443373517304, 16138676638004345802])
    False
    """
    if n in (2, 3, 5, 7):
        return True
    elif n < 10:
        return False
    elif math.gcd(n, 210) != 1:  # 210 = 2*3*5*7
        return False

    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d >> 1, s + 1
    # d*2^s == n-1

    def is_a_witness(a):
        a = a % n
        if a in (0, 1, n-1):
            return False
        if (m := pow(a, d, n)) == 1:
            return False
        for i in range(s):
            if m == n-1:  # here m = a^(d*2^i) mod n
                return False
            m = (m*m) % n
        return True

    if n.bit_length() <= 64:
        # deterministically correct up to 2^64 according to
        # https://miller-rabin.appspot.com/
        bases = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
    else:
        bases = (random.randrange(2, n-1) for _ in range(niter))

    for a in bases:
        if is_a_witness(a):
            return False
    return True


def primes_above(start):
    """
    Return a generator object to iterate over all primes larger or equal to
    start.

    >>> import itertools, random
    >>> first_n = lambda gen, n: list(itertools.islice(gen, n))

    >>> first_n(primes_above(0), 10)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    >>> first_n(primes_above(-100), 10)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    >>> first_n(primes_above(2), 10)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

    >>> first_n(primes_above(1024), 5)
    [1031, 1033, 1039, 1049, 1051]

    >>> first_n(primes_above(2**63), 5)  # doctest: +NORMALIZE_WHITESPACE
    [9223372036854775837, 9223372036854775907, 9223372036854775931, \
    9223372036854775939, 9223372036854775963]

    >>> m = random.randrange(2**100)
    >>> all(is_prime(p) for p in first_n(primes_above(m), 5)), m
    ... # doctest: +ELLIPSIS
    (True, ...)
    """
    n = start
    if n <= 2:
        yield 2
        n = 3
    else:
        n += (n+1) % 2
    # n is odd
    while True:
        if is_prime(n):
            yield n
        n += 2


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        import doctest
        doctest.testmod()
