import os
import re
import sys

from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.rings.number_field.number_field import NumberField
from .tools import get_verbose
from .tools import NOK
from .number_theory import CadoNumberTheory


class CadoPolyFile(object):
    """
    This class contains stuff to read a polynomial pair from a cado .poly
    file. Some effort goes into adapting to the various ways a poly file
    can be laid out, but these are certainly a few gaps.

    EXAMPLES:

    sage: from cado_sage import CadoPolyFile
    sage: data = '''
    ....: n: 123412341234618726123412341234611
    ....: skew: 418.982
    ....: c0: 1312885934640
    ....: c1: -181216146
    ....: c2: -19457299
    ....: c3: 280
    ....: Y0: -4970288929
    ....: Y1: 19570577
    ....: # MurphyE (Bf=1.311e+05,Bg=1.311e+05,area=6.554e+09) = 2.440e-04
    ....: # f(x) = 280*x^3-19457299*x^2-181216146*x+1312885934640
    ....: # g(x) = 19570577*x-4970288929
    ....: '''
    sage: p = CadoPolyFile()
    sage: p._CadoPolyFile__read(data.split('\n'))
    sage: print(p)
    NFS configuration for modulus 123412341234618726123412341234611
      K0: Number Field in alpha0 with defining polynomial
      19570577*x - 4970288929
      K1: Number Field in alpha1 with defining polynomial
      280*x^3 - 19457299*x^2 - 181216146*x + 1312885934640


    """
    def __init__(self,
                 filename=None,
                 wdir=None):

        self.filename = filename
        if wdir is None and filename is not None:
            self.wdir = os.path.dirname(self.filename)
        else:
            self.wdir = wdir

        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self.K = []
        self.f = []
        self.skewness = 0
        self.N = 0
        self.nt = None

    def read(self):
        try:
            if get_verbose():
                print(f"Reading {self.filename}")
            with open(self.filename, "r") as fm:
                self.__read(fm)

        except Exception as e:
            # We're really in bad shape if an exception occurs here.
            # We're not even trying to salvage the BwcMatrix object, as
            # the error is most probably obvious.
            print(f"Exception while reading {self.filename} {NOK}",
                  file=sys.stderr)
            raise e

    def __parseline(self, line):
        """

        EXAMPLES::

            sage: from cado_sage import CadoPolyFile
            sage: p = CadoPolyFile()
            sage: p._CadoPolyFile__parseline("n=-47")
            sage: p.N
            -47

            sage: p._CadoPolyFile__parseline("poly1=280*x^3-19*x^2-18*x+13")
            sage: p.f[1]
            280*x^3 - 19*x^2 - 18*x + 13

        """

        ZP = ZZ['x']
        x = ZP.gen()

        # regexps for integers, real numbers, equal sign, and lists of
        # coefficients.
        nb_nocapture = r"-?\d+"
        nbf_nocapture = r"-?\d+(?:\.\d*)?"
        nb = f"({nb_nocapture})"
        nbf = f"({nbf_nocapture})"
        eq = r"(?::\s*|\s*=\s*)"
        coeffs = rf"({nb_nocapture}(?:\s*,\s*{nb_nocapture})*)"

        if not line or line.startswith("#"):
            return
        elif m := re.match(rf"^[nN]{eq}{nb}$", line):
            self.N = ZZ(m.group(1))
        elif m := re.match(rf"^skew{eq}{nbf}$", line):
            self.skewness = RR(m.group(1))
        elif m := re.match(rf"^c(\d+){eq}{nb}$", line):
            if not self.f:
                self.f = [ZP(0)] * 2
            self.f[1] += ZZ(m.group(2)) * x**int(m.group(1))
        elif m := re.match(rf"^Y(\d+){eq}{nb}$", line):
            if not self.f:
                self.f = [ZP(0)] * 2
            self.f[0] += ZZ(m.group(2)) * x**int(m.group(1))
        elif m := re.match(rf"^poly(\d+){eq}{coeffs}$", line):
            i = int(m.group(1))
            coeffs = [ZZ(c) for c in m.group(2).split(',')]
            while i >= len(self.f):
                self.f.append(ZP(0))
            self.f[i] = ZP(coeffs)
        elif m := re.match(rf"^poly(\d+){eq}(.*)$", line):
            i = int(m.group(1))
            while i >= len(self.f):
                self.f.append(ZP(0))
            self.f[i] = ZP(m.group(2))
        else:
            raise ValueError(f"Cannot parse line in poly file: {line}")

    def __create_fields(self):
        self.K = [NumberField(f, names=f"alpha{i}")
                  for i, f in enumerate(self.f)]

    def __check_resultant(self):
        """
        Check that the resultant is a multiple of N. If we want to extend
        this to the multiple number field setting, there's a bit of a
        catch here because we really want _all_ polynomials to intersect
        at the same root above N (not just pairwise). E.g. x^2-3x+2,
        x^2-x-2, and x^2+N-1 should not pass this test, yet they
        currently do.
        """
        f = self.f[0]
        for i, g in enumerate(self.f[1:]):
            if f.resultant(g) % self.N != 0:
                raise ValueError("resultant of polynomials 0 and {i}" +
                                 " is not zero modulo N")

    def __read(self, contents):
        self.__clear_fields_for_read()

        for line in contents:
            self.__parseline(line)

        if len(self.f) < 2:
            raise ValueError("poly file must" +
                             " contain at least two polynomials")

        self.__create_fields()
        self.__check_resultant()
        self.nt = CadoNumberTheory(self)

    def __str__(self):
        if not self.N:
            if self.filename:
                return f"NFS configuration, to be read from {self.filename}"
            else:
                return "NFS configuration (empty placeholder)"
        else:
            rep = f"NFS configuration for modulus {self.N}\n"
            for i, K in enumerate(self.K):
                rep += f"  K{i}: {K}\n"
        return rep
