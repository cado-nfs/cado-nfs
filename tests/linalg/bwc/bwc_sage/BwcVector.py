import os
import re

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from .BwcParameters import BwcParameters
from .tools import NOK


class BwcVector(object):
    def __init__(self, params: BwcParameters, filename):
        self.params = params
        self.filename = filename
        vpat = r"^V(\d+)-(\d+).(\d+)$"
        if (m := re.match(vpat, os.path.basename(self.filename))) is None:
            raise ValueError("Wrong file name for vector: " + self.filename)
        self.j0 = int(m.groups()[0])
        self.j1 = int(m.groups()[1])
        self.iteration = int(m.groups()[2])
        sw = self.params.splitwidth
        if self.j1 - self.j0 != sw:
            what = f"does not match splitwidth={sw}"
            raise ValueError(f"{self.filename} {what} {NOK}")

        sz = os.stat(self.filename).st_size

        if sz % self.params.p_bytes != 0:
            expect = f"expected a multiple of {self.params.p_bytes}"
            what = f"has wrong size ({expect})"
            raise ValueError(f"{self.filename} {what} {NOK}")

        self.dimension = sz // self.params.p_bytes
        self.V = matrix(GF(self.params.p), self.dimension, sw)

    def __iter__(self):
        yield self.j0
        yield self.j1
        yield self.iteration
        yield self.filename

    def __str__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"block {self.j0}-{self.j1}"
        it = f"iteration {self.iteration}"
        return f"{dims} vector for {block} at {it}"

    def __repr__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"block {self.j0}-{self.j1}"
        it = f"iteration {self.iteration}"
        return f"{dims} vector for {block} at {it}"

    def read(self):
        print(f"Reading {self.filename}")
        if self.params.p == 2:
            i = 0
            with open(self.filename, "rb") as fv:
                while (b := bytearray(fv.read(self.params.p_bytes))) != b'':
                    for j in range(self.params.splitwidth):
                        if b[j // 8] & (1 << (j % 8)):
                            self.V[i, j] = 1
                    i += 1
        else:
            i = 0
            with open(self.filename, "rb") as fv:
                while (b := bytearray(fv.read(self.params.p_bytes))) != b'':
                    z = int.from_bytes(b, 'little')
                    assert self.params.splitwidth == 1
                    self.V[i, 0] = z
                    i += 1
