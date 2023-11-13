import os
import re

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from .BwcParameters import BwcParameters
from .tools import NOK


class BwcVectorBase(object):
    def __init__(self,
                 params: BwcParameters,
                 filename,
                 pattern=None,
                 _what=None):
        self.params = params
        self.filename = filename
        if pattern is not None:
            self.pattern = pattern
        if _what is not None:
            self._what = _what
        if (m := re.match(self.pattern, os.path.basename(filename))) is None:
            raise ValueError(f"{self.filename} : not a {self._what} {NOK}")
        self._filename_data = [int(c) for c in m.groups()]
        j0, j1 = self._filename_data[:2]
        sw = self.params.splitwidth
        if j1 - j0 != sw:
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
        for c in self._filename_data:
            yield c
        yield self.filename

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
