import os

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from .tools import NOK


class BwcXVector(object):
    """
    The X vector is somewhat peculiar, so it deserves its own class. It's
    also a very small file, so reading it multiple times is not really a
    problem.
    """
    def __init__(self, params, dims, dirname):
        """
        The dimension must be passed as a pairs (nrows, ncols).
        The expected size of the X
        vector is one of these.
        """
        self.params = params
        self.filename = os.path.join(dirname, "X")
        xdim, ydim = dims
        self.dim = xdim
        self.X = None

    def read(self):
        self.X = matrix(GF(self.params.p), self.dim, self.params.m)
        print(f"Reading {self.filename}")
        with open(self.filename, "r") as f:
            try:
                it = iter([int(x) for x in f.read().strip().split()])
                items_per_line = next(it)
                for j in range(self.params.m):
                    for k in range(items_per_line):
                        i = next(it)
                        self.X[i, j] += 1
            except StopIteration:
                raise ValueError(f"Short read in {self.filename} {NOK}")
