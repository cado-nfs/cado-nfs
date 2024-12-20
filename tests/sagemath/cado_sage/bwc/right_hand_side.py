from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from cado_sage.tools import OK


class BwcRightHandSide():
    def __init__(self, params, filename):
        self.params = params
        self.filename = filename
        self.R = None
        self.dimension = None

    def read(self):
        print(f"Reading {self.filename}")
        rhs = iter(open(self.filename))
        while line := next(rhs):
            if not line.startswith('#'):
                self.dimension, nrhs, p = (int(c) for c in line.strip().split())
                break
        assert p == self.params.p
        self.R = matrix(GF(p), self.dimension, nrhs)
        sw = self.params.splitwidth
        if self.params.p == 2:
            for i, line in zip(range(self.dimension), rhs):
                words = (int(c) for c in line.strip().split())
                assert len(words) == nrhs // sw
                for j in range(nrhs):
                    self.R[i,j] = (words[j // sw] >> (j % sw)) & 1
        else:
            assert sw == 1
            for i, line in zip(range(self.dimension), rhs):
                cc = line.strip().split()
                assert len(cc) == nrhs
                for j, c in zip(range(nrhs), cc):
                    self.R[i, j] = GF(p)(c)

    def check(self, M, Vs):
        r = self.R.ncols()
        if not r:
            return
        v = Vs.block_by_iteration(0)
        print("Checking that the RHS dimension is consistent ...")
        assert self.dimension == M.ncols_orig if self.params.nullspace == 0 else M.nrows_orig
        print("Checking that the RHS dimension is consistent ..." + OK)
        print("Checking that the RHS makes up the first rows of V ...")
        if self.params.is_nullspace_left():
            assert self.dimension == M.ncols_orig
            w = M.Q * v[:,:r]
            assert w[:self.dimension,:] == self.R
            assert 0 == v[self.dimension:, :r]
            assert 0 == w[self.dimension:, :]
        else:
            assert self.dimension == M.nrows_orig
            assert self.R == v[:self.dimension, :r]
            assert 0 == v[self.dimension:, :r]
        print("Checking that the RHS makes up the first rows of V ... " + OK)
