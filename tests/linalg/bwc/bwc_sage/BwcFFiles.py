import os
import re

from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.rings.finite_rings.finite_field_constructor import GF

from collections import defaultdict
from .BwcParameters import BwcParameters
from .tools import read_one_matrix, mcoeff, mdiv, mmod, mdivmod, mrev, mdeg
from .tools import OK, NOK, EXCL


class BwcFFiles(object):
    """
    This class is an abstraction for the collection of F files (and
    F*.rhs files) that are found in a bwc working directory. These files
    represent the linear generator that the lingen program has computed,
    based on a matrix sequence A.

    The collection of F files, together with the information on which of
    these have a .rhs companion file, carries the information of whether
    the linear system that we are currently solving is a homogenous or
    inhomogenous system (or something inbetween). In what follows, we
    often let r be the number of "right-hand-side columns". In
    particular, the first r columns of the degree 0 coefficient of A are
    target vectors. r is such that 0<=r<=n.

    Internally, data is kept here in a form that is close to what we have
    on disk. However, for verifications it might be easier to shuffle
    the data a little bit, as we do in the check() method.
    """

    def __init__(self, params: BwcParameters, dims, dirname):
        self.params = params
        self.dirname = dirname
        self.dims = dims
        self.F = matrix(GF(self.params.p)['x'], self.params.m, self.params.n)
        self.R = matrix(GF(self.params.p), self.params.m, self.params.n)
        self.KP = self.F.base_ring()
        self.ffiles = []
        self.ffiles_rhs = []
        self.rhs_columns = defaultdict(int)
        self._degrees = defaultdict(int)
        sw = self.params.splitwidth
        for filename in os.listdir(dirname):
            pat = r"^F.sols(\d+)-(\d+)\.(\d+)-(\d+)(\.rhs)?$"
            fn = os.path.join(dirname, filename)
            if (m := re.match(pat, filename)):
                tup = (fn, *[int(x) for x in m.groups()[:4]])
                s0, s1, j0, j1 = tup[1:]
                nj = j1 - j0
                ns = s1 - s0
                bytes_per_mat = ns * (nj // sw) * self.params.p_bytes
                if nj != sw or ns != sw:
                    why = f"does not match splitwidth={sw}"
                    raise ValueError(f"{fn} {why} {NOK}")
                st = os.stat(fn)
                nk = st.st_size // bytes_per_mat
                if not m.groups()[4]:
                    self._degrees[s0] = max(self._degrees[s0], nk - 1)
                    self.ffiles.append(tup)
                else:
                    if nk != 1:
                        raise ValueError(f"{fn} has wrong size")
                    self.ffiles_rhs.append(tup)
                    self.rhs_columns[j0] += 1
        self._degree = max(self._degrees.values())
        S = set(self.rhs_columns.values())
        self.rhs_columns = sorted(self.rhs_columns.keys())
        if not S:
            print("Homogenous system")
        elif S != {self.params.n // sw}:
            raise ValueError(f"Inconsistent counts wrt rhs ({S}) {NOK}")
        elif len(self.rhs_columns) == {self.params.n // sw}:
            print("Inhomogenous system")
        else:
            r = len(self.rhs_columns)
            print(f"Mixed (Homogenous / Inhomogenous) system (r={r})")

    def degree(self):
        return self._degree

    def degrees(self):
        return dict(self._degrees)

    def read(self):
        x = self.KP.gen()
        sw = self.params.splitwidth
        for filename, s0, s1, j0, j1 in self.ffiles:
            nj = j1 - j0
            ns = s1 - s0
            st = os.stat(filename)
            bytes_per_mat = ns * (nj // sw) * self.params.p_bytes
            nk = st.st_size // bytes_per_mat
            print(f"Reading {filename} (size: {nj}*{ns}, {nk} coeffs)")
            k = 0
            with open(filename, 'rb') as f:
                while (M := read_one_matrix(self.params, f, nj, ns)) is not None:  # noqa: E501
                    self.F[j0:j1, s0:s1] += x**k * M
                    k += 1

        for filename, s0, s1, j0, j1 in self.ffiles_rhs:
            nj = j1 - j0
            ns = s1 - s0
            st = os.stat(filename)
            print(f"Reading {filename} (size: {nj}*{ns})")
            with open(filename, 'rb') as f:
                # TODO: see above
                M = read_one_matrix(self.params, f, nj, ns)
                assert M is not None
                self.R[j0:j1, s0:s1] += M.transpose()

    def check(self, A):
        """
        This verifies that self is indeed a linear generator for the
        series in A.

        Note that because of the RHS setting, we know that what we
        computed is not exactly such that A*rev(self) is integral. So we
        shuffle data a little bit before we can do the verification
        """
        r = len(self.rhs_columns)
        n = self.params.n
        x = self.F.base_ring().gen()

        D = diagonal_matrix([x] * r + [1] * (n-r))

        q, re = mdivmod(A.A * D, 1)

        Fx = self.R + D * self.F

        d = self.degree()

        # mdeg(Fx) might be d+1
        Frx = mrev(Fx, mdeg(Fx))

        L = mdeg(A.A)

        print("Checking that F is a linear generator for A")
        if mdiv(mmod(q * Frx, L), mdeg(Fx)) != 0:
            raise ValueError("check failed " + NOK)
        print("Checking that F is a linear generator for A ... " + OK)

        # This is a consequence of the check above (but not equivalent).
        # When taking out the inherent left projection on x, the
        # corresponding vector is expected to also project to zero along
        # MQ.transpose()^i*x for all i>=0

        Fr = mrev(self.F, d)
        A0 = mcoeff(A.A, 0)
        A1 = mdiv(A.A, 1)
        R = self.R
        assert A0 * R + mcoeff(A1 * Fr, d) == 0

    def derive_solutions(self, A, Vs, MQ):
        """
        This returns a pair (U, V) such that each of the n columns (u in
        U and v in V) are such that  rhs * u + MQ * v = 0

        Therefore U has size r*n, and V has size N*n
        """

        # The magma code has complicated considerations about checking
        # for solutions that happen only after a few extra rounds of
        # applying M (in the homogenous case).

        r = len(self.rhs_columns)
        U = self.R[:r, :]
        v = Vs.block_by_iteration(0)
        rhs = v[:, :r]
        V = v.parent()()
        for k in range(self.degree()+1):
            V += v * mcoeff(self.F, k)
            v = MQ * v

        print("Checking solutions derived from the linear generator")

        should_be_zero = rhs * U + MQ * V
        if should_be_zero != 0:
            rk = should_be_zero.rank()
            event = f"rhs * U + MQ * V has rank {rk} (should be zero)"
            print(f"check failed: {event} {NOK}")
            print("This is typically _not_ recovered by the C code, " +
                  "but is otherwise rather harmless. Shit happens.")
            raise ValueError("check failed " + NOK)

        if self.params.is_nullspace_right():
            # see what happens with the original matrix.
            QV = MQ.parent.Q * V
            for j in range(QV.ncols()):
                if self.R[:, j] != 0:
                    continue
                warn = f"Warning (innocuous): solution {j} is zero"
                if QV[:, j] == 0:
                    print(f"{warn} {EXCL}")
                elif QV[:MQ.parent.ncols_orig, j] == 0:
                    print(f"{warn} on the interesting columns {EXCL}")

        print("Checking solutions derived from the linear generator " + OK)

        return (U, V)
