import os
import re
import copy

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF

from collections import defaultdict

from .BwcParameters import BwcParameters
from .tools import mcoeff, read_one_matrix
from .tools import OK, NOK, EXCL


class BwcAFiles(object):
    """
    This class gathers information on all the A files
    (transpose(X)*M^i*Y) that are computed in the Krylov stage of bwc.
    """
    def __init__(self, params: BwcParameters, dims, dirname):
        self.params = params
        self.A = matrix(GF(self.params.p)['x'], self.params.m, self.params.n)
        self.KP = self.A.base_ring()
        self.dims = dims
        self.dirname = dirname
        self.afiles = []
        sw = self.params.splitwidth
        for filename in os.listdir(dirname):
            fn = os.path.join(dirname, filename)
            apat = r"^A(\d+)-(\d+)\.(\d+)-(\d+)$"
            m = re.match(apat, filename)
            if not m:
                continue
            tup = (fn, *[int(x) for x in m.groups()])
            j0, j1 = tup[1:3]
            nj = j1 - j0
            if nj != sw and nj != self.params.n:
                raise ValueError(f"{fn} does not match splitwidth={sw} {NOK}")
            self.afiles.append(tup)

        self.__compute_occupancy_map()

    def __compute_occupancy_map(self):
        """
        Make self.occupancy a list indexed by the column indices from 0
        to n divided by splitwidth, and values being the list of ranges
        for which we happen to have some A files available
        """
        sw = self.params.splitwidth
        xn = self.params.n // sw
        self.occupancy = [[] for j in range(xn)]
        for tup in self.afiles:
            j0, j1, start, end = tup[1:]
            # note that we may have j1-j0 larger than splitwidth
            for j in range(j0, j1, sw):
                xj = j // sw
                self.occupancy[xj].append((start, end))

        # sort each list, and check for overlaps
        for xj in range(xn):
            j = xj * sw
            n = 0
            self.occupancy[xj] = sorted(self.occupancy[xj])
            for b, e in self.occupancy[xj]:
                if b < n:
                    where = f"position {0:self.params.m},{j}, degrees {b}:{n}"
                    raise ValueError(f"overlap in A files for {where}")
                n = e

    def __compute_check_map(self, vfiles):
        """
        Given a list of v files, determine what tests can be done. Warn
        if there are ranges for which verification is not possible.
        """

        check_map = []
        xn = len(self.occupancy)
        sw = self.params.splitwidth

        # First, expand the occupancy map to include which vectors can be
        # used to each range (possibly at the cost of doing many
        # multiplies, but we'll deal with that later on).
        candidates = [[(*x, []) for x in L] for L in self.occupancy]
        for v in vfiles:
            j0, j1, iteration, filename = v
            assert j1-j0 == sw
            xj = j0 // sw
            # Transtorm candidates[j] to a list where the vector v is
            # appended to [2] whenever it can be used for checking the
            # range
            iL = iter(candidates[xj])
            nL = []
            while (ell := next(iL, None)) is not None:
                s, e, usable = ell
                if iteration >= e:
                    # v cannot be used. ell carries over unchanged
                    nL.append(ell)
                elif iteration <= s:
                    # v can be used for the whole range
                    usable.append(v)
                    nL.append(ell)
                else:
                    # v can be used for a sub-range. Split in two
                    nL.append((s, iteration, usable))
                    nL.append((iteration, e, copy.copy(usable) + [v]))
            candidates[xj] = nL

        # This will be used as intermediary data.
        by_vector = defaultdict(lambda: [])

        # Do we have a way to check each range?
        for xj in range(xn):
            j = xj * sw
            winners = []
            # retain only one vector for each range.
            for b, e, usable in candidates[xj]:
                if not usable:
                    what = f"coefficients {b}:{e} in column {j}"
                    print(f"Warning: we have no way to check {what} {EXCL}")
                    continue
                winners.append((b, e, max(usable, key=lambda v: v.iteration)))

            for b, e, v in winners:
                by_vector[v].append((b, e))

        for k in sorted(by_vector.keys(), key=lambda x: tuple(x)):
            L = sorted(by_vector[k])
            b0 = None
            e0 = None
            for b, e in L:
                if b0 is None:
                    b0 = b
                    assert k.iteration == b
                elif b != e0:
                    # TODO: raise or just print a warning ?
                    fn = k.filename
                    where = f"degrees [{e0}:{b}]"
                    raise ValueError(f"Cannot check {where} with {fn} {NOK}")
                e0 = e
            check_map.append((k, e0))

        return check_map

    def read(self):
        x = self.KP.gen()
        sw = self.params.splitwidth
        for filename, j0, j1, start, end in self.afiles:
            ni = self.params.m
            nj = j1 - j0
            k = start
            st = os.stat(filename)
            bytes_per_mat = ni * (nj // sw) * self.params.p_bytes
            nk = st.st_size // bytes_per_mat
            print(f"Reading {filename} (size: {ni}*{nj}, {nk} coeffs)")
            if nk != end - start:
                expect = f"expected {nk*bytes_per_mat} for {nk} coefficients"
                raise ValueError(f"{filename} has wrong size ({expect}) {NOK}")
            with open(filename, 'rb') as f:
                while (M := read_one_matrix(self.params, f, ni, nj)) is not None:  # noqa: E501
                    self.A[:, j0:j1] += x**k * M
                    k += 1

    def check(self, M, x, Vs):
        """
        Takes a BwcMatrix (or DecorrelatedMatrix), a BwcXVector, and a
        list of BwcVector and check that the A data that we have is
        consistent
        """
        for v, k1 in self.__compute_check_map(Vs):
            what = f"coefficients [{v.iteration}:{k1}] of A"
            print(f"Check {what} using {v.filename}")
            w = v.V
            for k in range(v.iteration, k1):
                if x.X.transpose()*w != mcoeff(self.A, k)[:, v.j0:v.j1]:
                    where = f"at coefficient {k}"
                    raise ValueError(f"Inconsistency in A files {where}")
                w = M * w
            print(f"Check {what} using {v.filename} {OK}")
