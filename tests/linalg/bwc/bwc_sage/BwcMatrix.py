import os
import re
import sys
import math
import copy

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.special import block_matrix, zero_matrix
from sage.modules.free_module_element import vector

from .tools import OK, NOK, EXCL
from .tools import u32, s32
from .BwcParameters import BwcParameters
from .BwcBalancing import BwcBalancing, BwcShuffling


class BwcMatrix(object):
    """
    This class contains stuff to read a matrix in serialized format, but
    also to fetch information related to a given balancing of a matrix,
    or to provide an abstraction of a chain of matrices.

    Sanity checks that compare the on-disk representation with the
    original matrix are included.

    The left and right multiplication operators do not include the
    decorrelating permutation. This means that multplying by this matrix
    is really the same as multiplying by the original matrix that is
    pointed to by the filename, and eventually the output of the whole
    calculation can be checked against this algorithm.

    For cases where the decorrelated matrix M*Q is needed instead, a thin
    wrapper class is provided.
    """
    def __init__(self,
                 params: BwcParameters,
                 matrix=None,
                 wdir=None,
                 balancing_filename=None):
        self.params = params
        self.filename = matrix
        if wdir is None:
            self.wdir = os.path.dirname(self.filename)
        else:
            self.wdir = wdir

        if balancing_filename is not None:
            self.balancing = BwcBalancing(params, balancing_filename)
        else:
            self.balancing = None

        self.__clear_fields_for_read()
        self.__clear_fields_for_fetch_balancing()

    def dimensions(self):
        """
        This works only if the matrix has been read already. Returns a
        pair (d1,d2), where d1 is the number of coordinates of vectors X
        that go in a multiplication X*M, and d2 for vectors Y in M*Y.
        Therefore if nullspace=right this is (nrows, ncols) of the
        underlying matrix, and if nullspace=left it's (ncols,nrows).
        """
        assert self.nrows is not None
        if self.params.is_nullspace_right():
            return (self.nrows, self.ncols)
        else:
            return (self.ncols, self.nrows)

    def __clear_fields_for_read(self):
        self.M = None
        self.row_weights = None
        self.row_weights_filename = None
        self.nrows = None
        self.col_weights = None
        self.col_weights_filename = None
        self.ncols = None
        self.ncoeffs = 0
        # (nrows, ncols) are in most cases equal, they correspond to the
        # dimension of a square matrix. They're both equal to
        # max(nrows_orig, ncols_orig), which are the dimensions of the
        # matrix as it is found in the data files on disk.
        self.nrows_orig = None
        self.ncols_orig = None

    def __clear_fields_for_fetch_balancing(self):
        self.balancing = None
        self.S = None
        self.submatrices = None
        self.Mx = None
        self.sigma = None
        self.Q = None
        self.xQ = None
        self.P = None
        self.Mt = None

    def read(self, force_square=False):
        try:
            self.__read(force_square=force_square)
        except Exception as e:
            # We're really in bad shape if an exception occurs here.
            # We're not even trying to salvage the BwcMatrix object, as
            # the error is most probably obvious.
            print("Exception while reading {self.filename} {NOK}",
                  file=sys.stderr)
            raise e

    class DecorrelatedMatrix(object):
        def __init__(self, parent):
            self.parent = parent
            self.params = parent.params

        def __mul__(self, y):
            """
            multiplication operator for the decorrelated matrix
            (matrix times vector)
            """

            # TODO: This really tells us that we should not be using
            # operator overload blindly as we do here.
            if isinstance(y, BwcMatrix.DecorrelatedMatrix):
                raise ValueError("please do not do things like MQ*MQ")

            if self.params.is_nullspace_right():
                Qy = self.parent.Q * y
                return self.parent.M * Qy
            else:
                yM = y.transpose() * self.parent.M
                return (yM * self.parent.Q).transpose()

        def __rmul__(self, x):
            """
            multiplication operator for the decorrelated matrix
            (vector time matrixvector)
            """

            # TODO: This really tells us that we should not be using
            # operator overload blindly as we do here.
            if isinstance(x, BwcMatrix.DecorrelatedMatrix):
                raise ValueError("please do not do things like MQ*MQ")

            if self.params.is_nullspace_right():
                xM = x.transpose() * self.parent.M
                return (xM * self.parent.Q).transpose()
            else:
                Qx = self.parent.Q * x
                return self.parent.M * Qx

        def __pow__(self, exponent):
            return BwcMatrixPower(self, exponent)

    def decorrelated(self):
        return self.DecorrelatedMatrix(self)

    # Not absolutely sure we want to keep these two operators, in fact
    def __mul__(self, y):
        """
        multiplication operator for this matrix (matrix times vector).
        This multiplication does not include any decorrelated
        permutation.
        """
        if self.params.is_nullspace_right():
            return self.M * y
        else:
            return (y.transpose() * self.M).transpose()

    def __rmul__(self, x):
        """
        multiplication operator for this matrix (vector times matrix),
        This multiplication does not include any decorrelated
        permutation.
        """
        if self.params.is_nullspace_right():
            return (x.transpose() * self.M).transpose()
        else:
            return self.M * x

    def __read(self, nrows=None, ncols=None, force_square=False):
        """
        The parameters nrows and ncols are not meant for general use. We
        only use them when we know the balancing above, and the matrices
        that we're reading are just chunks of the bigger matrix
        """

        if nrows is not None:
            print(f"Reading {self.filename} (size: {nrows}*{ncols})")
        else:
            print(f"Reading {self.filename}")

        self.__clear_fields_for_read()

        self.row_weights = []
        try:
            fn = re.sub(r"\.bin$", ".rw.bin", self.filename)
            self.nrows_orig = os.stat(fn).st_size // 4
            self.row_weights_filename = fn
            with open(self.row_weights_filename, 'rb') as fm:
                while (w := u32(fm, may_fail=True)) is not None:
                    self.row_weights.append(w)
            assert len(self.row_weights) == self.nrows_orig
        except FileNotFoundError:
            self.row_weights_filename = None

        self.col_weights = []
        try:
            fn = re.sub(r"\.bin$", ".cw.bin", self.filename)
            self.ncols_orig = os.stat(fn).st_size // 4
            self.col_weights_filename = fn
            with open(self.col_weights_filename, 'rb') as fm:
                while (w := u32(fm, may_fail=True)) is not None:
                    self.col_weights.append(w)
            assert len(self.col_weights) == self.ncols_orig
        except FileNotFoundError:
            self.col_weights_filename = None

        inline_data = []
        inline_col_weights = []
        inline_row_weights = []
        if self.params.p == 2:
            with open(self.filename, "rb") as fm:
                i = 0
                j = 0
                try:
                    while (length := u32(fm, may_fail=True)) is not None:
                        for jj in range(length):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                pad = j + 1 - len(inline_col_weights)
                                inline_col_weights += [0] * pad
                            inline_col_weights[j] += 1
                            inline_data.append((i, j))
                        i += 1
                        inline_row_weights.append(length)
                except IndexError:
                    what = f"entry {i},{j} in matrix"
                    where = f"at byte 0x{fm.tell():04x} in matrix file"
                    raise ValueError(f"Cannot set {what} {where} {NOK}")
        else:
            with open(self.filename, "rb") as fm:
                i = 0
                j = 0
                v = 0
                try:
                    while (length := u32(fm, may_fail=True)) is not None:
                        for jj in range(length):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                pad = j + 1 - len(inline_col_weights)
                                inline_col_weights += [0] * pad
                            v = s32(fm)
                            inline_col_weights[j] += 1
                            inline_data.append((i, j, v))
                        i += 1
                        inline_row_weights.append(length)
                except IndexError:
                    what = f"entry {i},{j} in matrix"
                    where = f"at byte 0x{fm.tell():04x} in matrix file"
                    raise ValueError(f"Cannot set {what} {where} {NOK}")
        if self.row_weights:
            assert self.row_weights == inline_row_weights
        else:
            self.row_weights = inline_row_weights
            self.nrows_orig = len(self.row_weights)

        if self.col_weights:
            icw = inline_col_weights
            assert self.col_weights[:len(icw)] == icw
            assert vector(self.col_weights[len(icw):]) == 0
        else:
            self.col_weights = inline_col_weights
            self.ncols_orig = len(self.col_weights)

        if nrows is not None:
            assert self.nrows_orig <= nrows
            self.nrows = nrows
        else:
            self.nrows = self.nrows_orig

        if ncols is not None:
            assert self.ncols_orig <= ncols
            self.ncols = ncols
        else:
            # attention. We fill it for completeness, but it can be
            # that the matrix has some zero cols at the end. This is
            # the main use case for ncols, in fact.
            self.ncols = self.ncols_orig

        self.ncoeffs = len(inline_data)
        r2 = sum([float(x*x) for x in inline_row_weights]) / self.nrows_orig
        c2 = sum([float(x*x) for x in inline_col_weights]) / self.ncols_orig
        rmean = float(self.ncoeffs / self.nrows_orig)
        rsdev = math.sqrt(r2-rmean**2)
        cmean = float(self.ncoeffs / self.ncols_orig)
        csdev = math.sqrt(c2-cmean**2)

        rowscols = f"{self.nrows_orig} rows {self.ncols_orig} cols"
        coeffs = f"{self.ncoeffs} coefficients"
        stats = f"row: {rmean:.2f}~{rsdev:.2f}, col: {cmean:.2f}~{csdev:.2f}"
        print(f"{rowscols}, {coeffs} ({stats})")

        if force_square and self.nrows_orig != self.ncols_orig:
            print("Padding to a square matrix")
            self.nrows = max(self.nrows_orig, self.ncols_orig)
            self.ncols = max(self.nrows_orig, self.ncols_orig)

        self.M = matrix(GF(self.params.p), self.nrows, self.ncols, sparse=True)
        if self.params.p == 2:
            for i, j in inline_data:
                if self.M[i, j] != 0:
                    print(f"Warning: repeated entry in matrix data {EXCL}")
                self.M[i, j] += 1
        else:
            for i, j, v in inline_data:
                if self.M[i, j] != 0:
                    print(f"Warning: repeated entry in matrix data {EXCL}")
                self.M[i, j] += v
        if self.balancing is not None:
            self.balancing.read()
            self.S = BwcShuffling(self.params, self)

    def __subM(self, i, j):
        return self.submatrices[i][j].M

    def fetch_balancing(self, nh, nv):
        """
        Based on the filename of the matrix, try to see if a balancing is
        defined, read it, and read all submatrices. Of course this
        assumes that the submatrices were savec during the balancing
        operation done by the C++ code
        """

        self.__clear_fields_for_fetch_balancing()

        dir, base = os.path.split(self.filename)
        dir = self.wdir
        sdir = base.removesuffix(".bin") + f".{nh}x{nv}"
        bfile = os.path.join(dir, sdir, sdir + ".bin")
        if not os.path.exists(bfile):
            return FileNotFoundError(bfile)
        self.balancing = BwcBalancing(self.params, bfile)
        print("Reading balancing from " + bfile)
        self.balancing.read()
        self.S = BwcShuffling(self.params, self)
        self.submatrices = [[None for j in range(nv)] for i in range(nh)]
        bal = self.balancing
        for i in range(nh):
            for j in range(nv):
                # id = f"{self.balancing.checksum:08x}.h{i}.v{j}"
                id = f"h{i}.v{j}"
                fname = os.path.join(dir, sdir, f"{sdir}.{id}.bin")
                nrp = bal.tr // bal.nh
                ncp = bal.tc // bal.nv
                self.submatrices[i][j] = BwcMatrix(self.params, fname)
                self.submatrices[i][j].__read(nrows=nrp, ncols=ncp)
        # Now we want to check the consistency of the matrix that we just
        # read.

        # The submatrices may be stored in transposed order. This is done
        # because it is convenient for the matmul layers to read them so.
        # Sometimes. In truth, we should perhaps get rid of this
        # complication, I don't really know.
        t = lambda x:x
        if self.params.is_nullspace_left():
            t = lambda x:x.transpose()


        self.Mx = block_matrix(nh, nv,
                               [t(self.__subM(i, j))
                                for i in range(nh)
                                for j in range(nv)])
        self.sigma = self.S.sc.matrix()
        self.Q = self.S.shuf.matrix()
        self.xQ = self.S.xshuf.matrix()
        self.P = self.S.pr.matrix()

    def check_balancing_submatrices_are_well_formed(self):
        """
        This checks that given the balancing permutation sigma, the
        concatenation of the different submatrices at least gives a
        matrix that is consistent with the number of rows and columns of
        the matrix we supposedly started with.

        It is important to understand that this function does not check
        data that is read by self.read(), and it's on purpose.
        """

        bal = self.balancing

        A = (self.P*self.sigma).transpose() * self.Mx * self.sigma

        what = "Reconstructed matrix from submatrices"
        if A[bal.nr:bal.tr] != 0:
            raise ValueError(f"{what} has garbage trailing rows {NOK}")
        if A.transpose()[bal.nc:bal.tc] != 0:
            raise ValueError(f"{what} has garbage trailing columns {NOK}")

    def check_balancing_submatrices_consistency_with_M(self):
        """
        This does what the previous code does not do: check that the
        concatenation of the submatrices is consistent with the matrix we
        started with. Of course, to do so, we must have called
        self.read() first.
        """

        bal = self.balancing

        if self.M is None:
            raise ValueError(f".read() must be called first {EXCL}")

        if bal is None:
            raise ValueError(f".fetch_balancing() must be called first {EXCL}")

        if self.nrows_orig != bal.nr or self.ncols_orig != bal.nc:
            what = "inconsistent with number of rows/columns of the matrix"
            raise ValueError(f"{self.balancing.filename} is {what} {NOK}")

        A = (self.P*self.sigma).transpose() * self.Mx * self.sigma

        B = self.M * self.Q
        sub = lambda T: T.submatrix(0, 0, self.nrows, self.ncols)  # noqa: E731

        what = "Reconstructed matrix from submatrices"
        if sub(A) != B:
            raise ValueError(f"{what} does not match M " + NOK)

        self.Mt = zero_matrix(bal.tr, bal.tc)
        self.Mt[:self.nrows, :self.ncols] = self.M

        if A != self.Mt * self.xQ:
            # It's not really possible for this to happen if the above
            # check has passed, so let's just give the same error message
            # anyway.
            raise ValueError(f"{what} does not match M " + NOK)

    def check(self):
        nh = self.balancing.nh
        nv = self.balancing.nv

        what = f"that submatrices are well formed for {nh}*{nv} balancing"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_are_well_formed()
        print(f"Checking {what} ...  {OK}")

        what = "that submatrices are consistent with the matrix M"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_consistency_with_M()
        print(f"Checking {what} ... {OK}")

    def __pow__(self, exponent):
        return BwcMatrixPower(self, exponent)


class BwcMatrixPower(object):
    def __init__(self, M, exponent):
        self.M = M
        self.exponent = exponent

    def __mul__(self, y):
        v = copy.copy(y)
        for i in range(self.exponent):
            v = self.M * v
        return v

    def __rmul__(self, y):
        v = copy.copy(y)
        for i in range(self.exponent):
            v = v * self.M
        return v
