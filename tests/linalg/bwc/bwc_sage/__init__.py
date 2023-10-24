"""
This package is work in progress. As of late 2023, it is still far from
the functionality level of the original magma script (e.g. 

The work plan is:
    1. prime field case (nullspace=right), with arbitrary balancing
    2. we want to use this to debug the double matrix setting
    3. extend to the binary case (nullspace=left). We know that there are
       a few monsters down that path, like m4ri data being hard to
       initialize totally freely from within sage
    4. address the remaining cases like interleaving and so on.
"""

import os
import re
import sys
import itertools
import math

from sage.rings.integer import Integer
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.special import block_matrix,zero_matrix
from collections import defaultdict
import copy

FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES  = 64
MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES = 4

OK = "ok ✅"
NOK = "❌"
EXCL = "❗"
HURRAH = "🥳"

# class BwcDirectory(object):
#     def __init__(self, path):
#         self.path = os.path(path)
#     def detect_files(self):
#         BwcParallelismSettings
# 
# 
# class BwcParallelismSettings(object):
#     def __init__(self):
#         pass

def get_int(f, bytes=4, signed=False, repeat=None, may_fail=False):
    assert not (may_fail and repeat is not None)
    if repeat is not None:
        return [ get_int(f, bytes=bytes, signed=signed) for i in range(repeat) ]
    l = f.read(bytes)
    if not l and may_fail:
        return None
    assert l
    return int.from_bytes(l, 'little', signed=signed)

def u32(f, *args, **kwargs):
    return get_int(f, bytes=4, *args, **kwargs)

def u64(f, *args, **kwargs):
    return get_int(f, bytes=8, *args, **kwargs)

def s32(f, *args, **kwargs):
    return get_int(f, bytes=4, signed=True, *args, **kwargs)

def s64(f, *args, **kwargs):
    return get_int(f, bytes=8, signed=True, *args, **kwargs)


class BwcParameters(object):
    NULLSPACE_LEFT = 0
    NULLSPACE_RIGHT = 1
    def __init__(self, *args, **kwargs):
        try:
            self.m = kwargs['m']
            self.n = kwargs['n']
            self.p = kwargs['p']
        except KeyError:
            raise KeyError("Arguments m n p to BwcParameters are required")
        self.wordsize = kwargs.get('wordsize', 64)
        if self.p == 2:
            self.splitwidth = 64
            self.nullspace = self.NULLSPACE_LEFT
        else:
            self.splitwidth = 1
            self.nullspace = self.NULLSPACE_RIGHT
        self.p_words = len((self.p**self.splitwidth-1).digits(2**self.wordsize))
        ## self.p_bytes = len((self.p**self.splitwidth-1).digits(2**8))
        self.p_bytes = self.p_words * (self.wordsize // 8)
        n = kwargs.get('nullspace', None)
        if n == 'left':
            self.nullspace = self.NULLSPACE_LEFT
        elif n == 'right':
            self.nullspace = self.NULLSPACE_RIGHT
        else:
            # default behavior is to let the prime decide.
            pass
    def is_nullspace_left(self):
        return self.nullspace == self.NULLSPACE_LEFT
    def is_nullspace_right(self):
        # of course it's just for convenience matters that we have both.
        return self.nullspace == self.NULLSPACE_RIGHT

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
    def __init__(self, params : BwcParameters, matrix=None, wdir=None, balancing_filename=None):
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
            print("--- Exception while reading matrix {self.filename} --- {NOK}", file=sys.stderr)
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
            if self.params.is_nullspace_right():
                return self.parent.M * self.parent.Q * y
            else:
                return (y.transpose() * self.parent.M * self.parent.Q).transpose()

        def __rmul__(self, x):
            """
            multiplication operator for the decorrelated matrix
            (vector time matrixvector)
            """
            if self.params.is_nullspace_right():
                return (x.transpose() * self.parent.M * self.parent.Q).transpose()
            else:
                return self.parent.M * self.parent.Q * x

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

        self.row_weights=[]
        try:
            fn = re.sub(r"\.bin$", ".rw.bin", self.filename)
            self.nrows_orig = os.stat(fn).st_size // 4
            self.row_weights_filename = fn
            with open(self.row_weights_filename, 'rb') as fm:
                while (l := fm.read(4)):
                    self.row_weights.append(int.from_bytes(l, 'little'))
            assert len(self.row_weights) == self.nrows_orig
        except FileNotFoundError:
            self.row_weights_filename = None

        self.col_weights=[]
        try:
            fn = re.sub(r"\.bin$", ".cw.bin", self.filename)
            self.ncols_orig = os.stat(fn).st_size // 4
            self.col_weights_filename = fn
            with open(self.col_weights_filename, 'rb') as fm:
                while (l := fm.read(4)):
                    self.col_weights.append(int.from_bytes(l, 'little'))
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
                    while (l := u32(fm, may_fail=True)) is not None:
                        for jj in range(l):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                inline_col_weights += [0] * (j+1-len(inline_col_weights))
                            inline_col_weights[j] += 1
                            inline_data.append((i,j))
                        i += 1
                        inline_row_weights.append(l)
                except IndexError:
                    raise ValueError(f"Cannot set entry {i},{j} in matrix (at byte 0x{fm.tell():04x} in matrix file)")
        else:
            with open(self.filename, "rb") as fm:
                i = 0
                j = 0
                v = 0
                try:
                    while (l := u32(fm, may_fail=True)) is not None:
                        for jj in range(l):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                inline_col_weights += [0] * (j+1-len(inline_col_weights))
                            v = s32(fm)
                            inline_col_weights[j] += 1
                            inline_data.append((i,j,v))
                        i += 1
                        inline_row_weights.append(l)
                except IndexError:
                    raise ValueError(f"Cannot set entry {i},{j} to value 0x{v:x} in matrix (at byte 0x{fm.tell():04x} in matrix file)")
        if self.row_weights:
            assert self.row_weights == inline_row_weights
        else:
            self.row_weights = inline_row_weights
            self.nrows_orig = len(self.row_weights)

        if self.col_weights:
            assert self.col_weights == inline_col_weights
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
            for i,j in inline_data:
                self.M[i,j]=1
        else:
            for i,j,v in inline_data:
                self.M[i,j]=v
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

        dir,base = os.path.split(self.filename)
        dir = self.wdir
        sdir = base.removesuffix(".bin") + f".{nh}x{nv}"
        bfile = os.path.join(dir, sdir, sdir + ".bin")
        if not os.path.exists(bfile):
            return FileNotFoundError(bfile)
        self.balancing = BwcBalancing(self.params, bfile)
        print("Reading balancing from " + bfile)
        self.balancing.read()
        self.S = BwcShuffling(self.params, self)
        self.submatrices = [ [ None for j in range(nv) ] for i in range(nh) ]
        bal = self.balancing
        for i in range(nh):
            for j in range(nv):
                # id = f"{self.balancing.checksum:08x}.h{i}.v{j}"
                id = f"h{i}.v{j}"
                fname = os.path.join(dir,sdir,f"{sdir}.{id}.bin")
                nrp = bal.tr // bal.nh
                ncp = bal.tc // bal.nv
                self.submatrices[i][j] = BwcMatrix(self.params, fname)
                self.submatrices[i][j].__read(nrows = nrp, ncols = ncp)
        # Now we want to check the consistency of the matrix that we just
        # read.

        self.Mx = block_matrix(nh, nv, [ self.__subM(i,j) for i in range(nh) for j in range(nv) ])
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

        if A[bal.nr:bal.tr] != 0:
            raise ValueError("Reconstructed matrix from submatrices has garbage trailing rows " + NOK)
        if A.transpose()[bal.nc:bal.tc] != 0:
            raise ValueError("Reconstructed matrix from submatrices has garbage trailing columns " + NOK)


    def check_balancing_submatrices_consistency_with_M(self):
        """
        This does what the previous code does not do: check that the
        concatenation of the submatrices is consistent with the matrix we
        started with. Of course, to do so, we must have called
        self.read() first.
        """

        bal = self.balancing

        if self.M is None:
            raise ValueError(f"Cannot check consistency with matrix M, call .read() first {EXCL}")

        if bal is None:
            raise ValueError(f"Cannot check consistency with balancing that has not been read yet, call .fetch_balancing() first {EXCL}")

        if self.nrows_orig != bal.nr or self.ncols_orig != bal.nc:
            raise ValueError(f"balancing data from {self.balancing.filename} is inconsistent with number of rows/columns of the matrix {NOK}")

        A = (self.P*self.sigma).transpose() * self.Mx * self.sigma

        B = self.M * self.Q
        sub = lambda T: T.submatrix(0,0,self.nrows,self.ncols)

        if sub(A) != B:
            raise ValueError("Reconstructed matrix from submatrices does not match M " + NOK)

        self.Mt = zero_matrix(bal.tr, bal.tc)
        self.Mt[:self.nrows, :self.ncols] = self.M

        if A != self.Mt * self.xQ:
            # It's not really possible for this to happen if the above
            # check has passed, so let's just give the same error message
            # anyway.
            raise ValueError("Reconstructed matrix from submatrices does not match M " + NOK)

    def check(self):
        nh = self.balancing.nh
        nv = self.balancing.nv
        print(f"Checking that submatrices are well formed for {nh}*{nv} balancing ...")
        self.check_balancing_submatrices_are_well_formed()
        print(f"Checking that submatrices are well formed for {nh}*{nv} balancing ...  {OK}")
        print(f"Checking that submatrices are consistent with the matrix M ...")
        self.check_balancing_submatrices_consistency_with_M()
        print(f"Checking that submatrices are consistent with the matrix M ... {OK}")
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

class BwcBalancing(object):
    def __init__(self, params : BwcParameters, filename):
        self.params = params
        self.filename = filename

    def __flag_bytes_to_dict(self, flags):
        txflags = dict()
        if flags & 1:
            txflags['colperm']=1
        if flags & 2:
            txflags['rowperm']=1
        if flags & 8:
            txflags['replicate']=1
        return txflags

    def decorrelation_is_trivial(self):
        return self.pa == 1 and self.pb == 0

    def read(self):
        with open(self.filename, 'rb') as f:
            zero, magic = u32(f, repeat=2)
            assert zero == 0
            assert magic == 0xba1a0000
            self.nh, self.nv = u32(f, repeat=2)
            # (nr, nc) are the same as (nrows_orig, ncols_orig) in the
            # original matrix. However (nrows, ncols) in the BwcMatrix
            # object may be equal to the max of these two, and thus may
            # fail to match the (nr, nc) that we have here.
            self.nr, self.nc = u32(f, repeat=2)
            self.nzr, self.nzc = u32(f, repeat=2)
            self.ncoeffs = u64(f)
            self.checksum, flags = u32(f, repeat=2)

            self.pa, self.pb = u32(f, repeat=2)
            self.pai, self.pbi = u32(f, repeat=2)

            self.txflags = self.__flag_bytes_to_dict(flags)

            def pad(n, k, b=1):
                # return the next multiple of b*k
                return ((n + (b*k) - 1) // (b*k)) * (b*k)

            chunk = int(FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES/MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES)

            self.tr = pad(self.nr, self.nh * self.nv, chunk)
            self.tc = pad(self.nc, self.nh * self.nv, chunk)

            if 'replicate' in self.txflags:
                self.tr = max(self.tr, self.tc)
                self.tc = max(self.tr, self.tc)

            s = self.nh * self.nv

            self.tr = pad(self.tr, s)
            self.tc = pad(self.tc, s)

            print(f"{self.nr} rows {self.nc} cols, split {self.nh}x{self.nv}, checksum 0x{self.checksum:x} {list(self.txflags.keys())}")

            self.rowperm = None
            if 'rowperm' in self.txflags:
                self.rowperm = u32(f, repeat=self.tr)

            self.colperm = None
            if 'colperm' in self.txflags:
                self.colperm = u32(f, repeat=self.tc)

class BwcShuffling(object):
    def __init__(self, par:BwcParameters, mat:BwcMatrix):
        assert mat.balancing is not None
        self.params = par
        self.M = mat
        bal = self.M.balancing

        assert bal.nc == self.M.ncols_orig
        assert bal.nr == self.M.nrows_orig

        Sc = SymmetricGroup(range(bal.tc))
        Sr = SymmetricGroup(range(bal.tr))

        self.sc = None
        if bal.colperm is not None:
            self.sc = Sc(bal.colperm)
            self.scinv = self.sc**-1

        self.sr = None
        if bal.rowperm is not None:
            self.sr = Sr(bal.rowperm)
            self.scinv = self.sc**-1

        try:
            self.pr = Sr([self.__pr(Integer(x)) for x in range(bal.tr)])
            self.prinv = self.pr**-1
        except KeyError:
            print("Error in creating permutation pr " + NOK, file=sys.stderr)
            self.pr = None
            self.prinv = None

        # the minimal context in which the decorrelating permutation may
        # be defined is in dimension bal.nc == self.M.ncols_orig.
        # However this would be mostly useless, and we prefer to define
        # it minimally as a permutation of the indices of the (most often
        # square) matrix of dimension (nrows,ncols)
        # 
        # try:
        #     # the decorrelating permutation can be defined on the non-padded
        #     # indices too.
        #     Sc0 = SymmetricGroup(range(mat.ncols_orig))
        #     self.shuf = Sc0([self.__preshuf(Integer(x)) for x in range(mat.ncols_orig)])
        #     self.shufinv = self.shuf**-1
        # except KeyError:
        #     print("Error in creating permutation shuf " + NOK, file=sys.stderr)
        #     self.shuf = None
        #     self.shufinv = None

        try:
            Sc0 = SymmetricGroup(range(mat.ncols))
            self.shuf = Sc0([self.__preshuf(Integer(x)) for x in range(mat.ncols)])
            self.shufinv = self.shuf**-1
        except KeyError:
            print("Error in creating permutation shuf " + NOK, file=sys.stderr)
            self.shuf = None
            self.shufinv = None

        try:
            self.xshuf = Sc([self.__preshuf(Integer(x)) for x in range(bal.tc)])
            self.xshufinv = self.xshuf**-1
        except KeyError:
            print("Error in creating permutation xshuf " + NOK, file=sys.stderr)
            self.xshuf = None
            self.xshufinv = None

    def __preshuf(self, x):
        bal = self.M.balancing
        if x >= bal.nc:
            return x
        return (bal.pa * x + bal.pb) % bal.nc

    def __preshuf_inv(self, x):
        bal = self.M.balancing
        if x >= bal.nc:
            return x
        return (bal.pai * x + bal.pbi) % bal.nc


    def __pr(self, x):
        # nrp:=nv*(chunk*Ceiling(x/chunk)) where x is Ceiling (nr/(nh*nv));
        # ncp:=nh*(chunk*Ceiling(x/chunk)) where x is Ceiling (nc/(nh*nv));
        # (seems that nrp is tr/nh, so nrx=tr)

        bal = self.M.balancing
        nz = bal.tr // (bal.nh * bal.nv)
        q,r = x.quo_rem(nz)
        qq,qr = q.quo_rem(bal.nv)
        return (qr * bal.nh + qq) * nz + r

    def __prinv(self, x):
        bal = self.M.balancing
        nz = bal.tr // (bal.nh * bal.nv)
        q,r = x.quo_rem(nz)
        qq,qr = q.quo_rem(bal.nh)
        return (qr * bal.nv + qq) * nz + r

class BwcVector(object):
    def __init__(self, params : BwcParameters, filename):
        self.params = params
        self.filename = filename
        if (m := re.match(r"^V(\d+)-(\d+).(\d+)", os.path.basename(self.filename))) is None:
            raise ValueError("Wrong file name for vector: " + self.filename)
        self.j0 = int(m.groups()[0])
        self.j1 = int(m.groups()[1])
        self.iteration = int(m.groups()[2])
        if self.j1 - self.j0 != self.params.splitwidth:
            raise ValueError(f"Wrong file name for vector (does not match splitwidth): {self.filename} {NOK}")

        sz = os.stat(self.filename).st_size

        # note that p_bytes includes splitwidth

        if sz % self.params.p_bytes != 0:
            raise ValueError("Wrong file size for {} ({} is not a multiple of {})".format(self.filename, sz, self.params.p_bytes) + " " + NOK)

        self.dimension = sz // self.params.p_bytes

        self.V = matrix(GF(self.params.p), self.dimension, self.params.splitwidth)
    def __iter__(self):
        yield self.j0
        yield self.j1
        yield self.iteration
        yield self.filename

    def __str__(self):
        return f"{self.V.nrows()}x{self.V.ncols()} vector for block {self.j0}-{self.j1} at iteration {self.iteration}"

    def __repr__(self):
        return f"{self.V.nrows()}x{self.V.ncols()} vector for block {self.j0}-{self.j1} at iteration {self.iteration}"

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

def check_all_vectors(Vs, mat):
    for j0,j1 in sorted(list(set([(x.j0,x.j1) for x in Vs]))):
        isj0j1 = lambda x: x.j0==j0 and x.j1 == j1
        W = sorted([ x for x in Vs if isj0j1(x) ], key=lambda x:x.iteration)
        if len(W) == 1:
            print(f"Warning: only vector {W[0].filename} was found" +
                  " for block {j0}-{j1}: no check possible {EXCL}")
            continue
        print(f"Checking {len(W)} vectors for block {j0}-{j1}")
        print("Iterations: " + ", ".join([str(x.iteration) for x in W]))
        iW = iter(W)
        v = next(iW)
        i = v.iteration
        V = v.V
        print(f"Check base is {v.filename}")
        for t in iW:
            print(f"Checking {t.filename} ... ")
            V = mat**(t.iteration-i) * V
            i = t.iteration
            assert V == t.V
            print(f"Checking {t.filename} ... {OK}")

class BwcXVector(object):
    """
    The X vector is somewhat peculiar, so it deserves its own class. It's
    also a very small file, so reading it multiple times is not really a
    problem.
    """
    def __init__(self, params, dims, filename):
        """
        The dimension must be passed as a pairs (nrows, ncols).
        The expected size of the X
        vector is one of these.
        """
        self.params = params
        self.filename = filename
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
                        self.X[i,j] += 1
            except StopIteration:
                raise ValueError(f"Short read in {self.filename} {NOK}")




class BwcCheckData(object):
    """
    This class gathers information on all the C[vrdt] vectors in a given
    bwc directory.

    Ct0-$nchecks.0-$m (also referred to as T): a random matrix of size bw->m*nchecks (yes it's named like this because the data is written row-major, and the first interval in the name customarily denotes the number of items in a major division).
    Cr0-$nchecks.0-$nchecks (also referred to as R): a sequence of random matrices of size nchecks * nchecks
    Cv0-<splitwidth>.<j> == trsp(M)^j * X * Ct (See note (T))
    Cd0-<splitwidth>.<j> == \\sum_{0<=i<j} trsp(M)^i * X * Ct * Cr[i] (See note (T))

    """
    def __init__(self, params : BwcParameters, dims, dirname):
        """
        Out of convenience, we automatically search for check files in
        the given directory
        """
        self.params = params
        self.dims = dims
        self.dirname = dirname
        self.tfiles = []
        self.rfiles = []
        self.vfiles = []
        self.dfiles = []
        self.T = None
        self.R = []
        self.V = []
        self.D = []
        ls = os.listdir(dirname)
        for x in ls:
            if self.__try_rfile(x):
                pass
            elif self.__try_tfile(x):
                pass
            elif self.__try_vfile(x):
                pass
            elif self.__try_dfile(x):
                pass

        for a in [self.tfiles, self.rfiles, self.dfiles, self.vfiles]:
            if len(a) == 0:
                print(f"No complete set of check files found in {self.dirname}, cannot verify the check files {EXCL}")
                return

        if len(self.tfiles) != 1:
            raise ValueError(f"Conflicting sets of check files are found in {self.dirname}, this is not supported {NOK}")

        if len(self.rfiles) != 1:
            raise ValueError(f"Conflicting sets of check files are found in {self.dirname}, this is not supported {NOK}")

        self.nchecks = self.tfiles[0][1]

        if self.tfiles[0][2] != self.params.m:
            raise ValueError(f"A t file was found in {self.dirname}, but for a different m ({self.tfiles[0][2]} != {self.params.m}). Thiss is not supported {NOK}")

        for d in self.dfiles:
            if d[1] != self.params.splitwidth:
                raise ValueError(f"File {d[0]} does not refer to the correct splitwidth ({self.params.splitwidth})")

        for v in self.vfiles:
            if v[1] != self.params.splitwidth:
                raise ValueError(f"File {v[0]} does not refer to the correct splitwidth ({self.params.splitwidth})")

        self.x = BwcXVector(self.params, self.dims, os.path.join(self.dirname, "X"))

    def __try_tfile(self, x):
        tpat = r"^Ct0-(\d+)\.0-(\d+)$"
        if (m := re.match(tpat, x)):
            self.tfiles.append((x, int(m.groups()[0]), int(m.groups()[1])))

    def __try_rfile(self, x):
        rpat = r"^Cr0-(\d+)\.0-(\d+)$"
        if (m := re.match(rpat, x)):
            self.rfiles.append((x, int(m.groups()[0]), int(m.groups()[1])))

    def __try_vfile(self, x):
        vpat = r"^Cv0-(\d+)\.(\d+)$"
        if (m := re.match(vpat, x)):
            self.vfiles.append((x, int(m.groups()[0]), int(m.groups()[1])))

    def __try_dfile(self, x):
        dpat = r"^Cd0-(\d+)\.(\d+)$"
        if (m := re.match(dpat, x)):
            self.dfiles.append((x, int(m.groups()[0]), int(m.groups()[1])))

    def __read_tfile(self):
        """
        There is one single T file, a priori. It's a fixed size matrix.
        """
        filename = os.path.join(self.dirname, self.tfiles[0][0])
        st = os.stat(filename)
        s = self.nchecks * self.params.m * self.params.p_bytes
        if st.st_size != s:
            raise ValueError(f"{filename} has incorrect size (expected {s}) {NOK}")
        print(f"Reading {filename}")
        with open(filename, "rb") as fv:
            L = []
            for i in range(self.nchecks * self.params.m):
                b = bytearray(fv.read(self.params.p_bytes))
                if not b:
                    raise ValueError(f"{filename}: short read after {fv.tell()} bytes {NOK}")
                L.append(int.from_bytes(b, 'little'))
            self.T = matrix(GF(self.params.p), self.params.m, self.nchecks, L)

    def __read_rfile(self):
        """
        There is one single R file, but it's a sequence of arbitrary
        length
        """
        filename = os.path.join(self.dirname, self.rfiles[0][0])
        st = os.stat(filename)
        s = self.nchecks * self.nchecks * self.params.p_bytes
        if st.st_size % s != 0:
            raise ValueError(f"{filename} has incorrect size (expected a multiple of {s}) {NOK}")
        p = self.params.p
        print(f"Reading {filename}")
        with open(filename, "rb") as fv:
            while True:
                L = []
                for i in range(self.nchecks * self.nchecks):
                    b = bytearray(fv.read(self.params.p_bytes))
                    if not i and not b:
                        # we're done.
                        return
                    if not b:
                        raise ValueError(f"{filename}: short read after {fv.tell()} bytes {NOK}")
                    L.append(int.from_bytes(b, 'little'))
                
                self.R.append(matrix(GF(p), self.nchecks, self.nchecks, L))

    def __read_dv(self, files):
        """
        common code for D and V files, they have the exact same structure
        """
        ldim = self.dims[0]
        p = self.params.p
        ret = []
        for basename, sw, j in files:
            filename = os.path.join(self.dirname, basename)
            print(f"Reading {filename}")
            st = os.stat(filename)
            s = ldim * self.params.p_bytes
            if st.st_size != s:
                raise ValueError(f"{filename} has incorrect size (expected a multiple of {s}) {NOK}")
            with open(filename, "rb") as fv:
                L = []
                for i in range(ldim):
                    b = bytearray(fv.read(self.params.p_bytes))
                    if not b:
                        raise ValueError(f"{filename}: short read after {fv.tell()} bytes {NOK}")
                    L.append(int.from_bytes(b, 'little'))
                ### XXX this is clearly wrong if sw>1
                ### see email thread with malb, 19-23 Jan 2023
                assert sw==1
                ret.append((j, filename, matrix(GF(p), ldim, sw, L)))
        return sorted(ret)


    def __read_dfiles(self):
        self.D = self.__read_dv(self.dfiles)

    def __read_vfiles(self):
        self.V = self.__read_dv(self.vfiles)

    def check(self, M):
        """
        There's no such thing as checking the t and r files, since
        they're pure random data. However we do want to check the v and d
        files.
        Cv0-<splitwidth>.<j> == trsp(M)^j * X * Ct (See note (T))
        Cd0-<splitwidth>.<j> == \\sum_{0<=i<j} trsp(M)^i * X * Ct * Cr[i] (See note (T))
        """

        xt = self.x.X * self.T
        itv = iter(self.V)
        itd = iter(self.D)
        nvi, nvf, nvZ = next(itv, (None, None, None))
        ndi, ndf, ndZ = next(itd, (None, None, None))
        i = 0
        S = xt.parent()(0)
        errors = 0
        while nvi is not None and ndi is not None:
            ni = min(nvi, ndi)
            while i < ni:
                S += xt * self.R[i]
                i += 1
                xt = xt * M
            if ni == nvi:
                print(f"Checking {nvf} ... ")
                if nvZ != xt:
                    print(f"Checking {nvf} ... " + NOK)
                    errors += 1
                else:
                    print(f"Checking {nvf} ... " + OK)
                nvi, nvf, nvZ = next(itv, (None, None, None))
            if ni == ndi:
                print(f"Checking {ndf} ... ")
                if ndZ != S:
                    print(f"Checking {ndf} ... " + NOK)
                    errors += 1
                else:
                    print(f"Checking {ndf} ... " + OK)
                ndi, ndf, ndZ = next(itd, (None, None, None))
        if errors:
            raise ValueError(f"Errors found while checking V and D files {NOK}")


    def read(self):
        """
        This reads all check files that have been found. Consistency is
        not checked, but consistency of the data size *IS* checked. We
        also check for possible extraneous files in some cases.
        """
        self.__read_tfile()
        self.__read_rfile()
        self.x.read()
        self.__read_dfiles()
        self.__read_vfiles()



def scan_for_bwc_vectors(params, dirname, read=False):
    """
    This functions looks for files in a given directory that match the
    file name pattern of bwc vectors. The ones that are found are
    returned as BwcVector objects.

    By default, the .read() method is not called on the returned
    BwcVector objects, but that can be amended with the optional boolean
    read=True.
    """
    print(f"Scanning for bwc vectors in {dirname} ... ")
    l = lambda x: os.path.join(dirname, x)
    pat = r"^V(\d+)-(\d+).(\d+)$"
    L = [ BwcVector(params, l(x)) for x in os.listdir(dirname) if re.match(pat, x)]
    L = sorted(L, key=lambda x: (x.j0, x.j1, x.iteration))
    print(f"Scanning for bwc vectors in {dirname} ... {len(L)} vectors found")
    if read:
        for v in L:
            v.read()
    return L


class BwcAFiles(object):
    """
    This class gathers information on all the A files
    (transpose(X)*M^i*Y) that are computed in the Krylov stage of bwc.
    """
    def __init__(self, params : BwcParameters, dims, dirname):
        self.params = params
        self.A = matrix(GF(self.params.p)['x'], self.params.m, self.params.n)
        self.KP = self.A.base_ring()
        self.dims = dims
        self.dirname = dirname
        self.afiles = []
        ls = os.listdir(dirname)
        occupancy = [ [] for j in range(self.params.n) ]

        for filename in ls:
            apat = r"^A(\d+)-(\d+)\.(\d+)-(\d+)$"
            if (m := re.match(apat, filename)):
                tup = (os.path.join(dirname, filename), *[int(x) for x in m.groups()])
                j0, j1, start, end = tup[1:]
                nj = j1 - j0
                if nj != self.params.splitwidth and nj != self.params.n:
                    raise ValueError(f"Wrong file name for A file (does not match splitwidth): {filename} {NOK}")
                for j in range(j0, j1, self.params.splitwidth):
                    occupancy[j].append((start, end, []))
                self.afiles.append(tup)
                continue

        for j in range(0, self.params.n, self.params.splitwidth):
            n = 0
            occupancy[j] = sorted(occupancy[j])
            for b,e,L in occupancy[j]:
                if b < n:
                    raise ValueError(f"overlap in A files for position {0:self.params.m},{j}, degrees {b}:{n}")
                n = e
        self.vfiles = scan_for_bwc_vectors(self.params, dirname)

        # scan the v files, determine how each range is going to be
        # tested.
        for v in self.vfiles:
            j0,j1,iteration,filename = v
            assert j1-j0 == self.params.splitwidth
            iL = iter(occupancy[j0])
            nL = []
            while (ell := next(iL, None)) is not None:
                s, e, usable = ell
                if iteration >= e:
                    nL.append(ell)
                elif iteration <= s:
                    usable.append(v)
                    nL.append(ell)
                else:
                    nL.append(s, iteration, usable)
                    nL.append(iteration, e, copy(usable) + v)
            occupancy[j0] = nL
        for j in range(0, self.params.n, self.params.splitwidth):
            nL = []
            for b,e,usable in occupancy[j]:
                if not usable:
                    print(f"Warning: we have no way to check coefficients {b}:{e} in column {j} {EXCL}")
                    continue
                nL.append((b, e, max(usable, key=lambda v: v.iteration)))
            occupancy[j] = nL

        self.check_map=[]
        by_vector = defaultdict(lambda: [])
        for j in range(0, self.params.n, self.params.splitwidth):
            for x in occupancy[j]:
                by_vector[x[2]].append(x[:2])
        for k in sorted(by_vector.keys(), key=lambda x:tuple(x)):
            L = sorted(by_vector[k])
            b0 = None
            e0 = None
            for b,e in L:
                if b0 is None:
                    b0 = b
                    assert k.iteration == b
                else:
                    if b != e0:
                        raise ValueError(f"Found a gap in the set of values that can be checked using {k.filename}: nothing for degrees [{n}:{b}]")
                e0 = e
            self.check_map.append((k, e0))

        self.x = BwcXVector(self.params, self.dims, os.path.join(self.dirname, "X"))

    def read_one_matrix(self, f, ni, nj):
        M = matrix(self.KP, ni, nj)
        for i in range(ni):
            for j in range(nj):
                b = bytearray(f.read(self.params.p_bytes))
                if not b:
                    if i or j:
                        raise IOError(f"Short read, could not read a {ni}x{nj} matrix at offset {f.tell()}")
                    return None
                M[i,j] = int.from_bytes(b, 'little')
        return M


    def read(self):
        x = self.KP.gen()
        self.x.read()
        for v in self.vfiles:
            v.read()
        for filename, j0, j1, start, end in self.afiles:
            ni = self.params.m
            nj = j1 - j0
            k = start
            st = os.stat(filename)
            bytes_per_mat = ni * (nj // self.params.splitwidth) * self.params.p_bytes
            nk = st.st_size // bytes_per_mat
            print(f"Reading {filename} (size: {ni}*{nj}, {nk} coeffs)")
            if nk != end - start:
                raise ValueError(f"{filename} has incorrect size (expected {nk*bytes_per_mat} for {nk} coefficients) {NOK}")

            with open(filename, 'rb') as f:
                while (M := self.read_one_matrix(f, ni, nj)) is not None:
                    self.A[:,j0:j1] += x**k * M
                    k += 1


    def __getitem__(self, i):
        return matrix(self.KP.base_ring(),
                      self.A.nrows(),
                      self.A.ncols(),
                      [t[i] for t in self.A.coefficients()])
    def check(self, M):
        for v, k1 in self.check_map:
            print(f"Check iterations [{v.iteration}:{k1}] of A using {v.filename}")
            w = v.V
            for k in range(v.iteration, k1):
                if self.x.X.transpose()*w != self[k][:,v.j0:v.j1]:
                    raise ValueError(f"Inconsistency in A files at coefficient {k}")
                w = M * w
            print(f"Check iterations [{v.iteration}:{k1}] of A using {v.filename} {OK}")
