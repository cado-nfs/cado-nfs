import sys

from sage.rings.integer import Integer
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from .tools import u32, u64
from .tools import NOK
from .BwcParameters import BwcParameters


class BwcBalancing(object):
    def __init__(self, params: BwcParameters, filename):
        self.params = params
        self.filename = filename

    def __flag_bytes_to_dict(self, flags):
        txflags = dict()
        if flags & 1:
            txflags['colperm'] = 1
        if flags & 2:
            txflags['rowperm'] = 1
        if flags & 8:
            txflags['replicate'] = 1
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

            # copy b111d37a5: now the alignment is always on multiples of 8
            chunk = 8

            self.tr = pad(self.nr, self.nh * self.nv, chunk)
            self.tc = pad(self.nc, self.nh * self.nv, chunk)

            if 'replicate' in self.txflags:
                self.tr = max(self.tr, self.tc)
                self.tc = max(self.tr, self.tc)

            s = self.nh * self.nv

            self.tr = pad(self.tr, s)
            self.tc = pad(self.tc, s)

            rowcols = f"{self.nr} rows {self.nc} cols"
            split = f"split {self.nh}x{self.nv}"
            checksum = f"checksum 0x{self.checksum:x}"
            flags = f"{list(self.txflags.keys())}"
            print(", ".join([rowcols, split, checksum, flags]))

            self.rowperm = None
            if 'rowperm' in self.txflags:
                self.rowperm = u32(f, repeat=self.tr)

            self.colperm = None
            if 'colperm' in self.txflags:
                self.colperm = u32(f, repeat=self.tc)


class BwcShuffling(object):
    def __init__(self, par: BwcParameters, mat):
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

        try:
            nc = mat.ncols
            Sc0 = SymmetricGroup(range(nc))
            self.shuf = Sc0([self.__preshuf(Integer(x)) for x in range(nc)])
            self.shufinv = self.shuf**-1
        except KeyError:
            print("Error in creating permutation shuf " + NOK,
                  file=sys.stderr)
            self.shuf = None
            self.shufinv = None

        try:
            tc = bal.tc
            self.xshuf = Sc([self.__preshuf(Integer(x)) for x in range(tc)])
            self.xshufinv = self.xshuf**-1
        except KeyError:
            print("Error in creating permutation xshuf " + NOK,
                  file=sys.stderr)
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
        q, r = x.quo_rem(nz)
        qq, qr = q.quo_rem(bal.nv)
        return (qr * bal.nh + qq) * nz + r

    def __prinv(self, x):
        bal = self.M.balancing
        nz = bal.tr // (bal.nh * bal.nv)
        q, r = x.quo_rem(nz)
        qq, qr = q.quo_rem(bal.nh)
        return (qr * bal.nv + qq) * nz + r
