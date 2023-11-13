import os
import re

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer

from .BwcParameters import BwcParameters
from .BwcVectorBase import BwcVectorBase
from .tools import OK, NOK, EXCL


class BwcCheckData(object):
    """
    This class gathers information on all the C[vrdt] vectors in a given
    bwc directory.

    Ct0-$nchecks.0-$m (also referred to as T): a random matrix of size
    bw->m*nchecks (yes it's named like this because the data is written
    row-major, and the first interval in the name customarily denotes the
    number of items in a major division).

    Cr0-$nchecks.0-$nchecks (also referred to as R): a sequence of random
    matrices of size nchecks * nchecks

    Cv0-<splitwidth>.<j> == trsp(M)^j * X * Ct

    Cd0-<splitwidth>.<j> == \\sum_{0<=i<j} trsp(M)^i * X * Ct * Cr[i]
    """
    def __init__(self, params: BwcParameters, dims, dirname):
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

        what = "set of check files in " + self.dirname
        for a in [self.tfiles, self.rfiles, self.dfiles, self.vfiles]:
            if len(a) == 0:
                print(f"No complete {what}, cannot verify {EXCL}")
                return

        if len(self.tfiles) != 1:
            raise ValueError(f"Conflicting Ct files among {what} {NOK}")

        if len(self.rfiles) != 1:
            raise ValueError(f"Conflicting Cr files among {what} {NOK}")

        self.nchecks = self.tfiles[0][1]
        assert self.nchecks % self.params.splitwidth == 0

        if self.tfiles[0][2] != self.params.m:
            fn = self.tfiles[0]
            tm = self.tfiles[0][2]
            m = self.params.m
            raise ValueError(f"{fn} is for m={tm} (!= {m}). {NOK}")

        sw = self.params.splitwidth

        for d in self.dfiles:
            if d[1] != sw:
                raise ValueError(f"{d[0]} is not for splitwidth={sw} {NOK}")

        for v in self.vfiles:
            if v[1] != sw:
                raise ValueError(f"{v[0]} is not for splitwidth={sw} {NOK}")

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
        t = BwcVectorBase(self.params,
                          filename,
                          pattern = r"^Ct(\d+)-(\d+).0-(\d+)",
                          _what = "bwc T check vector")
        t.read()
        if t.V.nrows() != self.params.m:
            expect = f"expected {self.params.m} rows"
            raise ValueError(f"{filename} has wrong size ({expect}) {NOK}")
        self.T = t.V

    def __read_rfile(self):
        """
        There is one single R file, but it's a sequence of arbitrary
        length
        """
        filename = os.path.join(self.dirname, self.rfiles[0][0])
        r = BwcVectorBase(self.params,
                          filename,
                          pattern = r"^Cr(\d+)-(\d+).0-(\2)",
                          _what = "bwc R check vector")
        r.read()
        sd = list(range(self.params.splitwidth,
                        r.V.nrows(),
                        self.params.splitwidth))
        r.V.subdivide(sd)
        self.R = [ r.V.subdivision(i,0) for i in range(len(sd)+1) ]

    def __read_dfiles(self):
        ret = []
        for basename, sw, j in self.dfiles:
            filename = os.path.join(self.dirname, basename)
            d = BwcVectorBase(self.params,
                              filename,
                              pattern = r"^Cd(\d+)-(\d+).(\d+)$",
                              _what = "bwc D check vector")
            d.read()
            ret.append((j, filename, d.V))
        self.D = sorted(ret)

    def __read_vfiles(self):
        ret = []
        for basename, sw, j in self.vfiles:
            filename = os.path.join(self.dirname, basename)
            v = BwcVectorBase(self.params,
                              filename,
                              pattern = r"^Cv(\d+)-(\d+).(\d+)$",
                              _what = "bwc V check vector")
            v.read()
            ret.append((j, filename, v.V))
        self.V = sorted(ret)

    def check(self, M, x):
        """
        There's no such thing as checking the t and r files, since
        they're pure random data. However we do want to check the v and d
        files.

        Cv0-<splitwidth>.<j> == trsp(M)^j * X * Ct

        Cd0-<splitwidth>.<j> == \\sum_{0<=i<j} trsp(M)^i * X * Ct * Cr[i]
        """

        xt = x.X * self.T
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
            what = "checking V and D files"
            raise ValueError(f"Errors found while {what} {NOK}")

    def read(self):
        """
        This reads all check files that have been found. Consistency is
        not checked, but consistency of the data size *IS* checked. We
        also check for possible extraneous files in some cases.
        """
        self.__read_tfile()
        self.__read_rfile()
        self.__read_dfiles()
        self.__read_vfiles()
