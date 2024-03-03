import os
import gdb

from . import base
from . import display       # noqa: F401
from . import integers
from . import polynomials


class generic_matrix(object):
    def __init__(self, m, n, field_name):
        self.m = m
        self.n = n
        self.K = field_name

    def to_string(self):
        m = self.m
        n = self.n
        if not self.m or not self.n:
            return "<empty>"

        def S(x):
            if type(x) is str:
                return x
            else:
                return x.to_string()

        s = [[S(self[i, j]) for j in range(n)] for i in range(m)]
        # ell = max([len(x) for x in s])
        ell = [max([len(s[i][j]) for i in range(m)]) for j in range(n)]
        if bool(os.environ.get("CADO_NFS_GDB_SHOW_MATRIX_COMPACT")):
            return "".join([f"matrix({self.K}, {m}, {n}, [",
                            ", ".join([
                                f"{s[i][j]}"
                                for i in range(m)
                                for j in range(n)
                                ]), "])"])
        else:
            return "\n".join([f"matrix({self.K}, {m}, {n}, [",
                              ",\n".join([*sum([
                                  [", ".join([
                                      "{0:>{1}}".format(s[i][j], ell[j])
                                      for j in range(n)])]
                                  for i in range(m)], [])])
                              + "])"])


class mpz_mat_printer(generic_matrix):
    def __init__(self, match, val):
        self.match = match
        self.val = val
        super().__init__(int(val['m']), int(val['n']), 'ZZ')

    def __getitem__(self, ij):
        i, j = ij
        n = self.n
        return integers.mpz_printer(self.match, self.val['x'][i * n + j])


class cxx_mpz_mat_printer(mpz_mat_printer):
    def __init__(self, match, val):
        val = val['x']
        super().__init__(match, val)
        self.match = match


class mpq_mat_printer(generic_matrix):
    def __init__(self, match, val):
        self.val = val
        self.match = match
        super().__init__(int(val['m']), int(val['n']), 'QQ')

    def __getitem__(self, ij):
        i, j = ij
        n = self.n
        return integers.mpq_printer(self.match, self.val['x'][i * n + j])


class cxx_mpq_mat_printer(mpq_mat_printer):
    def __init__(self, match, val):
        val = val['x']
        super().__init__(match, val)
        self.match = match


class matpoly_binary_printer(generic_matrix):
    def __init__(self, match, val):
        self.val = val
        self.match = match
        super().__init__(int(val['m']), int(val['n']), 'KP')

    def __getitem__(self, ij):
        i, j = ij
        n = self.n
        W = int(self.val['alloc_words'])
        xij = self.val['x'] + (i * n + j) * W
        # read just one integer, and interpret it as a polynomial
        z = integers.limbs_printer(self.match, xij, W)
        return f"0x{z.to_int():x}"


class matpoly_nonbinary_printer(generic_matrix):
    def __init__(self, match, val):
        self.val = val
        self.match = match
        super().__init__(int(val['m']), int(val['n']), 'KP')

    def __getitem__(self, ij):
        i, j = ij
        n = self.n
        a = int(self.val['alloc'])
        K = self.val['ab'].dereference()
        W = polynomials.nlimbs(K, self.val['x'])
        raw_ptr_type = gdb.lookup_type('unsigned long').pointer()
        x0 = self.val['x'].cast(raw_ptr_type)
        xij = x0 + (i * n + j) * a * W
        # return f"{x0}+({i*n+j}*{a}) == {xij}"
        s = int(self.val['size'])
        return polynomials.flat_poly_printer(self.match, xij, K, W, s)


base.cado_nfs_printer.add('mpz_mat_s', mpz_mat_printer)
base.cado_nfs_printer.add('mpz_mat', mpz_mat_printer)
base.cado_nfs_printer.add('mpz_mat_ptr', mpz_mat_printer)
base.cado_nfs_printer.add('mpz_mat_srcptr', mpz_mat_printer)
base.cado_nfs_printer.add('cxx_mpz_mat', cxx_mpz_mat_printer)
base.cado_nfs_printer.add('mpq_mat_s', mpq_mat_printer)
base.cado_nfs_printer.add('mpq_mat', mpq_mat_printer)
base.cado_nfs_printer.add('mpq_mat_ptr', mpq_mat_printer)
base.cado_nfs_printer.add('mpq_mat_srcptr', mpq_mat_printer)
base.cado_nfs_printer.add('cxx_mpq_mat', cxx_mpq_mat_printer)
base.cado_nfs_printer.add('matpoly<true>', matpoly_binary_printer)
base.cado_nfs_printer.add('matpoly<false>', matpoly_nonbinary_printer)
