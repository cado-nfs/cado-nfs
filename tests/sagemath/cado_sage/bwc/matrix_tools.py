from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer
from cado_sage.tools import NOK


def read_one_matrix(params, f, ni, nj):
    M = matrix(GF(params.p), ni, nj)
    for i in range(ni):
        for j in range(0, nj, params.splitwidth):
            b = bytearray(f.read(params.p_bytes))
            if not b:
                if i or j:
                    what = f"could not read a {ni}x{nj} matrix"
                    where = f"at offset {f.tell()}"
                    raise IOError(f"Short read, {what} {where} {NOK}")
                return None
            z = int.from_bytes(b, 'little')
            if params.splitwidth == 1:
                M[i, j] = z
            elif params.splitwidth == 64 and params.p == 2:
                for k, bit in enumerate(Integer(z).bits()):
                    M[i, j + k] = bit

    return M


def mcoeff(M, i):
    KP = M.base_ring()
    m = M.nrows()
    n = M.ncols()
    return matrix(KP.base_ring(), m, n, [a[i] for a in M.list()])


def mrev(M, n):
    return M.parent()([a.reverse(degree=n) for a in M.list()])


def mdeg(M):
    return max([a.degree() for a in M.list()])

def mval(M):
    return min([a.valuation() for a in M.list()])


def mmod(M, n):
    KP = M.base_ring()
    t = KP.gen()
    return M.parent()([a.mod(t**n) for a in M.list()])


def mdiv(M, n):
    KP = M.base_ring()
    t = KP.gen()
    return M.parent()([a // t**n for a in M.list()])


def mdivmod(M, n):
    return (mdiv(M, n), mmod(M, n))
