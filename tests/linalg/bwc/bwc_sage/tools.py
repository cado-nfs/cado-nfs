
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer

OK = "ok ✅"
NOK = "❌"
EXCL = "❗"
NOTHING_TO_DO = "💤"
HURRAH = "🥳"


def get_int(f, bytes=4, signed=False, repeat=None, may_fail=False):
    assert not (may_fail and repeat is not None)
    if repeat is not None:
        return [get_int(f, bytes=bytes, signed=signed) for i in range(repeat)]
    data = f.read(bytes)
    if not data and may_fail:
        return None
    assert data
    return int.from_bytes(data, 'little', signed=signed)


def u32(f, *args, **kwargs):
    return get_int(f, bytes=4, *args, **kwargs)


def u64(f, *args, **kwargs):
    return get_int(f, bytes=8, *args, **kwargs)


def s32(f, *args, **kwargs):
    return get_int(f, bytes=4, signed=True, *args, **kwargs)


def s64(f, *args, **kwargs):
    return get_int(f, bytes=8, signed=True, *args, **kwargs)


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
    return matrix(KP.base_ring(), m, n, [a[i] for a in M.coefficients()])


def mrev(M, n):
    return M.parent()([a.reverse(degree=n) for a in M.coefficients()])


def mdeg(M):
    return max([a.degree() for a in M.coefficients()])


def mmod(M, n):
    KP = M.base_ring()
    t = KP.gen()
    return M.parent()([a.mod(t**n) for a in M.coefficients()])


def mdiv(M, n):
    KP = M.base_ring()
    t = KP.gen()
    return M.parent()([a // t**n for a in M.coefficients()])


def mdivmod(M, n):
    return (mdiv(M, n), mmod(M, n))
