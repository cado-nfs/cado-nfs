import gdb
import re
from . import base
from . import display
from . import integers


def nlimbs(K, x):
    """
    returns the limbs per coefficient, where x is a pointer to
    coefficients and K is the context in which this is to be interpreted
    """
    T = x.type
    while T.code == gdb.TYPE_CODE_PTR or T.code == gdb.TYPE_CODE_TYPEDEF:
        if T.code == gdb.TYPE_CODE_PTR:
            nT = T.target()
            # print(f"== {T} -> {nT}")
            T = nT
        if T.code == gdb.TYPE_CODE_TYPEDEF:
            nT = T.strip_typedefs()
            # print(f"== {T} -> {nT}")
            T = nT
    # It's ul (mp_limb_t == unsigned long) on 64-bit, but u (mp_limb_t ==
    # unsigned long == unsigned int) on 32-bit !
    rx = r"^arith_modp::details::gfp_base<(\d+)(?:ul?)?,.*$"
    if m := re.match(rx, str(T)):
        w = int(m.group(1))
        if w:
            return w
        else:
            # This is for variable width
            w = int(K['p']['x']['_mp_size'])
            return w
    elif str(T) == "unsigned long":
        w = int(K['p']['x']['_mp_size'])
        return w
    else:
        tag = base.gdb_explain_type_codes.get(x.type.code, "UNKNOWN")
        raise RuntimeError(f"Cannot get the type width for {x.type} ({tag})"
                           f"; the furthest we got is {T}")
        # return 0


class flat_poly_printer:
    def __init__(self, match, val, K, W, ncoeffs):
        self.match = match
        self.ncoeffs = ncoeffs
        self.nlimbs = W
        assert val.type.code == gdb.TYPE_CODE_PTR
        self.val = val
        self.K = K

    def to_string(self):
        raw_ptr_type = gdb.lookup_type('unsigned long').pointer()
        X = self.val.cast(raw_ptr_type)
        W = self.nlimbs
        C = [integers.limbs_printer(self.match, X + i * W, W).to_int()
             for i in range(self.ncoeffs)]
        # form a polynomial. It's not totally trivial.
        S = ""
        for i, z in enumerate(C):
            if not z:
                continue
            if z > 0:
                if S:
                    S += "+"
            if i and (z == 1 or z == -1):
                S += "x"
                if i > 1:
                    S += f"^{i}"
            else:
                S += str(z)
                if i:
                    S += "*x"
                    if i > 1:
                        S += f"^{i}"
        if not S:
            S = "0"
        return display.truncate_output(S)


# We're not exposing this printer directly.
# base.cado_nfs_printer.add('__mpz_struct', mpz_printer)
# base.cado_nfs_printer.add('mpz_t', mpz_printer)
