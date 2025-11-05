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


class mpz_poly_printer:
    def __init__(self, match, val, cxx=False):
        self.match = match
        self.val = val
        self.cxx = cxx
        if not self.cxx:
            t = "mpz_poly"
            tx = "cxx_" + t
            F = gdb.parse_and_eval(f"({tx}*) malloc(sizeof({tx}))")
            bt = val.type.strip_typedefs()
            gdb.execute(f"call (({tx}*){F})->{tx}()")
            if bt.code == gdb.TYPE_CODE_PTR:
                src = f"({t}_srcptr) {self.val.dereference().address}"
            else:
                src = f"({t}_srcptr) {self.val.address}"
            gdb.execute(f"call (void) {t}_set((({tx}*) {F})->x, {src})")
            self.cxx_val = gdb.parse_and_eval(f"* ({tx}*) {F}")
            self._cxx_temp = str(F)
        else:
            self.cxx_val = self.val

    def __del__(self):
        if not self.cxx:
            F = self._cxx_temp
            gdb.execute(f"call ((cxx_mpz_poly*) {F})->~cxx_mpz_poly()")
            gdb.execute(f"call free({F})")

    def to_string(self):
        """
        This implementation is just a test. We know that we can do it in
        a way that would be similar to the other pretty-printers, but
        this time we're going to try to trampoline via the code's own
        functions.
        """
        obj = f"(('cxx_mpz_poly' const *){self.cxx_val.address})"
        obj = f"{obj}->print_poly()"
        return eval(gdb.parse_and_eval(obj).format_string(max_elements=0))


class mpz_poly_bivariate_printer:
    def __init__(self, match, val):
        self.match = match
        self.val = val

    def to_string(self):
        """
        This implementation is just a test. We know that we can do it in
        a way that would be similar to the other pretty-printers, but
        this time we're going to try to trampoline via the code's own
        functions.
        """
        obj = f"(('cxx_mpz_poly_bivariate::self' const *){self.val.address})"
        obj = f"{obj}->print_poly()"
        return eval(gdb.parse_and_eval(obj).format_string(max_elements=0))


# it's really wreaking havoc with gdb.

if False:
    # We're not exposing flat_poly_printer directly.
    base.cado_nfs_printer.add('mpz_poly_ptr', mpz_poly_printer)
    base.cado_nfs_printer.add('mpz_poly_srcptr', mpz_poly_printer)
    base.cado_nfs_printer.add('mpz_poly_s', mpz_poly_printer)
    base.cado_nfs_printer.add('mpz_poly', mpz_poly_printer)
    base.cado_nfs_printer.add('cxx_mpz_poly', mpz_poly_printer, cxx=True)
    base.cado_nfs_printer.add('cxx_mpz_poly_bivariate',
                              mpz_poly_bivariate_printer)
