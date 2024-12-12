import gdb
from . import base
from . import display


class limbs_printer:
    def __init__(self, match, val, n):
        self.match = match
        assert val.type.code == gdb.TYPE_CODE_PTR
        self.val = val
        self.n = n

    def to_int(self):
        res = 0
        wordsize = self.val.type.target().sizeof * 8
        word = 1 << wordsize
        mp_limbs = self.val.reinterpret_cast(
            gdb.lookup_type("unsigned long").pointer())
        for i in range(0, self.n):
            a = int(mp_limbs[self.n-1-i])
            # python ints are signed, beware.
            if a < 0:
                a += word
            res = (res << wordsize) + a
        return res

    def to_string(self):
        return display.truncate_output(str(self.to_int()))


class mpz_printer:
    def __init__(self, match, val):
        self.match = match
        if val.type.code == gdb.TYPE_CODE_PTR:
            val = val.dereference()
        self.val = val

    def to_string(self):
        limbs = self.val['_mp_d']
        size = int(self.val['_mp_size'])
        mantissa = limbs_printer(self.match, limbs, abs(size)).to_int()
        if size < 0:
            mantissa = -mantissa
        # X = self.val
        # # There's apparently a bug in array member of structs. Their
        # # starting address seems to be constantly equal to the starting
        # # address of the struct itself...
        # # print "X at %s\n" % X.address
        # size = int(X['_mp_size'])
        # d = X['_mp_d']
        # nlimbs = size
        # if size < 0:
        #     nlimbs = -int(size)
        # # try:
        # mantissa = bigint_from_limbs(d, nlimbs)
        # # except RuntimeError:
        # # # it's not necessarily a good idea to do this.
        # # return "<error>"
        # if size < 0:y
        #     mantissa = -mantissa
        return display.truncate_output(str(mantissa))


class cxx_mpz_printer(mpz_printer):
    def __init__(self, match, val):
        self.match = match
        self.val = val['x']


class mpq_printer:
    def __init__(self, match, val):
        self.match = match
        self.num = mpz_printer(match, val['_mp_num'])
        self.den = mpz_printer(match, val['_mp_den'])

    def to_string(self):
        d = self.den.to_string()
        if d == "1":
            return self.num.to_string()
        else:
            return self.num.to_string() + "/" + d


class cxx_mpq_printer(mpq_printer):
    def __init__(self, match, val):
        val = val['x']
        super().__init__(match, val)
        self.match = match


base.cado_nfs_printer.add('__mpz_struct', mpz_printer)
base.cado_nfs_printer.add('mpz_t', mpz_printer)
base.cado_nfs_printer.add('mpz_ptr', mpz_printer)
base.cado_nfs_printer.add('mpz_srcptr', mpz_printer)
base.cado_nfs_printer.add('cxx_mpz', cxx_mpz_printer)
base.cado_nfs_printer.add('__mpq_struct', mpq_printer)
base.cado_nfs_printer.add('mpq_t', mpq_printer)
base.cado_nfs_printer.add('mpq_ptr', mpq_printer)
base.cado_nfs_printer.add('mpq_srcptr', mpq_printer)
base.cado_nfs_printer.add('cxx_mpq', cxx_mpq_printer)
