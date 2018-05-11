
# Simple python helper script to help printing gmp values within gdb.
#
# Emmanuel Thomé, 2013.

# This file is placed in the public domain

# Prerequisites:
#  - gdb 7.2 or later.
#  - gmpy2 bindings for python. On debian, this may mean either the
#  python-gmpy2 or python3-gmpy2 package. Your mileage may vary on other
#  systems. For a cross-system and realtively painless installation
#  procedure, you may consider:
#
#   easy_install --user gmpy2
#
# Usage: enable the pretty-printers of this file with the following gdb
# command:
#
#  source /path/to/this/file/gdb-gmp.py
#
# It is also possible to place this very command in your gdb startup file
# $HOME/.gdbinit.
#
# The code here currently reacts on the types mpz, mpfr, mpc, mpfrx
# (defined in their respective libraries, but it does not hurt to
# pretty-print them all). Please do not hesitate to expand to other
# types, and contribute to this file.
#
# Big integers being sometimes big, there is a natural risk of filling up
# terminal space. Therefore, this script has a limit on the size of
# numbers displayed, set by default to 40 characters. Numbers above this
# limit are printed as follows:
#
# Breakpoint 1, foo (x=12, z=68812...<6642>...9205452) at toto.c:6
#
# The display limit may be updated from within gdb as follows:
#
# (gdb) python update_display_limit(200)


# Referene on python scripting:
# http://sourceware.org/gdb/current/onlinedocs/gdb/Python.html

# Prefer gdb 7.2 or later for this, because of the following bug (which
# bites you if you have mpz_t's in structs):
# http://comments.gmane.org/gmane.comp.gdb.patches/59238
# Otherwise gdb 7.1 will work.

import gdb

import string
from gmpy2 import mpz,mpfr


__all__=['update_display_limit', 'make_mp_printer_objects', 'hook_mp_printers', 'unhook_mp_printers' ]

try:
    gdb.pretty_printers.remove(make_mp_printer_objects)
except NameError:
    pass
except ValueError:
    print("Caught ValueError, ignoring")
    pass


# arrange so that at most 40 characters of each field are printed.
display_limit=40

def gmpy_from_mpn(d, nlimbs, wordsize):
    res=mpz(0)
    word=pow(mpz(2),wordsize)
    for i in range(0,nlimbs):
        a=mpz(int(d[nlimbs-1-i]))
        # python ints are signed, beware.
        if a < 0: a+=word
        res=res*word+a
    return res

def truncate_output(r):
    # not completely satisfactory
    if display_limit == 0:
        return r
    l=len(r)
    if l >= display_limit:
        r=r[0:display_limit-15] +"...<%d>..."%(l-(display_limit-12)) +r[l-7:l]
    return r
    
def update_display_limit(l):
    global display_limit
    display_limit=l

class mpfr_printer:
    def __init__(self, val):
        self.val = val
    def mantissa_exponent(self):
        X = self.val
        sign = int(X['_mpfr_sign'])
        exp =  int(X['_mpfr_exp'])
        prec = int(X['_mpfr_prec'])
        d =    X['_mpfr_d']
        wordsize=gdb.lookup_type("unsigned long").sizeof * 8
        nlimbs=(prec+wordsize-1)/wordsize
        # print "(via %s) Have prec %d, %d limbs\n" % (self.vt, prec,nlimbs)
        # try:
        mantissa=gmpy_from_mpn(d, nlimbs, wordsize)
        mantissa*=sign
        e=exp-nlimbs*wordsize
        return mantissa, e
    def to_string(self):
        X = self.val
        mantissa,e = self.mantissa_exponent()
        prec = int(X['_mpfr_prec'])
        exp =  int(X['_mpfr_exp'])
        wordsize=gdb.lookup_type("unsigned long").sizeof * 8
        special=-pow(2,wordsize-1)
        sign = int(X['_mpfr_sign'])
        if exp == special+2:
            return "NaN"
        if exp == special+1:
            if sign < 0:
                return "-0"
            else:
                return "0"
        if exp == special+3:
            if sign < 0:
                return "-inf"
            else:
                return "+inf"
        res=mpfr(mantissa, prec)
        if e>0: res*=pow(mpz(2),e)
        else: res/=pow(mpz(2),-e)
        return truncate_output(str(res))

class mpz_printer:
    def __init__(self, val):
        self.val = val
    def to_string(self):
        X = self.val
        # There's apparently a bug in array member of structs. Their
        # starting address seems to be constantly equal to the starting
        # address of the struct itself...
        # print "X at %s\n" % X.address
        size = int(X['_mp_size'])
        d =    X['_mp_d']
        wordsize=gdb.lookup_type("unsigned long").sizeof * 8
        nlimbs=size
        if size<0: nlimbs=-int(size)
        # try:
        mantissa=gmpy_from_mpn(d, nlimbs, wordsize)
        # except RuntimeError:
        # # it's not necessarily a good idea to do this.
        # return "<error>"
        if size<0: mantissa=-mantissa
        return truncate_output(str(mantissa))

class mpq_printer:
    def __init__(self, val):
        self.num = mpz_printer(val['_mp_num'])
        self.den = mpz_printer(val['_mp_den'])
    def to_string(self):
        d = self.den.to_string()
        if d == "1":
            return self.num.to_string()
        else:
            return self.num.to_string() + "/" + d

class mpz_mat_printer:
    def __init__(self, val):
        self.val = val
        self.m = int(val['m'])
        self.n = int(val['n'])
    def to_string(self):
        s = []
        for i in range(self.m):
            for j in range(self.n):
                foo = mpz_printer(self.val['x'][i * self.n + j])
                s.append(foo.to_string())
        l=max([len(x) for x in s])
        fmt = "%%-%ds" % l
        out="\n"
        for i in range(self.m):
            for j in range(self.n):
                if j > 0:
                    out += " "
                out += (fmt % s[i * self.n + j])
            out += "\n"
        return out

class mpq_mat_printer:
    def __init__(self, val):
        self.val = val
        self.m = int(val['m'])
        self.n = int(val['n'])
    def to_string(self):
        s = []
        for i in range(self.m):
            for j in range(self.n):
                foo = mpq_printer(self.val['x'][i * self.n + j])
                s.append(foo.to_string())
        l=max([len(x) for x in s])
        fmt = "%%-%ds" % l
        out="\n"
        for i in range(self.m):
            for j in range(self.n):
                if j > 0:
                    out += " "
                out += (fmt % s[i * self.n + j])
            out += "\n"
        return out

class mpc_printer:
    def __init__(self, val):
        # print val['re'].address
        self.re = mpfr_printer(val['re'].dereference())
        self.im = mpfr_printer(val['im'].dereference())
    def to_string(self):
        return self.re.to_string() + "+i*" + self.im.to_string()

class mpfrx_printer:
    def __init__(self, val):
        self.val = val
    def to_string(self):
        X = self.val
        n=int(X['deg'])+1
        coeffs=[mpfr_printer(X['coeff'][i]).to_string() for i in range(n)]
        res="(%s %s)" % (n-1, " ".join(coeffs))
        return res

def make_mp_printer_objects(val):
    try:
        # beware. If we hook on mpfr_t to remove the extra braces, then
        # some bugs pop up -- apparently dereferencing this pointer does
        # not work as it should in arrays of mpfr_t's. Haven't checked,
        # but I assume it's similar for other types.
        t=str(val.type)
        # print("[request for %s]\n" % t)
        if (t == 'mpz_ptr' or t == 'mpz_srcptr'):
            return mpz_printer(val.dereference())
        if (t == 'mpz_t' or t == '__mpz_struct [1]'):
            return mpz_printer(val.dereference())
        if (t == '__mpz_struct'):
            return mpz_printer(val)

        if (t == 'mpq_ptr' or t == 'mpq_srcptr'):
            return mpq_printer(val.dereference())
        if (t == 'mpq_t' or t == '__mpq_struct [1]'):
            return mpq_printer(val.dereference())
        if (t == '__mpq_struct'):
            return mpq_printer(val)

        if (t == 'mpfr_ptr' or t == 'mpfr_srcptr'):
            return mpfr_printer(val.dereference())
        if (t == 'mpfr_t' or t == '__mpfr_struct [1]'):
            return mpfr_printer(val.dereference())
        if (t == '__mpfr_struct'):
            return mpfr_printer(val)

        if (t == 'mpc_ptr' or t == 'mpc_srcptr'):
            return mpc_printer(val.dereference())
        if (t == 'mpc_t' or t == '__mpc_struct [1]'):
            return mpc_printer(val.dereference())
        if (t == '__mpc_struct'):
            return mpc_printer(val)

        if (t == 'mpfrx_ptr' or t == 'mpfrx_srcptr'):
            return mpfrx_printer(val.dereference())
        if (t == 'mpfrx_t' or t == '__mpfrx_struct [1]'):
            return mpfrx_printer(val.dereference())
        if (t == '__mpfrx_struct'):
            return mpfrx_printer(val)

        if (t == 'mpz_mat_ptr' or t == 'mpz_mat_srcptr'):
            return mpz_mat_printer(val.dereference())
        if (t == 'mpz_mat'):
            return mpz_mat_printer(val.dereference())
        if (t == 'mpz_mat_s'):
            return mpz_mat_printer(val)

        if (t == 'mpq_mat_ptr' or t == 'mpq_mat_srcptr'):
            return mpq_mat_printer(val.dereference())
        if (t == 'mpq_mat'):
            return mpq_mat_printer(val.dereference())
        if (t == 'mpq_mat_s'):
            return mpq_mat_printer(val)



    except RuntimeError:
    # constructors may abandon building if the object looks too complicated.
        return None
    return None

def remove_all_printers():
    while len(gdb.pretty_printers):
        gdb.pretty_printers.pop(0)

def hook_mp_printers():
    gdb.pretty_printers.append(make_mp_printer_objects)

def unhook_mp_printers():
    gdb.pretty_printers.remove(make_mp_printer_objects)

# this is just for easily sourcing this again and again while debugging.
remove_all_printers()
hook_mp_printers()

