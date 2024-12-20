from . import base      # noqa: F401
from . import display   # noqa: F401
from . import integers  # noqa: F401

# TODO

# class mpfr_printer:
#     def __init__(self, val):
#         self.val = val
#
#     def mantissa_exponent(self):
#         X = self.val
#         sign = int(X['_mpfr_sign'])
#         exp = int(X['_mpfr_exp'])
#         prec = int(X['_mpfr_prec'])
#         d = X['_mpfr_d']
#         wordsize = gdb.lookup_type("unsigned long").sizeof * 8
#         nlimbs = (prec+wordsize-1)/wordsize
#         # print "(via %s) Have prec %d, %d limbs\n" % (self.vt, prec,nlimbs)
#         # try:
#         mantissa = gmpy_from_mpn(d, nlimbs, wordsize)
#         mantissa *= sign
#         e = exp-nlimbs * wordsize
#         return mantissa, e
#
#     def to_string(self):
#         X = self.val
#         mantissa, e = self.mantissa_exponent()
#         prec = int(X['_mpfr_prec'])
#         exp = int(X['_mpfr_exp'])
#         wordsize = gdb.lookup_type("unsigned long").sizeof * 8
#         special = -pow(2, wordsize-1)
#         sign = int(X['_mpfr_sign'])
#         if exp == special+2:
#             return "NaN"
#         if exp == special+1:
#             if sign < 0:
#                 return "-0"
#             else:
#                 return "0"
#         if exp == special+3:
#             if sign < 0:
#                 return "-inf"
#             else:
#                 return "+inf"
#         res = mpfr(mantissa, prec)
#         if e > 0:
#             res *= pow(mpz(2), e)
#         else:
#             res /= pow(mpz(2), -e)
#         return truncate_output(str(res))
#
#
# class mpc_printer:
#     def __init__(self, val):
#         # print val['re'].address
#         self.re = mpfr_printer(val['re'].dereference())
#         self.im = mpfr_printer(val['im'].dereference())
#
#     def to_string(self):
#         return self.re.to_string() + "+i*" + self.im.to_string()
#
