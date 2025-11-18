import gdb
import gdb.printing


debug = False


gdb_explain_type_codes = dict()

for c in ("PTR ARRAY STRUCT UNION ENUM FLAGS FUNC INT FLT VOID"
          " SET RANGE STRING BITSTRING ERROR METHOD METHODPTR MEMBERPTR"
          " REF RVALUE_REF CHAR BOOL COMPLEX TYPEDEF NAMESPACE DECFLOAT"
          " INTERNAL_FUNCTION XMETHOD FIXED_POINT NAMESPACE").split():
    # gdb_type_codes[c] = eval(f"gdb.TYPE_CODE_{c}")
    try:
        gdb_explain_type_codes[eval(f"gdb.TYPE_CODE_{c}")] = c
    except Exception as e:
        print("Exception", type(e), e)


class CadoNFSSubPrinter(object):
    def __init__(self, name, function, **kwargs):
        self.name = name
        self._function = function
        self.kwargs = kwargs

    def invoke(self, value):
        if value.type.code == gdb.TYPE_CODE_REF:
            value = value.referenced_value()

        return self._function(self.name, value, **self.kwargs)


class CadoNFSPrinter(object):
    def __init__(self, name):
        super().__init__()
        self.name = name
        self._lookup = {}
        self.enabled = True

    def add(self, name, function, **kwargs):
        if debug:
            print(f"# registered pretty printer for {name}")
        self._lookup[name] = CadoNFSSubPrinter(name, function, **kwargs)

    @staticmethod
    def get_basic_type(type):
        # If it points to a reference, get the reference.
        if type.code == gdb.TYPE_CODE_REF:
            type = type.target()

        # Get the unqualified type, stripped of typedefs.
        type = type.unqualified().strip_typedefs()

        return type.tag

    def __call__(self, val):
        tag = gdb_explain_type_codes.get(val.type.code, "UNKNOWN")
        bt = gdb.types.get_basic_type(val.type)
        # print(gdb.types.get_type_recognizers())
        if debug:
            print(f"[request for {val.type} ({tag}) {bt}]")

        # for C types, it's convenient to be able to match on the
        # typedef.
        if (s := self._lookup.get(str(val.type))) is not None:
            return s.invoke(val)

        typename = self.get_basic_type(val.type)
        if not typename:
            return None

        if debug:
            print(f"[dereferenced to {typename}]")

        if val.type.code == gdb.TYPE_CODE_REF:
            if hasattr(gdb.Value, "referenced_value"):
                val = val.referenced_value()

        if (s := self._lookup.get(typename)) is not None:
            return s.invoke(val)

        return None


cado_nfs_printer = CadoNFSPrinter("cado-nfs")


# class mpfrx_printer:
#     def __init__(self, val):
#         self.val = val
#
#     def to_string(self):
#         X = self.val
#         n = int(X['deg']) + 1
#         coeffs = [mpfr_printer(X['coeff'][i]).to_string() for i in range(n)]
#         res = "(%s %s)" % (n - 1, " ".join(coeffs))
#         return res
#

# class abpelt_printer:
#     def __init__(self, val, width):
#         self.val = val
#         self.width = width
#     def to_string(self):
#         mantissa=gmpy_from_mpn(self.val, self.width, 64)
#         return str(mantissa)
#
#
# class polymat_printer:
#     def __init__(self, val, width):
#         self.val = val
#         self.width = width
#         self.m = val['m']
#         self.n = val['n']
#     def to_string(self):
#         X = self.val
#         m = self.m
#         n = self.n
#         s = ""
#         if X['size'] == 0:
#             return "0"
#         for k in range(X['size']):
#             z = True
#             t = ""
#             t = t + "Matrix(KP,%d,%d,[" % (m,n)
#             ptr = X['x']
#             for i in range(m):
#                 for j in range(n):
#                     if str(ptr.type) == "abase_pz_vec":
#                         pp = ptr + (k*m*n+i*n+j) * self.width
#                     else:
#                         pp = ptr[k*m*n+i*n+j]
#                     nb = abpelt_printer(pp, self.width).to_string()
#                     if (i,j) != (0,0):
#                         t = t + ", "
#                     t = t + nb
#                 t = t + "\n"
#             t = t + "])"
#             if k==1:
#                 s = s + "+X*"
#             elif k>1:
#                 s = s + "+X^%d*" % k
#             s = s + t
#         return s
#
#
# last_fti=None
#
#
# class fft_transform_info_printer:
#     def __init__(self, val):
#         self.fti_n = 2**val['depth']
#         self.fti_w = val['w']
#         self.val=val
#         global last_fti
#         last_fti = self
#     def to_string(self):
#         return "depth-%d FFT description object in Z/2^(%d*%d)+1" % (
#             self.val['depth'], 2**self.val['depth'], self.val['w'])
#     pass
#
#
# def ssfft_coeff(d, rsize0, wordsize):
#     res=mpz(0)
#     word=pow(mpz(2),wordsize)
#     for i in range(0,rsize0+1):
#         a=mpz(int(d[rsize0+1-1-i]))
#         # python ints are signed, beware.
#         if a < 0: a+=word
#         res=res*word+a
#     if res >= pow(mpz(2),wordsize*rsize0-1):
#         res -= (pow(word,rsize0)+1)
#     return res
#
#
# class ss_transform_printer:
#     def __init__(self, val):
#         self.val = val
#     def to_string(self):
#         global last_fti
#         if last_fti is None or last_fti.fti_n == 0 or last_fti.fti_w == 0:
#             return "<undef>"
#         ptr = self.val
#         wordsize = 64
#         sizeof_limbs = wordsize // 8
#         rsize0 = (last_fti.fti_n * last_fti.fti_w + wordsize - 1) // wordsize
#         rsize0 = int(rsize0)
#         N = last_fti.fti_n
#         nptrs = 4*N+2
#         entrysize = (nptrs + nptrs * (rsize0 + 1)) * sizeof_limbs
#         mat=[]
#         base = ptr
#         ptrpool = base.cast(gdb.lookup_type("mp_limb_t").pointer().pointer())
#         coeffs=[]
#         for k in range(4*N):
#             q = ptrpool[k]
#             coeffs.append(str(ssfft_coeff(q, rsize0, 64)))
#         return "[" + ", ".join(coeffs) + "]"
#
#
# class pft(gdb.Function):
#   def __init__(self):
#     super(pft, self).__init__("pft")
#   def invoke(self, var):
#     return ss_transform_printer(var).to_string()
#
# pft()
#
#
# class matpoly_ft_printer:
#     def __init__(self, val):
#         self.val = val
#         self.m = val['m']
#         self.n = val['n']
#     def to_string(self):
#         global last_fti
#         if last_fti is None or last_fti.fti_n == 0 or last_fti.fti_w == 0:
#             return "<undef>"
#         X = self.val
#         ptr = X['data']
#         m = self.m
#         n = self.n
#         wordsize = 64
#         sizeof_limbs = wordsize // 8
#         rsize0 = (last_fti.fti_n * last_fti.fti_w + wordsize - 1) // wordsize
#         rsize0 = int(rsize0)
#         N = last_fti.fti_n
#         nptrs = 4*N+2
#         entrysize = (nptrs + nptrs * (rsize0 + 1)) * sizeof_limbs
#         mat=[]
#         print("---- %d ----\n" % rsize0)
#         for i in range(m):
#             row=[]
#             for j in range(n):
#                 base = ptr + entrysize*(i*n+j)
#                 ptrpool = base.cast(
#                     gdb.lookup_type("mp_limb_t").pointer().pointer())
#                 coeffs=[]
#                 for k in range(4*N):
#                     q = ptrpool[k]
#                     coeffs.append(str(ssfft_coeff(q, rsize0, 64)))
#                 row.append("[" + ", ".join(coeffs) + "]")
#             mat.append("[" + ", ".join(row) + "]")
#         return "[" + ",\n".join(mat) + "]"
#
#
# global_nlimbs=None
# global_nlimbs_default=1
#
# def set_global_nlimbs(x):
#     global global_nlimbs
#     global_nlimbs=x
#
#
#
#
# def make_mp_printer_objects(val):
#     global global_nlimbs
#     if global_nlimbs is None:
#         print("#### warning: setting nlimbs to %d #####"
#             % global_nlimbs_default)
#         print("#### use \"python set_global_nlimbs(xxx)\" to change it ####")
#         set_global_nlimbs(global_nlimbs_default)
#     nlimbs = global_nlimbs
#
#     try:
#         # beware. If we hook on mpfr_t to remove the extra braces, then
#         # some bugs pop up -- apparently dereferencing this pointer does
#         # not work as it should in arrays of mpfr_t's. Haven't checked,
#         # but I assume it's similar for other types.
#         t = str(val.type)
#         tag = gdb_explain_type_codes.get(val.type.code, "UNKNOWN")
#         bt = gdb.types.get_basic_type(val.type)
#         print(gdb.types.get_type_recognizers())
#         print(f"[request for {t} ({tag}) {bt}]")
#
#         # if t == 'abase_p_{nlimbs}_dst_elt':
#         #     return abpelt_printer(val, nlimbs)
#         # if t == 'abase_p_{nlimbs}_src_elt':
#         #     return abpelt_printer(val, nlimbs)
#         # if (t == 'abase_p_%d_elt'%nlimbs):
#         #     return abpelt_printer(val, nlimbs)
#         # if (t == 'polymat'):
#         #     return polymat_printer(val.dereference(), nlimbs)
#         # if (t == 'polymat_ptr'):
#         #     return polymat_printer(val.dereference(), nlimbs)
#         # if (t == 'struct polymat_s *'):
#         #     return polymat_printer(val.dereference(), nlimbs)
#         # if (t == 'polymat_ur'):
#         #     return polymat_printer(val.dereference(), 2*nlimbs)
#         # if (t == 'struct polymat_ur_s *'):
#         #     return polymat_printer(val.dereference(), 2*nlimbs)
#         # if (t == 'matpoly'):
#         #     return matpoly_printer(val, nlimbs)
#         # if (t == 'matpoly &'):
#         #     return matpoly_printer(val.dereference(), nlimbs)
#         # if (t == 'matpoly_ptr'):
#         #     return matpoly_printer(val.dereference(), nlimbs)
#         # if (t == 'struct matpoly_s *'):
#         #     return matpoly_printer(val.dereference(), nlimbs)
#         # if (t == 'matpoly_ur'):
#         #     return matpoly_printer(val.dereference(), 2*nlimbs)
#         # if (t == 'struct matpoly_ur_s *'):
#         #     return matpoly_printer(val.dereference(), 2*nlimbs)
#         # if (t == 'matpoly_ft'):
#         #     return matpoly_ft_printer(val.dereference())
#         # if (t == 'matpoly_ft_ptr'):
#         #     return matpoly_ft_printer(val.dereference())
#         # if (t == 'struct matpoly_ft_s *'):
#         #     return matpoly_ft_printer(val.dereference())
#         # if (t == 'matpoly_ft_ur'):
#         #     return matpoly_ft_printer(val.dereference())
#         # if (t == 'struct matpoly_ft_ur_s *'):
#         #     return matpoly_ft_printer(val.dereference())
#         # if (t == 'struct fft_transform_info *'):
#         #     return fft_transform_info_printer(val.dereference())
#         # if (t == 'fft_transform_info *'):
#         #     return fft_transform_info_printer(val.dereference())
#
#         # if (t == 'mpz_ptr' or t == 'mpz_srcptr'):
#         #     return mpz_printer(val.dereference())
#         # if (t == 'mpz_t' or t == '__mpz_struct [1]'):
#         #     return mpz_printer(val.dereference())
#         # if (t == '__mpz_struct'):
#         #     return mpz_printer(val)
#         # if (t == 'cxx_mpz'):
#         #     return mpz_printer(val['x'])
#
#         # if (t == 'mpq_ptr' or t == 'mpq_srcptr'):
#         #     return mpq_printer(val.dereference())
#         # if (t == 'mpq_t' or t == '__mpq_struct [1]'):
#         #     return mpq_printer(val.dereference())
#         # if (t == '__mpq_struct'):
#         #     return mpq_printer(val)
#         # if (t == 'cxx_mpq'):
#         #     return mpq_printer(val['x'])
#
#         # if (t == 'mpfr_ptr' or t == 'mpfr_srcptr'):
#         #     return mpfr_printer(val.dereference())
#         # if (t == 'mpfr_t' or t == '__mpfr_struct [1]'):
#         #     return mpfr_printer(val.dereference())
#         # if (t == '__mpfr_struct'):
#         #     return mpfr_printer(val)
#
#         # if (t == 'mpc_ptr' or t == 'mpc_srcptr'):
#         #     return mpc_printer(val.dereference())
#         # if (t == 'mpc_t' or t == '__mpc_struct [1]'):
#         #     return mpc_printer(val.dereference())
#         # if (t == '__mpc_struct'):
#         #     return mpc_printer(val)
#
#         # if (t == 'mpfrx_ptr' or t == 'mpfrx_srcptr'):
#         #     return mpfrx_printer(val.dereference())
#         # if (t == 'mpfrx_t' or t == '__mpfrx_struct [1]'):
#         #     return mpfrx_printer(val.dereference())
#         # if (t == '__mpfrx_struct'):
#         #     return mpfrx_printer(val)
#
#         # if (t == 'mpz_mat_ptr' or t == 'mpz_mat_srcptr'):
#         #     return mpz_mat_printer(val.dereference())
#         # if (t == 'mpz_mat'):
#         #     return mpz_mat_printer(val.dereference())
#         # if (t == 'mpz_mat_s'):
#         #     return mpz_mat_printer(val)
#         # if (t == 'cxx_mpz_mat'):
#         #     return mpz_mat_printer(val['x'])
#
#         # if (t == 'mpq_mat_ptr' or t == 'mpq_mat_srcptr'):
#         #     return mpq_mat_printer(val.dereference())
#         # if (t == 'mpq_mat'):
#         #     return mpq_mat_printer(val.dereference())
#         # if (t == 'mpq_mat_s'):
#         #     return mpq_mat_printer(val)
#         # if (t == 'cxx_mpq_mat'):
#         #     return mpq_mat_printer(val['x'])
#
#     except RuntimeError:
#         # constructors may abandon building if the object looks too
#         # complicated.
#         return None
#     return None
#
#
# def remove_all_printers():
#     while len(gdb.pretty_printers):
#         gdb.pretty_printers.pop(0)
#
#
# def hook_mp_printers():
#     gdb.pretty_printers.append(make_mp_printer_objects)
#
#
# def unhook_mp_printers():
#     gdb.pretty_printers.remove(make_mp_printer_objects)
#
#
# # this is just for easily sourcing this again and again while debugging.
# remove_all_printers()
# hook_mp_printers()
