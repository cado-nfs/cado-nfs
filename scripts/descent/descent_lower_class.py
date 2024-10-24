import os
import re
import subprocess

from .ideals_above_p import ideals_above_p
from .descent_utils import a_over_b_mod_p
from .descent_utils import check_result


class DescentLowerClass(object):
    def declare_args(parser):
        pass

    def __init__(self, general, args):
        self.general = general
        self.args = args

    def __count_multiplicites(self, L):
        LL = []
        prev_p = L[0]
        m = 1
        for i in range(1, len(L)):
            p = L[i]
            if p == prev_p:
                m += 1
            else:
                LL.append([prev_p, m])
                prev_p = p
                m = 1
        LL.append([prev_p, m])
        return LL

    def do_descent(self, relsfile, initial_split):
        general = self.general
        # args = parser.parse_args()
        tmpdir = general.tmpdir()
        prefix = f"{general.prefix()}.descent." \
                 f"{general.short_target()}.lower."
        relsforSM = os.path.join(tmpdir, prefix + "relsforSM")
        SMfile = os.path.join(tmpdir, prefix + "SM")

        # Read descent relations
        descrels = []
        for rfile in relsfile:
            with open(rfile, 'r') as file:
                with open(relsforSM, 'a') as fileSM:
                    for line in file:
                        foo = re.match(r"^Taken: (-?\d+),(-?\d+):", line)
                        if foo:
                            r = line.split(':')[1:]
                            r[0] = r[0].lstrip()
                            fileSM.write(r[0] + ":" + r[1] + ":" + r[2])
                            a, b = r[0].split(',')
                            a = int(a)
                            b = int(b)
                            list_p = [ [], [] ]
                            for side in range(2):
                                for p in r[side + 1].strip().split(','):
                                    list_p[side].append(int(p, 16))
                            list_p = [self.__count_multiplicites(list_p[0]),
                                      self.__count_multiplicites(list_p[1])]
                            descrels.append(([a, b], list_p))
        nrels = len(descrels)
        print("--- Final reconstruction (from %d relations) ---" % nrels)

        # Compute SM
        call_that = [general.sm_simple_bin(),
                     "-poly", general.poly(),
                     "-inp", relsforSM,
                     "-out", SMfile,
                     "-ell", general.ell()
                     ]
        if self.args.sm_mode is not None:
            call_that += [ "-sm-mode", self.args.sm_mode ]
        call_that = [str(x) for x in call_that]
        print("command line:\n" + " ".join(call_that))
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(call_that, stderr=devnull)

        SM = []
        with open(SMfile, 'r') as file:
            for line in file:
                r = line.split()
                sm = [ int(x) for x in r ]
                SM.append(sm)
        assert len(SM) == nrels

        # Reverse the order of relations to get only one unknown log
        # per relation while processing them
        descrels.reverse()
        SM.reverse()

        # Fill-in the log database
        logDB = general.logDB
        irel = 0
        for rel in descrels:
            unk = None
            a, b = rel[0]
            list_p = rel[1]
            acc_log = 0
            if logDB.fullcolumn is not None:
                acc_log += logDB.fullcolumn
            sm = SM[irel]
            for i in range(len(sm)):
                acc_log += logDB.allSM(i) * sm[i]
            for side in range(2):
                for p, k in list_p[side]:
                    ideal = ideals_above_p(p, k, a, b, side, general)
                    log = ideal.get_log()
                    if log is None:
                        if unk is not None:
                            raise ValueError(
                                "Two unknown ideals in relation"
                                " a,b=%d,%d: %d (side %d)"
                                " and %d (side %d)"
                                % (a, b, unk[0], unk[2], p, side))
                        else:
                            unk = [p, a_over_b_mod_p(a, b, p), side]
                    else:
                        acc_log += ideal.get_log()
            acc_log = acc_log % general.ell()
            if unk is None:
                assert acc_log == 0
            else:
                log = general.ell() - acc_log
                print("Deduced log of (%d, %d, %d) from rel: %d"
                      % (unk[0], unk[1], unk[2], log))
                logDB.add_log(unk[0], unk[1], unk[2], log)
            irel += 1

        if general.has_rational_side():
            # Deduce the log of the target
            Num, Den, factNum, factDen = initial_split
            log_target = 0
            errors = []
            for p in factNum:
                lp = logDB.get_log(p, -1, 0)
                if lp is None:
                    errors.append(p)
                else:
                    log_target = log_target + lp
            for p in factDen:
                lp = logDB.get_log(p, -1, 0)
                if lp is None:
                    errors.append(p)
                else:
                    log_target = log_target - lp
            if len(errors):
                msg = "Some logarithms missing:\n"
                msg += "\n".join(["\t" + str(x) for x in errors])
                raise RuntimeError(msg)
            p = general.p()
            ell = general.ell()
            log_target = log_target % ell
            if general.initrandomizer != 1:
                # divide result by randomizer modulo ell
                multiplier = pow(general.initrandomizer, ell - 2, ell)
                log_target = (log_target * multiplier) % ell
            print("# p=%d" % p)
            print("# ell=%d" % ell)
            print("log(2)=%d" % logDB.get_log(2, -1, 0))
            print("log(3)=%d" % logDB.get_log(3, -1, 0))
            print("# target=%s" % self.args.target)
            print("log(target)=%d" % log_target)
            check_result(2,
                         logDB.get_log(2, -1, 0),
                         int(self.args.target),
                         log_target, p, ell)
        else:
            # No rational side; more complicated.
            # We need to compute the SMs for U and V.
            polyforSM = os.path.join(tmpdir, prefix + "polyforSM")
            SM2file = os.path.join(tmpdir, prefix + "SM2")

            with open(polyforSM, 'w') as f:
                for poly in [ general.initU, general.initV ]:
                    f.write("p %d" % (len(poly) - 1))
                    for c in poly:
                        f.write(" %d" % c)
                    f.write("\n")

            call_that = [general.sm_simple_bin(),
                         "-poly", general.poly(),
                         "-inp", polyforSM,
                         "-out", SM2file,
                         "-ell", general.ell()
                         ]
            call_that = [str(x) for x in call_that]
            print("command line:\n" + " ".join(call_that))
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(call_that, stderr=devnull)

            SM2 = []
            with open(SM2file, 'r') as file:
                for line in file:
                    r = line.split()
                    sm = [ int(x) for x in r ]
                    SM2.append(sm)
            assert len(SM2) == 2

            ell = general.ell()
            vlog = [0, 0]
            factored = [ general.initfacu, general.initfacv ]
            for i in range(0, 2):
                for xx in factored[i]:
                    vlog[i] += logDB.get_log(xx[0], xx[1], general.init_side)
                ind_shift = 0
                if general.init_side == 1:
                    ind_shift = logDB.nSM(0)
                for j in range(logDB.nSM(general.init_side)):
                    logsm = logDB.SM(general.init_side, j)
                    sm = SM2[i][ind_shift + j]
                    vlog[i] += logsm * sm
                vlog[i] = vlog[i] % ell

            log_target = (vlog[0] - vlog[1]) % ell
            multiplier = pow(general.initrandomizer, ell - 2, ell)
            log_target = (log_target * multiplier) % ell
            print("# p=%d" % general.p())
            print("# ell=%d" % ell)
            print("# target=%s" % self.args.target)
            print("log(target)=%d" % log_target)
