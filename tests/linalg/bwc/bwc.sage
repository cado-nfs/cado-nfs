import re
import os
import sys
from bwc_sage import *
from bwc_sage.tools import HURRAH, OK, NOK

if __name__ == '__main__':
    args = dict()
    for v in sys.argv[1:]:
        if (m := re.match(r"^(p)(?:rime)?=(\d+)$", v)):
            args[m.groups()[0]] = Integer(m.groups()[1])
        elif (m := re.match(r"^(m|n|nh|nv|wordsize)=(\d+)$", v)):
            args[m.groups()[0]] = int(m.groups()[1])
        elif (m := re.match(r"^(wdir|matrix)=(.*)$", v)):
            args[m.groups()[0]] = m.groups()[1]
        else:
            raise ValueError(f"Unknown parameter {v}")
    if "wdir" not in args:
        args["wdir"] = os.path.dirname(args["matrix"])

    def filter_dict(D, pat):
        return dict([kv for kv in D.items() if re.match(pat, kv[0])])

    par = BwcParameters(**filter_dict(args, r"^([mnp]|wordsize)$"))

    M = BwcMatrix(par, **filter_dict(args, r"^(matrix|wdir)$"))
    M.read(force_square=True)
    M.fetch_balancing(**filter_dict(args, r"^n[hv]$"))
    M.check()
    MQ = M.decorrelated()

    Vs = BwcVectorSet(par, args["wdir"])
    Vs.read()
    Vs.check(MQ)

    x = BwcXVector(par, M.dimensions(), args["wdir"])
    x.read()

    c = BwcCheckData(par, M.dimensions(), args["wdir"])
    c.read()
    c.check(MQ, x)

    a = BwcAFiles(par, M.dimensions(), args["wdir"])
    a.read()
    a.check(MQ, x, Vs)

    f = BwcFFiles(par, M.dimensions(), args["wdir"])
    f.read()
    f.check(a)
    U, V = f.derive_solutions(a, Vs, MQ)

    s = BwcSVectorSet(par, args["wdir"])
    s.read()
    s.check(Vs, MQ, f, known_solutions=(U, V))

    print("All checks passed " + HURRAH)
