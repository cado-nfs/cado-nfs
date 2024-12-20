from .descent_utils import a_over_b_mod_p
from .descent_utils import is_a_over_b_equal_r_mod_pk


class ideals_above_p(object):
    def __is_badideal(self, p, r, side, general):
        return (p, r, side) in general.list_badideals

    def __handle_badideal(self, p, k, a, b, r, side, general):
        baddata = general.badidealdata
        expo = []
        badid = None
        for X in baddata:
            if side != X[3] or p != X[0]:
                continue
            pk = pow(p, X[1])
            rk = X[2]
            if not is_a_over_b_equal_r_mod_pk(a, b, rk, p, pk):
                continue
            vals = X[4]
            badid = (p, r, side)
            for v in vals:
                if v > 0:
                    exp = v
                else:
                    assert k >= -v
                    exp = k + v
                expo.append(exp)
            break
        if badid is None:
            raise ValueError("Error while handling badideal"
                             " p=%d side=%d a=%d b=%d" % (p, side, a, b))
        return {"badid": badid, "exp": expo}

    def __init__(self, p, k, a, b, side, general):
        self.general = general
        self.logDB = general.logDB
        self.p = p
        self.k = k
        self.side = side
        if general.has_rational_side() and side == 0:
            self.r = -1
        else:
            self.r = a_over_b_mod_p(a, b, p)
        self.isbad = self.__is_badideal(p, self.r, side, general)
        if self.isbad:
            self.bads = self.__handle_badideal(p, k, a, b,
                                               self.r, side, general)

    # return an unreduced virtual log or None if unknown.
    def get_log(self):
        general = self.general
        if not self.isbad:
            l = self.logDB.get_log(self.p, self.r, self.side)
            if l is None:
                return None
            else:
                return self.k * l
        else:
            ind_b = 0
            ind_l = 0
            while general.list_badideals[ind_b] != (self.p, self.r, self.side):
                ind_l += general.list_ncols[ind_b]
                ind_b += 1
            logs = [self.logDB.bad_ideal(x)
                    for x in range(ind_l, ind_l + general.list_ncols[ind_b]) ]
            assert len(logs) == len(self.bads["exp"])
            log = 0
            for i in range(len(logs)):
                log = log + logs[i] * self.bads["exp"][i]
            return log
