from sage.rings.integer_ring import ZZ
from cado_sage.tools import get_verbose, cat_or_zcat


class CadoIdealsDebugFile(object):
    def __init__(self, poly, filename):
        self.poly = poly
        self.filename = filename
        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self._ideals = []

    def __repr__(self):
        return ("CadoIdealsDebugFile("
                + f"CadoPolyFile(\"{self.poly.filename}\")"
                + f", \"{self.filename}\")")

    def __str__(self):
        rep = "cado-nfs ideals debug file"
        if self.filename:
            rep += f" (path={self.filename})"
        else:
            rep += " (transient)"
        if not self.relsets:
            rep += ", no data read yet"
            return rep
        else:
            rep += f", {len(self.ideals)} ideals"
            return rep

    def read(self):
        self.__clear_fields_for_read()
        if get_verbose():
            print(f"Reading {self.filename}")

        K = self.poly.K
        J = self.poly.nt.J()
        OK = self.poly.nt.maximal_orders()

        for t in cat_or_zcat(self.filename):
            if t.startswith('#'):
                continue

            parser, *data = t.split(' ')

            if parser == 'J':
                if len(data) == 1:
                    I = J[int(data[0])]
                else:
                    raise ValueError("How can I deal with this consistently?")
            elif parser == 'rat':
                side, p = data
                side = int(side)
                p = ZZ(p)
                I = OK[side].fractional_ideal(p)
            elif parser == 'proj':
                side, p = data
                side = int(side)
                p = ZZ(p)
                I = OK[side].fractional_ideal(p) + J[side]
            elif parser == 'easy':
                side, p, r = data
                side = int(side)
                p = ZZ(p)
                r = ZZ(r)
                I = OK[side].fractional_ideal(p, K[side].gen() - r) * J[side]
            elif parser == 'generic':
                side, p, denom, *coeffs = data
                side = int(side)
                p = ZZ(p)
                denom = ZZ(denom)
                theta = K[side]([ZZ(c) for c in coeffs]) / denom
                I = OK[side].fractional_ideal(p, theta)

            self._ideals.append(I)

    def __iter__(self):
        return iter(self._ideals)

    def __len__(self):
        return len(self._ideals)

    def __getitem__(self, i):
        return self._ideals[i]
