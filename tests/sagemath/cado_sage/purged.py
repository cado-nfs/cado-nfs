from sage.rings.integer_ring import ZZ
from cado_sage.tools import get_verbose, cat_or_zcat


class CadoPurgedFile(object):
    """
    The cado-nfs purged file is documented in filter/README.fileformat
    """

    def __init__(self, filename):
        self.filename = filename
        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self.ab = []
        self.relations = []

    def __repr__(self):
        return f"CadoPurgedFile(\"{self.filename}\")"

    def __str__(self):
        rep = "cado-nfs purged file"
        if self.filename:
            rep += f" (path={self.filename})"
        else:
            rep += " (transient)"
        if not self.ab:
            rep += ", no data read yet"
            return rep
        elif not self.relations:
            rep += f", {len(self.ab)} pairs, relations not read yet"
            return rep
        else:
            rep += f", {len(self.ab)} pairs with relations"
            return rep

    def read_abpairs(self):
        self.__clear_fields_for_read()
        if get_verbose():
            print(f"Reading {self.filename} (only (a,b) pairs)")
        for t in cat_or_zcat(self.filename):
            if t.startswith(b'#'):
                continue
            ab, rel = t.decode('ascii').split(':')
            a, b = [ZZ(c, 16) for c in ab.split(',')]
            self.ab.append((a, b))

    def read_relations(self):
        self.__clear_fields_for_read()
        if get_verbose():
            print(f"Reading {self.filename}")
        for t in cat_or_zcat(self.filename):
            if t.startswith(b'#'):
                continue
            ab, rel = t.decode('ascii').split(':')
            a, b = [ZZ(c, 16) for c in ab.split(',')]
            rel = [ZZ(c, 16) for c in rel.split(',')]
            self.ab.append((a, b))
            self.relations.append(rel)

    def __iter__(self):
        if self.relations:
            for i, ab in enumerate(self.ab):
                yield *ab, self.relations[i]
        else:
            for ab in self.ab:
                yield *ab, None

    def __getitem__(self, i):
        if type(i) is slice:
            start, stop, step = i.indices(len(self.ab))
            return [ self[i] for i in range(start, stop, step) ]
        else:
            if self.relations:
                return *self.ab[i], self.relations[i]
            else:
                return *self.ab[i], None

    def __len__(self):
        return len(self.ab)
