from sage.rings.integer_ring import ZZ
from cado_sage.tools import get_verbose, cat_or_zcat
from collections import defaultdict


class CadoIdealsMapFile(object):
    """
    The cado-nfs ideals map file is documented in
    filter/README.fileformat (the file name is XXX.ideal)
    """

    def __init__(self, filename):
        self.filename = filename
        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self.ideals_map = []
        self.reverse_ideals_map = []

    def __repr__(self):
        return f"CadoIdealsMapFile(\"{self.filename}\")"

    def __str__(self):
        rep = "cado-nfs ideals map file"
        if self.filename:
            rep += f" (path={self.filename})"
        else:
            rep += " (transient)"
        if not self.ideals_map:
            rep += ", no data read yet"
            return rep
        else:
            rep += f", {len(self.ideals_map)} ideal mappings"
            return rep

    def read(self):
        self.__clear_fields_for_read()
        if get_verbose():
            print(f"Reading {self.filename}")

        contents = cat_or_zcat(self.filename)

        ideals_map_length = ZZ(next(contents)[1:])
        self.reverse_ideals_map = defaultdict(lambda:None)
        for ii, t in enumerate(contents):
            i, j = t.decode('ascii').split(' ')
            assert ii == ZZ(i)
            j = ZZ(j, 16)
            self.ideals_map.append(j)
            self.reverse_ideals_map[j] = ii
        assert len(self.ideals_map) == ideals_map_length

        self.reverse_ideals_map = dict(self.reverse_ideals_map)

    def __iter__(self):
        return iter(self.ideals_map)

    def __getitem__(self, i):
        return self.ideals_map[i]

    def __len__(self):
        return len(self.ideals_map)

    def reverse(self):
        return self.reverse_ideals_map
