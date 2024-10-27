from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from cado_sage.tools import get_verbose, cat_or_zcat


class CadoIndexFile(object):
    """
    The cado-nfs index file is documented in filter/README.fileformat
    """

    def __init__(self, filename):
        self.filename = filename
        self.__clear_fields_for_read()

    def __clear_fields_for_read(self):
        self.relsets = []

    def __repr__(self):
        return f"CadoIndexFile(\"{self.filename}\")"

    def __str__(self):
        rep = "cado-nfs index file"
        if self.filename:
            rep += f" (path={self.filename})"
        else:
            rep += " (transient)"
        if not self.relsets:
            rep += ", no data read yet"
            return rep
        else:
            rep += f", {len(self.relsets)} relation-sets"
            return rep

    def read(self):
        self.__clear_fields_for_read()
        if get_verbose():
            print(f"Reading {self.filename}")

        contents = cat_or_zcat(self.filename)

        relsets_length = ZZ(next(contents))

        for t in contents:
            if t.startswith(b'#'):
                continue
            line_length, *entries = t.decode('ascii').split(' ')
            assert len(entries) == ZZ(line_length)
            self.relsets.append([tuple([ZZ(col, 16), ZZ(val)])
                                 for col, val in
                                 [e.split(':') for e in entries]])
        assert len(self.relsets) == relsets_length

    def matrix(self):
        """
        makes an integer matrix from the relsets
        """
        if not self.relsets:
            raise ValueError("No data read from {self.filename}," +
                             " call .read() first")
        return matrix(dict(sum([
            [((i, j), v) for j, v in e]
            for i, e in enumerate(self.relsets)], [])))

    def __iter__(self):
        return iter(self.relsets)

    def __getitem__(self, i):
        return self.relsets[i]

    def __len__(self):
        return len(self.relsets)
