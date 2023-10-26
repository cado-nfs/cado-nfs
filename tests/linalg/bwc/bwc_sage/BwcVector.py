from .BwcParameters import BwcParameters
from .BwcVectorBase import BwcVectorBase


class BwcVector(BwcVectorBase):
    pattern = r"^V(\d+)-(\d+).(\d+)$"
    _what = "bwc vector"
    _whats = "bwc vectors"

    def __init__(self, params: BwcParameters, filename):
        super().__init__(params, filename)
        self.j0, self.j1, self.iteration = self._filename_data

    def __str__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"block {self.j0}-{self.j1}"
        it = f"iteration {self.iteration}"
        return f"{dims} {self._what} for {block} at {it}"

    def __repr__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"block {self.j0}-{self.j1}"
        it = f"iteration {self.iteration}"
        return f"{dims} {self._what} for {block} at {it}"
