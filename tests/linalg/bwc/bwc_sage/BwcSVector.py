from .BwcParameters import BwcParameters
from .BwcVectorBase import BwcVectorBase
from .tools import OK, NOK
from .tools import mcoeff


class BwcSVector(BwcVectorBase):
    pattern = r"^S.sols(\d+)-(\d+).(\d+)-(\d+)$"
    _what = "partial sum"
    _whats = "partial sums"

    def __init__(self, params: BwcParameters, filename):
        super().__init__(params, filename)
        self.j0, self.j1, self.start, self.end = self._filename_data

    def __str__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"solutions {self.j0}-{self.j1}"
        it = f"iterations {self.start}..{self.end}"
        return f"{dims} {self._what} for {block} at {it}"

    def __repr__(self):
        dims = f"{self.V.nrows()}x{self.V.ncols()}"
        block = f"solutions {self.j0}-{self.j1}"
        it = f"iterations {self.start}..{self.end}"
        return f"{dims} {self._what} for {block} at {it}"

    def check(self, vs, MQ, F, return_vector_out=False):
        """
        This is not an easy check, since it hauls lots of stuff under the
        hood. Recompute and verify the partial sum from the following data
         - vs is a block matrix, typically made from the collection of
           vectors at a given iteration that matches the starting point
           for this partial sum
         - F is the matrix inside the BwcFFiles data
        (you most probably want to use the check method from
        BwcSVectorSet instead)
        """

        it = f"iterations {self.start}..{self.end}"
        block = f"solutions {self.j0}-{self.j1}"
        what = f"that the {self._what} for {block} at {it} is correct"
        print(f"Checking {what}")
        s = self.V.parent()()
        for k in range(self.start, self.end):
            s += vs * mcoeff(F, k)[:, self.j0:self.j1]
            vs = MQ * vs

        if self.V != s:
            raise ValueError(f"check failed {NOK}")
        print(f"Checking {what} ... {OK}")

        if return_vector_out:
            return vs
