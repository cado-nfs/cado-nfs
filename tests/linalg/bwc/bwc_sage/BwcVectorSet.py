from sage.matrix.special import block_matrix
from .BwcVector import BwcVector
from .BwcVectorSetBase import BwcVectorSetBase
from .tools import OK, EXCL


class BwcVectorSet(BwcVectorSetBase):
    def __init__(self, params, dirname):
        super().__init__(BwcVector, params, dirname)

    def block_by_iteration(self, it):
        """
        Returns the n-column matrix formed by the concatenation of the
        vectors at this iteration (assuming all are available on disk)
        """
        sw = self.params.splitwidth
        xn = self.params.n // sw
        blocks = sorted([v for v in self.Vs if v.iteration == it],
                        key=lambda v: v.j0)
        for xj in range(xn):
            assert blocks[xj].j0 == xj * sw

        assert len(blocks) == xn

        return block_matrix(1, xn, [v.V for v in blocks])

    def blocks_by_iteration(self):
        """
        This returns a dictionary of (block) matrices, each of them
        collecting the vectors at a given iteration (which is the key to
        that value)
        """

        V_blocks = dict()

        for it in set([v.iteration for v in self.Vs]):
            V_blocks[it] = self.block_by_iteration(it)

        return V_blocks

    def check(self, mat):
        for j0, j1 in sorted(list(set([(x.j0, x.j1) for x in self.Vs]))):
            isj0j1 = lambda x: x.j0 == j0 and x.j1 == j1  # noqa: E731
            W = sorted([x for x in self.Vs if isj0j1(x)],
                       key=lambda x: x.iteration)
            if len(W) == 1:
                print(f"Warning: only vector {W[0].filename} was found" +
                      f" for block {j0}-{j1}: no check possible {EXCL}")
                continue
            print(f"Checking {len(W)} vectors for block {j0}-{j1}")
            print("Iterations: " + ", ".join([str(x.iteration) for x in W]))
            iW = iter(W)
            v = next(iW)
            i = v.iteration
            V = v.V
            print(f"Check base is {v.filename}")
            for t in iW:
                doing = f"Checking {t.filename}"
                print(doing)
                V = mat**(t.iteration-i) * V
                i = t.iteration
                assert V == t.V
                print(f"{doing} ... {OK}")
