from .BwcVector import BwcVector
from .BwcVectorSetBase import BwcVectorSetBase
from .tools import OK, EXCL


class BwcVectorSet(BwcVectorSetBase):
    def __init__(self, params, dirname):
        super().__init__(BwcVector, params, dirname)

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
