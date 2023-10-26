from .BwcSVector import BwcSVector
from .BwcVectorSetBase import BwcVectorSetBase
from .tools import NOTHING_TO_DO, OK, NOK


class BwcSVectorSet(BwcVectorSetBase):
    def __init__(self, params, dirname):
        super().__init__(BwcSVector, params, dirname)

    def check(self, Vs, MQ, F, known_solutions=None):
        """
        This is not an easy check, since it hauls lots of stuff under the
        hood. Recompute and verify all partial sum from the following data
         - vs is the BwcVectorSet data
         - F is the BwcFFiles data
        """

        # TODO (maybe): compute with a Horner scheme instead (or maybe do
        # both and check that we're doing it right).

        sw = self.params.splitwidth
        xn = self.params.n // sw

        V_blocks = Vs.blocks_by_iteration()

        # It's a bit different here, because we want to check solution
        # blocks separately

        S_blocks = [[] for xj in range(xn)]
        for s in self.Vs:
            S_blocks[s.j0 // sw].append(s)
        for xj in range(xn):
            S_blocks[xj] = sorted(S_blocks[xj], key=lambda x: tuple(x))

        for xj in range(xn):
            j0 = xj * sw
            j1 = j0 + sw
            doing = f"Checking solutions {j0}-{j1}"

            # We're a bit lazy here. In most occasions, when we're doing
            # very simple checks, all intermediary vectors will still be
            # there, so there's no problem in relying on them being
            # present on disk. If some vectors happen to be missing, we
            # may conceivably chain the computations one after another
            # with the return_vector_out argument of BwcSVector.check

            if S_blocks[xj]:
                print(doing)
                for s in S_blocks[xj]:
                    s.check(V_blocks[s.start], MQ, F.F)

                if known_solutions is not None:
                    doing += " against known data"
                    U, V = known_solutions
                    print(doing)
                    if sum([s.V for s in S_blocks[xj]]) != V[:, j0:j1]:
                        raise ValueError("check failed " + NOK)
                    print(f"{doing} ... {OK}")
            else:
                print(f"{doing}: no partial sum found {NOTHING_TO_DO}")
