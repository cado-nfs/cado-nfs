# You probably want to run this with export PYTHONPATH=./tests/sagemath

%load_ext autoreload
%autoreload 2

import subprocess
import cado_sage
from cado_sage import CadoPolyFile
from cado_sage import CadoPurgedFile
from cado_sage import CadoIndexFile
from cado_sage import CadoIdealsMapFile
from cado_sage import CadoIdealsDebugFile
from cado_sage import CadoNumberTheory
from cado_sage import CadoMontgomeryReductionProcess
from cado_sage.bwc import BwcParameters, BwcMatrix
from cado_sage.tools import cat_or_zcat
from collections import defaultdict
import os.path

cado_sage.set_verbose(True)

exponent=79
wdir="/tmp/blah"
name="c30"
number_of_solutions=3
par = BwcParameters(m=8,n=4,p=exponent,wordsize=64)
extra_checks = False

# not sure we'll ever need more.
prec = 53


matrixfile = f"{wdir}/{name}.sparse.bin"
solution_files = [f"{wdir}/{name}.bwc.{exponent}/K.sols{i}-{i+1}.0.txt"
                  for i in range(number_of_solutions)]
purgedfile = f"{wdir}/{name}.purged.gz"
purgedfile_withsm = f"{wdir}/{name}.purged_withsm{exponent}.gz"
indexfile = f"{wdir}/{name}.index.gz"
polyfile = f"{wdir}/{name}.poly"
idealsmapfile = f"{wdir}/{name}.ideals.gz"
# This is created by explain_indexed_relation -python -all
idealsdebugfile = f"{wdir}/{name}.debug-ideals.txt"


M = BwcMatrix(par, matrix=matrixfile, wdir=wdir)
M.read()
def get_solution_file(f):
    print(f"Reading {f}")
    return [ZZ(c) for c in open(f).readlines()]

ker = matrix([get_solution_file(f) for f in solution_files])

assert ker * M.M == 0


abpairs = CadoPurgedFile(purgedfile)

# read_abpairs() is sufficient for most of what we do, but for additional
# triple-checking, we're interested in the full relations as well.
abpairs.read_abpairs()


     


relsets = CadoIndexFile(indexfile)
relsets.read()
R = relsets.matrix()

poly = CadoPolyFile(polyfile)
poly.read()

K0 = poly.K[0]; alpha0 = K0.gen()
K1 = poly.K[1]; alpha1 = K1.gen()

unit_rank = [sum(K.signature()) - 1 for K in poly.K]

if sum(unit_rank) > number_of_solutions:
    raise ValueError(f"number_of_solutions={number_of_solutions} will not" +
                     f" disambiguate units (ranks: {unit_rank})")


nt = CadoNumberTheory(poly)

# It's perhaps a bit too much to ask, because sagemath will factor the
# discriminant, here.
OK0, OK1 = nt.maximal_orders()
J0, J1 = nt.J()


ideals_map = CadoIdealsMapFile(idealsmapfile)
ideals_map.read()
ideals_map_rev = ideals_map.reverse()

# read the ideals in sage format as well. This takes about 4 seconds for
# a {name}...
all_ideals = CadoIdealsDebugFile(poly, idealsdebugfile)
all_ideals.read()


# The matrix that the linear algebra step considers is the result of
# merge-dl. It is built from a subset of the columns (a subset of the
# ideals). More precisely, column i in the linear algebra matrix is
# all_ideals[ideals_map[i]], and if an ideal index j in
# range(len(all_ideals)) is equal to ideals_map[i] for some i, then that
# i is ideals_map_rev[j].

# We'll need to do some gymnastics between the different representation
# of the factored part of the relations. A relation involves some columns
# that end up in the sparse matrix, and some that do not. Our approach
# will be to construct two sparse matrices MB and ML. Filtering 
# amounts to left multiplication by the matrix R (the relsets matrix) and
# the linear algebra matrix will be R*MB, while R*ML will be zero.

# We also want to incorporate a fix so that the valuations at J1 (which
# is not a prime ideal) do not appear. This amounts to right multiplication
# by a square matrix. By doing so, we want to reach the point where each
# row is exactly the factorization of a-b*alpha_i, or the corresponding
# relset.

ncB = len(ideals_map)
ncL = len(all_ideals)

print("Fixing the contributions of J everywhere")
fixup_matrix = matrix(ZZ, ncB + ncL, ncB + ncL, sparse=True)
fixup_matrix.subdivide([],[ncB])

if all_ideals.index_of_J:
    def store_factors(j, J):
        for I,e in J.factor():
            cI = all_ideals.index(I)
            if (jI := ideals_map_rev.get(cI, None)) is not None:
                fixup_matrix[j, jI] = e
            else:
                # It may happen that one of the prime ideals dividing
                # J was filtered out during merge. In which case
                # we'll put it in the second part of the matrix.
                fixup_matrix[j, ncB + cI] = e

    if not all_ideals.has_merged_J:
        assert False
        # untested, but maybe it works.
        for c in all_ideals.index_of_J:
            J = all_ideals[c]
            if (j := ideals_map_rev.get(c, None)) is None:
                # seems very weird that we never encountered this J,
                # right?
                continue
            fixup_matrix[j, j] = 1
            store_factors(j, J)
    else:
        # This is the "merged J" case.
        c, = all_ideals.index_of_J
        j = ideals_map_rev.get(c, None)
        assert j is not None
        fixup_matrix[j, j] = 1
        # all_ideals[c] is a tuple
        for J in all_ideals[c]:
            store_factors(j, J)

MB_and_ML = block_matrix(1,2,[M.MZ, matrix(ZZ, M.MZ.nrows(), ncL, sparse=True)])
MB_and_ML.subdivide([],[ncB])
MB = MB_and_ML.subdivision(0,0)
ML = MB_and_ML.subdivision(0,1)
assert ML == 0 # by construction.

MB_and_ML_noJ = MB_and_ML - MB_and_ML * fixup_matrix
MB_and_ML_noJ.subdivide([],[ncB])
MB_noJ = MB_and_ML_noJ.subdivision(0,0)
ML_noJ = MB_and_ML_noJ.subdivision(0,1)



if extra_checks:
    print("Doing some extra checking that merge-dl did its job correctly")
    abpairs.read_relations()
    
    def apply_merge_renumbering(a,b,rel):
        """
        return a pair of sparse vectors, compatible with MB and ML
        """
        v = vector(ZZ, ncB + ncL, sparse=True)
        F = [ n.K.maximal_order().fractional_ideal(a-b*n.K.gen()) for n in nt ]
        for c in rel:
            if c in all_ideals.index_of_J:
                v[ideals_map_rev[c]] += 1
                if all_ideals.has_merged_J:
                    # it's a special case, really. It's done in a fairly
                    # messy way, with this idea of the combined column.
                    F[0] *= J0
                    F[1] *= J1
                else:
                    J = all_ideals[c]
                    side = 0 if I.number_field() == nt[0].K else 1
                    F[side] *= J
                    assert c in ideals_map_rev
                    v[ideals_map_rev[c]] += 1
                continue
            I = all_ideals[c]
            side = 0 if I.number_field() == nt[0].K else 1
            # only J is not prime, and it's certainly in ideals_map_rev
            assert I.is_prime()
            if c in ideals_map_rev:
                v[ideals_map_rev[c]] += 1
            else:
                v[ncB + c] += 1
            F[side] /= I
        for ff in F:
            assert ff == 1
        return v

    print("Re-checking all relations")
    # in essence, it's a matrix, of course. It's just not presented so.
    expanded_relations = [(a,b,apply_merge_renumbering(a,b,rel)) for a,b,rel in abpairs]

    # the relsets * the expanded relations should give us MZ

    print("Checking that the sparse matrix is indeed what it should be")
    e_cut = matrix([e[2][:ncB] for e in expanded_relations])
    e = matrix([v for a,b,v in expanded_relations])
    e_noJ = e - e * fixup_matrix

    assert M.M == R * e_cut
    # This happens mod p

    assert MB_and_ML == R * e
    # note that since ML == 0, this also ensures that all large primes
    # cancel.

    assert MB_and_ML_noJ == R * e_noJ
    # ML_noJ is not necessarily zero, but the kernel vector will cancel
    # it.

    # Note that if we had used M.M.lift_centered() instead of M.MZ in the
    # definition of MB_and_ML, we would have had some wraparound every
    # now and then.


print("Computing the column renumbering")

# We're only really interested in the valuations per number field. Since
# we now have all the information available, we can split the matrix in
# two blocks, one per number field.
#
# Note that because of the fixup operation, this renumbering will not
# necessarily preserve the number of columns: we'll remove J ideals, and
# add their factors that might be missing in the matrix MB currently.

matrix_renumbering = []

MLJ_cols = [[j for j in {j for (i,j),v in ML_noJ.items()}
             if all_ideals[j].number_field() == K]
            for K in poly.K]

# what are the indices within range(ncB) + range(ncB, ncB+ncL) that
# correspond to the given side ?
indices_per_field = [[i for i,j in enumerate(ideals_map)
                      if j in all_ideals.index_of_J or
                      all_ideals[j].number_field() == K]
                     + [ncB + j for j in MLJ_cols[side] ]
                     for side,K in enumerate(poly.K)]

ideals_per_field = [[ all_ideals[j] for j in list(ideals_map) + MLJ_cols[side]
                     if j in all_ideals.index_of_J or
                     all_ideals[j].number_field() == K] for side,K in
                    enumerate(poly.K)]

matrix_renumbering = block_matrix(1,2,
                                  [matrix(ZZ,
                                          ncB + ncL,
                                          len(indices),
                                          entries={(j,i):1
                                                   for i,j in
                                                   enumerate(indices)})
                                   for indices in indices_per_field])

if all_ideals.has_merged_J:
    assert len(list(matrix_renumbering[ideals_map_rev[all_ideals.index_of_J[0]]].items())) == 2

split_valuations_matrix = MB_and_ML_noJ * matrix_renumbering
split_valuations_matrix.subdivide(*matrix_renumbering.subdivisions())

# cut = split_valuations_matrix.subdivisions()[1][0]
# split_valuations_matrix += split_valuations_matrix[:,cut+ideals_per_field[1].index(J1)] * ambiguous
# split_valuations_matrix.subdivide(*matrix_renumbering.subdivisions())
# 

if extra_checks:
    print("Checking the ideal factorization of each relation-set (LONG)")
    for i,r in enumerate(R):
        big = ZZ['x'](1)
        x = ZZ['x'].gen()
        sval = vector(ZZ, ncB + ncL)
        for j,v in r.items():
            a,b,val = expanded_relations[j]
            big  *= (a-b*x)^v
            sval += v * val
        sval -= sval * fixup_matrix
        sval *= matrix_renumbering
        cut = matrix_renumbering.subdivisions()[1][0]
        svals = sval[:cut], sval[cut:]
        for side in [0,1]:
            K = nt[side].K
            alpha = K.gen()
            OK = K.maximal_order()
            evalbig = big(alpha)
            I = OK.fractional_ideal(evalbig)
            factored = Factorization([(ideals_per_field[side][j],v) for j,v in svals[side].items()])
            assert I == factored.prod()
            # print({all_ideals.index(ii):v for ii,v in (I/factored.prod()).factor()})
    print("")


# Use Schirokauer maps to guarantee an ell-th power. We do this modulo
# both fields, but it can well be that it's useless modulo one of them.
# It's not much of an issue, since in that case we're going to have a
# zero-rank block and that's all

if True:
    # not os.path.exists(purgedfile_withsm):
    print("Computing the Schirokauer maps with sage")
    if os.path.exists(purgedfile_withsm):
        print("Note: {purgedfile_withsm} is ignored because cado-nfs doesn't like small e yet")

    # we need to precompute the maps and not create them multiple times.
    sm = [ n.schirokauer_maps(exponent) for n in nt ]
    def sm_apply(m, x):
        L = [f(x) for f in m]
        return flatten(L, ltypes=(type(L[0]),))

    schirokauer_block = block_matrix([[matrix([sm_apply(sm[i],a - b * n.K.gen())
                                               for a,b,rel in abpairs])
                                       for i,n in enumerate(nt)]])

    S = ker * R * schirokauer_block
    assert S.subdivision(0,0).rank() <= unit_rank[0]
    assert S.subdivision(0,1).rank() <= unit_rank[1]

else:
    print("Fetching the Schirokauer maps from file")
    # Alternatively, we can get the SMs from the sm_append output.
    def get_schirokauer_block_from_file(f):
        rows = []
        i = 0
        print(f"Reading {f}")
        for t in cat_or_zcat(f):
            if t.startswith(b'#'):
                continue
            ab, rel, smdata = t.decode('ascii').split(':')
            a, b = [ZZ(c, 16) for c in ab.split(',')]
            smdata = [ZZ(c) for c in smdata.split(',')]
            assert (a,b) == abpairs[i][:2]
            rows.append(vector(GF(exponent), smdata))
            i += 1
        return matrix(rows)

    schirokauer_block = get_schirokauer_block_from_file(purgedfile_withsm)
    S = ker * R * schirokauer_block


if False:
    # The code here is much simpler, but it is also a lot slower. The
    # good thing is that we do compute the same thing! The
    # code here is not smart enough to cope with e dividing the norm of
    # a-b*alpha, so it clearly does not have the level of generality that
    # we get from n.schirokauer_map() -- it will probably crash on your
    # example if you try to run it.L
    sm2 = [n.schirokauer_map_simple(exponent) for n in nt]
    schirokauer_block2 = matrix([vector(sum([list(sm2[i](a - b * n.K.gen()))
                                        for i,n in enumerate(nt)
                                        ],[]))
                                 for a,b,rel in abpairs])
    S2 = ker * R * schirokauer_block2

    assert S.left_kernel() == S2.left_kernel()


# This is such that the product of (a,b) pairs with these exponents is
# now guaranteed to be an ell-th power. It's important to take a centered
# lift, as it minimizes the valuations that we get later on!
ker_power = (S.left_kernel().basis()[0] * ker).lift_centered()

# Now let's focus on one side in particular.
side = 1

K = poly.K[side]
OK = K.maximal_order()
d = K.degree()

ideals = ideals_per_field[side]

# all the valuations that we have for our current e-th power.
valuations = ker_power * split_valuations_matrix.subdivision(0,side) / exponent


L = nt[side].LogMap(prec)

E = matrix([ L(a - b * K.gen()) for a,b,rel in abpairs ])

# these are the log embeddings of our root
embeddings = ker_power * R * E / exponent

big_power = None

def compute_big_power():
    print("computing the absurdly big algebraic number (LONG!)")
    return prod([ (ab[0] - ab[1] * alpha1)**(ker_power * R)[i] for i,ab in enumerate(abpairs) ])

if extra_checks:
    # even for exponent=23, this takes an absurdly long time. In production,
    # we'll never do this.
    big_power = compute_big_power()

def check_invariant(MM, ideals, big_power):
    print("Doing consistency check on MM (LONG!)")
    for i,I in enumerate(ideals):
        if type(I) is tuple:
            continue
        else:
            vb = big_power.valuation(I)
            vm = MM.current.get(I,0)
            vg = MM.OK.fractional_ideal(MM.accumulated.prod()).valuation(I)
            if vb == (vm+vg)*exponent:
                continue
            print(f"Weird valuation at ideal {i} [[{I}]]: {vb/exponent},{vm}")
    print("Doing consistency check on MM: OK")

MM = CadoMontgomeryReductionProcess(poly, side, ideals, valuations, embeddings)
MM.status()

if extra_checks:
    check_invariant(MM, ideals, big_power)

while (MM.nplus + MM.nminus)[1] > 100:
    MM.one_reduction_step(256)
    MM.status()

# one last round.
MM.one_reduction_step(256)
MM.status()

gamma = MM.accumulated.prod()

if big_power is None:
    big_power = compute_big_power()

# Now that it seems that we do it correctly, we're apparently landing on
# 1 as a remaining e-th power, which is kinda cool.
if MM.bound_on_coefficients_of_remaining_part() == ([1] + [0] * (d-1)):
    print("Remaining part is 1, we don't have to compute it")
    delta = 1
    # it can also be -1. Which of the two it is can be decided by
    # computing modulo an inert p, but it's not clear that we want to
    # bother with doing so. After all it's only a pair of possible
    # solutions.
else:

    print("computing the remaining bit after all cancellations")
    rest = big_power / gamma**exponent

    print("computing the final e-th root")
    delta = rest.nth_root(exponent)

# then (gamma * delta) is an e-th root of big_power
print("Final check: ", big_power / (gamma*delta)**exponent, " (both +1 and -1 are fine)")
