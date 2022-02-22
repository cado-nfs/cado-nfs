

def cado_encoding_of_projective_line_element(a, b, q):
    """
    returns an integer in [0,2q[ that encodes (a:b) in P1(Z/q)
    """
    g = gcd(b, q)
    if g == 1:
        a = Integers()(GF(q)(a/b))
        return a
    assert gcd(a, q) == 1
    return q + Integers()(GF(q)(b/a))


def primary_basis(r, q):
    """
    given r, which is the cado encoding of some element of q(Z/q),
    return the "straightforward" basis of the corresponding lattice

    sage: primary_basis(49,49)
    [ 1  0]
    [ 0 49]
    sage: primary_basis(56,49)
    [7 0]
    [1 7]
    sage: primary_basis(768,512)
    [  2   0]
    [  1 256]
    sage: primary_basis(13,37)
    [37  0]
    [13  1]
    """

    if r < q:
        return matrix([(q,0), (r, 1)])

    # otherwise it's significantly more complicated. we want a basis that
    # fits well our reduction routine... In particular we want L[(0,0)]
    # to be at least q
    gs = r - q
    L0 = matrix([(1, gs), (0,q)])
    g, t, h = xgcd(gs, q)
    s = gs / g
    if t < 0:
        t += q//g
        h -= s
    assert s*t+h*q//g == 1
    L = matrix([(q//g, 0), (t, g)])
    assert L == matrix([(q//g, -s), (t, h)]) * L0
    return L

def reduce_plattice(L, I):
    """
    This reduces the input lattice (which has to meet certain conditions,
    see below) following the procedure of FrKl05. We also cover special
    cases such that _all_ primes, prime powers, projective primes and so
    on can be covered.

    sage: reduce_plattice(primary_basis(0,881),512)
    [  0   1]
    [881   0]
    sage: reduce_plattice(primary_basis(5632,16384),512)
    [  0  32]
    [512   3]
    sage: reduce_plattice(primary_basis(3844,31^3),512)
    [  0  31]
    [961   8]
    sage: reduce_plattice(primary_basis(233,2187),512)
    [-323    8]
    [ 233    1]
    sage: reduce_plattice(primary_basis(231,2187),512)
    [-339    8]
    [ 231    1]
    sage: reduce_plattice(primary_basis(10092,29^3),512)
    [  0  29]
    [841  17]

    """


    assert L[0,1] == 0
    assert L[0,0] > 0
    assert L[1,1] > 0
    assert L[1,0] >= 0
    assert L.determinant() >= I

    needs_special_treatment = L[1,0] == 0 or L[1,1] > 1 and L[0,0] < I
    # note that L[0,0] < I implies L[1,1] >= I/L[0,0] > 1
    # L[1,0] is for the case of vertical lattices, which as a matter of
    # fact we handle the exact same way

    def handle_orthogonal_lattice(L):
        # in this case, the "good" vector is clear, but it does come
        # from the normal reduction procedure
        u0, u1 = L
        # let k be the smallest integer such that k * u0[0] >= I
        # let k be the largest integer such that k1 * u0[0] < I
        # let k be the largest integer such that k1 * u0[0] <= I-1
        #  k1 = (I-1) // u0[0]
        #  u2 = vector([-k1 * u0[0], 0]) 
        a = (I-1) % u0[0]
        u2 = vector([a-(I-1),0])    # this is a multiple of u0
        u2 += u1
        # This is a bit of a kludge, but I can't see how this test can be
        # avoided.
        if u0[0] - u2[0] < I:
            u2 -= u0
        Lx = matrix([u2,u0])
        assert (Lx*L^-1).determinant()^2 == 1
        return Lx
    if needs_special_treatment:
        return handle_orthogonal_lattice(L)
    else:
        assert L[0, 0] >= I
        mi0 =  L[0,0] ; j0 = 0       # encodes vector ( i0, j0) with i0=-mi0
        i1  =  L[1,0] ; j1 = L[1,1]  # encodes vector ( i1, j1)
        while True:
            if i1 < I:
                assert mi0 >= I
                assert i1 >= 0
                if i1 == 0:
                    # then it's really a straight vertical lattice, with
                    # no return vector at all. Note that we only have one
                    # vertical vector a priori, and no horizontal one.
                    # But we don't really care.
                    Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                    return handle_orthogonal_lattice(Lo)
                # let a be the smallest integer such that mi0 - a * i1 < I
                # i.e. 1 + the largest integer such that mi0 - a * i1 >= I
                # i.e. 1 + (mi0 - I) // i1
                # a = (mi0 - I) // i1
                # assert mi0 - a * i1 >= I
                # a += 1
                # assert mi0 - a * i1 < I
                a = (mi0 + i1 - I) // i1
                Lx = matrix([ ((-mi0) + a * i1, j0 + a * j1), (i1, j1) ])
                break
            k,mi0 = mi0.quo_rem(i1)
            j0 += k * j1
            if mi0 < I:
                assert i1 >= I
                assert mi0 >= 0
                if mi0 == 0:
                    Lo=matrix([ (i1, j1), (-mi0, j0)])
                    return handle_orthogonal_lattice(Lo)
                # let a be the smallest integer such that i1 - a * mi0 < I
                # i.e. 1 + the largest integer such that i1 - a * mi0 >= I
                # i.e. 1 + (i1 - I) // mi0
                # a = (i1 - I) // mi0
                # assert i1 - a * mi0 >= I
                # a += 1
                # assert i1 - a * mi0 < I
                a = (mi0 + i1 - I) // mi0
                Lx = matrix([ (-mi0, j0), (i1 - a * mi0, j1 + a * j0) ])
                break
            k,i1 = i1.quo_rem(mi0)
            j1 += k * j0

    return Lx

def plattice_post(Lx, I):
    """
    This verifies that reduce_plattice obeys its output conditions.

    sage: plattice_post(reduce_plattice(primary_basis(0,881),512),512)
    [  0   1]
    [881   0]
    sage: plattice_post(reduce_plattice(primary_basis(5632,16384),512),512)
    [  0  32]
    [512   3]
    sage: plattice_post(reduce_plattice(primary_basis(3844,31^3),512),512)
    [  0  31]
    [961   8]
    sage: plattice_post(reduce_plattice(primary_basis(233,2187),512),512)
    [-323    8]
    [ 233    1]
    sage: plattice_post(reduce_plattice(primary_basis(231,2187),512),512)
    [-339    8]
    [ 231    1]
    sage: plattice_post(reduce_plattice(primary_basis(10092,29^3),512),512)
    [  0  29]
    [841  17]

    """

    assert -I < Lx[0,0] <= 0
    assert Lx[0,1] > 0
    assert Lx[1,1] >= 0 # because we want to handle the projective case
    assert  0 <= Lx[1,0]
    # This assertion is possibly violated for vertical lattices
    # assert Lx[1,0] < I
    assert Lx[1,0] < I or Lx[0,0] == 0
    assert Lx[1,0] - Lx[0,0] >= I
    return Lx

def test_reduce_plattice(I, ntests=10, quiet=False):
    """
    This is just doing a few randomized tests
    
    sage: test_reduce_plattice(512, ntests=10^3, quiet=True)
    """
    B=I^2
    NN=[]
    while True:
        k=len(NN)+1
        y=float(log_integral(B^(1/k)))
        x=float(log_integral(I^(1/k)))
        if y < 1:
            break
        NN.append(ZZ(round(y-x)))
    sum_NN=sum(NN)
    def random_exponent():
        i=ZZ.random_element(sum_NN)
        for k in range(len(NN)):
            if i < NN[k]:
                return k+1
            else:
                i -= NN[k]
        assert False

    stats={}
    for t in range(ntests):
        k=random_exponent()
        p=random_prime(floor(B^(1/k)), lbound=ceil(I^(1/k)))
        q=p^k
        r=ZZ.random_element(q+p^(k-1))
        if k not in stats:
            stats[k] = {'affine':0, 'proj':0, 'affine+even':0, 'proj+even': 0 }
        desc='affine'
        if r >= q:
            r = q + p * (r - q)
            desc='proj'
        if p == 2:
            desc += "+even"
        try:
            Lx=plattice_post(reduce_plattice(primary_basis(r, q), I), I)
            stats[k][desc] += 1
        except Exception as e:
            print(f"Failed test (#{t}) for p^k={p}^{k}={q}, r={r}: {e}")
    if quiet:
        return
    for k in sorted(stats.keys()):
        s=stats[k]
        print(f'k={k}: ' + ', '.join([ f'{u}:{s[u]}' for u in sorted(s.keys()) ]))
