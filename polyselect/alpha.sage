
# Code for computing alpha, and some variations around it.
# The main function exported here is alpha(f,p)

# We also provide splits alpha_affine and alpha_projective.
# alpha_projective, roughly said, depends only on the top coefficients
# (except for very special cases where touching a low degree coefficients
# might change the value).

def number_of_roots(f,p):
    """
    Counts the roots of f mod p, without multiplicities. Projective roots
    are also counted (without multiplicities).
    """
    fp= GF(p)['x'](f)
    if fp != 0:
        s=len([r[1] for r in fp.roots()])
        if (f.degree()>fp.degree()):
            s+=1
    else:
        # f reducing to zero mod p is a degenerate case. Not clear what
        # we should return...
        print "Warning, counting roots of zero polynomial\n"
        s=f.degree()
    return s

def number_of_roots_sub(f,p,r,s):
    """
    Counts the roots of f equal to r/s mod p, without multiplicities.
    Projective roots are also counted (without multiplicities).
    """
    fp = GF(p)['x'](f)
    t = 0
    if s != 0: # r/s is non-projective
        r = (r/s) % p
        for x,_ in fp.roots():
            if x == r:
                return 1
        return 0
    else: # we are looking for projective roots
        if fp.degree() < f.degree():
            return 1
        else:
            return 0

def average_valuation_affine(f,p):
    """
    returns the average p-valuation of the polynomial f. Works recursively.
    """
    v = valuation (f.content(), p)
    ZP= f.parent()
    x = ZP.gen()
    fv = ZP(f/p^v)
    Q = GF(p)['x'](fv.derivative())
    for r in fv.roots(GF(p)):
        if Q(r[0]) != 0:
            # we count an extra 1 when p divides (with probability 1/p)
            # + 1 when p^2 divides (with probability 1/p^2), and so on,
            # thus 1/p + 1/p^2 + ... = 1/(p-1)
            v += 1/(p-1)
        else:
            # we expand fv(r0+p*x) and divide by p since here we consider
            # only one of the p classes r0+p*x for 0 <= r0 < p
            f2 = fv(Integers()(r[0])+p*x)
            v += average_valuation_affine(f2, p)/p
    return v

# same as average_valuation_affine, but for the class a/b = r mod p
def average_valuation_affine_sub(f,p,r):
    v = valuation (f.content(), p)
    ZP= f.parent()
    x = ZP.gen()
    fv = ZP(f/p^v)
    Q = GF(p)['x'](fv.derivative())
    for t,_ in fv.roots(GF(p)):
        t = ZZ(t)
        if t != r:
           continue
        if Q(t) != 0: # single root
            # we count 1 for this root (with probability 1 since we are
            # in the residue class of this root), plus 1 when p^2 divides
            # (with probability 1/p), ..., thus 1 + 1/p + ... = p/(p-1)
            v += p/(p-1)
        else:
            f2 = fv(t+p*x)
            # since we are in this case in the general case (no restriction
            # any more on x), we call average_valuation_affine() [no _sub]
            v += average_valuation_affine(f2, p)
    return v

#def ava_deviation(f,p):
#    """
#    discrepancy wrt the simplistic estimate (for the simplistic
#    estimates, roots are counted only once)
#    """
#    v = valuation (f.content(), p)
#    ZP= f.parent()
#    x = ZP.gen()
#    fv = ZP(f/p^v)
#    Q = GF(p)['x'](fv.derivative())
#    for r in fv.roots(GF(p)):
#        if Q(r[0]) == 0:
#            f2 = fv(Integers()(r[0])+p*x)
#            v -= 1/(p-1)
#            v += average_valuation_affine(f2, p)/p
#    return v
#
def average_valuation_homogeneous_coprime(f,p):
    """
    returns the average p-valuation of f(a,b) for a,b coprime. Projective
    roots are counted as well.
    """
    if disc(f) % p > 0:
        # Then we know that the average valuation is
        # number_of_roots/(p-1), multiplied by p/(p+1) so as to take into
        # account the non-coprime pairs.
        return number_of_roots(f,p)/(p-1)*p/(p+1)
    # modulo p^n, roots touch one class having p^n*(1-1/p) representatives,
    # amongst the p^(2n)*(1-1/p^2) coprime representative pairs for
    # P^1(Z/p^n). So the contribution is, for p^n, p/(p+1) * p^-n
    r = average_valuation_affine(f, p) * p
    ZP= f.parent()
    x = ZP.gen()
    # Infinity expands as 0 when swapped. So we already know the
    # first-order term. For this reason, we we may drop one p, since it's
    # already counting one order ahead.
    r += average_valuation_affine((f.reverse())(p*x), p)
    r /= p+1
    return r

# same as average_valuation_homogeneous_coprime, but for a/b = r/s mod p
def average_valuation_homogeneous_coprime_sub(f,p,r,s):
    if disc(f) % p > 0:
        # Then we know that the average valuation is
        # number_of_roots/(p-1), multiplied by p/(p+1) so as to take into
        # account the non-coprime pairs.
        return number_of_roots_sub(f,p,r,s)/(p-1)*p
    # modulo p^n, roots touch one class having p^n*(1-1/p) representatives,
    # amongst the p^2n*(1-1/p^2) coprime representative pairs for
    # P^1(Z/p^n). So the contribution is, for p^n, p/(p+1) * p^-n
    if s <> 0:
        # we consider the class a/b = r/s mod p
        return average_valuation_affine_sub(f, p, (r/s) % p)
    else: # we consider the projective class
       ZP = f.parent()
       x = ZP.gen()
       return average_valuation_affine_sub((f.reverse())(p*x), p, 0)

def alpha_p(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f
    """
    return float((1/(p-1)-average_valuation_homogeneous_coprime(f,p))*log(p))

def alpha_p_nodisc(f,p):
    """
    Computes the contribution at p of the alpha value of f
    """
    return alpha_p(f,f.discriminant(),p)

#def alpha_p_affine_deviation(f,disc,p):
#    """
#    Computes the contribution at p of the alpha value of f, as a
#    deviation in comparison to the simplistic estimate.
#    """
#    if disc % p == 0:
#        s = float(-ava_deviation(f,p)*p/(p+1)*log(p))
#    else:
#        s = 0
#    return s
#
def alpha_p_affine(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f (affine part
    only)
    """
    if disc % p == 0:
        e = average_valuation_affine(f, p) * p / (p+1)
        s = float((1/(p-1)-e)*log(p))
    else:
        s = float((1-len(f.roots(GF(p)))*p/(p+1))*log(p)/(p-1))
    return s

def alpha_p_projective(f,disc,p):
    """
    Computes the contribution at p of the alpha value of f (projective part
    only). Always <= 0.
    """
    x=f.parent().gen()
    q=f.reverse()(p*x)
    if disc % p == 0:
        e = average_valuation_affine(q,p)/(p+1)
        s = float((-e)*log(p))
    else:
        fp= GF(p)['x'](f)
        assert fp != 0, "fp <> 0"
        assert(f.degree()-fp.degree() <= 1)
        s = float((-(f.degree()-fp.degree())*p/(p+1))*log(p)/(p-1))
    return s

#def alpha_p_projective_deviation(f,disc,p):
#    x=f.parent().gen()
#    q=f.reverse()(p*x)
#    if disc % p == 0:
#        s = float((-ava_deviation(f,p)/(p+1))*log(p))
#    else:
#        s = 0
#    return s
#

def alpha_p_affine_nodisc(f,p):
    return alpha_p_affine(f,f.discriminant(),p)

def alpha_p_projective_nodisc(f,p):
    return alpha_p_projective(f,f.discriminant(),p)

def alpha(f,B):
    """
    Computes the alpha value of f, up to prime bound B
    """
    x = f.variables()[0]
    R.<x> = PolynomialRing(ZZ)
    f = R(f)
    disc = f.discriminant()
    return sum([alpha_p(f, disc, p) for p in prime_range(2,B+1)])

def alpha_affine(f,B):
    """
    Computes the affine part of the alpha value of f, up to prime bound B
    """
    disc = f.discriminant()
    return sum([alpha_p_affine(f, disc, p) for p in prime_range(2,B+1)])

# If we consider 10000 random irreducible polynomials (of degree 5) with
# leading coefficient ad, the records for the average projective alpha are:
# ad alpha
# 1 0
# 2 -0.360690687857
# 3 -0.369037600417
# 4 -0.509931774087
# 6 -0.730588833099
# 12 -0.878720091232
# 24 -0.946481583101
# 30 -1.04940896209
# 60 -1.19986394549
# 120 -1.26667891471
# 180 -1.31231327574
# 210 -1.32294380312
# 360 -1.38818816912
# 420 -1.47641320724
# 840 -1.54315029826
# 1260 -1.5916463863
# 2520 -1.66107247254
# 4620 -1.69245754414
# 5040 -1.70643311481
# 9240 -1.76238361373
# 13860 -1.8092879226
# 27720 -1.87554651074
# 55440 -1.92230149269
# 110880 -1.9409975078
# 120120 -1.95548978794
# 138600 -1.9579224421
# 180180 -2.00300045817
# 360360 -2.07759306437
# 720720 -2.11776935996
# 1441440 -2.12997717254
# 1801800 -2.14915176062
# 2162160 -2.15017287971
# 3063060 -2.17025643763
# 3603600 -2.19357626159
# 6126120 -2.23852322631
# 12252240 -2.28536035419
# 24504480 -2.29913075626
# 30630600 -2.32088761142
# 58198140 -2.32116657149
# 61261200 -2.35980719547
# 116396280 -2.39334666315
# 183783600 -2.39692462272
# 232792560 -2.43637124145
# 465585120 -2.45262230569
def alpha_projective(f,B):
    """
    Computes the projective part of the alpha value of f, up to prime bound B.
    Always <= 0.
    """
    disc = f.discriminant()
    return sum([alpha_p_projective(f, disc, p) for p in prime_range(2,B+1)])

# estimate the average projective alpha of polynomials of degree d
# with leading coefficient ad divisible by 'incr', for incr <= ad < K*incr
def estimate_alpha_projective(d,incr,K,B):
    ad = incr
    n = s = 0
    R.<x> = ZZ[]
    for ad in range(incr,(K+1)*incr,incr):
       f = ad*x^d + R.random_element(degree=d-1)
       while not f.is_irreducible():
          f = ad*x^d + R.random_element(degree=d-1)
       s += alpha_projective (f, B)
       n += 1
       print ad, s/n

# -*-*- debug -*-*-
def estimate_average_valuation_homogeneous_coprime(f, p, x0, x1, y0, y1):
    """
    this computes the same thing as average_valuation_homogeneous_coprime(f,p),
    but does so experimentally.  Only for debugging.
    """
    s=0
    n=0
    x=f.parent().gen()
    for a in range(x0,x0+x1):
        for b in range(y0,y0+y1):
            if gcd(a,b) == 1:
                s+=valuation(f.resultant(a*x-b), p)
                n+=1
    return float(s/n) if n > 0 else Infinity

def estimate_alpha_p(f, p, nt):
    """
    Should compute the same thing as alpha_p, but experimentally
    """
    s=0
    n=0
    x=f.parent().gen()
    l=log(p)
    for i in range(nt):
        while True:
           a=randrange(0,nt^2)
           b=randrange(0,nt^2)
           c=randrange(0,nt^2)
           if gcd(a,b)==1 and c<>0:
              break
        s+=valuation(c,p)-valuation(f.resultant(a*x-b), p)
        n+=1
        # sys.stdout.write("%f\r" % float(s*l/n))
    return float(s*l/n) if n > 0 else Infinity

# same as estimate_alpha_p, but for (a,b) = (r,s) mod p
def estimate_alpha_p_2(f, p, nt, r, s):
    S=0
    n=0
    x = f.variables()[0]
    var('y')
    F = (f(x=x/y)*y^f.degree()).expand()
    l=log(p)
    for i in range(nt):
        while True:
           a=randrange(r,nt^2,p)
           b=randrange(s,nt^2,p)
           c=randrange(0,nt^2)
           if gcd(a,b)==1 and c<>0:
              break
        S+=valuation(c,p)-valuation(F(x=a,y=b), p)
        n+=1
    return float(S*l/n)

# same as above, but only counts the average valuation at p
def estimate_average_valuation_p_2(f, p, nt, r, s):
    S=0
    n=0
    x = f.variables()[0]
    var('y')
    F = (f(x=x/y)*y^f.degree()).expand()
    for i in range(nt):
        while True:
           a=randrange(r,nt^2,p)
           b=randrange(s,nt^2,p)
           if gcd(a,b)==1:
              break
        S+=valuation(F(x=a,y=b), p)
        n+=1
    return float(S/n)

def special_valuation_2(F, p, max_depth):
   c = F.content()
   e = c.valuation(p)
   F = F//c
   Rp.<xx,yy> = GF(p)[]
   Fp = Rp(F)
   if Fp.degree(xx)==0 and Fp.degree(yy)==0 and Fp<>0:
      # Fp is a non-zero constant polynomial: valuation is zero
      return e
   if max_depth==0:
      return e
   x, y = F.variables()
   for r in range(p):
      for s in range(p):
         e += special_valuation_2(F(x=p*x+r,y=p*y+s), p, max_depth-1)/p^2
   return e

def estimate_average_valuation_p_3(f, p, r, s, max_depth=2):
   x = f.variables()[0]
   y = var('y')
   F = (f(x=x/y)*y^f.degree()).expand()
   F = F(x=p*x+r,y=p*y+s).expand()
   R.<x,y> = ZZ[]
   F = R(F)
   return special_valuation_2(F, p, max_depth)

# auxiliary
def alpha_p_simplistic(f,p):
    """
    This ignores the multiple roots altogether.
    """
    return float((1-number_of_roots(f,p)*p/(p+1))*log(p)/(p-1))
# auxiliary
def alpha_simplistic(f,B):
    """
    This ignores the multiple roots altogether.
    """
    return sum([alpha_p_simplistic(f,p) for p in prime_range(2,B+1)])

# this function loops for f=144*x^4 + 576*x^3 + 864*x^2 + 576*x + 144
# = 144*(x+1)^4 and p=2, or more simply with f=(x+1)^2 and p=2
def special_val0 (f, p):
    c = f.content ()
    v = 0.0
    while c % p == 0:
        v += 1
        c = c // p
    if v <> 0:
        g = f // p^v
    else:
        g = f
    # g(x) = f(x)/p^v
    h = g(p*x)
    roots = g.roots(GF(p))
    nroots = len(roots)
    r0 = 0
    for i in range(nroots):
        r = roots[i][0]
        fp = f.derivative()
        c = fp(r)
        if c % p <> 0:
            v += 1.0 / (p - 1)
        else:
            h = h(x+ZZ(r)/p)
            h = PolynomialRing(ZZ,'x')(h)
            r0 = r
	    # here we have h(x) = f(p*x+r)/p^v
	    # we can have h(x) = f(x) only when v=d, and then the roots x of f
	    # are invariant under x -> p*x+r, which means that x = r/(1-p),
	    # thus f has one single root of multiplicity d, i.e.,
	    # f = lc(f) * (x + r/(p-1))^d
	    if h <> f: # avoid infinite loop
               v += special_val0 (h, p) / p
    return v

def special_valuation (f, p):
    disc = f.discriminant()
    pvaluation_disc = 0
    if disc % p == 0:
        pvaluation_disc += 1
        t = disc // p
        if t % p == 0:
            pvaluation_disc += 1
    if pvaluation_disc == 0:
        e = len(f.roots(GF(p)))
        if f.leading_coefficient() % p == 0:
            e += 1
        return float((p * e) / (p^2 - 1))
    elif pvaluation_disc == 1:
        e = len(f.roots(GF(p)))
        if f.leading_coefficient() % p == 0:
            e += 1
        return float((p * e - 1) / (p^2 - 1))
    else:
        d = f.degree()
        v = float(special_val0 (f, p) * p)
        if f.leading_coefficient() % p == 0:
            g = sum([f.coeffs()[d-i]*p^i*x^i for i in range (0,d+1)])
            v += float(special_val0 (g, p))
        return v / (p + 1)

def alpha_p_projective2(f,p):
    e = special_valuation (f, p)
    return float((-e)*log(p))

def check_alpha_projective(f,B):
    s = 0
    s2 = 0
    disc = f.discriminant()
    for p in prime_range(2,B+1):
        a = alpha_p_projective (f,disc,p)
        a2 = alpha_p_projective2 (f, p)
        s += a
        s2 += a2
        print p, a, a2, s, s2

# given a rootsieve space of n points, estimate the best alpha value
# see http://maths-people.anu.edu.au/~brent/pd/Bai-thesis.pdf 
# (page 40, formula 3.6) and
# A note on the first moment of extreme order statistics from the normal
# distribution, Max Petzold, https://gupea.ub.gu.se/handle/2077/3092
def expected_alpha(n):
    mu = 0
    sigma = 0.824
    return mu - sigma*(sqrt(2*log(n))-(log(log(n))+1.3766)/(2*sqrt(2*log(n))))

# print expected_alpha(2^k) for 0 <= k < bound (in 5 columns)
def expected_alpha_est(bound=150):
    c = 0
    print('%8.3f,' % 0),
    for k in range(1,bound):
        c = c + 1
        E = expected_alpha(2^k)
        if (c % 5 == 0):
            print
        print('%8.3f,' % E),

# compute some kind of "combined" alpha-value for two polynomials:
# for each prime p < B, and for each common root r mod p, add e*log(p)
# where e is the least multiplicity of the root for f and g
def combined_alpha(f,g,B):
   s = 0
   for p in prime_range(B):
      lf = f.roots(ring=GF(p))
      lg = g.roots(ring=GF(p))
      i = j = 0
      while i < len(lf) and j < len(lg):
         if lf[i][0] == lg[j][0]:
            s += min(lf[i][1],lg[j][1])*log(1.0*p)
            i += 1
            j += 1
         elif lf[i][0] < lg[j][0]:
            i += 1
         else:
            j += 1
   # roots at infinity
   lf = f.reverse().roots(ring=GF(p))
   lg = g.reverse().roots(ring=GF(p))
   i = j = 0
   while i < len(lf) and lf[i][0] != 0:
      i += 1
   while j < len(lg) and lg[j][0] != 0:
      j += 1
   if i < len(lf) and j < len(lg):
      assert lf[i][0] == 0
      assert lg[i][0] == 0
      s += min(lf[i][1],lg[i][1])*log(1.0*p)
   return s
