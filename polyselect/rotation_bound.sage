
# Here follows some code for estimating an intelligent rotation bound for
# a given polynomial pair, by making a compromise between:
# - possibly increasing the norm
# - getting better alpha values

# use as follows:
# attach rotation_bound.sage
# define f and g
# Try the ``Try this'' examples below.

import os,time

attach alpha.sage
attach integral.sage
attach E.sage

# optional
attach gnuplot_stuff.sage

# Experimentally, for degree 5, the *affine* contribution to alpha has
# approximately this shape. The projective contribution is disregarded
# because rotation leaves it unaffected.
# Note that the table has indices scaled by 10.
#reduced_alpha_affine_table=build_reduced_minimum_normal_table(0.571,0.851,0.05)
reduced_alpha_affine_table=build_reduced_minimum_normal_table(1.14,0.75,0.005)

# Note however that the affine contribution is not independent from the
# projective contribution. Hence for rotation starting from a given
# polynomial f, the answer is not at all clear.


# Indicates how far we hope to be from the expected value for the best
# alpha. 0 means exactly on the mean, hence 50% chances. +1 means 1
# standard deviation above, hence 84% chances. -1 means 16% chances.
deviation=0

##### Norm computations.

PRECISION=53

def deskew_polynomial(f,s):
    """Utility function. Returns a polynomial over the reals."""
    d=f.degree()
    g=[]
    ss=1
    for i in range(d):
        g.append(float(f[i])*ss)
        ss*=float(s)
    g.append(f[d]*ss)
    return RealField(PRECISION)['x'](g), ss

def supnorm(f,s):
    g,ss=deskew_polynomial(f,s)
    return max([abs(c) for c in g.coefficients()])/sqrt(ss)

def l1norm(f,s):
    g,ss=deskew_polynomial(f,s)
    return sum([abs(c) for c in g.coefficients()])/sqrt(ss)

def l2norm(f,s):
    g,ss=deskew_polynomial(f,s)
    return sqrt(sum([c^2 for c in g.coefficients()])/ss)

def square_evenpart(f):
    """Utility function. Computes the even part of the square of the
    polynomial f """
    d=f.degree()
    g=f.parent()(([sum([f[2*i-j]*f[j] for j in [0..2*i]]) for i in [0..d]]))
    return g

def l2norm_tk(f,s):
    """
    This norm gives the square-root of the integral of f^2 over the
    square [-1,1]^2, taking into account the given skewness.
    Note: We no longer take the half-sqrt here, which is what TK
    does (out of a mistake, apparently).
    """
    g,ss=deskew_polynomial(f,s)
    g2=square_evenpart(g)
    d=f.degree()
    coeffs=[4/(2*i+1)/(2*(d-i)+1) for i in [0..d]]
    return sqrt(vector(g2.coefficients())*vector(coeffs)/ss)

# see symbolic_l2 to generate formulas for larger degrees
def l2norm_tk_circular(f,s):
   if f.degree()==8: # symbolic_l2(8,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      a7 = f[7] * s^7
      a8 = f[8] * s^8
      n = 6435*a0^2 + 429*a1^2 + 858*a0*a2 + 99*a2^2 + 198*a1*a3 + 45*a3^2 + 198*a0*a4 + 90*a2*a4 + 35*a4^2 + 90*a1*a5 + 70*a3*a5 + 45*a5^2 + 90*a0*a6 + 70*a2*a6 + 90*a4*a6 + 99*a6^2 + 70*a1*a7 + 90*a3*a7 + 198*a5*a7 + 429*a7^2 + 70*a0*a8 + 90*a2*a8 + 198*a4*a8 + 858*a6*a8 + 6435*a8^2
      n = n * pi / 294912
      return RealField(PRECISION)(1/2 * log(n / s^8))
   elif f.degree()==7: # symbolic_l2(7,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      a7 = f[7] * s^7
      n = 429 * (a0^2 + a7^2) + 33 * (a1^2 + a6^2) + 66 * (a0*a2 + a5*a7) + 9 * (a2^2 + a5^2) + 18 * (a1*a3 + a0*a4 + a4*a6 + a3*a7) + 5 * (a3^2 + a4^2) + 10 * (a2*a4 + a1*a5 + a3*a5 + a0*a6 + a2*a6 + a1*a7)
      n = n * pi / 16384
      return RealField(PRECISION)(1/2 * log(n / s^7))
   elif f.degree()==6: # symbolic_l2(6,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      n = 231 * (a6 * a6 + a0 * a0) + 42 * (a6 * a4 + a2 * a0) + 21 * (a5 * a5 + a1 * a1) + 7 * (a4 * a4 + a2 * a2) + 14 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1) + 10 * (a6 * a0 + a5 * a1 + a4 * a2) + 5 * a3 * a3
      n = n * pi / 7168
      # return float(1/2 * log(n / (s * s * s * s * s * s)))
      return RealField(PRECISION)(1/2 * log(n / (s * s * s * s * s * s)))
   elif f.degree()==5: # symbolic_l2(5,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      n = 6 * (a3 * a1 + a1 * a5 + a4 * a2 + a0 * a4) + 14 * (a0 * a2 + a3 * a5) + 63.0 * (a0 * a0 + a5 * a5) + 7 * (a4 * a4 + a1 * a1) + 3 * (a3 * a3 + a2 * a2)
      n = n * pi / 1536
      return RealField(PRECISION)(1/2 * log(n / (s * s * s * s * s)))
   elif f.degree()==4: # symbolic_l2(4,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      n = 35 * (a4^2 + a0^2) + 10 * (a4*a2 + a2*a0) + 5 * (a3^2 + a1^2) + 6 * (a4*a0 + a3*a1) + 3 * a2^2
      n = n * pi / 640
      return RealField(PRECISION)(1/2 * log(n / (s * s * s * s)))
   elif f.degree()==3: # symbolic_l2(3,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      n = 5 * (a3^2 + a0^2) + 2 * (a3*a1 + a0*a2) + a1^2 + a2^2
      n = n * pi / 64
      return RealField(PRECISION)(1/2 * log(n / (s * s * s)))
   elif f.degree()==2: # symbolic_l2(2,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      n = 3 * (a2 * a2 + a0 * a0) + 2 * a0 * a2 + a1 * a1
      n = n * pi / 24
      return RealField(PRECISION)(1/2 * log(n / (s * s)))
   elif f.degree()==1: # symbolic_l2(1,0)[0](s).collect(pi).normalize()
      a0 = f[0]
      a1 = f[1] * s
      n = a0 * a0 + a1 * a1
      n = n * pi / 4
      return RealField(PRECISION)(1/2 * log(n / s))
   else:
      raise (ValueError, "circular norm not yet implemented for degree %d" % \
            f.degree())

##### Optimizing norms.

# For the sup norm, it's rather special: there's a fairly easy way to
# compute the best skewness from the graph log(norm) as a function of
# log(s)

def supnorm_hull_inner(f,logs):
    d=f.degree()
    i=0
    while f[i] == 0:
        i+=1
    l=[(-Infinity,i)]
    i+=1
    while i<=d:
        while True:
            t=l.pop()
            sx=float((logs[i]-logs[t[1]])/(t[1]-i))
            if sx >= t[0]:
                # This one was good, put it back in
                l.append(t)
                break
        l.append((sx,i))
        i+=1
    return l

def supnorm_hull(f):
    """
    Return the sup of the lines giving log(supnorm) as a function of
    log(s)
    """
    d=f.degree()
    logs=[float(log(abs(a))) if a != 0 else -Infinity for a in f.coefficients()]
    l=supnorm_hull_inner(f,logs)
    if l[0][0] == -Infinity:
        l.pop(0)
    # returns index, logskew, skew, lognorm
    return [ (v[1],v[0],exp(v[0]),logs[v[1]]+(v[1]-d/2)*v[0]) for v in l ]

def skew_supnorm(f):
    h=[(v[3],v[2]) for v in supnorm_hull(f)]
    m=min(h)
    return m[1]

def unique_positive_real_root(f,e):
    """Utility function. Amongst the positive real roots r of f (presumably
    there's only one), return the one giving the smallest value e(r)"""
    r=f.real_roots()
    root_pos=[s for s in r if s > 0]
    if len(root_pos)==1:
        return root_pos[0]
    else:
        ev=[(e(r),r) for r in root_pos]
        return min(ev)[1]

def skew_l1norm(f):
    d=f.degree()
    ZP=f.parent()
    # This is not the derivative, but rather 2s^(d/2+1) times the
    # derivative.
    g=ZP([abs(f[i])*(2*i-d) for i in [0..d]])
    return unique_positive_real_root(g,lambda s:l1norm(f,s))

def skew_l2norm(f):
    d=f.degree()
    ZP=f.parent()
    # This is not the derivative. It's rather, once evaluated at s^2,
    # s^(d+1) times the derivative of the square.
    g=ZP([f[i]^2*(2*i-d) for i in [0..d]])
    return sqrt(unique_positive_real_root(g,lambda s:l2norm(f,sqrt(s))))

def skew_l2norm_tk(f):
    """Gives the optimal skew for the L2 norm of a polynomial f"""
    d=f.degree()
    ZP=f.parent()
    x=ZP.gen()
    coeffs=[4*(2*i-d)*x^i/(2*i+1)/(2*(d-i)+1) for i in [0..d]]
    dd=vector(square_evenpart(f).coefficients())*vector(coeffs)
    return sqrt(unique_positive_real_root(dd,lambda s:l2norm_tk(f,sqrt(s))))

# Return the best chosen norm
def best_supnorm(f):
    return supnorm(f, skew_supnorm(f))

def best_l1norm(f):
    return l1norm(f, skew_l1norm(f))

def best_l2norm(f):
    return l2norm(f, skew_l2norm(f))

def best_l2norm_tk(f):
    return l2norm_tk(f, skew_l2norm_tk(f))

def best_l2norm_tk_circular(f):
    return l2norm_tk_circular(f, skew_l2norm_tk_circular(f))

def skew_l2norm_tk_circular(f, verbose=false):
   if f.degree()==8: # symbolic_l2(8,0)[0](s).collect(pi).normalize()
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      a7 = f[7] * s^7
      a8 = f[8] * s^8
      d = (6435*a0^2 + 429*a1^2 + 858*a0*a2 + 99*a2^2 + 198*a1*a3 + 45*a3^2 + 198*a0*a4 + 90*a2*a4 + 35*a4^2 + 90*a1*a5 + 70*a3*a5 + 45*a5^2 + 90*a0*a6 + 70*a2*a6 + 90*a4*a6 + 99*a6^2 + 70*a1*a7 + 90*a3*a7 + 198*a5*a7 + 429*a7^2 + 70*a0*a8 + 90*a2*a8 + 198*a4*a8 + 858*a6*a8 + 6435*a8^2)/s^8
      if verbose:
         print ("trying to minimize", d)
      # derivative of d wrt s (numerator, divided by 18)
      e = -2860*a0^2 - 143*a1^2 - 286*a0*a2 - 22*a2^2 - 44*a1*a3 - 5*a3^2 - 44*a0*a4 - 10*a2*a4 - 10*a1*a5 + 5*a5^2 - 10*a0*a6 + 10*a4*a6 + 22*a6^2 + 10*a3*a7 + 44*a5*a7 + 143*a7^2 + 10*a2*a8 + 44*a4*a8 + 286*a6*a8 + 2860*a8^2
      if verbose:
         print ("derivative is", e)
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      if verbose:
         print ("positive real roots of derivative:", r)
      best_norm = infinity
      for r in root_pos:
         no = d(r)
         if verbose:
            print ("value at r=", r, "is", no)
         if no < best_norm:
            best_norm = no
            best_root = r
      return best_root
   elif f.degree()==7: # symbolic_l2(7,0)[0](s).collect(pi).normalize()
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      a7 = f[7] * s^7
      d = (429*a0^2 + 33*a1^2 + 66*a0*a2 + 9*a2^2 + 18*a1*a3 + 5*a3^2 + 18*a0*a4 + 10*a2*a4 + 5*a4^2 + 10*a1*a5 + 10*a3*a5 + 9*a5^2 + 10*a0*a6 + 10*a2*a6 + 18*a4*a6 + 33*a6^2 + 10*a1*a7 + 18*a3*a7 + 66*a5*a7 + 429*a7^2)/s^7
      if verbose:
         print ("trying to minimize", d)
      # derivative of d wrt s (numerator)
      e = -3003*a0^2 - 165*a1^2 - 330*a0*a2 - 27*a2^2 - 54*a1*a3 - 5*a3^2 - 54*a0*a4 - 10*a2*a4 + 5*a4^2 - 10*a1*a5 + 10*a3*a5 + 27*a5^2 - 10*a0*a6 + 10*a2*a6 + 54*a4*a6 + 165*a6^2 + 10*a1*a7 + 54*a3*a7 + 330*a5*a7 + 3003*a7^2
      if verbose:
         print ("derivative is", e)
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      if verbose:
         print ("positive real roots of derivative:", r)
      best_norm = infinity
      for r in root_pos:
         no = d(r)
         if verbose:
            print ("value at r=", r, "is", no)
         if no < best_norm:
            best_norm = no
            best_root = r
      return best_root
   elif f.degree()==6:
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      a6 = f[6] * s^6
      d = (231*a0^2+42*a0*a2+14*a0*a4+10*a0*a6+21*a1^2+14*a1*a3+10*a1*a5+7*a2^2+10*a2*a4+14*a2*a6+5*a3^2+14*a3*a5+7*a4^2+42*a4*a6+21*a5^2+231*a6^2)/s^6
      if verbose:
         print ("trying to minimize", d)
      # derivative of d wrt s (numerator divided by 14)
      e = -99*a0^2+12*a6*a4+2*a6*a2-2*a4*a0-12*a2*a0+6*a5^2-6*a1^2+99*a6^2+a4^2-a2^2+2*a5*a3-2*a3*a1
      if verbose:
         print ("derivative is", e)
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      if verbose:
         print ("positive real roots of derivative:", r)
      best_norm = infinity
      for r in root_pos:
         no = d(r)
         if verbose:
            print ("value at r=", r, "is", no)
         if no < best_norm:
            best_norm = no
            best_root = r
      return best_root
   elif f.degree()==5:
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      a5 = f[5] * s^5
      e = 105*a5^2+7*(2*a3*a5+a4^2)+(2*a1*a5+2*a2*a4+a3^2)-(2*a0*a4+2*a1*a3+a2^2)-7*(2*a0*a2+a1^2)-105*a0^2
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      if len(root_pos) != 1:
         raise (ValueError, "number of positive roots != 1")
      return root_pos[0]
   elif f.degree()==4:
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      a4 = f[4] * s^4
      e = 14*a4^2 + (a3^2 + 2*a2*a4) - (a1^2 + 2*a0*a2) - 14*a0^2
      n = 35 * (a4^2 + a0^2) + 10 * (a4*a2 + a2*a0) + 5 * (a3^2 + a1^2) + 6 * (a4*a0 + a3*a1) + 3 * a2^2
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      i0 = 0
      v0 = n(root_pos[0])
      for i in range(1,len(root_pos)):
         v = n(root_pos[i])
         if v < v0:
           i0, v0 = i, v
      return root_pos[i0]
   elif f.degree()==3:
      R.<s> = RealField(PRECISION)[]
      a0 = f[0]
      a1 = f[1] * s
      a2 = f[2] * s^2
      a3 = f[3] * s^3
      e = 15*a3^2 + (a2^2 + 2*a1*a3) - (a1^2 + 2*a0*a2) - 15*a0^2
      r = e.real_roots()
      root_pos=[s for s in r if s > 0]
      if len(root_pos) != 1:
         raise (ValueError, "number of positive roots != 1")
      return root_pos[0]
   elif f.degree()==2:
      return RR(sqrt(abs(f[0]/f[2])))
   elif f.degree()==1:
      return RR(abs(f[0]/f[1]))
   else:
      raise (ValueError, "skew_l2norm_tk_circular not yet implemented for degree %d" % f.degree())

# return the skewness optimizing the Murphy-E value
def L2_combined_skewness (f, g, Bf, Bg, area, verbose=false):
   a = skew_l2norm_tk_circular (f)
   if verbose:
      print ("optimal skewness(f) = ", a)
   b = skew_l2norm_tk_circular (g)
   if verbose:
      print ("optimal skewness(g) = ", b)
   # for SNFS polynomials, b might be huge, which will take
   # many steps to converge, thus we search between a and 2^k*a
   a = min(a,b)
   b = 2*a
   # now a <= b
   va = MurphyE (f, g, a, Bf, Bg, area)
   vb = MurphyE (f, g, b, Bf, Bg, area)
   while vb > va:
      b = 2*b
      vb = MurphyE (f, g, b, Bf, Bg, area)
   while b - a > a * 0.001:
      if verbose:
         print ("a=", a, "va=", va, "b=", b, "vb=", vb)
      c = (2 * a + b) / 3.0
      vc = MurphyE (f, g, c, Bf, Bg, area)
      d = (a + 2 * b) / 3.0
      vd = MurphyE (f, g, d, Bf, Bg, area)
      max_left = max(va, vc)
      max_right = max(vd, vb)
      if max_left > max_right: # maximum in a or c
         b, vb = d, vd
      else:
         a, va = c, vc
   return (a + b) * 0.5

# return the skewness optimizing the lognorm sum
def L2_combined_skewness2 (f, g, verbose=false):
   a = skew_l2norm_tk_circular (f)
   b = skew_l2norm_tk_circular (g)
   if b < a:
      a, b = b, a
   # now a <= b
   va = l2norm_tk_circular (f, a) + l2norm_tk_circular (g, a)
   if verbose:
      print ("a=", a, "va=", va)
   vb = l2norm_tk_circular (f, b) + l2norm_tk_circular (g, b)
   if verbose:
      print ("b=", b, "vb=", vb)
   while b - a > a * 0.001:
      c = (2 * a + b) / 3.0
      vc = l2norm_tk_circular (f, c) + l2norm_tk_circular (g, c)
      if verbose:
      	 print ("c=", c, "vc=", vc)
      d = (a + 2 * b) / 3.0
      vd = l2norm_tk_circular (f, d) + l2norm_tk_circular (g, d)
      if verbose:
      	 print ("d=", d, "vd=", vd)
      max_left = max(va, vc)
      max_right = max(vd, vb)
      if max_left < max_right: # maximum in a or c
         b, vb = d, vd
      else:
         a, va = c, vc
      if verbose:
      	 print ("a=", a, "b=", b)
   return (a + b) * 0.5

def square_l2norm_tk_sym(f,s):
    d=f.degree()
    g=[]
    ss=1
    for i in range(d):
        g.append(f[i]*ss)
        ss*=float(s)
    g.append(f[d]*ss)
    g2=square_evenpart(f.parent()(g))
    d=f.degree()
    coeffs=[4/(2*i+1)/(2*(d-i)+1) for i in [0..d]]
    return 1/4*vector(g2.coefficients())*vector(coeffs)/ss


def bounds_l2norm_tk(f,g,multiplier):
    """Given a polynomial pair, computes eight points on the boundary of
    the area given by Norm(f+lambda*g)/Norm(f)<multiplier, where Norm is
    the L2 integral norm, and lambda is at most linear. Then return
    the largest enclosing rectangle."""
    s=skew_l2norm_tk(f)
    RP2.<x,u>=RR['x','u']
    RP1.<uu>=RR['u']
    r=RP1(square_l2norm_tk_sym((f(x)+u*g(x)).polynomial(x),s))
    r=r/r(0)-multiplier^2
    haxis=r.roots()
    r=RP1(square_l2norm_tk_sym((f(x)+u*x*g(x)).polynomial(x),s))
    r=r/r(0)-multiplier^2
    vaxis=r.roots()
    r=RP1(square_l2norm_tk_sym((f(x)+(u*(x+s))*g(x)).polynomial(x),s))
    r=r/r(0)-multiplier^2
    d1axis=r.roots()
    r=RP1(square_l2norm_tk_sym((f(x)+(u*(x-s))*g(x)).polynomial(x),s))
    r=r/r(0)-multiplier^2
    d2axis=r.roots()
    haxis=[(ra[0],0) for ra in haxis]
    vaxis=[(0,ra[0]) for ra in vaxis]
    d1axis=[(ra[0]*s,ra[0]) for ra in d1axis]
    d2axis=[(ra[0]*s,ra[0]) for ra in d2axis]
    # So what's the largest enclosing rectangle ?
    all=haxis+vaxis+d1axis+d2axis
    vmin=min([x[0] for x in all])
    vmax=max([x[0] for x in all])
    umin=min([x[1] for x in all])
    umax=max([x[1] for x in all])
    res=(vmin,vmax,umin,umax)
    volume=(umax-umin)*(vmax-vmin)
    # print "[%.2f..%.2f]x[%.2f..%.2f]"%res
    # print "vol=%.2f"%volume
    return res




### We can thus plot the norm with respect to the chosen skew:

def norms_from_s_for_plotting(f,g):
    def mine(c,t):
        return ((lambda s: flog(c(f, fexp10(s)))),t)
    l=[]
    l.append(mine(l1norm,"L1"))
    l.append(mine(l2norm,"L2"))
    l.append(mine(supnorm,"sup"))
    l.append(mine(l2norm_tk,"L2-int"))
    return l
# Try this:
#       myplot(norms_from_s_for_plotting(f,g), 3,6, 100)

# This shows (using for example the RSA768 polynomial) that the L2
# integral norm takes well into account the presence of real roots.



### We now turn to the sensibility of the norm(s) with respect to a0.

# utility stuff
def flog(x): return float(log(abs(float(x))))
def fexp(x): return float(exp(float(x)))
def flog10(x): return float(log(abs(float(x)))/log(10))
def fexp10(x): return float(exp(float(x*log(10))))

def change_ccoeff(f,a0):
    return PolynomialRing(RealField(PRECISION),'x')(f-f[0])+a0



# For the supnorm, the algorithm is slightly special
def supnorm_hull_extra(f):
    """
    Return the hull of the lines giving log(supnorm) as a function of
    log(s). Also gives extra information utilized in order to give the
    graph log(best supnorm) as a function of log(a0) and log(a1)
    """
    d=f.degree()
    c=1
    logs=[float(log(abs(a))) if a != 0 else -Infinity for a in f.coefficients()]
    f1=f-f[0]
    l=[(v[0],v[1],logs[v[1]]+(v[1]-d/2)*v[0])           \
                for v in supnorm_hull_inner(f1,logs)]
    while len(l) > 1 and l[0][2] > l[1][2]:
        l.pop(0)
    cuts=[]
    for v in l:
        # Which value of a0 matches this ?
        # We have to solve
        # log(a0)-d/2*log(s) == v[2]
        log_a0_cutoff = v[2]+d/2*v[0]
        log_a0_slope = (v[1]-d/2)/v[1]
        cuts.append((log_a0_cutoff,log_a0_slope,v[0],v[2]))
    return cuts;

def log_best_supnorm_from_log_a_quick(h,la):
    if la < h[0][0]:
        r=h[0][3]
    else:
        k=-1
        while la < h[k][0]:
            k-=1
        r=h[k][3]+h[k][1]*(la-h[k][0])
    return float(r)*1.0

def log_best_supnorm_from_log_a(f,la):
    h=supnorm_hull_extra(f)
    return log_best_supnorm_from_log_a_quick(h,la)

# Apart from the supnorm, other norms are accounted for in a standard
# manner.

def norms_from_a0_for_plotting(f,g):
    def mine(c,t,m):
        return ((lambda k: flog(c(change_ccoeff(f, m*fexp10(k))))),t)
    l=[]
    l.append(mine(best_l1norm,"L1",1))
    l.append(mine(best_l2norm,"L2",1))
    l.append(mine(best_supnorm,"sup",1))
    l.append(mine(best_l2norm_tk,"L2-int",1))
    l.append(mine(best_l2norm_tk,"L2-int",-1))
    return l
# Try this:
#       myplot(norms_from_a0_for_plotting(f,g), 38, 47, 100)
# Notice how the behaviour is a lot more intricate for the L2 integral norm.
# Our two plots give respecively f+k*g and f-k*g, which turn out to be
# somewhat different.

# From now on, we disregard norms other that the L2 integral norm.

def norms_from_a0_for_plotting_sconstant(f,g):
    def smine(c,s,t,m):
        return ((lambda k: flog(c(change_ccoeff(f, m*fexp10(k)), \
                s))),t)
    l=[]
    #l.append(smine(l1norm,skew_l1norm(f),"L1 (constant s)",1))
    #l.append(smine(l2norm,skew_l2norm(f),"L2 (constant s)",1))
    #l.append(smine(supnorm,skew_supnorm(f),"sup (constant s)",1))
    l.append(smine(l2norm_tk,skew_l2norm_tk(f),"L2-int (constant s)",1))
    l.append(smine(l2norm_tk,skew_l2norm_tk(f),"L2-int (constant s)",-1))
    return l

def sconstant_or_optimal(f,g):
    def mine(c,t,m):
        return ((lambda k: flog(c(change_ccoeff(f, m*fexp10(k))))),t)
    def smine(c,s,t,m):
        return ((lambda k: flog(c(change_ccoeff(f, m*fexp10(k)), \
                s))),t)
    l=[]
    l.append(mine(best_l2norm_tk,"L2-int",1))
    l.append(mine(best_l2norm_tk,"L2-int",-1))
    l.append(smine(l2norm_tk,skew_l2norm_tk(f),"L2-int (constant s)",1))
    l.append(smine(l2norm_tk,skew_l2norm_tk(f),"L2-int (constant s)",-1))
    return l
# Try this:
#       myplot(sconstant_or_optimal(f,g), 40, 45, 100)
# This indicates how a sensible choice for the skewness impacts the norm
# in the end.


# Now we focus on the rotation f+k*g. Given a bound 10^w on the size of
# the rotation area, try to have a guess on what is attainable in terms
# of log(norm)+bbest alpha.


# again, special case for the sup norm
def logsupnorm_plus_alpha_rot(f,g,w,h):
    """
    Gives the expected log(sup norm) + alpha value with rotation bounded
    by 10^w. The result takes into account an alpha value that is
    separated from the mean by as many standard deviations as given by
    the global ``deviation'' parameter.
    """
    lf0=flog(f[0])
    lm=flog(g[0])
    lognorm=log_best_supnorm_from_log_a_quick(h,max(lf0,w*flog(10)+lm))
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*w)
    return float(lognorm+alpha[0]+deviation*alpha[1])

def lognorm_plus_alpha_rot(f,g,normfunc,w):
    """
    Gives the expected log(chosen norm) + alpha value with rotation bounded
    by 10^w. The result takes into account an alpha value that is
    separated from the mean by as many standard deviations as given by
    the global ``deviation'' parameter.
    """
    RP=PolynomialRing(RealField(PRECISION),'z');
    E=[fexp10(w), -fexp10(w)]
    lognorm=min([flog(normfunc(RP(f)+e*RP(g))) for e in E])
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*(w))
    return float(lognorm+alpha[0]+deviation*alpha[1])

def lognorm_plus_alpha_rot_scons(f,g,normfunc,skew,w):
    """
    Does the same as lognorm_plus_alpha_rot, except that we assume that
    the skewness does not change.
    Gives the expected log(chosen norm) + alpha value with rotation bounded
    by 10^w. The result takes into account an alpha value that is
    separated from the mean by as many standard deviations as given by
    the global ``deviation'' parameter.
    """
    RP=PolynomialRing(RealField(PRECISION),'z');
    E=[fexp10(w), -fexp10(w)]
    lognorm=min([flog(normfunc(RP(f)+e*RP(g), skew)) for e in E])
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*(w))
    return float(lognorm+alpha[0]+deviation*alpha[1])

def lognorm_plus_alpha_rot_scons_linear(f,g,normfunc,skew,w):
    """
    Gives the expected log(chosen norm) + alpha value with degree-1
    rotation bounded by 10^w. The result takes into account an alpha
    value that is separated from the mean by as many standard deviations
    as given by the global ``deviation'' parameter.
    This code assumes that the skewness does not change.
    """
    s=flog(skew)
    logb=(w*flog(10)+s)/2
    loga=(w*flog(10)-s)/2
    RP=PolynomialRing(RealField(PRECISION),'z');
    z=RP.gen()
    Ea=fexp(loga); Eb=fexp(logb)
    E=[ Ea*z+Eb, Ea*z-Eb, -Ea*z+Eb, -Ea*z-Eb ];
    lognorm=min([flog(normfunc(RP(f)+e*RP(g), skew)) for e in E])
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*(w))
    return float(lognorm+alpha[0]+deviation*alpha[1])

def lognorm_plus_alpha_rot_linear(f,g,normfunc,skew,w):
    """
    Gives the expected log(chosen norm) + alpha value with degree-1
    rotation bounded by 10^w. The result takes into account an alpha
    value that is separated from the mean by as many standard deviations
    as given by the global ``deviation'' parameter.
    This code re-computes the optimal skewness from the corner polynomial
    pairs. However, the constant intitial skewness is used for the linear
    polynomial coefficients.
    """
    s=flog(skew)
    logb=(w*flog(10)+s)/2
    loga=(w*flog(10)-s)/2
    RP=PolynomialRing(RealField(PRECISION),'z');
    z=RP.gen()
    Ea=fexp(loga); Eb=fexp(logb)
    E=[ Ea*z+Eb, Ea*z-Eb, -Ea*z+Eb, -Ea*z-Eb ];
    lognorm=min([flog(normfunc(RP(f)+e*RP(g))) for e in E])
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*(w))
    return float(lognorm+alpha[0]+deviation*alpha[1])

def lognorm_plus_alpha_rot_linear_sopt(f,g,normfunc,w):
    """
    Gives the expected log(chosen norm) + alpha value with degree-1
    rotation bounded by 10^w. The result takes into account an alpha
    value that is separated from the mean by as many standard deviations
    as given by the global ``deviation'' parameter.
    This code re-computes the optimal skewness from the corner polynomial
    pairs. The ratio of the coefficients of the linear multiplier is
    chosen dynamically.
    This is D-O-G S-L-O-W !!!
    """
    cands=[]
    nsteps=20
    RP=PolynomialRing(RealField(PRECISION),'z');
    z=RP.gen()
    for s in [(k+1)/nsteps for k in range(nsteps)]:
        #+ [ flog(skew_l2norm_tk(f)) ]:
        # Try skewness A0/A1 = exp(1/s)
        # That's kind of silly, since we should probably try some
        # variation around the original skewness.
        loga=(w*flog(10)-1/s)/2; Ea=fexp(loga)
        logb=(w*flog(10)+1/s)/2; Eb=fexp(logb)
        E=[ Ea*z+Eb, Ea*z-Eb, -Ea*z+Eb, -Ea*z-Eb ];
        lognorm=min([flog(normfunc(RP(f)+e*RP(g))) for e in E])
        cands.append(lognorm)
    alpha=eval_dichotomy(reduced_alpha_affine_table,10.0*(w))
    return float(min(cands)+alpha[0]+deviation*alpha[1])


##########################################################################
# PLOTS

def rotatebound_dichotomy_sbest(f,g):
    """ this can be used to plot optimal-skew choices, for degree 0 rotation"""
    def mine(c,t):
        return (lambda v: lognorm_plus_alpha_rot(f,g,c,v), t)
    l=[]
    #l.append(mine(best_l2norm, "L2"))
    #l.append(mine(best_l1norm, "L1"))
    ## For the sup norm, we may use the faster algorithm.
    #h=supnorm_hull_extra(f)
    #l.append((lambda v: logsupnorm_plus_alpha_rot(f,g,v,h), "sup"))
    # l.append(mine(best_supnorm, "sup"))
    l.append(mine(best_l2norm_tk, "L2-int"))
    return l

def rotatebound_dichotomy_sconst(f,g):
    """ this can be used to plot constant-skew choices, for degree 0 rotation"""
    def mine(c,s,t):
        tx = t + " (constant skew)"
        return (lambda v: lognorm_plus_alpha_rot_scons(f,g,c,s,v), tx)
    l=[]
    #l.append(mine(l2norm, skew_l2norm(f), "L2"))
    #l.append(mine(l1norm, skew_l1norm(f), "L1"))
    #l.append(mine(supnorm, skew_supnorm(f), "sup"))
    l.append(mine(l2norm_tk, skew_l2norm_tk(f), "L2-int"))
    return l

def rotatebound_choices(f,g):
    def mine(c,t):
        return (lambda v: lognorm_plus_alpha_rot(f,g,c,v), t)
    def scmine(c,s,t):
        tx = t + " (constant skew)"
        return (lambda v: lognorm_plus_alpha_rot_scons(f,g,c,s,v), tx)
    def sclmine(c,s,t):
        tx = t + " (constant skew, degree 1)"
        return (lambda v: lognorm_plus_alpha_rot_scons_linear(f,g,c,s,v), tx)
    def lmine(c,s,t):
        tx = t + " (degree 1)"
        return (lambda v: lognorm_plus_alpha_rot_linear(f,g,c,s,v), tx)
    l=[]
    l.append(mine(best_l2norm_tk, "L2-int"))
    l.append(scmine(l2norm_tk, skew_l2norm_tk(f), "L2-int"))
    l.append(sclmine(l2norm_tk, skew_l2norm_tk(f), "L2-int"))
    l.append(lmine(best_l2norm_tk, skew_l2norm_tk(f), "L2-int"))
    return l

def rotatebound_choices_all(f,g):
    # Includes expensive computations.
    l=rotatebound_choices(f,g)
    tx="L2-int (degree 1, optimal skew)"
    c=best_l2norm_tk
    l.append((lambda v: lognorm_plus_alpha_rot_linear_sopt(f,g,c,v), tx))

    return l

# Try this:
#               myplot(rotatebound_choices(f,g), 5, 10, 100)
#               myplot(rotatebound_choices_all(f,g), 6, 9, 20)

# This shows that degree-1 rotation, as soon as some work margin is
# available, give a very clear win over constant rotation, as expected.

# optimize the norm of a given polynomial-pair (f,g) by applying successively
# translation, rotation by x*g, and rotation by g, until we find a local
# minimum
def optimize(f,g,maxsteps=200,eps=0.001):
    R = f.parent()
    x = f.parent().gen()
    logmu0 = best_l2norm_tk_circular(f)
    count = 0
    kt = kr0 = kr1 = kr2 = 1
    while best_l2norm_tk_circular(R(f+kr0*g)) <= logmu0 + eps:
       kr0 *= 2
    while best_l2norm_tk_circular(R(f+kr1*x*g)) <= logmu0 + eps:
       kr1 *= 2
    while best_l2norm_tk_circular(R(f+kr2*x^2*g)) <= logmu0 + eps:
       kr2 *= 2
    while best_l2norm_tk_circular(R(f(x=x+kt))) <= logmu0 + eps:
       kt *= 2
    while count < maxsteps:
        count += 1
        changedt = changedr2 = changedr1 = changedr0 = False

        # try translation by k
        f = f(x+kt)
        g = g(x+kt)
        logmu = best_l2norm_tk_circular(R(f))
        if logmu < logmu0:
            changedt = True
            logmu0 = logmu
        else:
            f = f(x-2*kt)
            g = g(x-2*kt)
            logmu = best_l2norm_tk_circular(R(f))
            if logmu < logmu0:
               changedt = True
               logmu0 = logmu
            else:
               f = f(x+kt)
               g = g(x+kt)

        # try rotation by x^2*g
        f = f + kr2*x^2*g
        logmu = best_l2norm_tk_circular(R(f))
        if logmu < logmu0:
            changedr2 = True
            logmu0 = logmu
        else:
            f = f - 2*kr2*x^2*g
            logmu = best_l2norm_tk_circular(R(f))
            if logmu < logmu0:
                changedr2 = True
                logmu0 = logmu
            else:
                f = f + kr2*x^2*g

        # try rotation by x*g
        f = f + kr1*x*g
        logmu = best_l2norm_tk_circular(R(f))
        if logmu < logmu0:
            changedr1 = True
            logmu0 = logmu
        else:
            f = f - 2*kr1*x*g
            logmu = best_l2norm_tk_circular(R(f))
            if logmu < logmu0:
                changedr1 = True
                logmu0 = logmu
            else:
                f = f + kr1*x*g

        # try rotation by g
        f = f + kr0*g
        logmu = best_l2norm_tk_circular(R(f))
        if logmu < logmu0:
            changedr0 = True
            logmu0 = logmu
        else:
            f = f - 2*kr0*g
            logmu = best_l2norm_tk_circular(R(f))
            if logmu < logmu0:
                changedr0 = True
                logmu0 = logmu
            else:
                f = f + kr0*g

        if changedt == False and changedr2 == False and changedr1 == False and changedr0 == False and kt == 1 and kr2 == 1 and kr1 == 1 and kr0 == 1:
            break
        if changedt == True:
            kt = 2*kt
        elif kt > 1:
            kt = kt//2
        if changedr2 == True:
            kr2 = 2*kr2
        elif kr2 > 1:
            kr2 = kr2//2
        if changedr1 == True:
            kr1 = 2*kr1
        elif kr1 > 1:
            kr1 = kr1//2
        if changedr0 == True:
            kr0 = 2*kr0
        elif kr0 > 1:
            kr0 = kr0//2
    return f, g


###################################################
# 1. symbolically generate l2norm
# usage sage: [f, fd, fd2] = symbolic_l2norm(6, 0)
#             6 is the degree, 0 means ellipse region (circular);
#                              1 is rectangular region.
###################################################
def symbolic_l2(d, method):
    if ((d < 1) | (d > 8)):
        return "Not Implemented"
    a8, a7, a6, a5, a4, a3, a2, a1, a0 = var('a8 a7 a6 a5 a4 a3 a2 a1 a0')
    x, y, s, r, t = var('x, y, s, r, t')
    if d <= 1:
        a2 = 0
    if d <= 2:
        a3 = 0
    if d <= 3:
        a4 = 0
    if d <= 4:
        a5 = 0
    if d <= 5:
        a6 = 0
    if d <= 6:
        a7 = 0
    if d <= 7:
        a8 = 0
    f = a8*x^8+a7*x^7+a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
    g(x, y) = f(x=x/y)*y^d

    # 0 means ellipse, others means rectangular
    if method == 0:
        g(x=r*cos(t), y=r*sin(t))
        h(r,t,s)= g(x=s^(1/2)*r*cos(t), y=r/s^(1/2)*sin(t))
        l(t,s) = integrate(h^2*r, r, 0, 1)
        F(s) = integrate(l, t, 0, 2*pi)
    else:
        l(x,s) = integrate(g^2/s^d, y, 0, 1)
        F(s) = integrate(l, x, 0, 1)

    return [F, F.derivative(s), F.derivative(s).derivative(s)]

# Usage: define the polynomial f and g, then use
# symbolic_l2_f(f, method) to return  [normf, normfd, normfd2]
def symbolic_l2_f(f, method):
    d = f.degree()
    if ((d < 3) | (d > 6)):
        return "Not Implemented"
    a6, a5, a4, a3, a2, a1, a0 = var('a6 a5 a4 a3 a2 a1 a0')
    x, y, s, r, t = var('x, y, s, r, t')
    if (d == 3):
        a6 = 0
        a5 = 0
        a4 = 0
    if (d == 4):
        a6 = 0
        a5 = 0
        a4 = f[4]
    if (d == 5):
        a6 = 0
        a5 = f[5]
        a4 = f[4]
    if (d == 6):
        a6 = f[6]
        a5 = f[5]
        a4 = f[4]
    a3 = f[3]
    a2 = f[2]
    a1 = f[1]
    a0 = f[0]
    f = a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
    g(x, y) = f(x=x/y)*y^d
    # 0 means ellilpse, others means rectangular
    if (method == 0):
        g(x=r*cos(t), y=r*sin(t))
        h(r,t,s)= g(x=s^(1/2)*r*cos(t), y=r/s^(1/2)*sin(t))
        l(t,s) = integrate(h^2*r, r, 0, 1)
        F(s) = integrate(l, t, 0, 2*pi)
    else:
        l(x,s) = integrate(g^2/s^d, y, 0, 1)
        F(s) = integrate(l, x, 0, 1)

    return [F, F.derivative(s), F.derivative(s).derivative(s)]

###################################
# 2. symbolically do translation
###################################
def symbolic_translate(d):
    if ((d < 3) | (d > 6)):
        return "Not Implemented"
    a6, a5, a4, a3, a2, a1, a0 = var('a6 a5 a4 a3 a2 a1 a0')
    x, y, l, m, k = var('x, y, l, m, k')
    if (d == 3):
        a6 = 0
        a5 = 0
        a4 = 0
    if (d == 4):
        a6 = 0
        a5 = 0
    if (d == 5):
        a6 = 0
    f = a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0
    g = l*x - m
    f = f(x = x + k)
    g = g + k*l
    return [f, g]
