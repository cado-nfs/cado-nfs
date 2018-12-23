# cf https://hal.inria.fr/hal-01273045v2

# estimate average p-valuation of f
# by trying all degree-2 polynomials with coefficients in [0..k-1]
# f is the polynomial
# p is the prime
# Example:
# R.<x> = ZZ[]
# f = x^9+9*x^8+30*x^7+62*x^6+87*x^5+87*x^4+62*x^3+30*x^2+9*x+3
# estimate_valuation_p(f,3,81) gives 6.968/13
def estimate_valuation_p(f,p,k):
   assert f.is_irreducible()
   s = n = 0
   R = f.parent()
   x = f.variables()[0]
   for a2 in range(1,k): # to ensure degree 2, we start at a2=1
      for a1 in range(k):
         for a0 in range(k):
            a = a2*x^2+a1*x+a0
            if a.content()==1 and a.is_irreducible():
               s += f.resultant(a).valuation(p)
               n += 1
   return s/n
   
   
