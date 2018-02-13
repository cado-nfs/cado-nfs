# number of special-q in [q0,q1] with factors in [qfac_min,qfac_max]
def nb_special_q(q0,q1,qfac_min=None,qfac_max=infinity,verbose=false):
   if qfac_min == None: # we consider only primes
      return prime_pi(q1)-prime_pi(q0)
   else:
      # determine maximum number of prime factors
      kmax = floor(log(q1)/log(qfac_min))
      assert kmax <= 3, "number of special-q factors is assumed <= 3"
      pi1 = prime_pi(min(q1,qfac_max))-prime_pi(max(q0,qfac_min))
      # count number of composites p*q with q0 <= p*q < q1 and p <= q
      # we thus have p^2 < q1, thus p < sqrt(q1)
      # for each p, we have q0/p <= q < q1/p
      pi2 = 0
      pmax = min(qfac_max,floor(sqrt(q1-1)))
      for p in prime_range(qfac_min,pmax+1):
         qmin = max(ceil(q0/p),p)
         qmax = floor((q1-1)/p)
         assert qmin <= qmax
         pi2 += prime_pi(qmax) - prime_pi(qmin-1)
      # q0 <= p*q*r < q1 with p <= q <= r
      pi3 = 0
      pmax = min(qfac_max,floor((q1-1)^(1/3)))
      for p in prime_range(qfac_min,pmax+1):
         # q0/p <= q*r < q1/p
         qmin = p
         qmax = min(qfac_max,floor(sqrt((q1-1)//p)))
         for q in prime_range(p,qmax+1):
            rmin = max(ceil(q0/(p*q)),q)
            rmax = floor((q1-1)/(p*q))
            assert rmin <= rmax
            pi3 += prime_pi(rmax) - prime_pi(rmin-1)
      if verbose:
         print pi1, pi2, pi3
      return pi1 + pi2 + pi3

# exact count for k factors (slow)
def nb_special_q_k(q0,q1,k,qfac_min=None,qfac_max=infinity):
   q = q0
   count = 0
   while q < q1:
      l = factor(q)
      nb_factors = sum(l[i][1] for i in range(len(l)))
      if nb_factors == k:
         qmin = min(l[i][0] for i in range(len(l)))
         qmax = max(l[i][0] for i in range(len(l)))
         if qfac_min<=qmin and qmax<=qfac_max:
            count += 1
      q = q + 1
   return count
