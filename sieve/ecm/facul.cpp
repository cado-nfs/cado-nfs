#include "cado.h" // IWYU pragma: keep

#include <algorithm>    // for sort
#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "modset.hpp"   // for FaculModulusBase

struct cxx_mpz_cmp {
    inline bool operator()(cxx_mpz const& a, cxx_mpz const & b) {
        return mpz_cmp(a, b) < 0;
    }
};

int
facul (std::vector<cxx_mpz> & factors, cxx_mpz const & N, facul_strategy_oneside const & strategy)
{
    int found = 0;

    /* XXX ATTENTION: This function may be called recursively. In
     * particular it may happen that the factors[] vector is not empty. */

    size_t factors_previous_size = factors.size();

#ifdef PARI
    gmp_fprintf (stderr, "%Zd", N);
#endif

    if (mpz_sgn (N) <= 0)
        return -1;
    if (mpz_cmp_ui (N, 1UL) == 0)
        return 0;

    /* Use the fastest modular arithmetic that's large enough for this input */
    const FaculModulusBase *m = FaculModulusBase::init_mpz(N);
    /* If the composite does not fit into our modular arithmetic, return
       no factor */
    if (m == NULL)
        return 0;

    found = m->facul_doit(factors, strategy, 0);

    delete m;
    m = NULL;

    if (found > 1)
    {
        /* Sort the factors we found */
        std::sort(factors.begin() + factors_previous_size, factors.end(), cxx_mpz_cmp());
    }

    return found;
}

/*
 * This is our auxiliary factorization.
 * It applies a bunch of ECM curves with larger bounds to find
 * a factor with high probability. It returns -1 if the factor
 * is not smooth, otherwise the number of
 * factors.
 * 
 * It tries the factoring methods in "methods", starting at the
 * "method_start"-th one. If any factor is found, then we try to factor any
 * composite factor or cofactor by recursively calling facul_aux.
 */
static int
facul_aux (std::vector<cxx_mpz> & factors, const FaculModulusBase *m,
	   facul_strategies const & strategies,
           std::vector<facul_method_side> const & methods,
	   size_t method_start, int side)
{
  int found = 0;

  if (method_start >= methods.size())
    return found;

  for (size_t i = method_start ; i < methods.size() ; ++i) {
      auto const & meth(methods[i]);
      if (meth.side != side)
	continue; /* this method is not for this side */

#ifdef ENABLE_UNSAFE_FACUL_STATS
      if (i < STATS_LEN)
	stats_called_aux[i]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
      const FaculModulusBase *fm = NULL, *cfm = NULL;

      int res_fac = m->facul_doit_onefm(factors,
            *meth.method, fm, cfm,
            strategies.lpb[side],
            strategies.BB[side],
            strategies.BBB[side]);
      // check our result
      // res_fac contains the number of factors found
      if (res_fac == FACUL_NOT_SMOOTH)
	{
	  /*
	    The cofactor m is not smooth. So, one stops the
	    cofactorization.
	  */
	  found = FACUL_NOT_SMOOTH;
	  break;
	}
      if (res_fac == 0)
	{
	  /* Zero factors found. If it was the last method for this
	     side, then one stops the cofactorization. Otherwise, one
	     tries with an other method */
	    continue;
	}

      found += res_fac;
      if (res_fac == 2)
	break;

      /*
	res_fac == 1  Only one factor has been found. Hence, our
	factorization is not finished.
      */
      if (fm != NULL)
	{
	  int found2 = facul_aux (factors, fm, strategies,
				  methods, i+1, side);
          if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
	    {
	      found = FACUL_NOT_SMOOTH;
	      delete cfm;
              cfm = NULL;
	      delete fm;
              fm = NULL;
	      break;
	    }
	  else
	    found += found2;
	  delete fm;
          fm = NULL;
	}
      if (cfm != NULL)
	{
	  int found2 = facul_aux (factors, cfm, strategies,
				  methods, i+1, side);
          if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
          {
              found = FACUL_NOT_SMOOTH;
          }
	  else
	    found += found2;
	  delete cfm;
          cfm = NULL;
	  break;
	}
      break;
    }
  return found;
}





/*
  This function tries to factor a pair of cofactors (m[0], m[1]) from
  strategies. It returns the number of factors found on each side, or
  -1 if the factor is not smooth.
  Remarks: - the values of factors found are stored in 'factors'.
           - the variable 'is_smooth' allows to know if a cofactor
             is already factored.

 XXX this is a mess. Cleanup needed.
*/

static std::vector<int>
facul_both_src (std::vector<std::vector<cxx_mpz>> & factors, const FaculModulusBase** m,
		facul_strategies const & strategies, int* cof,
		int* is_smooth)
{
    int nsides = factors.size();
    std::vector<int> found(nsides, 0);

    std::vector<facul_method_side> const & methods = strategies(cof[0],cof[1]);

    if (methods.empty())
        return found;

    const FaculModulusBase *f[2][2] = {{NULL, NULL}, {NULL, NULL}};
#ifdef ENABLE_UNSAFE_FACUL_STATS
    int stats_nb_side = 0, stats_index_transition = 0;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
    int last_i[2]; /* last_i[s] is the index of last method tried on side s */
    for (int i = 0; i < (int) methods.size() ; i++)
    {
#ifdef ENABLE_UNSAFE_FACUL_STATS
        stats_current_index = i - stats_nb_side * stats_index_transition;
        if (methods[i].is_last)
        {
            stats_nb_side = 1;
            stats_index_transition = i+1;
        }
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
        int side = methods[i].side;
        if (is_smooth[side] != FACUL_MAYBE)
        {
            /* If both sides are smooth, we can exit the loop,
               otherwise we must continue with the next methods,
               since methods might be interleaved between side 0 and 1,
               thus we don't have an easy way to skip all methods for this side.
               We could do this with another representation, say methods[0][i]
               for side 0, 0 <= i < m, methods[1][j] for side 1, 0 <= j < n,
               and which_method[k] = {0, 1} for 0 <= k < m+n. */
            if (is_smooth[side] == FACUL_SMOOTH &&
                    is_smooth[1-side] == FACUL_SMOOTH)
                break;
            continue;
        }

#ifdef ENABLE_UNSAFE_FACUL_STATS
        if (stats_current_index < STATS_LEN)
            stats_called[stats_current_index]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
        int res_fac = 0;
        last_i[side] = i;
        res_fac = m[side]->facul_doit_onefm (factors[side], *methods[i].method,
                f[side][0], f[side][1], strategies.lpb[side],
                strategies.BB[side], strategies.BBB[side]);
        // check our result
        // res_fac contains the number of factors found, or -1 if not smooth
        if (res_fac == -1)
        {
            /*
               The cofactor m[side] is not smooth. So, one stops the
               cofactorization.
               */
            found[side] = -1;
            break;
        }
        if (res_fac == 0)
        {
            /* No factor found. If it was the last method for this
               side, then one stops the cofactorization. Otherwise, one
               tries with an other method.
               */
            if (methods[i].is_last)
                break;
            else
                continue;
        }
        found[side] = res_fac;

        if (res_fac == 2)
        {
            /*
               Indeed, if using only one factoring method we found two
               prime factors of m (f, m/f) then m is factored and work is
               finished for this cofactor.
               */
            is_smooth[side] = FACUL_SMOOTH;
            continue;
        }
        /*
           res_fac == 1. Only one factor has been found. Hence, an
           auxiliary factorization will be necessary.
           */
        is_smooth[side] = FACUL_AUX;
    }
    // begin the auxiliary factorization
    if (is_smooth[0] >= 1 && is_smooth[1] >= 1) {
        for (int side = 0; side < 2; side++) {
            if (is_smooth[side] == FACUL_AUX) {
                for (int ind_cof = 0; ind_cof < 2; ind_cof++) {
                    // factor f[side][0] or/and f[side][1]
                    if (f[side][ind_cof] != NULL)
                    {
                        // **IF** we reach here, then some is_smooth[side] was
                        // set to FACUL_AUX somehow, and this can only happen if
                        // we passed through last_i = side in the loop above
                        // (because we *never* set to FACUL_AUX elsewhere).
                        //
                        // XXX honestly, this can be understood as a sign that
                        // this code deserves some long overdue cleanup.
                        //
                        // coverity[uninit_use]
                        int found2 = facul_aux (factors[side],
                                f[side][ind_cof], strategies,
                                methods, last_i[side] + 1, side);
                        if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
                        {
                            is_smooth[side] = FACUL_NOT_SMOOTH;
                            found[side] = found2;// FACUL_NOT_SMOOTH or FACUL_MAYBE
                            goto clean_up;
                        } else {
                            is_smooth[side] = FACUL_SMOOTH;
                            found[side] += found2;
                        }
                    }
                }
            }
        }
    }

clean_up:
    delete f[0][0];
    delete f[0][1];
    delete f[1][0];
    delete f[1][1];
    return found;
}


/*
  This function is like facul, but we will work with both norms
  together.  It returns the number of factors for each side.
*/
std::vector<int>
facul_both (std::vector<std::vector<cxx_mpz>> & factors,
            std::vector<cxx_mpz> & N,
	    facul_strategies const & strategies, int* is_smooth)
{
  int nsides = factors.size();
  ASSERT_ALWAYS(factors.size() == (size_t) nsides);
  ASSERT_ALWAYS(N.size() ==(size_t)  nsides);
  int cof[2];
  size_t bits;
  std::vector<int> found(nsides, 0);

  const FaculModulusBase *n[2];

#ifdef PARI
  gmp_fprintf (stderr, "(%Zd %Zd)", N[0], N[1]);
#endif

  /* cofactors should be positive */
  ASSERT (mpz_sgn (N[0]) > 0 && mpz_sgn (N[1]) > 0);

  if (mpz_cmp_ui (N[0], 1UL) == 0)
    is_smooth[0] = FACUL_SMOOTH;
  if (mpz_cmp_ui (N[1], 1UL) == 0)
    is_smooth[1] = FACUL_SMOOTH;

  for (int side = 0; side < 2; side++)
    {
      /* If the composite does not fit into our modular arithmetic, return
	 no factor */
      bits = mpz_sizeinbase (N[side], 2);
      cof[side] = bits;
      if (bits > MODMPZ_MAXBITS)
	return found;

      /* Use the fastest modular arithmetic that's large enough for
	 this input */
      n[side] = FaculModulusBase::init_mpz(N[side]);
    }

  found = facul_both_src (factors, n, strategies, cof, is_smooth);
  for (int side = 0; side < 2; side++)
    {
      if (found[side] > 1)
	{
	  /* Sort the factors we found */
          ASSERT_ALWAYS(factors[side].size() == (size_t) found[side]);
          std::sort(factors[side].begin(), factors[side].end(), cxx_mpz_cmp());
	}
    }

  // Free
  delete n[0];
  delete n[1];

  return found;
}
