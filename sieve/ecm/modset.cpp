#include "cado.h" // IWYU pragma: keep

#include <limits.h>        // for LONG_BIT
#include <stdlib.h>        // for abort, size_t

#include "modset.hpp"

#include "macros.h"        // for ASSERT_ALWAYS

#include "facul.hpp"       // for facul_strategy_t
#include "facul_doit.hpp"  // for facul_doit, facul_doit_onefm
struct cxx_mpz;


/* This should be overloaded, but all the mod_int types are passed as
 * pointers to unsigned long, so they are indistinguishable in the
 * function signature. Thus we cannot use a function template, either,
 * as all the instances of the template would have the same signature. */

const FaculModulusBase *
FaculModulusBase::init_ul (const modintredcul_t n)
{
  return new FaculModulusUl(n);
}

const FaculModulusBase *
FaculModulusBase::init_15ul (const modintredc15ul_t n)
{
  const size_t bits = modredc15ul_intbits (n);
  unsigned long t1[2];
  size_t nr_words = modredc15ul_intget_uls(t1, n);
  ASSERT_ALWAYS(nr_words <= 2);

  if (bits <= MODREDCUL_MAXBITS) {
    modintredcul_t i;
    modredcul_intinit(i);
    modredcul_intset_uls(i, t1, nr_words);
    FaculModulusUl *m = new FaculModulusUl(i);
    modredcul_intclear(i);
    return m;
  }
  else if (bits <= MODREDC15UL_MAXBITS) {
      return new FaculModulus15Ul (n);
  } else
    abort();
}

const FaculModulusBase *
FaculModulusBase::init_2ul2 (const modintredc2ul2_t n)
{
  const size_t bits = modredc2ul2_intbits (n);
  unsigned long t1[2];
  size_t nr_words = modredc2ul2_intget_uls(t1, n);
  ASSERT_ALWAYS(nr_words <= 2);

  if (bits <= MODREDCUL_MAXBITS) {
    modintredcul_t i;
    modredcul_intinit(i);
    modredcul_intset_uls(i, t1, nr_words);
    FaculModulusUl *m = new FaculModulusUl(i);
    modredcul_intclear(i);
    return m;
  } else if (bits <= MODREDC15UL_MAXBITS) {
      modintredc15ul_t t2;
      modredc15ul_intinit (t2);
      modredc15ul_intset_uls (t2, t1, nr_words);
      FaculModulus15Ul *m = new FaculModulus15Ul (t2);
      modredc15ul_intclear (t2);
      return m;
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      return new FaculModulus2Ul2(n);
    }
  else
      abort();
}

const FaculModulusBase *
FaculModulusBase::init_mpz (const modintmpz_t n)
{
  const size_t bits = modmpz_intbits (n);
  FaculModulusBase *m;
  if (bits <= MODREDCUL_MAXBITS) {
      unsigned long t1[1];
      size_t nr_words = modmpz_intget_uls(t1, n);
      ASSERT_ALWAYS(nr_words <= 1);
      modintredcul_t i;
      modredcul_intinit(i);
      modredcul_intset_uls(i, t1, nr_words);
      m = new FaculModulusUl(i);
      modredcul_intclear(i);
      return m;
  } else if (bits <= MODREDC15UL_MAXBITS) {
      unsigned long t1[2];
      modintredc15ul_t t2;
      size_t nr_words = modmpz_intget_uls(t1, n);
      ASSERT_ALWAYS(nr_words <= 2);
      modredc15ul_intinit (t2);
      modredc15ul_intset_uls (t2, t1, nr_words);
      m = new FaculModulus15Ul (t2);
      modredc15ul_intclear (t2);
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      unsigned long t1[2];
      modintredc2ul2_t t2;
      size_t nr_words = modmpz_intget_uls(t1, n);
      ASSERT_ALWAYS(nr_words <= 2);
      modredc2ul2_intinit (t2);
      modredc2ul2_intset_uls (t2, t1, nr_words);
      m = new FaculModulus2Ul2 (t2);
      modredc2ul2_intclear (t2);
    }
  else if (bits <= MODMPZ_MAXBITS)
    {
      /* We assume for now that m is a modintmpz_t */
      m = new FaculModulusMpz(n);
    }
  else
      abort();
  return m;
}

int FaculModulusUl::facul_doit(std::vector<cxx_mpz> & factors,
    facul_strategy_oneside const & strategy, const int method_start) const
{
    return ::facul_doit (factors, m, strategy, method_start);
}

int FaculModulusUl::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method, const FaculModulusBase * &fm,
    const FaculModulusBase * &cfm, unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, fm, cfm, lpb,
        assume_prime_thresh, BBB);
}

void
FaculModulusUl::get_z (mpz_t z) const
{
    mpz_set_ui (z, m->m);
}

int FaculModulus15Ul::facul_doit(std::vector<cxx_mpz> & factors, 
    facul_strategy_oneside const & strategy, const int method_start) const
{
    return ::facul_doit (factors, m, strategy, method_start);
}

int FaculModulus15Ul::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method, const FaculModulusBase * &fm,
    const FaculModulusBase * &cfm, unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, fm, cfm, lpb,
        assume_prime_thresh, BBB);
}

void
FaculModulus15Ul::get_z (mpz_t z) const
{
    mpz_set_ui (z, m->m[1]);
    mpz_mul_2exp (z, z, LONG_BIT);
    mpz_add_ui (z, z, m->m[0]);
}

int FaculModulus2Ul2::facul_doit(std::vector<cxx_mpz> & factors, 
    facul_strategy_oneside const & strategy, const int method_start) const
{
    return ::facul_doit (factors, m, strategy, method_start);
}

int FaculModulus2Ul2::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method, const FaculModulusBase * &fm,
    const FaculModulusBase * &cfm, unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, fm, cfm, lpb,
        assume_prime_thresh, BBB);
}

void
FaculModulus2Ul2::get_z (mpz_t z) const
{
    mpz_set_ui (z, m->m[1]);
    mpz_mul_2exp (z, z, LONG_BIT);
    mpz_add_ui (z, z, m->m[0]);
}

int FaculModulusMpz::facul_doit(std::vector<cxx_mpz> & factors, 
    facul_strategy_oneside const & strategy, const int method_start) const
{
    return ::facul_doit (factors, m, strategy, method_start);
}

int FaculModulusMpz::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method, const FaculModulusBase * &fm,
    const FaculModulusBase * &cfm, unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, fm, cfm, lpb,
        assume_prime_thresh, BBB);
}

void
FaculModulusMpz::get_z (mpz_t z) const
{
    mpz_set (z, m);
}
