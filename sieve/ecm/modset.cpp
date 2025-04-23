#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdlib>

#include <vector>
#include <memory>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "modset.hpp"
#include "macros.h"
#include "facul.hpp"
#include "facul_doit.hpp"
#include "facul_method.hpp"
#include "arith/mod_mpz.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/modredc_ul.h"

struct cxx_mpz;


/* This should be overloaded, but all the mod_int types are passed as
 * pointers to unsigned long, so they are indistinguishable in the
 * function signature. Thus we cannot use a function template, either,
 * as all the instances of the template would have the same signature. */

FaculModulusBase *
FaculModulusBase::init_ul (const modintredcul_t n)
{
  return new FaculModulusUl(n);
}

FaculModulusBase *
FaculModulusBase::init_15ul (const modintredc15ul_t n)
{
  const size_t bits = modredc15ul_intbits (n);
  unsigned long t1[2];
  size_t const nr_words = modredc15ul_intget_uls(t1, n);
  ASSERT_ALWAYS(nr_words <= 2);

  if (bits <= MODREDCUL_MAXBITS) {
    modintredcul_t i;
    modredcul_intinit(i);
    modredcul_intset_uls(i, t1, nr_words);
    auto *m = new FaculModulusUl(i);
    modredcul_intclear(i);
    return m;
  }
  else if (bits <= MODREDC15UL_MAXBITS) {
      return new FaculModulus15Ul (n);
  } else
    abort();
}

FaculModulusBase *
FaculModulusBase::init_2ul2 (const modintredc2ul2_t n)
{
  const size_t bits = modredc2ul2_intbits (n);
  unsigned long t1[2];
  size_t const nr_words = modredc2ul2_intget_uls(t1, n);
  ASSERT_ALWAYS(nr_words <= 2);

  if (bits <= MODREDCUL_MAXBITS) {
    modintredcul_t i;
    modredcul_intinit(i);
    modredcul_intset_uls(i, t1, nr_words);
    auto *m = new FaculModulusUl(i);
    modredcul_intclear(i);
    return m;
  } else if (bits <= MODREDC15UL_MAXBITS) {
      modintredc15ul_t t2;
      modredc15ul_intinit (t2);
      modredc15ul_intset_uls (t2, t1, nr_words);
      auto *m = new FaculModulus15Ul (t2);
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

FaculModulusBase *
FaculModulusBase::init_mpz (const modintmpz_t n)
{
  const size_t bits = modmpz_intbits (n);
  FaculModulusBase *m;
  if (bits <= MODREDCUL_MAXBITS) {
      unsigned long t1[1];
      size_t const nr_words = modmpz_intget_uls(t1, n);
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
      size_t const nr_words = modmpz_intget_uls(t1, n);
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
      size_t const nr_words = modmpz_intget_uls(t1, n);
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
      return nullptr;
  return m;
}

facul_status FaculModulusUl::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method,
    std::vector<std::unique_ptr<FaculModulusBase>> & composites,
    unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, composites, lpb,
        assume_prime_thresh, BBB);
}

cxx_mpz
FaculModulusUl::get_z () const
{
    return { m->m };
}

facul_status FaculModulus15Ul::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method,
    std::vector<std::unique_ptr<FaculModulusBase>> & composites,
    unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, composites, lpb,
        assume_prime_thresh, BBB);
}

cxx_mpz
FaculModulus15Ul::get_z () const
{
    cxx_mpz z;
    mpz_set_ui (z, m->m[1]);
    mpz_mul_2exp (z, z, ULONG_BITS);
    mpz_add_ui (z, z, m->m[0]);
    return z;
}

facul_status FaculModulus2Ul2::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method,
    std::vector<std::unique_ptr<FaculModulusBase>> & composites,
    unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, composites, lpb,
        assume_prime_thresh, BBB);
}

cxx_mpz
FaculModulus2Ul2::get_z () const
{
    cxx_mpz z;
    mpz_set_ui (z, m->m[1]);
    mpz_mul_2exp (z, z, ULONG_BITS);
    mpz_add_ui (z, z, m->m[0]);
    return z;
}

facul_status FaculModulusMpz::facul_doit_onefm (std::vector<cxx_mpz> & factors,
    facul_method const & method,
    std::vector<std::unique_ptr<FaculModulusBase>> & composites,
    unsigned long lpb, double assume_prime_thresh,
    double BBB) const
{
    return ::facul_doit_onefm (factors, m, method, composites, lpb,
        assume_prime_thresh, BBB);
}

cxx_mpz
FaculModulusMpz::get_z () const
{
    return { m };
}
