#include "cado.h"
#include "modset.hpp"
#include "facul_doit.hpp"

void
modset_t::clear ()
{
  switch (arith) {
      case CHOOSE_NONE: /* already clear */
    break;
      case CHOOSE_UL:
    modredcul_clearmod (m_ul);
    break;
      case CHOOSE_15UL:
    modredc15ul_clearmod (m_15ul);
    break;
      case CHOOSE_2UL2:
    modredc2ul2_clearmod (m_2ul2);
    break;
      case CHOOSE_MPZ:
    modmpz_clearmod (m_mpz);
    break;
  default:
    ASSERT_ALWAYS(0);
  }
  arith = CHOOSE_NONE;
}

/* This should be overloaded, but all the mod_int types are passed as
 * pointers to unsigned long, so they are indistinguishable in the
 * function signature. Thus we cannot use a function template, either,
 * as all the instances of the template would have the same signature. */

void
modset_t::init_ul (modintredcul_t m)
{
  ASSERT_ALWAYS(arith == CHOOSE_NONE);
  arith = CHOOSE_UL;
  modredcul_initmod_int (m_ul, m);
}

void
modset_t::init_15ul (modintredc15ul_t m)
{
  const size_t bits = modredc15ul_intbits (m);
  ASSERT_ALWAYS(arith == CHOOSE_NONE);
  if (bits <= MODREDCUL_MAXBITS)
    {
      unsigned long t1[1];
      size_t nr_words = modredc15ul_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 1);
      arith = modset_t::CHOOSE_UL;
      modredcul_initmod_ul (m_ul, t1[0]);
    }
  else if (bits <= MODREDC15UL_MAXBITS)
    {
      arith = CHOOSE_15UL;
      modredc15ul_initmod_int (m_15ul, m);
    }
  else
      abort();
}

void
modset_t::init_2ul2 (modintredc2ul2_t m)
{
  const size_t bits = modredc2ul2_intbits (m);
  ASSERT_ALWAYS(arith == CHOOSE_NONE);
  if (bits <= MODREDCUL_MAXBITS)
    {
      unsigned long t1[1];
      size_t nr_words = modredc2ul2_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 1);
      arith = CHOOSE_UL;
      modredcul_initmod_ul (m_ul, t1[0]);
    }
  else if (bits <= MODREDC15UL_MAXBITS)
    {
      unsigned long t1[2];
      modintredc15ul_t t2;
      size_t nr_words = modredc2ul2_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 2);
      modredc15ul_intset_uls (t2, t1, nr_words);
      arith = CHOOSE_15UL;
      modredc15ul_initmod_int (m_15ul, t2);
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      arith = CHOOSE_2UL2;
      modredc2ul2_initmod_int (m_2ul2, m);
    }
  else
      abort();
}

void
modset_t::init_mpz (modintmpz_t m)
{
  const size_t bits = modmpz_intbits (m);
  ASSERT_ALWAYS(arith == CHOOSE_NONE);
  if (bits <= MODREDCUL_MAXBITS)
    {
      unsigned long t1[1];
      size_t nr_words = modmpz_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 1);
      arith = CHOOSE_UL;
      modredcul_initmod_ul (m_ul, t1[0]);
    }
  else if (bits <= MODREDC15UL_MAXBITS)
    {
      unsigned long t1[2];
      modintredc15ul_t t2;
      size_t nr_words = modmpz_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 2);
      modredc15ul_intset_uls (t2, t1, nr_words);
      arith = CHOOSE_15UL;
      modredc15ul_initmod_int (m_15ul, t2);
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      unsigned long t1[2];
      modintredc2ul2_t t2;
      size_t nr_words = modmpz_intget_uls(t1, m);
      ASSERT_ALWAYS(nr_words <= 2);
      modredc2ul2_intset_uls (t2, t1, nr_words);
      arith = CHOOSE_2UL2;
      modredc2ul2_initmod_int (m_2ul2, t2);
    }
  else if (bits <= MODMPZ_MAXBITS)
    {
      /* We assume for now that m is a modintmpz_t */
      arith = CHOOSE_MPZ;
      modmpz_initmod_int (m_mpz, m);
    }
  else
      abort();
}

int 
modset_t::call_facul(std::vector<cxx_mpz> & factors, 
    const facul_strategy_t *strategy, const int method_start) const
{
  switch (arith) {
      case CHOOSE_UL:
          return facul_doit_ul (factors, m_ul, strategy, method_start);
      case CHOOSE_15UL:
          return facul_doit_15ul (factors, m_15ul, strategy, method_start);
          break;
      case CHOOSE_2UL2:
          return facul_doit_2ul2 (factors, m_2ul2, strategy, method_start);
          break;
      case CHOOSE_MPZ:
          return facul_doit_mpz (factors, m_mpz, strategy, method_start);
          break;
      default: abort();
  }
}

void
modset_t::get_z (mpz_t z) const
{
  ASSERT_ALWAYS(arith != CHOOSE_NONE);
  switch (arith)
    {
    case CHOOSE_UL:
      mpz_set_ui (z, m_ul->m);
      break;
    case CHOOSE_15UL:
      mpz_set_ui (z, m_15ul->m[1]);
      mpz_mul_2exp (z, z, LONG_BIT);
      mpz_add_ui (z, z, m_15ul->m[0]);
      break;
    case CHOOSE_2UL2:
      mpz_set_ui (z, m_2ul2->m[1]);
      mpz_mul_2exp (z, z, LONG_BIT);
      mpz_add_ui (z, z, m_2ul2->m[0]);
      break;
    case CHOOSE_MPZ:
      mpz_set (z, m_mpz);
      break;
    default:
      ASSERT_ALWAYS(0);
    }
}
