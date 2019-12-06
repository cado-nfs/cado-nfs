#ifndef MODSET_HPP_
#define MODSET_HPP_

#include <vector>
#include "cxx_mpz.hpp"
#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "facul_fwd.hpp"

/* This is kind of a poor man's runtime polymorphism. It's like a base class
 * of Modulus classes with a factory method that initialises the most
 * efficient modular arithmetic for a given modulus.
*/


struct modset_t {
  /* The arith variable tells which modulus type has been initialised for 
     arithmetic. It has a value of CHOOSE_NONE if no modulus currently 
     initialised. */
  enum {
      CHOOSE_NONE,
      CHOOSE_UL,
      CHOOSE_15UL,
      CHOOSE_2UL2,
      CHOOSE_MPZ,
  } arith;

  modulusredcul_t m_ul;
  modulusredc15ul_t m_15ul;
  modulusredc2ul2_t m_2ul2;
  modulusmpz_t m_mpz;
  
  modset_t() : arith(CHOOSE_NONE) {}
  ~modset_t() {clear();}
  
  void clear ();
  void init_ul (modintredcul_t m);
  void init_15ul (modintredc15ul_t m);
  void init_2ul2 (modintredc2ul2_t m);
  void init_mpz (modintmpz_t m);

  void get_z (mpz_t) const;
  int call_facul(std::vector<cxx_mpz> & factors,
    const facul_strategy_t *strategy, const int method_start) const;

  /* Run the relevant mod_isprime() function, using the arithmetic selected in the modset */
  int isprime () const
  {
    switch (arith) {
    case CHOOSE_UL:
      return modredcul_isprime (m_ul);
    case CHOOSE_15UL:
      return modredc15ul_isprime (m_15ul);
    case CHOOSE_2UL2:
      return modredc2ul2_isprime (m_2ul2);
    case CHOOSE_MPZ:
      return modmpz_isprime (m_mpz);
    default:
      abort();
    }
  }
};

#endif
