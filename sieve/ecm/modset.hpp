#ifndef MODSET_HPP_
#define MODSET_HPP_

#include <vector>
#include "cxx_mpz.hpp"
#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "facul_fwd.hpp"


class FaculModulusBase {
public:
  FaculModulusBase() {}
  virtual ~FaculModulusBase() {}

  /* Factory methods that initialise the subclass with the fastest modular
   * arithmetic that's large enough for the integer n*/
  static const FaculModulusBase *init_ul (const modintredcul_t n);
  static const FaculModulusBase *init_15ul (const modintredc15ul_t n);
  static const FaculModulusBase *init_2ul2 (const modintredc2ul2_t n);
  static const FaculModulusBase *init_mpz (const modintmpz_t n);

  virtual void get_z (mpz_t) const = 0;
  virtual int call_facul(std::vector<cxx_mpz> &, const facul_strategy_t *,
    const int) const = 0;
  virtual int facul_doit_onefm (std::vector<cxx_mpz> &,
    const facul_method_t, const FaculModulusBase * &,
    const FaculModulusBase * &, unsigned long, double, double) const = 0;
  virtual int isprime () const = 0;
};

class FaculModulusUl : public FaculModulusBase {
    modulusredcul_t m;
public:
    FaculModulusUl(const modintredcul_t n) {modredcul_initmod_int(m, n);}
    ~FaculModulusUl() {modredcul_clearmod(m);}
    int isprime () const {return modredcul_isprime (m);}
    void get_z (mpz_t) const;
    int call_facul(std::vector<cxx_mpz> & factors,
        const facul_strategy_t *strategy, const int method_start) const;
    int facul_doit_onefm (std::vector<cxx_mpz> &,
        const facul_method_t, const FaculModulusBase * &,
        const FaculModulusBase * &, unsigned long, double, double) const;
};

class FaculModulus15Ul : public FaculModulusBase {
    modulusredc15ul_t m;
public:
    FaculModulus15Ul(const modintredc15ul_t n) {modredc15ul_initmod_int(m, n);}
    ~FaculModulus15Ul() {modredc15ul_clearmod(m);}
    int isprime () const {return modredc15ul_isprime (m);}
    void get_z (mpz_t) const;
    int call_facul(std::vector<cxx_mpz> & factors,
        const facul_strategy_t *strategy, const int method_start) const;
    int facul_doit_onefm (std::vector<cxx_mpz> &,
        const facul_method_t, const FaculModulusBase * &,
        const FaculModulusBase * &, unsigned long, double, double) const;
};

class FaculModulus2Ul2 : public FaculModulusBase {
    modulusredc2ul2_t m;
public:
    FaculModulus2Ul2(const modintredc2ul2_t n) {modredc2ul2_initmod_int(m, n);}
    ~FaculModulus2Ul2() {modredc2ul2_clearmod(m);}
    int isprime () const {return modredc2ul2_isprime (m);}
    void get_z (mpz_t) const;
    int call_facul(std::vector<cxx_mpz> & factors,
        const facul_strategy_t *strategy, const int method_start) const;
    int facul_doit_onefm (std::vector<cxx_mpz> &,
        const facul_method_t, const FaculModulusBase * &,
        const FaculModulusBase * &, unsigned long, double, double)const ;
};

class FaculModulusMpz : public FaculModulusBase {
    modulusmpz_t m;
public:
    FaculModulusMpz(const modintmpz_t n) {modmpz_initmod_int(m, n);}
    ~FaculModulusMpz() {modmpz_clearmod(m);}
    int isprime () const {return modmpz_isprime (m);}
    void get_z (mpz_t) const;
    int call_facul(std::vector<cxx_mpz> & factors,
        const facul_strategy_t *strategy, const int method_start) const;
    int facul_doit_onefm (std::vector<cxx_mpz> &,
        const facul_method_t, const FaculModulusBase * &,
        const FaculModulusBase * &, unsigned long, double, double) const;
};

#endif
