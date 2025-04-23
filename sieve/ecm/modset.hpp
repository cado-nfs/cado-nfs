#ifndef CADO_MODSET_HPP
#define CADO_MODSET_HPP

#include <memory>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_method.hpp"
#include "arith/mod_mpz.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/modredc_ul.h"

struct cxx_mpz;

class FaculModulusBase
{
  public:
    FaculModulusBase() = default;
    virtual ~FaculModulusBase() = default;

    /* Factory methods that initialise the subclass with the fastest modular
     * arithmetic that's large enough for the integer n*/
    static FaculModulusBase * init_ul(modintredcul_t const n);
    static FaculModulusBase * init_15ul(modintredc15ul_t const n);
    static FaculModulusBase * init_2ul2(modintredc2ul2_t const n);
    static FaculModulusBase * init_mpz(modintmpz_t const n);

    virtual cxx_mpz get_z() const = 0;
    virtual facul_status facul_doit_onefm(
            std::vector<cxx_mpz> &,
            facul_method const &,
            std::vector<std::unique_ptr<FaculModulusBase>> &,
            unsigned long,
            double, double) const = 0;
    virtual int isprime() const = 0;
};

class FaculModulusUl : public FaculModulusBase
{
    modulusredcul_t m;

  public:
    FaculModulusUl(modintredcul_t const n) { modredcul_initmod_int(m, n); }
    ~FaculModulusUl() override { modredcul_clearmod(m); }
    int isprime() const override { return modredcul_isprime(m); }
    cxx_mpz get_z() const override;
    facul_status facul_doit_onefm(std::vector<cxx_mpz> &, facul_method const &,
            std::vector<std::unique_ptr<FaculModulusBase>> &,
                         unsigned long, double, double) const override;
};

class FaculModulus15Ul : public FaculModulusBase
{
    modulusredc15ul_t m;

  public:
    FaculModulus15Ul(modintredc15ul_t const n)
    {
        modredc15ul_initmod_int(m, n);
    }
    ~FaculModulus15Ul() override { modredc15ul_clearmod(m); }
    int isprime() const override { return modredc15ul_isprime(m); }
    cxx_mpz get_z() const override;
    facul_status facul_doit_onefm(std::vector<cxx_mpz> &, facul_method const &,
            std::vector<std::unique_ptr<FaculModulusBase>> &,
                         unsigned long, double, double) const override;
};

class FaculModulus2Ul2 : public FaculModulusBase
{
    modulusredc2ul2_t m;

  public:
    FaculModulus2Ul2(modintredc2ul2_t const n)
    {
        modredc2ul2_initmod_int(m, n);
    }
    ~FaculModulus2Ul2() override { modredc2ul2_clearmod(m); }
    int isprime() const override { return modredc2ul2_isprime(m); }
    cxx_mpz get_z() const override;
    facul_status facul_doit_onefm(std::vector<cxx_mpz> &, facul_method const &,
                         std::vector<std::unique_ptr<FaculModulusBase>> &,
                         unsigned long, double, double) const override;
};

class FaculModulusMpz : public FaculModulusBase
{
    modulusmpz_t m;

  public:
    FaculModulusMpz(modintmpz_t const n) { modmpz_initmod_int(m, n); }
    ~FaculModulusMpz() override { modmpz_clearmod(m); }
    int isprime() const override { return modmpz_isprime(m); }
    cxx_mpz get_z() const override;
    facul_status facul_doit_onefm(std::vector<cxx_mpz> &, facul_method const &,
                         std::vector<std::unique_ptr<FaculModulusBase>> &,
                         unsigned long, double, double) const override;
};

#endif
