#ifndef CADO_MODSET_HPP
#define CADO_MODSET_HPP

#include <memory>
#include <vector>

#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_method.hpp"

/* A number to be factored, with the arithxx modulus that suits its size.
 * The modulus is set up once, and serves for all the factoring methods
 * that are tried on this number.
 */
class FaculModulusBase
{
  public:
    FaculModulusBase() = default;
    virtual ~FaculModulusBase() = default;
    FaculModulusBase(FaculModulusBase const &) = delete;
    FaculModulusBase & operator=(FaculModulusBase const &) = delete;

    /* Factory method that uses the fastest modular arithmetic that's
     * large enough for n, which must be odd. */
    static FaculModulusBase * init_mpz(cxx_mpz const & n);

    virtual cxx_mpz get_z() const = 0;
    virtual facul_status facul_doit_onefm(
            std::vector<cxx_mpz> & factors,
            facul_method const & method,
            std::vector<std::unique_ptr<FaculModulusBase>> & composites,
            unsigned long lpb,
            double BB, double BBB) const = 0;
    virtual int isprime() const = 0;
};

/* modset.cpp instantiates this for modredc64, modredc96, modredc126 and
 * mod_mpz_new. */
template<typename layer>
class FaculModulus : public FaculModulusBase
{
    typename layer::Modulus m;

  public:
    explicit FaculModulus(typename layer::Integer const & n) : m(n) {}
    cxx_mpz get_z() const override;
    int isprime() const override;
    facul_status facul_doit_onefm(std::vector<cxx_mpz> &,
            facul_method const &,
            std::vector<std::unique_ptr<FaculModulusBase>> &,
            unsigned long, double, double) const override;
};

#endif
