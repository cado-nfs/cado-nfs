#ifndef CADO_FACUL_DOIT_HPP
#define CADO_FACUL_DOIT_HPP

#include <memory>
#include <vector>

#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_method.hpp"
#include "modset.hpp"

/* Apply one factoring method to the number n, given by its modulus m
 * over an arithxx layer.
 *
 * The prime factors found are appended to [factors]. The [composites]
 * argument must be empty on input, and all composites that are found are
 * stored there. Note that it is possible to find composites even if no
 * prime factor is found!
 *
 * This returns FACUL_NOT_SMOOTH, FACUL_MAYBE, or FACUL_SMOOTH.
 *
 * facul_doit.cpp instantiates this for modredc64, modredc96, modredc126
 * and mod_mpz_new.
 */
template<typename layer>
facul_status facul_doit_onefm(std::vector<cxx_mpz> & factors,
        typename layer::Modulus const & m,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long lpb, double BB, double BBB);

/* The same, for a number n given as an mpz (with the arithxx layer that
 * suits its size). */
facul_status facul_doit_onefm(std::vector<cxx_mpz> & factors,
        cxx_mpz const & n,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long lpb, double BB, double BBB);

#endif /* CADO_FACUL_DOIT_HPP */
