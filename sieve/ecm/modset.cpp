#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <memory>
#include <vector>

#include <gmp.h>

#include "arithxx/mod_mpz_new.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc96.hpp"
#include "cxx_mpz.hpp"
#include "facul.hpp"
#include "facul_doit.hpp"
#include "facul_method.hpp"
#include "gmp_aux.h"
#include "modset.hpp"

FaculModulusBase * FaculModulusBase::init_mpz(cxx_mpz const & n)
{
    size_t const bits = mpz_sizeinbase(n, 2);
    /* The Montgomery representation needs an odd modulus */
    if (mpz_odd_p(n)) {
        if (bits <= 64)
            return new FaculModulus<arithxx_modredc64>(Integer64(mpz_getlimbn_uint64(n, 0)));
        if (bits <= 96)
            return new FaculModulus<arithxx_modredc96>(Integer128(mpz_getlimbn_uint64(n, 0), mpz_getlimbn_uint64(n, 1)));
        if (bits <= 126)
            return new FaculModulus<arithxx_modredc126>(Integer128(mpz_getlimbn_uint64(n, 0), mpz_getlimbn_uint64(n, 1)));
    }
    return new FaculModulus<arithxx_mod_mpz_new>(n);
}

template<typename layer>
cxx_mpz FaculModulus<layer>::get_z() const
{
    return cxx_mpz(m.getmod());
}

template<typename layer>
int FaculModulus<layer>::isprime() const
{
    return m.is_prime();
}

template<typename layer>
facul_status FaculModulus<layer>::facul_doit_onefm(
        std::vector<cxx_mpz> & factors,
        facul_method const & method,
        std::vector<std::unique_ptr<FaculModulusBase>> & composites,
        unsigned long const lpb, double const BB, double const BBB) const
{
    return ::facul_doit_onefm<layer>(factors, m, method, composites, lpb, BB, BBB);
}

template class FaculModulus<arithxx_modredc64>;
template class FaculModulus<arithxx_modredc96>;
template class FaculModulus<arithxx_modredc126>;
template class FaculModulus<arithxx_mod_mpz_new>;
