#ifndef CADO_SIEVE_ECM_FACUL_ARITHXX_HPP
#define CADO_SIEVE_ECM_FACUL_ARITHXX_HPP

/* Running a factoring method written for the arithxx layers on the
 * moduli of the old arith layer, which facul_doit still uses.
 *
 * facul_arithxx_run(f, m, method) views the old modulus m as an mpz
 * without copying it, picks the arithxx layer that suits its size, calls
 * method.template operator()<layer>(g, mm) with an arithxx modulus mm
 * and an Integer g of that layer, and stores the factor g in f.
 */

#include <cstddef>
#include <cstdint>
#include <climits>

#include <type_traits>

#include <gmp.h>

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc96.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/mod_mpz_new.hpp"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"

namespace facul_arithxx_details {

/* Run the method on the modulus N with the smallest arithxx layer that
 * fits, and give the factor found to set_f, in the Integer type of that
 * layer. */
template<typename Method, typename SetF>
int dispatch(mpz_srcptr N, Method && method, SetF && set_f)
{
    size_t const bits = mpz_sizeinbase(N, 2);
    auto run = [&]<typename layer>(typename layer::Integer const & n) {
        typename layer::Modulus const m(n);
        typename layer::Integer f;
        int const bt = method.template operator()<layer>(f, m);
        set_f(f);
        return bt;
    };
    if (bits <= 64)
        return run.template operator()<arithxx_modredc64>(
                Integer64(mpz_getlimbn_uint64(N, 0)));
    if (bits <= 96)
        return run.template operator()<arithxx_modredc96>(
                Integer128(mpz_getlimbn_uint64(N, 0), mpz_getlimbn_uint64(N, 1)));
    if (bits <= 126)
        return run.template operator()<arithxx_modredc126>(
                Integer128(mpz_getlimbn_uint64(N, 0), mpz_getlimbn_uint64(N, 1)));
    return run.template operator()<arithxx_mod_mpz_new>(cxx_mpz(N));
}

/* an mpz view of the words of an old modulus */
inline void mpz_view(mpz_t z, unsigned long const * w, size_t n)
{
    static_assert(sizeof(mp_limb_t) == sizeof(unsigned long));
    mpz_roinit_n(z, reinterpret_cast<mp_limb_t const *>(w), mp_size_t(n));
}

/* the unsigned long words of an arithxx Integer, least significant
 * first, without leading zeros (but at least one word) */
template<typename I>
size_t ulong_words(unsigned long * w, size_t n, I const & x)
{
    unsigned long t[8];
    size_t k = 0;
    if constexpr (std::is_same_v<I, cxx_mpz>) {
        ASSERT_ALWAYS(mpz_sizeinbase(x, 2) <= n * ULONG_BITS);
        mpz_export(t, &k, -1, sizeof(unsigned long), 0, 0, x);
    } else {
        for (size_t i = 0; i < I::max_size_in_words; i++) {
            uint64_t const v = x[i];
#if ULONG_BITS == 64
            t[k++] = v;
#else
            t[k++] = (unsigned long) v;
            t[k++] = (unsigned long) (v >> 32);
#endif
        }
        for (; k > 0 && t[k - 1] == 0; k--);
    }
    ASSERT_ALWAYS(k <= n);
    if (k == 0)
        t[k++] = 0;
    for (size_t i = 0; i < k; i++)
        w[i] = t[i];
    return k;
}

} /* namespace facul_arithxx_details */

template<typename Method>
int facul_arithxx_run(modintredcul_t f, const modulusredcul_t m, Method && method)
{
    using namespace facul_arithxx_details;
    mpz_t N;
    mpz_view(N, &m[0].m, 1);
    return dispatch(N, method, [&](auto const & g) {
            unsigned long w[MODREDCUL_SIZE];
            modredcul_intset_uls(f, w, ulong_words(w, MODREDCUL_SIZE, g));
            });
}

template<typename Method>
int facul_arithxx_run(modintredc15ul_t f, const modulusredc15ul_t m, Method && method)
{
    using namespace facul_arithxx_details;
    mpz_t N;
    mpz_view(N, m[0].m, MODREDC15UL_SIZE);
    return dispatch(N, method, [&](auto const & g) {
            unsigned long w[MODREDC15UL_SIZE];
            modredc15ul_intset_uls(f, w, ulong_words(w, MODREDC15UL_SIZE, g));
            });
}

template<typename Method>
int facul_arithxx_run(modintredc2ul2_t f, const modulusredc2ul2_t m, Method && method)
{
    using namespace facul_arithxx_details;
    mpz_t N;
    mpz_view(N, m[0].m, MODREDC2UL2_SIZE);
    return dispatch(N, method, [&](auto const & g) {
            unsigned long w[MODREDC2UL2_SIZE];
            modredc2ul2_intset_uls(f, w, ulong_words(w, MODREDC2UL2_SIZE, g));
            });
}

template<typename Method>
int facul_arithxx_run(modintmpz_t f, const modulusmpz_t m, Method && method)
{
    using namespace facul_arithxx_details;
    return dispatch(m, method, [&](auto const & g) {
            if constexpr (std::is_same_v<std::decay_t<decltype(g)>, cxx_mpz>) {
                mpz_set(f, g);
            } else {
                unsigned long w[4];
                size_t const n = ulong_words(w, 4, g);
                mpz_import(f, n, -1, sizeof(unsigned long), 0, 0, w);
            }
            });
}

#endif /* CADO_SIEVE_ECM_FACUL_ARITHXX_HPP */
