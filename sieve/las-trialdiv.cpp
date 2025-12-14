#include "cado.h" // IWYU pragma: keep
#include <cstddef>                    // for NULL, size_t
#include <algorithm>                   // for is_sorted, min
#include <array>                       // for array
#include <iterator>                    // for back_insert_iterator, back_ins...
#include <map>                         // for map, operator!=, _Rb_tree_iter...
#include <mutex>                       // for lock_guard, mutex
#include <utility>                     // for swap, pair
#include <vector>                      // for vector, vector<>::iterator
#include <gmp.h>                       // for mpz_divisible_ui_p, mpz_get_ui

#include "fb.hpp"                      // for fb_factorbase::key_type, fb_en...
#include "cxx_mpz.hpp"   // cxx_mpz
#include "getprime.h"   // prime_info getprime_mt
#include "las-sieve-shared-data.hpp"   // for sieve_shared_data::side_data
#include "lock_guarded_container.hpp"  // for lock_guarded_container
#include "macros.h"                    // for MIN, ASSERT, ASSERT_ALWAYS
#include "mpz_poly.h"
#include "rootfinder.h" // mpz_poly_roots_ulong
#include "trialdiv.hpp"                // for trialdiv_data, trialdiv_data::...

namespace {
template<typename T>
unsigned long
append_prime_list (T inserter, prime_info pi, unsigned long pmax, cxx_mpz_poly const & f, gmp_randstate_ptr rstate, unsigned int minroots = 1)
{
    unsigned long p;
    ASSERT_ALWAYS(f->deg >= 0 && minroots <= (unsigned int) f->deg);
    if (f->deg == 1) {
        for (; (p = getprime_mt (pi)) < pmax; )
            *inserter++ = p;
    } else {
        for (; (p = getprime_mt (pi)) < pmax; )
            if (mpz_divisible_ui_p (mpz_poly_lc(f), p) ||
                    mpz_poly_roots_ulong (nullptr, f, p, rstate) >= minroots)
                *inserter++ = p;
    }
    return p;
}
}

trialdiv_data const * sieve_shared_data::side_data::get_trialdiv_data(fb_factorbase::key_type fbK, fb_factorbase::slicing const * fbs)
{
    const std::lock_guard<std::mutex> foo(trialdiv_data_cache.mutex());
    auto it = trialdiv_data_cache.find(fbK);
    if (it != trialdiv_data_cache.end()) {
        return &it->second;
    }

    ASSERT_ALWAYS(fbs);

    /* Note that since we have trialdiv_data_cache.mutex() unlock, we may
     * safely access the random state in this->rstate
     */

    /* Now compute the trialdiv data for these thresholds. */

    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    const unsigned long pmax = std::min((unsigned long) fbK.thresholds[0],
                             trialdiv_data::max_p);

    std::vector<unsigned long> trialdiv_primes = fbs->small_sieve_entries.skipped;

    /* Maybe we can use the factor base. If we have one, of course ! */
    unsigned long pmax_sofar = 0;

    for(auto const & pp : fbs->small_sieve_entries.rest) {
        if (pp.k > 1) continue;
        trialdiv_primes.push_back(pp.p);
    }
    if (!trialdiv_primes.empty()) {
        cxx_mpz zz(trialdiv_primes.back());
        mpz_nextprime(zz, zz);
        pmax_sofar = MIN(pmax, mpz_get_ui(zz));
    }

    if (pmax_sofar < pmax) {
        /* we need some more. */
        prime_info pi;
        prime_info_init(pi);

        prime_info_seek(pi, pmax_sofar);

        for(int minroots = 1 ; minroots <= f->deg ; minroots++) {
            append_prime_list(std::back_inserter(trialdiv_primes),
                    pi, MIN(pmax, minroots * fbK.td_thresh), f, rstate, minroots);
        }
        prime_info_clear (pi);
    }

    ASSERT(std::ranges::is_sorted(trialdiv_primes));
    // std::ranges::sort(trialdiv_primes);
    
    /* note that we might have several "2"'s in the factor base because
     * of powers: when bucket-sieving powers, we separate factor base
     * entries that lead to distinct exponents */
    size_t skip2 = 0;
    for( ; skip2 < trialdiv_primes.size() && trialdiv_primes[skip2] == 2 ; skip2++);

    trialdiv_data td(trialdiv_primes, skip2);
    trialdiv_data_cache[fbK];
    std::swap(trialdiv_data_cache[fbK], td);

    return &trialdiv_data_cache[fbK];
}
