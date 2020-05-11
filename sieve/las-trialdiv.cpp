#include "cado.h"
#include <iterator>
#include "las-info.hpp"

template<typename T>
static unsigned long
append_prime_list (T inserter, prime_info pi, unsigned long pmax, cxx_mpz_poly const & f, int minroots = 1)
{
    unsigned long p;
    ASSERT_ALWAYS(minroots <= f->deg);
    if (f->deg == 1) {
        for (; (p = getprime_mt (pi)) < pmax; )
            *inserter++ = p;
    } else {
        for (; (p = getprime_mt (pi)) < pmax; )
            if (mpz_divisible_ui_p (f->coeff[f->deg], p) ||
                    mpz_poly_roots_ulong (NULL, f, p) >= minroots)
                *inserter++ = p;
    }
    return p;
}


trialdiv_data const * sieve_shared_data::side_data::get_trialdiv_data(fb_factorbase::key_type fbK, fb_factorbase::slicing const * fbs)
{
    std::lock_guard<std::mutex> foo(trialdiv_data_cache.mutex());
    auto it = trialdiv_data_cache.find(fbK);
    if (it != trialdiv_data_cache.end()) {
        return &it->second;
    }

    /* Now compute the trialdiv data for these thresholds. */

    /* Our trial division needs odd divisors, 2 is handled by mpz_even_p().
       If the FB primes to trial divide contain 2, we skip over it.
       We assume that if 2 is in the list, it is the first list entry,
       and that it appears at most once. */

    unsigned long pmax = std::min((unsigned long) fbK.thresholds[0],
                             trialdiv_data::max_p);

    std::vector<unsigned long> trialdiv_primes = fbs->small_sieve_entries.skipped;

    /* Maybe we can use the factor base. If we have one, of course ! */
    unsigned long pmax_sofar = 0;
    if (fbs) {
        for(auto const & pp : fbs->small_sieve_entries.rest) {
            if (pp.k > 1) continue;
            trialdiv_primes.push_back(pp.p);
        }
        if (!trialdiv_primes.empty()) {
            cxx_mpz zz(trialdiv_primes.back());
            mpz_nextprime(zz, zz);
            pmax_sofar = MIN(pmax, mpz_get_ui(zz));
        }
    }
    if (pmax_sofar < pmax) {
        /* we need some more. */
        prime_info pi;
        prime_info_init(pi);
        unsigned long p;
        /* first seek to the end of the fb. */
        for ( ; (p = getprime_mt (pi)) < pmax_sofar ; );

        for(int minroots = 1 ; minroots <= f->deg ; minroots++) {
            p = append_prime_list(std::back_inserter(trialdiv_primes),
                    pi, MIN(pmax, minroots * fbK.td_thresh), f, minroots);
        }
        prime_info_clear (pi);
    }

    ASSERT(std::is_sorted(trialdiv_primes.begin(), trialdiv_primes.end()));
    // std::sort(trialdiv_primes.begin(), trialdiv_primes.end());
    
    size_t skip2 = !trialdiv_primes.empty() && trialdiv_primes[0] == 2;

    trialdiv_data td(trialdiv_primes, skip2);
    trialdiv_data_cache[fbK];
    std::swap(trialdiv_data_cache[fbK], td);

    return &trialdiv_data_cache[fbK];
}
