#include "cado.h"
#include "las-todo-entry.hpp"
#include "getprime.h"

void las_todo_entry::set(const mpz_t _p, const mpz_t _r, const int _side, const int _depth, const int _iteration) {
    mpz_set(p, _p);
    mpz_set(r, _r);
    side = _side;
    depth = _depth;
    iteration = _iteration;
    prime_sq = mpz_probab_prime_p(p, 2);
    if (prime_sq) return;

    // Need to pre-compute the prime factors of the special-q in order to
    // skip them while sieving.
    prime_info pi;
    prime_info_init (pi);
    unsigned long f = 2;
    mpz_t B;
    mpz_init(B);
    mpz_sqrt(B, p);
    unsigned long bound = mpz_get_ui(B);
    mpz_clear(B);
    mpz_t pp;
    mpz_init_set(pp, p);
    while ((f <= bound) && (mpz_cmp_ui(pp, 1) > 0)) {
        if (mpz_divisible_ui_p(pp, f)) {
            mpz_divexact_ui(pp, pp, f);
            // Powers are not allowed in special-q
            ASSERT_ALWAYS(!mpz_divisible_ui_p(pp, f));
            prime_factors.push_back(f);
            if (mpz_probab_prime_p(pp, 2)) {
                ASSERT_ALWAYS(mpz_fits_ulong_p(pp));
                prime_factors.push_back(mpz_get_ui(pp));
                mpz_set_ui(pp, 1);
            }
        }
        f = getprime_mt(pi);
    }
    prime_info_clear (pi);
    mpz_clear(pp);
}
