#include "cado.h" // IWYU pragma: keep

#include <inttypes.h> // for PRId64, PRIu64 // IWYU pragma: keep
#include <array>                       // for array, array<>::value_type
#include <cstdint>                     // for uint8_t
#include <mutex>                       // for lock_guard, mutex
#include <vector>                      // for vector

#include <gmp.h>                       // for mpz_srcptr, gmp_fprintf, mpz_c...

#include "las-cofac-standalone.hpp"    // for cofac_standalone
#include "gcd.h"       // for bin_gcd_int64_safe // IWYU pragma: keep
#include "ecm/batch.hpp"                   // for cofac_list
#include "las-cofactor.hpp"            // for factor_both_leftover_norms
#include "las-coordinates.hpp"         // for convert_Nx_to_ab
#include "las-divide-primes.hpp"       // for factor_list_t
#include "las-siever-config.hpp"       // for siever_config::side_config
#include "las-threads-work-data.hpp"   // for nfs_work_cofac
#include "las-todo-entry.hpp"          // for las_todo_entry
#include "las-where-am-i-proxy.hpp"    // for extern_trace_on_spot_ab
#include "lock_guarded_container.hpp"  // for lock_guarded_container
#include "relation.hpp"                // for relation

struct qlattice_basis; // IWYU pragma: keep

/* This is one input to the late cofactoring process (aka ECM). Here, we
 * mean the stuff that is done detached from the rest of the siever
 * stuff: we no longer care about purging buckets and so on, these may
 * safely be used for later work.
 */
cofac_standalone::cofac_standalone() : a(0), b(0) {/*{{{*/
#ifdef SUPPORT_LARGE_Q
    mpz_set_ui(az, 0);
    mpz_set_ui(bz, 0);
#endif
}/*}}}*/
cofac_standalone::cofac_standalone(int nsides, int N, size_t x, int logI, qlattice_basis const & Q)
    : S(nsides, 0)
    , norm(nsides, 0)
    , factors(nsides)
    , lps(nsides)
{/*{{{*/
    convert_Nx_to_ab (a, b, N, x, logI, Q);
#ifdef SUPPORT_LARGE_Q
    convert_Nx_to_abmpz (az, bz, N, x, logI, Q);
#endif
}/*}}}*/
bool cofac_standalone::trace_on_spot() const {/*{{{*/
    return extern_trace_on_spot_ab(a, b);
}/*}}}*/
bool cofac_standalone::gcd_coprime_with_q(las_todo_entry const & E) {/*{{{*/
    /* Since the q-lattice is exactly those (a, b) with
       a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
    /* In case of composite sq, have to check all factors... */
    /* FIXME: fast divisibility test here! */
    if (E.is_prime()) {
#ifndef SUPPORT_LARGE_Q
        if (b == 0 || (mpz_cmp_ui(E.p, b) <= 0 && b % mpz_get_ui(E.p) == 0))
#else
            if ((mpz_cmp_ui(bz, 0) == 0) || 
                    (mpz_cmp(E.p, bz) <= 0 &&
                     mpz_divisible_p(bz, E.p)))
#endif
                return false;
    } else {
#ifdef SUPPORT_LARGE_Q
        if (mpz_cmp_ui(bz, 0) == 0)
            return false;
        for (auto const& facq : E.prime_factors) {
            if ((mpz_cmp_ui(bz, facq) >= 0) && (mpz_divisible_ui_p(bz, facq))) {
                return false;
            }
        }
#else
        if (b == 0)
            return false;
        for (auto const& facq : E.prime_factors) {
            if (facq <= b && b % facq == 0) {
                return false;
            }
        }
#endif
    }
    return true;
}/*}}}*/
bool cofac_standalone::ab_coprime() const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
    return bin_gcd_int64_safe(a,b) == 1;
#else
    cxx_mpz g;
    mpz_gcd(g, az, bz);
    return mpz_cmp_ui(g, 1) == 0;
#endif
}/*}}}*/
void cofac_standalone::print_as_survivor(FILE * f) {/*{{{*/
#ifndef SUPPORT_LARGE_Q
    gmp_fprintf(f, "%" PRId64 " %" PRIu64 " %Zd %Zd\n", a, b,
            (mpz_srcptr) norm[0],
            (mpz_srcptr) norm[1]);
#else
    gmp_fprintf(f, "%Zd %Zd %Zd %Zd\n",
            (mpz_srcptr) az,
            (mpz_srcptr) bz,
            (mpz_srcptr) norm[0],
            (mpz_srcptr) norm[1]);
#endif
}/*}}}*/
relation cofac_standalone::get_relation(las_todo_entry const & doing) {/*{{{*/
#ifndef SUPPORT_LARGE_Q
    relation rel(a, b);
#else
    relation rel(az, bz);
#endif

    /* Note that we explicitly do not bother about storing r in
     * the relations below */
    for (int side = 0; side < 2; side++) {
        for (auto const& z : factors[side])
            rel.add(side, z, 0);
        for (auto const& z : lps[side])
            rel.add(side, z, 0);
    }
    if (doing.is_prime()) {
        rel.add(doing.side, doing.p, 0);
    } else {
        for (auto const& facq : doing.prime_factors)
            rel.add(doing.side, facq, 0);
    }

    rel.compress();
    return rel;
}/*}}}*/
void cofac_standalone::transfer_to_cofac_list(lock_guarded_container<cofac_list> & L, las_todo_entry const & doing) {/*{{{*/
    std::lock_guard<std::mutex> foo(L.mutex());
    /* "doing" must be an object that lives in the main todo list,
     * and will stay alive until the end of the program. Yes, it's a
     * bit dangerous. */
    L.emplace_back(a, b, norm, &doing);
#if 0
    /* make sure threads don't write the cofactor list at the
     * same time !!! */
#warning "possible contention"
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    cofac_list_add (L, a, b, norm, doing.side, doing.p);
    pthread_mutex_unlock(&lock);
#endif
}/*}}}*/
int cofac_standalone::factor_both_leftover_norms(nfs_work_cofac & wc) {/*{{{*/
    /* This proxies to las-cofactor.cpp */
    return ::factor_both_leftover_norms(norm,
            lps,
            {{ wc.sc.sides[0].lim, wc.sc.sides[1].lim }},
            *wc.strategies);
}/*}}}*/

