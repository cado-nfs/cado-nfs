#include "cado.h" // IWYU pragma: keep

#ifndef SUPPORT_LARGE_Q
#include <cinttypes>
#endif
#include <cstddef>

#include <array>
#include <mutex>
#include <vector>

#include <gmp.h>

#include "ecm/batch.hpp"
#ifndef SUPPORT_LARGE_Q
#include "gcd.h"
#endif
#ifdef SUPPORT_LARGE_Q
#include "cxx_mpz.hpp"
#endif
#include "las-cofac-standalone.hpp"
#include "las-cofactor.hpp"
#include "las-coordinates.hpp"
#include "las-siever-config.hpp"
#include "las-threads-work-data.hpp"
#include "special-q.hpp"
#include "las-where-am-i-proxy.hpp"
#include "lock_guarded_container.hpp"
#include "relation.hpp"

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
bool cofac_standalone::gcd_coprime_with_q(special_q const & E) const {/*{{{*/
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
    gmp_fprintf(f, "%" PRId64 " %" PRIu64, a, b);
#else
    gmp_fprintf(f, "%Zd %Zd", (mpz_srcptr) az, (mpz_srcptr) bz);
#endif
    for (auto const & n: norm)
        gmp_fprintf(f, " %Zd", (mpz_srcptr) n);
    fprintf(f, "\n");
}/*}}}*/
relation cofac_standalone::get_relation(special_q const & doing) const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
    relation rel(a, b);
#else
    relation rel(az, bz);
#endif

    /* Note that we explicitly do not bother about storing r in
     * the relations below */
    for (unsigned int side = 0; side < std::min(rel.sides.size(), factors.size()); side++) { // FIXME workaround for HARDCODED 2
        for (auto const& z : factors[side])
            rel.add(side, z, 0);
        for (auto const& z : lps[side])
            rel.add(side, z, cxx_mpz(0));
    }
    if (doing.is_prime()) {
        rel.add(doing.side, doing.p, cxx_mpz(0));
    } else {
        for (auto const& facq : doing.prime_factors)
            rel.add(doing.side, facq, cxx_mpz(0));
    }

    rel.compress();
    return rel;
}/*}}}*/
void cofac_standalone::transfer_to_cofac_list(lock_guarded_container<std::list<cofac_candidate>> & L) {/*{{{*/
    std::lock_guard<std::mutex> const foo(L.mutex());
    L.emplace_back(a, b, norm);
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
int cofac_standalone::factor_leftover_norms(nfs_work_cofac & wc) {/*{{{*/
    /* This proxies to las-cofactor.cpp */
    std::vector<unsigned long> Bs;
    Bs.reserve(wc.sc.sides.size());
    for (auto const & s: wc.sc.sides)
        Bs.push_back(s.lim);
    return ::factor_leftover_norms(norm,
            lps,
            Bs,
            *wc.strategies);
}/*}}}*/

