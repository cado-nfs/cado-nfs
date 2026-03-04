#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <limits>

#include "las-arith.hpp"
#include "las-smallsieve.hpp"
#include "macros.h"
#include "bucket-push-update.hpp"
#include "bucket.hpp"
#include "las-where-am-i.hpp"
#include "las-fbroot-qlattice.hpp"
#include "las-qlattice.hpp"
#include "las-sieve2357.hpp"
#include "siqs-smallsieve.hpp"
#include "siqs-smallsieve-glue.hpp"
#include "las-where-am-i-proxy.hpp"
#include "verbose.h"

using preferred_simd_type = sieve2357base::preferred_simd_type;
#if GNUC_VERSION_ATLEAST(6,1,0)
/* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69884 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
using preferred_sieve2357 = sieve2357<preferred_simd_type, uint8_t>;
#if GNUC_VERSION_ATLEAST(6,1,0)
#pragma GCC diagnostic pop
#endif


siqs_ssp_pdata_t::siqs_ssp_pdata_t(
        fb_entry_general const & e,
        siqs_special_q_data const & Q)
    : p(e.q)
    , pmask((p & 1u) ? std::numeric_limits<fbprime_t>::max() : (p-1U))
    , offset(((uint64_t (p-1u)) << LOG_BUCKET_REGION) % p)
{
    /* call to e.compute_crt_data_modp sets invq and crt_data_modp */
    rq0 = e.compute_crt_data_modp(invq, crt_data_modp, Q, true);
}

fbroot_t
siqs_ssp_pdata_t::transform_root(fbroot_t rp, redc_invp_t invp) const
{
    if (UNLIKELY(!(p & 1))) { /* p power of 2 */
        return (rp * invq - rq0) & pmask;
    } else {
        fbroot_t t = mulmodredc_u32<true>(rp, invq, p, invp);
        return submod_u32(t, rq0, p);
    }
}

auto
fmt::formatter<siqs_ssp_simple_t>::format(
        siqs_ssp_simple_t const & a,
        format_context & ctx) const -> format_context::iterator
{
    format_to(ctx.out(), "# p = {}, r = {}, logp = {}",
            a.get_p(), a.get_r(), a.get_logp());
    return ctx.out();
}

auto
fmt::formatter<siqs_ssp_t>::format(
        siqs_ssp_t const & a,
        format_context & ctx) const -> fmt::format_context::iterator
{
    format_to(ctx.out(), "# p = {}, r = {}, logp = {}",
            a.get_p(), a.get_r(), a.get_logp());
    if (a.is_pow2())
        format_to(ctx.out(), " (power of 2)");
    if (a.is_pattern_sieved())
        format_to(ctx.out(), " (pattern-sieved)");
    return ctx.out();
}

siqs_ssp_t::siqs_ssp_t(
        std::shared_ptr<siqs_ssp_pdata_t> const & pdata,
        fbprime_t r,
        unsigned char logp)
    : siqs_ssp_simple_t(pdata, r, logp)
{
    if (pdata->get_p() % 2 == 0) {
        flags |= SSP_POW2;
    }
    if (preferred_sieve2357::can_sieve(this->get_p())) {
        flags |= SSP_PATTERN_SIEVED;
    }
}

/* Functor for sorting siqs_ssp_t in the order in which sieve2357 expects them
 */
struct order_ssp_t {
    bool operator()(siqs_ssp_t const & ssp1, siqs_ssp_t const & ssp2) const {
        return sieve2357base::order_lt(ssp1.get_p(), ssp2.get_p());
    }
};
static order_ssp_t order_ssp;


void
siqs_small_sieve_data::small_sieve_info(const char * what, int side) const
{
    int nnice = ssps.size();
    int npow2 = 0;
    int npattern = 0;

    for(auto const & sp : ssp) {
        npow2 += sp.is_pow2();
        npattern += sp.is_pattern_sieved();
        nnice += sp.is_nice();
    }

    verbose_output_start_batch();
    verbose_fmt_print(0, 3, "# {}(side {}): {} nice primes", what, side, nnice);
    if (npow2) verbose_fmt_print(0, 3, ", {} powers of 2", npow2);
    if (npattern) verbose_fmt_print(0, 3, ", {} pattern-sieved", npattern);
    verbose_fmt_print(0, 3, ".\n");
    verbose_fmt_print (0, 4, "# Dump of small sieve data:\n{}\n{}\n",
            join(ssp, "\n"),
            join(ssps, "\n"));
    verbose_output_end_batch();
}

void
siqs_small_sieve_data::small_sieve_clear()
{
    ssps.clear();
    ssp.clear();
    ssdpos_many.clear();
    ssdpos_many_next.clear();
}

void
siqs_small_sieve_data::small_sieve_init(
        std::vector<fb_entry_general> const & resieved,
        std::vector<fb_entry_general> const & rest,
        int /* logI */,
        int /* side */,
        fb_factorbase::key_type const & factorbaseK,
        siqs_special_q_data const & Q,
        double scale)
{
    const unsigned int thresh = factorbaseK.thresholds[0];
    where_am_I w MAYBE_UNUSED;

    fbK = factorbaseK;

    // This zeroes out all vectors, but keeps storage around nevertheless
    small_sieve_clear();

    ssps.reserve(resieved.size() + rest.size());

    for (auto const & c : { &resieved, &rest }) {
        if (c == &rest)
            resieve_end_offset = ssps.size();
        for (auto const & e : *c) {
            /* p=pp^k, the prime or prime power in this entry, and pp is prime*/
            const fbprime_t p = e.q, pp = e.p;
            WHERE_AM_I_UPDATE(w, p, p);

            ASSERT_ALWAYS(p <= thresh);

            if (!Q.doing.is_coprime_to(pp)) {
                continue;
            }

            auto pdata = std::make_shared<siqs_ssp_pdata_t>(e, Q);

            for (int nr = 0; nr < e.nr_roots; nr++) {
                const fb_general_root &root = e.roots[nr];
                fb_root_p1 const Rab { root.r, root.proj };
                WHERE_AM_I_UPDATE(w, r, Rab.to_old_format(p));


                const unsigned char logp = fb_log_delta (pp, root.exp,
                                                             root.oldexp,
                                                             scale);
                fbroot_t const rp = pdata->transform_root(root.r, e.invq);
                siqs_ssp_t new_ssp(pdata, rp, logp);
                if (new_ssp.is_pow2() || new_ssp.is_pattern_sieved()) {
                    ssp.push_back(new_ssp);
                } else {
                    ssps.push_back(new_ssp);
                }
            }
        }
    }

    /* arrange so that the small_sieve() ctor is happy */
    /* I _think_ that normally, if the new code does its job correctly,
     * then this should be already sorted */
    ASSERT(std::is_sorted(ssps.begin(), ssps.begin() + resieve_end_offset));
    ASSERT(std::is_sorted(ssps.begin() + resieve_end_offset, ssps.end()));

    /* Sort general ssp vector in the order in which sieve2357::sieve expects
       them. siqs_small_sieve::do_pattern_sieve may drop some of these entries
       but preserves the ordering. */
    std::sort(ssp.begin(), ssp.end(), order_ssp);
}

void
siqs_small_sieve_data::small_sieve_activate_many_start_positions()
{
    std::swap(ssdpos_many, ssdpos_many_next);
}

void
siqs_small_sieve_data::small_sieve_prepare_many_start_positions(
        unsigned int first_region_index,
        int nregions,
        int logI,
        sublat_t const &)
{
    /* nregions should be a positive power of 2 */
    ASSERT_ALWAYS(nregions > 0 && !(nregions & (nregions - 1u)));
    siqs_small_sieve_base const C(logI, first_region_index);

    auto & res(ssdpos_many_next);
    res.clear();

    /* 1 if logI <= logB else 2^(logI-logB) */
    unsigned int regions_per_line = 1u << C.log_regions_per_line;

    res.assign(nregions+regions_per_line,
               std::vector<siqs_pos_t>(ssps.size(), 0));

    if (ssdpos_many.empty()) {
        auto & ssdpos = res.front();
        ASSERT_ALWAYS(C.j0 == 0);
        ASSERT_ALWAYS(C.i0 == -(1 << (logI-1)));
        for(size_t s = 0; auto const & ssp: ssps) {
            ssdpos[s] = ssp.first_position_first_line(logI);
            s++;
        }
        /* only for regions_per_line > 1, i.e., logI > LOG_BUCKET_REGION */
        for(unsigned idx = 1; idx < regions_per_line; ++idx) {
            /* complete this row */
            auto const & prev = res[idx-1];
            auto & cur = res[idx];
            for(size_t s = 0; auto const & ssp: ssps) {
                cur[s] = addmod_u32(prev[s], ssp.get_offset(), ssp.get_p());
                s++;
            }
        }
    } else {
        ASSERT_ALWAYS(ssdpos_many.size() >= regions_per_line);
        for(unsigned idx = 0; idx < regions_per_line; ++idx) {
            std::swap(ssdpos_many[ssdpos_many.size()-regions_per_line+idx],
                      res[idx]);
        }
    }
    ASSERT(res.front().size() == ssps.size());

    for(unsigned int k = regions_per_line; k < nregions + regions_per_line;
                                                k += regions_per_line) {
        /* infer from previous bucket region */
        ASSERT(res[k].size() == ssps.size());
        auto & prev_ssdpos = res[k-regions_per_line];
        auto & ssdpos = res[k];
        //XXX beware that k can be J here
        siqs_small_sieve_base Ct(logI, first_region_index+k);
        for(size_t s = 0; auto const & ssp: ssps) {
            ssdpos[s] = Ct.first_position_in_region(ssp, prev_ssdpos[s]);
            s++;
        }

        /* only for regions_per_line > 1, i.e., logI > LOG_BUCKET_REGION */
        for(unsigned idx = 1; idx < regions_per_line; ++idx) {
            /* complete this row */
            auto const & prev = res[k+idx-1];
            auto & cur = res[k+idx];
            for(size_t s = 0; auto const & ssp: ssps) {
                cur[s] = addmod_u32(prev[s], ssp.get_offset(), ssp.get_p());
                s++;
            }
        }
    }

#ifndef NDEBUG
    for(unsigned int k = 0; k < nregions + regions_per_line; k++) {
        ASSERT(res[k].size() == ssps.size());
        int N = first_region_index + k;
        siqs_small_sieve_base Ct(logI, N);
        for(size_t s = 0 ; s < ssps.size(); ++s) {
            ASSERT(res[k][s] == Ct.first_position_in_region(ssps[s]));
        }
    }
#endif
}

void
siqs_small_sieve_data::sieve_small_bucket_region(
        unsigned char *S,
        unsigned int N,
        int bucket_relative_index,
        int logI,
        sublat_t const &,
        where_am_I & w) const
{
    std::vector<siqs_pos_t> const & ssdpos = ssdpos_many[bucket_relative_index];
    siqs_small_sieve SS(ssdpos, ssps, ssp, S, logI, N);
    SS.do_pattern_sieve(w);
    SS.exceptional_sieve(w);
    SS.normal_sieve(w);
}

void
siqs_small_sieve_data::resieve_small_bucket_region(
        bucket_primes_t *BP,
        unsigned char *S,
        unsigned int N,
        int bucket_relative_index,
        int logI,
        sublat_t const &,
        where_am_I & w MAYBE_UNUSED)
{
    auto const & ssdpos = ssdpos_many[bucket_relative_index];
    siqs_small_sieve_base const C(logI, N);

    unsigned char *S_ptr;
    const int resieve_very_verbose = 0;

    for(size_t index = 0 ; index < resieve_end_offset ; index++) {
        auto const & ssps_cur(ssps[index]);

        const fbprime_t p = ssps_cur.get_p();
        fbprime_t const r = ssps_cur.get_r();
        WHERE_AM_I_UPDATE(w, p, p);
        siqs_pos_t pos = ssdpos[index];
        S_ptr = S;
        ASSERT(pos < p);

        for (unsigned int j = C.j0; j < C.j1; ++j) {
            WHERE_AM_I_UPDATE(w, j, j);
            pos = C.first_position_in_line(ssps_cur, pos, j);
            for (siqs_pos_t i = pos; i < C.F ; i += p) {
                if (LIKELY(S_ptr[i] == 255)) continue;
                bucket_update_t<1, primehint_t> prime;
                unsigned int const x = ((size_t) (j - C.j0) << logI) + i;
                if (resieve_very_verbose) {
                    verbose_fmt_print(0, 1, "resieve_small_bucket_region:"
                            " root {},{} divides at x = " "{} = {} * {} + {}",
                            p, r, x, j, C.I, i);
                }
                prime.p = p;
                prime.x = x;
                ASSERT(prime.p >= fbK.td_thresh);
                BP->push_update(prime);
            }
            S_ptr += C.I;
        }
    }

    /* "not nice" cases are power of two. As we obviously won't resieve
     * powers of two, nothing to be done for ssp.
     */
}

/* Pattern-sieve primes with the is_pattern_sieved flag */
void siqs_small_sieve::do_pattern_sieve(where_am_I & w MAYBE_UNUSED)
{
    sieve2357base::prime_t psp[not_nice_primes.size() + 1];

    for (unsigned int j = j0; j < j1; ++j) {
        const size_t x0 = (size_t) (j - j0) << logI;
        WHERE_AM_I_UPDATE(w, j, (j - j0));
#ifdef TRACE_K
        unsigned char orig_Sx = 0;
        if (trace_on_range_Nx(N, x0, x0 + F)) {
            orig_Sx = S[trace_Nx.x];
        }
#endif
        size_t l = 0;
        for (auto const & ssp : not_nice_primes) {
            if (!ssp.is_pattern_sieved()) {
                continue; /* Nothing to do here */
            }
            WHERE_AM_I_UPDATE(w, p, ssp.get_p());
            const fbprime_t pos =
                (j == j0) ? first_position_in_region(ssp)
                          : first_position_in_line(ssp, psp[l].idx, j);
            psp[l++] = {ssp.get_p(), pos, ssp.get_logp()};

#ifdef TRACE_K
            if (trace_on_range_Nx(N, x0, x0 + F)) {
                /* We are in the correct line (fragment). */
                const fbprime_t q = psp[l - 1].q, pos = psp[l - 1].idx;
                const unsigned char logp = psp[l - 1].logp;
                if (0) {
                    fmt::print("# Pattern sieve side {}, line {}"
                            " (N={}, x0={}, trace_Nx.x={}):"
                            " Adding psp[{}] = {{{}, {}, {}}}, from  {}\n",
                        w->side, j, N, x0, trace_Nx.x, l - 1, q, pos, logp, ssp);
                }
                if ((trace_Nx.x - x0) % q == pos) {
                    WHERE_AM_I_UPDATE(w, x, trace_Nx.x);
                    sieve_increase(S + trace_Nx.x, logp, w);
                }
            }
#endif
        }
#ifdef TRACE_K
        unsigned char new_Sx = 0;
        if (trace_on_range_Nx(N, x0, x0 + F)) {
            new_Sx = S[trace_Nx.x];
            S[trace_Nx.x] = orig_Sx;
        }
#endif
        psp[l++] = {0, 0, 0};
        preferred_sieve2357::sieve(
            (preferred_simd_type *) (S + x0), F, psp,
            0, sieve2357base::update_add, w);
#ifdef TRACE_K
        if (trace_on_range_Nx(N, x0, x0 + F)) {
            ASSERT_ALWAYS(new_Sx == S[trace_Nx.x]);
        }
#endif
    }
}
