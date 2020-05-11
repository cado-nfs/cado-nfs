#include "cado.h"

#include <cstdio>
#include <cstdarg>
#include <gmp.h>

#include "las-detached-cofac.hpp"
#include "las-coordinates.hpp"
#include "las-auxiliary-data.hpp"
#include "las-duplicate.hpp"
#include "las-globals.hpp"
#include "las-galois.hpp"


/*{{{ asynchronous cofactorization */
/* This is one input to the late cofactoring process (aka ECM). Here, we
 * mean the stuff that is done detached from the rest of the siever
 * stuff: we no longer care about purging buckets and so on, these may
 * safely be used for later work.
 */
cofac_standalone::cofac_standalone() : a(0), b(0) {/*{{{*/
    S[0] = S[1] = 0;
#ifdef SUPPORT_LARGE_Q
    mpz_set_ui(az, 0);
    mpz_set_ui(bz, 0);
#endif
}/*}}}*/
cofac_standalone::cofac_standalone(int N, size_t x, int logI, qlattice_basis const & Q) {/*{{{*/
    S[0] = S[1] = 0;
    NxToAB (a, b, N, x, logI, Q);
#ifdef SUPPORT_LARGE_Q
    NxToABmpz (az, bz, N, x, logI, Q);
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
            wc.strategies);
}/*}}}*/

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int) /* {{{ */
{
    auto clean_param = call_dtor([_param]() { delete _param; });

    /* We must exit by cleaning the param structure we've been given. But
     * everything we do with the objects whose life is dependent on our
     * param structure must of course be completed at this point. This
     * holds as well for the timer. Yet ACTIVATE_TIMER below registers
     * some stuff to be done at dtor time, so it's important that we
     * clean up the parameters *after* the timer cleans up.
     */
    detached_cofac_parameters *param = static_cast<detached_cofac_parameters *>(_param);

    /* Import some contextual stuff. */
    int id = worker->rank();
    nfs_work_cofac & wc(*param->wc_p);
    nfs_aux & aux(*param->aux_p);
    nfs_aux::thread_data & taux(aux.th[id]);
    las_info const & las(wc.las);
    las_report & rep(taux.rep);
    timetree_t & timer(aux.get_timer(worker));
    /* The timer is normally not running, as we're in a thread task.
     * However, in descent mode, this is called synchronously, and then
     * the situation is different since the timer has already been
     * activated above.
     */
#ifndef DLP_DESCENT
    ENTER_THREAD_TIMER(timer);
#else
    CHILD_TIMER(timer, __func__);
#endif
    nfs_aux::rel_hash_t& rel_hash(aux.get_rel_hash());

    cofac_standalone & cur(*param);

    std::array<int, 2> cof_bitsize {{ 0,0 }}; /* placate compiler */
    las.cofac_stats.call(cur.norm, cof_bitsize);

    SIBLING_TIMER(timer, "cofactoring"); // aka factor_both_leftover_norms
    TIMER_CATEGORY(timer, cofactoring_mixed());

    int pass = cur.factor_both_leftover_norms(wc);
    rep.survivors.cofactored += (pass != 0);

    auto res = new detached_cofac_result;

    if (cur.trace_on_spot() && pass == 0) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# factor_both_leftover_norm failed for (%" PRId64 ",%" PRIu64 "), ", cur.a, cur.b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "remains %Zd, %Zd unfactored\n",
                (mpz_srcptr) cur.norm[0],
                (mpz_srcptr) cur.norm[1]);
    }
    if (pass <= 0) {
        /* a factor was > 2^lpb, or some
           factorization was incomplete */
        return res;
    }

    rep.survivors.smooth++;

    /* yippee: we found a relation! */
    SIBLING_TIMER(timer, "print relations");
    TIMER_CATEGORY(timer, bookkeeping());

    las.cofac_stats.success(cof_bitsize);

    relation rel = cur.get_relation(aux.doing);

    if (cur.trace_on_spot()) {
        verbose_output_print(TRACE_CHANNEL, 0, "# Relation for (%"
                PRId64 ",%" PRIu64 ") printed\n", cur.a, cur.b);
    }

    {
        int do_check = las.suppress_duplicates;

        /* note that if we have large primes which don't fit in
         * an unsigned long, then the duplicate check will
         * quickly return "no".
         */

        nfs_aux::abpair_t ab(cur.a, cur.b);
        bool is_new_rel;
        {
            std::lock_guard<std::mutex> foo(rel_hash.mutex());
            is_new_rel = rel_hash.insert(ab).second;
        }

        const char * dup_comment = NULL;

        if (do_check && relation_is_duplicate(rel, wc.doing, wc.las)) {
            dup_comment = "# DUPE ";
        } else {
            if (!is_new_rel) {
                /* the relation was already printed in a prior attempt,
                 * that was aborted because of an exception. */
                dup_comment = "# DUP ";
            }
            /* Even if the relation was already printed, the las_report
             * object is (now) specific to the current attempt, and has
             * no memory of the number of reports for the failed
             * attempt. Therefore we need to count the relation as a
             * report no matter what.
             */
            rep.reports ++;
            /* Not clear what gives when we have Galois relations.  */
        }

        /* In some rare cases, the norm on one side is exactly 1, which
         * creates undefined behaviour later on. (bug # 21707) */
        if (rel.nb_polys > 2)
            for (int i = 0; i < rel.nb_polys; ++i)
                if (rel.sides[i].size() == 0)
                    dup_comment = "# NORM1 ";

        if (!dup_comment) dup_comment = "";

        std::ostringstream os;

        if (prepend_relation_time)
            os << "(" << seconds() - tt_qstart << ") ";

        // verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n", i, j, cur.S[0], cur.S[1]);

        os << dup_comment << rel << "\n";

        if(las.galois != NULL) {
            // adding relations on the fly in Galois cases
            // once filtering is ok for all Galois cases, 
            // this entire block would have to disappear
            add_relations_with_galois(las.galois, os, dup_comment,
                    &rep.reports, rel);
        }

        /* print all in one go */
        verbose_output_start_batch();     /* unlock I/O */
        verbose_output_print(0, 1, "%s", os.str().c_str());
        verbose_output_end_batch();     /* unlock I/O */
#ifdef DLP_DESCENT
        res->rel_p = std::make_shared<relation>(std::move(rel));
#endif
    }

    /* Build histogram of lucky S[x] values */
    rep.mark_report(cur.S[0], cur.S[1]);

    return (task_result*) res;
}

/* }}} */
/*}}}*/

