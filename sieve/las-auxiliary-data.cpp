#include "cado.h"
#include "las-auxiliary-data.hpp"
#include "las-threads.hpp"
#include "las-info.hpp"
#include <iomanip>

/* in las.cpp */
extern int sync_at_special_q;

void nfs_aux::thread_data::update_checksums(nfs_work::thread_data & tws)
{
    for(int side = 0 ; side < 2 ; side++)
        checksum_post_sieve[side].update(tws.sides[side].bucket_region, BUCKET_REGION);
}

void
sieve_checksum::update(const unsigned int other)
{
    unsigned long r;
    ularith_addmod_ul_ul(&r, checksum, other, checksum_prime);
    checksum = r;
}

void
sieve_checksum::update(const unsigned char *data, const size_t len)
{
    cxx_mpz mb;
    unsigned int new_checksum;

    mpz_import(mb, len, -1, sizeof(unsigned char), -1, 0, data);
    new_checksum = mpz_tdiv_ui(mb, checksum_prime);

    this->update(new_checksum);
}

nfs_aux::~nfs_aux()
{
    ASSERT_ALWAYS(!timer_special_q.running());

    if (!complete)
        return;

    for (auto & T : th) {
        rep.accumulate_and_clear(std::move(T.rep));
        timer_special_q += T.timer;
        for (int side = 0; side < 2; side++)
            checksum_post_sieve[side].update(T.checksum_post_sieve[side]);
    }

    verbose_output_start_batch();

    if (tdict::global_enable >= 2) {
        verbose_output_print (0, 1, "%s", timer_special_q.display().c_str());

        timetree_t::timer_data_type t = 0;
        for(auto const &c : timer_special_q.filter_by_category()) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << c.second;
            verbose_output_print (0, 1, "# %s: %s\n",
                    coarse_las_timers::explain(c.first).c_str(),
                    os.str().c_str());
            t += c.second;
        }
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << t;
        verbose_output_print (0, 1, "# total counted time: %s\n", os.str().c_str());
    }

    rep.display_survivor_counters();

    verbose_output_print(0, 2,
            "# Checksums over sieve region: "
            "after all sieving: %u, %u\n",
            checksum_post_sieve[0].get_checksum(),
            checksum_post_sieve[1].get_checksum());

    verbose_output_vfprint(0, 1, gmp_vfprintf,
            "# %lu %s\n",
            rep.reports,
            las.batch ? "survivor(s) saved" : "relation(s)"
            );

    if (las.suppress_duplicates)
        verbose_output_print(0, 1, "# number of eliminated duplicates: %lu\n", rep.duplicates);

    qt0 = seconds() - qt0;
    wct_qt0 = wct_seconds() - wct_qt0;

    /* qt0 also aggregates time that was spent while we were processing
     * stuff asynchronously. Other thread might have been processing
     * other special-qs during that time.
     */
    double qtts = qt0 - rep.tn[0] - rep.tn[1] - rep.ttf - rep.ttcof;
    if (rep.survivors.after_sieve != rep.survivors.not_both_even) {
        verbose_output_print(0, 1, "# Warning: found %ld hits with i,j both even (not a bug, but should be very rare)\n", rep.survivors.after_sieve - rep.survivors.not_both_even);
    }
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (dont_print_tally && las.number_of_threads_total() > 1) {
        verbose_output_print(0, 2, "# Time for this special-q: %1.4fs [tally available only in mono-thread] in %1.4f elapsed s\n", qt0, wct_qt0);
    } else {
        const char * fmt_always[2] = {
            /* for default, special-q-overlapping mode: */
            "# Time for this special-q: *%1.4fs",
            /* for legacy synchronous mode: */
            "# Time for this special-q: %1.4fs"
        };
        const char * fmt[2] = {
            /* for default, special-q-overlapping mode: */
                " [norm %1.4f+%1.4f, sieving *%1.4f"
                " (%1.4f + %1.4f + *%1.4f),"
                " factor %1.4f (%1.4f + %1.4f)]"
                " in %1.4f elapsed s"
                " [*: incl overlap time]",
            /* for legacy synchronous mode: */
                " [norm %1.4f+%1.4f, sieving %1.4f"
                " (%1.4f + %1.4f + %1.4f),"
                " factor %1.4f (%1.4f + %1.4f)]"
                " in %1.4f elapsed s"
        };
        verbose_output_print(0, 1, fmt_always[sync_at_special_q != 0], qt0);
        verbose_output_print(0, 2, fmt[sync_at_special_q != 0],
                rep.tn[0],
                rep.tn[1],
                qtts,
                rep.ttbuckets_fill,
                rep.ttbuckets_apply,
                qtts-rep.ttbuckets_fill-rep.ttbuckets_apply,
                rep.ttf + rep.ttcof, rep.ttf, rep.ttcof,
                wct_qt0);
        verbose_output_print(0, 1, "\n");
    }

    verbose_output_end_batch();
}

#ifndef DISABLE_TIMINGS
/* we really wish to have a single timing slot for all the instantiations
 * of fill_in_buckets_toplevel_wrapper */
tdict::slot tdict_slot_for_fibt("fill_in_buckets_toplevel");
tdict::slot tdict_slot_for_alloc_buckets("allocate_buckets");
tdict::slot tdict_slot_for_threads("multithreaded tasks");
tdict::slot_parametric tdict_slot_for_side("side ", "");
#endif

