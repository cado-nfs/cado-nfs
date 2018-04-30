#include "cado.h"
#include "las-auxiliary-data.hpp"
#include "las-threads.hpp"
#include "las-types.hpp"

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


/*  las_report_accumulate_threads_and_display
 * This function does three distinct things.
 *  - accumulates the timing reports for all threads into a collated report
 *  - display the per-sq timing relative to this report, and the given
 *    timing argument (in seconds).
 */
nfs_aux::~nfs_aux()
{
    /* If complete==false, then the processing of this special-q was
     * interrupted by an exception (buckets_are_full). We won't print
     * stuff in this case. Maybe we should still make some effort so that
     * we can print the lost time ?
     */

    timer_special_q.stop();

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

        double t = 0;
        for(auto const &c : timer_special_q.filter_by_category()) {
            verbose_output_print (0, 1, "# %s: %.2f\n",
                    coarse_las_timers::explain(c.first).c_str(),
                    c.second);
            t += c.second;
        }
        verbose_output_print (0, 1, "# total counted time: %.2f\n", t);
    }

    rep.display_survivor_counters();

    verbose_output_print(0, 2,
            "# Checksums over sieve region: "
            "after all sieving: %u, %u\n",
            checksum_post_sieve[0].get_checksum(),
            checksum_post_sieve[1].get_checksum());
    verbose_output_vfprint(0, 1, gmp_vfprintf,
            "# %lu %s for side-%d (%Zd,%Zd)\n",
            rep.reports,
            las.batch ? "survivor(s) saved" : "relation(s)",
            doing.side,
            (mpz_srcptr) doing.p,
            (mpz_srcptr) doing.r);

    if (las.suppress_duplicates)
        verbose_output_print(0, 1, "# number of eliminated duplicates: %lu\n", rep.duplicates);

    qt0 = seconds() - qt0;
    double qtts = qt0 - rep.tn[0] - rep.tn[1] - rep.ttf;
    if (rep.survivors.after_sieve != rep.survivors.not_both_even) {
        verbose_output_print(0, 1, "# Warning: found %ld hits with i,j both even (not a bug, but should be very rare)\n", rep.survivors.after_sieve - rep.survivors.not_both_even);
    }
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (dont_print_tally && las.nb_threads > 1) {
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [tally available only in mono-thread]\n", qt0);
    } else {
        //convert ttcof from microseconds to seconds
        rep.ttcof *= 0.000001;
        verbose_output_print(0, 1, "# Time for this special-q: %1.4fs [norm %1.4f+%1.4f, sieving %1.4f"
                " (%1.4f + %1.4f + %1.4f),"
                " factor %1.4f (%1.4f + %1.4f)]\n", qt0,
                rep.tn[0],
                rep.tn[1],
                qtts,
                rep.ttbuckets_fill,
                rep.ttbuckets_apply,
                qtts-rep.ttbuckets_fill-rep.ttbuckets_apply,
                rep.ttf, rep.ttf - rep.ttcof, rep.ttcof);
    }

    verbose_output_end_batch();
}

